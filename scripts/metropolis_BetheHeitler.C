#include <TH2D.h> 
#include <TH1D.h> 
#include <random> 
#include <cmath> 
#include <iostream> 
#include <functional> 
#include <stdio.h> 
#include <TCanvas.h> 
#include <TStyle.h> 
#include <TRandom3.h> 
#include <TVector3.h> 
#include <ApexOptics.h>
#include <TBox.h>  

using ApexOptics::Trajectory_t; 

using namespace std; 

/// @brief Basic Lortenz 4-vector. metric is (+---). 
struct LortenzVector_t {

    double t, x,y,z; 

    inline double Mag2() const { return t*t - (x*x + y*y* + z*z); }

    inline double& operator[](int i) { 
        switch (i) {
            case 0 : return t; 
            case 1 : return x; 
            case 2 : return y; 
            case 3 : return z; 
            default : std::fprintf(stderr, "in <LortenzVector_t::%s>: invalid index given: %i\n", __func__, i); return t; 
        }
    }

    double operator*(const LortenzVector_t& r) const { return t*r.t - (x*r.x + y*r.y + z*r.z); } 
};

class BetheHeitlerGenerator {
private: 
    
    double fBeam_energy = 2200.;                 /// beam energy, MeV
    static constexpr double fMass_e  = 0.501;   /// e+/e- rest mass, MeV/c^2 
    double fE_electron;                         /// electron/positrion energy in A' rest frame, MeV
    
    /// ---- geometrical / momentum acceptance of both detectors
    double fRange_L_x_sv[2] = {-85.e-3, +85.e-3}; 
    double fRange_L_y_sv[2] = {-70.e-3, +80.e-3}; 

    double fRange_R_x_sv[2] = {-85.e-3, +85.e-3}; 
    double fRange_R_y_sv[2] = {-80.e-3, +70.e-3}; 

    double fRange_dp  = +0.06; 

    static constexpr double pi = 3.14159265359; 
    static constexpr double sqrt_2 = 1.41421356237; 

    static constexpr int fN_foils = 10;         /// number of production foils
    double fSpectrometer_p0 = 1104.;            /// Central momentum of spectrometer, MeV/c 

    const double x_sweep_range      = 1.e-6;   //maximum change in x = |K|/E0, where |K| is photon energy    
    const double theta_sweep_range  = 0.05e-3;  //maximum change in outgoing photon theta
    const double phi_sweep_range    = pi/4.;    //maximum change in angle between electron theta & photon theta
    const double foil_jump_prob     = 0.20;      //probablility to 'jump' from the current foil to another foil
    
    const double rotation_RMS_Pbeam = 0.100/std::sqrt(3.);  /// proportional to beam electron's angle-change in radians. 
    const double rotation_RMS       = 0.020/std::sqrt(3.);  /// proportional to positron/electron's angle-change in radians. 
    
    //minimum value of x_photon
    const double fMin_x_photon = 0.90;
    
    //max angle between brem. photon and beam
    const double fMax_theta_photon = 5.*(pi/180.);  
    
    //uniform random number generator
    std::uniform_real_distribution<double> fRandGen_uniform; 

    //normally-distributed number generator
    std::normal_distribution<double> fRandGen_gaus; 
    
    //pseudorandom engine
    std::mt19937 fTwister; 

    /// @return Uniformly-distributed random number on the interval [0, 1)
    inline double Rand() { return fRandGen_uniform(fTwister); }

    /// @return Normally-distributed random number with mean=0 and RMS=1 
    inline double Gaus() { return fRandGen_gaus(fTwister); }

    struct State_t {

        TVector3 vertex;              /// react vertex, HCS (m)
        double x_photon=1.;           /// E_photon / E_beam, for the bremstrahlung photon (1st B-H process)  
        double theta_photon=0.;       /// Angle between brem. photon and beam, lab frame (rad)
        double phi_photon=0.;         /// azimuthal angle of brem. photon about beam (rad)
        
        TVector3 P1;                  /// 3-momentum of electron after bremstrahlung (lab-frame)
        
        TVector3 Pp;                  /// 3-momentum of the positron, HCS (MeV)
        TVector3 Pe;                  /// 3-momentum of the electron, HCS (MeV)
        
        int foil_num=4;               /// the number of the production foil we're generating from 

        double amplitude=0.;          /// 'Amplitude' of the PDF(x,theta,P1) at current values of x,theta
        bool in_acceptance=false;     /// is this phase-point in the acceptance?  
    }; 

    State_t fState; /// current state of the system

    /// @param  i_foil the index of production foil to use
    /// @return z-position of production foil
    double Get_foil_z(int i_foil) { return -238.8e-3  + ((double)i_foil)*55e-3; }

    /// perform random rotation of given TVector3, with given RMS rotation magnitude. preserve magnitude of vector.  
    TVector3 Random_rotation(const TVector3& v, double RMS_rotation_mag); 

public: 

    BetheHeitlerGenerator(double _beam_E=2200., double _spec_p0=1104., double _x=0.999); 

    ~BetheHeitlerGenerator() {};

    void SetRange_x_sv(bool is_RHRS, double min, double max);
    void SetRange_y_sv(bool is_RHRS, double min, double max);
    
    void SetRange_p0(double _p0) { fRange_dp=_p0; }

    /// perform metropolis update 
    void Update(); 

    /// @return true if current state is inside acceptance
    inline bool Inside_Acceptance() const { return fState.in_acceptance; }

    /// @return event vertex, in HCS (m)
    inline TVector3 GetVertex() const { return fState.vertex; }; 

    /// @return electron 3-momentum after bremstrahlung (lab frame)
    inline TVector3 Get_P1() const { return fState.P1; }

    /// @return photon 3-momentum (lab frame)
    TVector3 Get_K() const; 
    

    /// @brief Compute the differential Bethe-Heitler production amplitude
    /// @param P 3-momentum of electorn, post bremstrahlung (lab frame)
    /// @param P0 3-momentum of electron, pre-bremstrahlung (lab frame)
    /// @param K 3-momentum of produced photon (lab frame)
    double BetheHeitler_amplitude(const TVector3& P, const TVector3& P0, const TVector3& K) const; 
}; 

namespace {
    constexpr double beam_energy = 2200.;    /// beam energy, GeV
    constexpr double m_A         =  320.;    /// Rest mass of A', GeV/c2
    constexpr double m_e         =    0.501; /// e+/e- rest mass, GeV/c2 
    constexpr double E_electron = m_A/2.;    /// electron/positrion energy in A' rest frame, GeV
        
    constexpr double beam_E2 = beam_energy*beam_energy; 
    constexpr double m_A2 = m_A*m_A; 
    constexpr double m_e2 = m_e*m_e; 

    constexpr double pi = 3.14159265359; 
    constexpr double sqrt_2 = 1.41421356237; 

    constexpr int N_foils = 10;     /// number of production foils

    constexpr double spectrometer_p0 = 1104.; //GeV/c 
    
    uniform_real_distribution<double> rand_generator(0., 1.); 

    //normally-distributed number generator
    normal_distribution<double> gaus_generator(0., 1.); 
    
    //used to pick a random foil
    uniform_int_distribution<int> random_foil(0, N_foils-1);

    mt19937 twister; 

    /// @return Uniformly-distributed random number on the interval [0, 1)
    inline double Rand() { return rand_generator(twister); }

    /// @return Normally-distributed random number with mean=0 and RMS=1 
    inline double Gaus() { return gaus_generator(twister); }

    /// @return The z-position of a randomly-chosen production foil (in meters) 
    double random_prod_foil_z() {
        return -238.8e-3  + ((double)random_foil(twister))*55e-3; 
    }

    /// @param  i_foil the index of production foil to use
    /// @return z-position of production foil
    double Get_foil_z(int i_foil) { return -238.8e-3  + ((double)i_foil)*55e-3; }

    //perform random rotation of given TVector3, with given RMS rotation magnitude. preserve magnitude of vector.  
    TVector3 Random_rotation(const TVector3& v, double RMS_rotation_mag) {

        //implementing: vR = v + v x R 
        //              vR_i = \delta_ij v_j   +  \epsilon_ijk R_j v_k    
        // where: 
        //      \delta_ij is the kroneker delta tensor, and 
        //      \epsilon_ijk is the antisymmetric tensor. 
        // and we use the fact that 'R' should be a small vector, so this is a first-order approximation w/r/t mag(R). 
        double vx(v[0]), vy(v[1]), vz(v[2]); 
        TVector3 vR( vx, vy, vz );  

        double mag = v.Mag(); 

        double RX = RMS_rotation_mag*Gaus();     
        double RY = RMS_rotation_mag*Gaus();
        double RZ = RMS_rotation_mag*Gaus();
        
        vR[0] +=  RY*vz - RZ*vy;  
        vR[1] +=  RZ*vx - RX*vz; 
        vR[2] +=  RX*vy - RY*vx; 

        return vR.Unit() * mag; 
    }   
}; 

#define VERBOSE false

int metropolis_BetheHeitler(const bool is_RHRS)
{
    double x_start = 1. - 1.e-5; 

    BetheHeitlerGenerator BH_generator(beam_energy, spectrometer_p0, x_start); 

    TH2D* x_and_theta  = new TH2D("h", "A' produced;x;#theta_{A}", 200, 1. - 16e-5, 1., 200, -1.5, +1.5); 

    TH2D* h_x_theta = new TH2D("h_momentum", ";E_{photon} / E_{beam};#theta_{photon}", 200, 1. - 2.5e-4, 1., 200, -1e-5, +0.001); 

    TH1D* h_x = new TH1D("h_theta", "E_{photon} / E_{beam};;", 200, 1. - 1e-4, 1.); 

    TH2D* h_dx_dy_hall = new TH2D("h_dx_dy", "Photon angle / energy;energy;angle", 200, -200,+200, 200, -200,+200); 

    TH2D* h_x_y = new TH2D(
        "h_x_y_sv", 
        Form("electrons on LHRS sieve face - m_{A} = %.1f MeV/c^{2};x_{sv} mm (lab-vertical);y_{sv} mm (lab-horizontal)",m_A), 
        200,-120.,+120., 
        200,-120.,+120.
    ); 


    const long int n_stats = 5e6;
    long int n_generated=0; 

    const long int max_starting_attempts = 1e6; 
    long int i_start=0; 

    printf("generating %li events...", n_stats); cout << flush; 
    while (n_generated < n_stats) {

        BH_generator.Update(); 

        if (!BH_generator.Inside_Acceptance()) { 
            if (++i_start > max_starting_attempts) break; 
            continue; 
        }
        
        auto K = BH_generator.Get_K(); 

        h_x_theta->Fill( K.Mag()/beam_energy, std::sqrt(K[0]*K[0] + K[1]*K[1])/K.z() );

        h_x->Fill( K.Mag()/beam_energy ); 

        TVector3 vertex = BH_generator.GetVertex(); 

        Trajectory_t traj_HCS{
            vertex.x() + (K.x()/K.z())*(0. - vertex.z()), 
            vertex.y() + (K.y()/K.z())*(0. - vertex.z()), 
            K.x()/K.z(), 
            K.y()/K.z(), 
            K.Mag()
        }; 

        auto traj_SCS = ApexOptics::HCS_to_SCS(is_RHRS, traj_HCS);

        h_x_y->Fill( traj_SCS.x*1e3, traj_SCS.y*1e3 ); 

        n_generated++; 
    }
    cout << "done" << endl; 

    new TCanvas; 
    h_x->Draw();
    //h_x_theta->Draw("colz");
    
    return 0; 

    new TCanvas; 
    gPad->SetLeftMargin(0.15); 
    gStyle->SetOptStat(0); 

    h_x_y->Draw("col"); 

    auto box = new TBox(-50,-30, +50,+30); 
    box->SetFillStyle(0); 
    box->Draw(); 
    
    return 0; 
}

//_________________________________________________________________________________________________
BetheHeitlerGenerator::BetheHeitlerGenerator(double _beam_E, double _spec_p0, double _x)
    : 
    fBeam_energy{_beam_E},
    fSpectrometer_p0{_spec_p0}
{
    std::random_device rd;

    fTwister = std::mt19937(rd()); 

    fState = State_t{
        .vertex             = TVector3(0., 0., Get_foil_z(4)), 
        .x_photon           = _x,
        .theta_photon       = 0.100e-3,
        .phi_photon         = 0.
    }; 

    fState.P1               = TVector3(0., 0., (1.-_x)*fBeam_energy);
    fState.Pp               = TVector3(_x*fBeam_energy*std::sin(-5.*3.14159256/180.)/2., 0., _x*fBeam_energy*std::cos(5.*3.14159256/180.)/2.); 
    fState.Pe               = TVector3(_x*fBeam_energy*std::sin(+5.*3.14159256/180.)/2., 0., _x*fBeam_energy*std::cos(5.*3.14159256/180.)/2.); 
    fState.foil_num         = 4; 
    fState.amplitude        = 0.; /// a negative amplitude indicates that this state is invalid  
    fState.in_acceptance    = true;
    

    //photon energy s
    fState.amplitude = BetheHeitler_amplitude(
        fState.P1, 
        TVector3(0., 0., fBeam_energy), 
        Get_K()
    ); 

    std::printf("amplitude: %f\n", fState.amplitude); 
}
//_________________________________________________________________________________________________
TVector3 BetheHeitlerGenerator::Random_rotation(const TVector3& v, double RMS_rotation_mag) 
{
    //implementing: vR = v + v x R 
    //              vR_i = \delta_ij v_j   +  \epsilon_ijk R_j v_k    
    // where: 
    //      \delta_ij is the kroneker delta tensor, and 
    //      \epsilon_ijk is the antisymmetric tensor. 
    // and we use the fact that 'R' should be a small vector, so this is a first-order approximation w/r/t mag(R). 
    double vx(v[0]), vy(v[1]), vz(v[2]); 
    TVector3 vR( vx, vy, vz );  

    double mag = v.Mag(); 

    double RX = RMS_rotation_mag*Gaus();     
    double RY = RMS_rotation_mag*Gaus();
    double RZ = RMS_rotation_mag*Gaus();
    
    vR[0] +=  RY*vz - RZ*vy;  
    vR[1] +=  RZ*vx - RX*vz; 
    vR[2] +=  RX*vy - RY*vx; 

    return vR.Unit() * mag; 
}   
//_________________________________________________________________________________________________
double BetheHeitlerGenerator::BetheHeitler_amplitude(const TVector3& P1, const TVector3& P0, const TVector3& K) const
{
    //compute some quantities we will need.
    // we assume that the incoming electron is relativistic, while the outgoing electron is not. 
    
    //magnitude of vectors
    double p1 = P1.Mag(); 
    double p0 = P0.Mag(); 
    double k  = K.Mag(); 

    // energy of particles 
    double e0 = std::sqrt( fMass_e*fMass_e + p0*p0 ); 
    double e1 = std::sqrt( fMass_e*fMass_e + p1*p1 ); 
    
    //four-vector dot product of K and P1,P0 
    double K_dot_P1 = k*e1 - (K*P1); 
    double K_dot_P0 = k*e0 - (K*P0); 

    //compute cos_phi (angle between (K,P0) and (K,P1) plane)
    TVector3 P0_x_K = P0.Cross(K); 
    TVector3 P1_x_K = P1.Cross(K); 
    double cos_phi  = (P1_x_K * P0_x_K) / ( P0_x_K.Mag() * P1_x_K.Mag() ); 

    //momentum transfer 
    double q2 = (P0 - P1 - K).Mag2(); 

    double cos_theta1 = K * P1 / (k*p1);
    double cos_theta0 = K * P0 / (k*p0);
    
    double sin_theta1 = std::sqrt( 1. - cos_theta1*cos_theta1 ); 
    double sin_theta0 = std::sqrt( 1. - cos_theta0*cos_theta0 ); 
    
    double amplitude=0.;  
        
    amplitude += (p1*p1*sin_theta1*sin_theta1)*( 4.*e0*e0 - q2 )/(K_dot_P1*K_dot_P1); 

    amplitude += (p0*p0*sin_theta0*sin_theta0)*( 4.*e1*e1 - q2 )/(K_dot_P0*K_dot_P0); 

    amplitude += ( 2.*k*k*( p1*p1*sin_theta1*sin_theta1 + p0*p0*sin_theta0*sin_theta0 ) 
            - 2.*p0*sin_theta0*p1*sin_theta1*cos_phi*( 4.*e1*e0 + 2*k*k - q2 ) )/(K_dot_P0*K_dot_P1); 

    amplitude *= p1*k/(p0 * q2*q2); 

    return amplitude; 
}
//_________________________________________________________________________________________________
void BetheHeitlerGenerator::SetRange_x_sv(bool is_RHRS, double min, double max)
{
    if (is_RHRS) { fRange_R_x_sv[0] =min; fRange_R_x_sv[1] =max; }
    else         { fRange_L_x_sv[0] =min; fRange_L_x_sv[1] =max; }
}
//_________________________________________________________________________________________________
void BetheHeitlerGenerator::SetRange_y_sv(bool is_RHRS, double min, double max)
{
    if (is_RHRS) { fRange_R_y_sv[0] =min; fRange_R_y_sv[1] =max; }
    else         { fRange_L_y_sv[0] =min; fRange_L_y_sv[1] =max; }
}
//_________________________________________________________________________________________________
void BetheHeitlerGenerator::Update()
{
    auto pt = fState; 

    //set the 'x' value; x = E_photon / E_beam
    pt.x_photon     = min( pt.x_photon + (1. - 2.*Rand())*x_sweep_range, 1. ); 
    pt.x_photon     = max( pt.x_photon, fMin_x_photon ); 

    //set the angle between the photon and the beam
    pt.theta_photon += (1. - 2.*Rand())*theta_sweep_range; 

    if (fabs(pt.theta_photon) > fMax_theta_photon) {
        pt.theta_photon = fMax_theta_photon * pt.theta_photon/fabs(pt.theta_photon); 
    }

    //this is the azimuthal angle of the photon
    pt.phi_photon += (1. - 2.*Rand())*phi_sweep_range; 

    double jump_r = Rand(); 
    if (jump_r < foil_jump_prob) {
        //decide whether to jump to the next or previous foil
        if (jump_r / foil_jump_prob < 0.5) { pt.foil_num += +1; } else { pt.foil_num += -1; } 
        //fix the new foil num if its out-of-range
        if (pt.foil_num < 0) pt.foil_num = 0; 
        if (pt.foil_num >= N_foils) pt.foil_num = N_foils-1; 
    }

    //rotate the direction of the beam electron
    pt.P1 = Random_rotation(pt.P1, rotation_RMS_Pbeam); 
    pt.P1 = pt.P1.Unit() * fBeam_energy * (1. - pt.x_photon); 

    pt.amplitude = BetheHeitler_amplitude(
        fState.P1, 
        TVector3(0., 0., fBeam_energy), 
        Get_K()
    );
    
    pt.vertex = TVector3( 
        2.e-3*( 1. - 2.*Rand() ), 
        2.e-3*( 1. - 2.*Rand() ), 
        Get_foil_z(pt.foil_num)
    );
    
    TVector3 K(     
        0., 
        fBeam_energy*pt.x_photon*std::sin(pt.theta_photon), 
        fBeam_energy*pt.x_photon*std::cos(pt.theta_photon) 
    ); 
    K.RotateZ( pt.phi_photon );
    pt.amplitude = BetheHeitler_amplitude( fState.P1, TVector3(0., 0., fBeam_energy), K ); 
    
    bool accepted=false; 
    if ( pt.amplitude > fState.amplitude || pt.amplitude/fState.amplitude > Rand() ) {

        fState = pt; 
        accepted=true; 
    }

#if VERBOSE

    printf(
        "~~~~~~~~~~~~~~~~~~~~~ phase space ~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "Photon:\n"
        "   x (k/E0):   1. - %.2e\n"
        "   theta_0     %+.3f mrad\n"
        "   phi_0       %+.3f rad\n"
        "   K           (%+.2f, %+.2f, %.2f) MeV/c\n"
        "   amplitude:  %.3e\n"
        "\n"
        "Beam electron\n"
        "   P           (%+.3f, %+.3f, %+.3f) MeV/c\n"
        "accepted? %s\n"
        "\n",
        1. - pt.x_photon,
        pt.theta_photon*1.e3, 
        pt.phi_photon, 
        K[0], K[1], K[2],
        pt.amplitude,
        pt.P1[0], pt.P1[1], pt.P1[2],
        accepted ? "yes" : "no"
    ); 

#endif 
    
    /*  
    //get a new electorn direction (in A' rest frame)
    pt.P_electron = Random_rotation(pt.P_electron, rotation_RMS);  

    double p_electron_mag = std::sqrt( fE_electron*fE_electron - fMass_e*fMass_e ); 

    pt.P_electron = pt.P_electron.Unit() * p_electron_mag; 

    double E_a = pt.x * fBeam_energy; 
    double p_A_mag = sqrt( (E_a*E_a) - fMass_Aprime*fMass_Aprime ); 

    TVector3 p_A(
        p_A_mag * cos(pt.theta) * cos(pt.phi), 
        p_A_mag * sin(pt.theta) * sin(pt.phi), 
        p_A_mag * cos(pt.theta)
    ); 

    //now, in the rest-frame of the A', we generate the directions of the positron / electron randomly
    double gamma_A = E_a / fMass_Aprime; 
    double beta_A  = sqrt( 1. - 1./(gamma_A*gamma_A) ); 

    //_________________________________________________________________________________________
    //take a vector in the rest-frame of the A', and boost it to the lab frame
    auto boost_and_rotate = [gamma_A, beta_A, &pt](const TVector3 &p, double E)
    {
        TVector3 p_boost(
            p[0], 
            p[1], 
            gamma_A * (p[2] + beta_A*E)
        ); 
        
        //so, now that we've boosed the momentum back to the lab, 
        //let's rotate this vector, so that it's boost is in line with the momentum of the A' 
        p_boost.RotateX( pt.theta );
        p_boost.RotateZ( pt.phi );
        
        return p_boost; 
    };
    //_________________________________________________________________________________________
    
    //randomly generate the vertex, in line with the production target geometry        
    pt.vertex = TVector3( 
        2.e-3*( 1. - 2.*Rand() ), 
        2.e-3*( 1. - 2.*Rand() ), 
        Get_foil_z(pt.foil_num)
    );
    
    //check the acceptances of the electron and positron
    pt.Pe_boosted = boost_and_rotate(pt.P_electron, fE_electron); 
    pt.Pp_boosted = boost_and_rotate(-1.*pt.P_electron, fE_electron); 

    pt.in_acceptance = true; 
    
    //reject update if the new electron is outside the momentum range
    if (fabs(pt.Pe_boosted.Mag() - fSpectrometer_p0)/fSpectrometer_p0 > fRange_dp) { pt.in_acceptance=false; } 
    if (fabs(pt.Pp_boosted.Mag() - fSpectrometer_p0)/fSpectrometer_p0 > fRange_dp) { pt.in_acceptance=false; } 

    //___________________________________________________________________________________________________
    auto check_acceptance = [&](const TVector3& p, bool is_RHRS)
    {   
        Trajectory_t traj_HCS{
            pt.vertex.x() + (p.x()/p.z())*(0. - pt.vertex.z()), 
            pt.vertex.y() + (p.y()/p.z())*(0. - pt.vertex.z()), 
            p.x()/p.z(), 
            p.y()/p.z(), 
            p.Mag()
        }; 

        auto traj_SCS = ApexOptics::HCS_to_SCS(is_RHRS, traj_HCS);
        
        double xmin = is_RHRS ? fRange_R_x_sv[0] : fRange_L_x_sv[0]; 
        double xmax = is_RHRS ? fRange_R_x_sv[1] : fRange_L_x_sv[1]; 

        if (traj_SCS.x > xmax || traj_SCS.x < xmin) return false;   
        
        double ymin = is_RHRS ? fRange_R_y_sv[0] : fRange_L_y_sv[0]; 
        double ymax = is_RHRS ? fRange_R_y_sv[1] : fRange_L_y_sv[1];     

        if (traj_SCS.y > ymax || traj_SCS.y < ymin) return false; 

        return true; 
    };
    //___________________________________________________________________________________________________

    if (check_acceptance(pt.Pp_boosted, true)  == false) { pt.in_acceptance=false; } 
    if (check_acceptance(pt.Pe_boosted, false) == false) { pt.in_acceptance=false; } 
    
#if VERBOSE
    Trajectory_t tj_e_HCS{
        pt.vertex.x() + (pt.Pe_boosted.x()/pt.Pe_boosted.z())*(0. - pt.vertex.z()), 
        pt.vertex.y() + (pt.Pe_boosted.y()/pt.Pe_boosted.z())*(0. - pt.vertex.z()), 
        pt.Pe_boosted.x()/pt.Pe_boosted.z(), 
        pt.Pe_boosted.y()/pt.Pe_boosted.z(), 
        pt.Pe_boosted.Mag()
    }; 

    auto tj_e_SCS = ApexOptics::HCS_to_SCS(false, tj_e_HCS);

    Trajectory_t tj_p_HCS{
        pt.vertex.x() + (pt.Pp_boosted.x()/pt.Pp_boosted.z())*(0. - pt.vertex.z()), 
        pt.vertex.y() + (pt.Pp_boosted.y()/pt.Pp_boosted.z())*(0. - pt.vertex.z()), 
        pt.Pp_boosted.x()/pt.Pp_boosted.z(), 
        pt.Pp_boosted.y()/pt.Pp_boosted.z(), 
        pt.Pp_boosted.Mag()
    }; 

    auto tj_p_SCS = ApexOptics::HCS_to_SCS(true, tj_p_HCS);

    printf(
        " ~~~~~~~~~~~~~~~~ phase-space: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "   x_hcs       %-3.4f mm\n"
        "   y_hcs       %-3.4f mm\n"
        "   foil #      %i\n"
        "   P[electron] - A' rest frame:  (%-.1f, %-.1f, %-.1f) MeV/c\n"
        "   A':\n"
        "       EA/E0   = 1.0 - %.2e\n"
        "       theta_A = %-.3f mrad\n"
        "       phi_A   = %-3.1f deg\n"
        "   Electron:\n"
        "       x_sv = %-.3f mm\n"
        "       y_sv = %-.3f mm\n"
        "   Positron:\n"
        "       x_sv = %-.3f mm\n"
        "       y_sv = %-.3f mm\n"
        " inside acceptance? %s\n", 
        pt.vertex.x()*1e3,
        pt.vertex.y()*1e3, 
        pt.foil_num,
        pt.P_electron.x(), pt.P_electron.y(), pt.P_electron.z(), 
        1. - pt.x, 
        pt.theta * 1.e3, 
        pt.phi * (180./pi), 
        tj_e_SCS.x*1e3, 
        tj_e_SCS.y*1e3, 
        tj_p_SCS.x*1e3, 
        tj_p_SCS.y*1e3, 
        pt.in_acceptance ? "yes" : "no"
    );
#endif 

    //reject update
    if ((!pt.in_acceptance) && (fState.in_acceptance)) { return; }
    
    pt.amplitude = A_production_amplitude(pt.x, pt.theta);  

    //if we've gotten here, we can consider an update. 
    if ( pt.amplitude > fState.amplitude || pt.amplitude/fState.amplitude > Rand() ) {

        fState = pt; 
    }
    return; */  
}; 
//_________________________________________________________________________________________________
TVector3 BetheHeitlerGenerator::Get_K() const 
{
    TVector3 K( 
        0., 
        fBeam_energy*fState.x_photon*std::sin(fState.theta_photon), 
        fBeam_energy*fState.x_photon*std::cos(fState.theta_photon) 
    ); 
    K.RotateZ( fState.phi_photon ); 

    return K; 
}
//_________________________________________________________________________________________________
//_________________________________________________________________________________________________
//_________________________________________________________________________________________________


