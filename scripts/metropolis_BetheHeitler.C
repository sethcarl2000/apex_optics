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

    const double foil_jump_prob     = 0.10;      //probablility to 'jump' from the current foil to another foil    
    const double rotation_RMS       = 0.020/sqrt(3.);    /// proportional to positron/electron's angle-change in radians. 
    const double inv_mass_sweep_range = 5.; /// sweep range of the invariant mass, in MeV
    const double inv_mass_range[2] = {110., 270}; 
    
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
        
        TVector3 Pp;                  /// 3-momentum of the positron, HCS (MeV)
        TVector3 Pe;                  /// 3-momentum of the electron, HCS (MeV)

        int foil_num=4;               /// the number of the production foil we're generating from 

        double amplitude=0.;          /// 'Amplitude' of the PDF(x,theta,P1) at current values of x,theta
        bool in_acceptance=false;     /// is this phase-point in the acceptance?  

        TVector3 P_electron;          /// 3-momentum of electron in rest frame of beam electron
        double   inv_mass;            /// invariant mass of virtual photon 
    }; 

    State_t fState; /// current state of the system

    /// @param  i_foil the index of production foil to use
    /// @return z-position of production foil
    double Get_foil_z(int i_foil) { return -238.8e-3  + ((double)i_foil)*55e-3; }

    /// perform random rotation of given TVector3, with given RMS rotation magnitude. preserve magnitude of vector.  
    TVector3 Random_rotation(const TVector3& v, double RMS_rotation_mag); 

public: 

    BetheHeitlerGenerator(double _beam_E=2200., double _spec_p0=1104.); 

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
    
    /// @return 3-momentum of positron (MeV/c)
    inline TVector3 Get_Pp() const { return fState.Pp; }

    /// @return 3-momentum of electron (MeV/c)
    inline TVector3 Get_Pe() const { return fState.Pe; }

    /// @return Invariant mass of e+e- pair 
    inline double Get_invariantMass() const { return fState.inv_mass; }

    /// @brief Compute the differential Bethe-Heitler production amplitude
    /// @param P 3-momentum of electorn, post bremstrahlung (lab frame)
    /// @param P0 3-momentum of electron, pre-bremstrahlung (lab frame)
    /// @param K 3-momentum of produced photon (lab frame)
    double BetheHeitler_bremsstrahlung_amplitude(const TVector3& P, const TVector3& P0, const TVector3& K) const; 

    /// @brief Compute the differential Bethe-Heitler production amplitude
    /// @param Pp 3-momentum of positron, (lab frame)
    /// @param Pe 3-momentum of electron, (lab frame)
    /// @param k energy of photon (lab frame)
    double BetheHeitler_pairprod_amplitude(const TVector3& Pp, const TVector3& Pm, const double k) const; 
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

    BetheHeitlerGenerator BH_generator(beam_energy, spectrometer_p0); 

    TH2D* x_and_theta  = new TH2D("h", "A' produced;x;#theta_{A}", 200, 1. - 16e-5, 1., 200, -1.5, +1.5); 

    TH2D* h_x_theta = new TH2D("h_momentum", ";E_{photon} / E_{beam};#theta_{photon}", 200, 1. - 2.5e-4, 1., 200, -1e-5, +0.001); 

    TH1D* h_dxdz = new TH1D("h_theta", "dx/dz;;", 200, -0.1, 0.1); 

    TH2D* h_dx_dy_hall = new TH2D("h_dx_dy", "Photon angle / energy;energy;angle", 200, -200,+200, 200, -200,+200); 

    TH1D* h_m = new TH1D("h_m", "Invariant mass;MeV/c^{2};", 200, 105., 275.); 

    TH2D* h_x_y = new TH2D(
        "h_x_y_sv", 
        "electron on sv-face;dx/dz_{sv};dy/dz_{sv}", 
        200,-0.1,+0.1, 
        200,-0.1,+0.1
    ); 


    const long int n_stats = 1e6;
    long int n_generated=0; 

    const long int max_starting_attempts = 10; 
    long int i_start=0; 

    printf("generating %li events...", n_stats); cout << flush; 
    while (n_generated < n_stats) {

        for (int i=0; i<10; i++) BH_generator.Update(); 

        if (!BH_generator.Inside_Acceptance()) { 
            if (++i_start > max_starting_attempts) break; 
            continue; 
        }
        
        auto Pe = BH_generator.Get_Pe(); 

        TVector3 vertex = BH_generator.GetVertex(); 
           
        Trajectory_t traj_HCS{
            vertex.x() + (Pe.x()/Pe.z())*(0. - vertex.z()), 
            vertex.y() + (Pe.y()/Pe.z())*(0. - vertex.z()), 
            Pe.x()/Pe.z(), 
            Pe.y()/Pe.z()
        };

        auto traj_SCS = ApexOptics::HCS_to_SCS(false, traj_HCS);
        h_x_y->Fill( traj_SCS.dxdz, traj_SCS.dydz ); 

        h_dxdz->Fill( traj_SCS.dxdz ); 

        h_m->Fill( BH_generator.Get_invariantMass() ); 

        n_generated++; 
    }
    cout << "done" << endl; 
    
    new TCanvas; 
    gStyle->SetPalette(kSunset); 
    gPad->SetLeftMargin(0.15); 
    gStyle->SetOptStat(0); 

    h_m->Draw(); 


    new TCanvas; 
    
    h_x_y->Draw("col"); 

    return 0; 

    auto box = new TBox(-50,-30, +50,+30); 
    box->SetFillStyle(0); 
    box->Draw(); 
    
    return 0; 
}

//_________________________________________________________________________________________________
BetheHeitlerGenerator::BetheHeitlerGenerator(double _beam_E, double _spec_p0)
    : 
    fBeam_energy{_beam_E},
    fSpectrometer_p0{_spec_p0}
{
    std::random_device rd;

    fTwister = std::mt19937(rd()); 

    fState = State_t{
        .vertex = TVector3(0., 0., Get_foil_z(4)), 
    }; 

    fState.Pp = TVector3(
        fBeam_energy*std::sin(-5.*3.14159256/180.)/2., 
        0., 
        fBeam_energy*std::cos(5.*3.14159256/180.)/2.
    ); 
    
    fState.Pe = TVector3(
        fBeam_energy*std::sin(+5.*3.14159256/180.)/2., 
        0., 
        fBeam_energy*std::cos(5.*3.14159256/180.)/2.
    ); 
    
    fState.foil_num         = 4; 
    fState.amplitude        = 0.; /// a negative amplitude indicates that this state is invalid  
    fState.in_acceptance    = true;

    fState.inv_mass = 200.;
    fState.P_electron = TVector3(fState.inv_mass/2., 0., 0.); 

    //photon energy
    fState.amplitude = BetheHeitler_pairprod_amplitude(
        fState.Pp, 
        fState.Pe, 
        fBeam_energy 
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
double BetheHeitlerGenerator::BetheHeitler_bremsstrahlung_amplitude(const TVector3& P1, const TVector3& P0, const TVector3& K) const
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

    amplitude *= sin_theta0*sin_theta1 * p1*k/(p0 * q2*q2); 

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

    double jump_r = Rand(); 
    if (jump_r < foil_jump_prob) {
        //decide whether to jump to the next or previous foil
        if (jump_r / foil_jump_prob < 0.5) { pt.foil_num += +1; } else { pt.foil_num += -1; } 
        //fix the new foil num if its out-of-range
        if (pt.foil_num < 0) pt.foil_num = 0; 
        if (pt.foil_num >= N_foils) pt.foil_num = N_foils-1; 
    }

    
    pt.vertex = TVector3( 
        2.e-3*( 1. - 2.*Rand() ), 
        2.e-3*( 1. - 2.*Rand() ), 
        Get_foil_z(pt.foil_num)
    );
    
    pt.inv_mass += inv_mass_sweep_range*(1. - 2.*Rand());
    if (pt.inv_mass < inv_mass_range[0]) pt.inv_mass = inv_mass_range[0];
    if (pt.inv_mass > inv_mass_range[1]) pt.inv_mass = inv_mass_range[1];
    
    pt.in_acceptance =false; 

    long long int max_tries = 2e3; 

    long long int i_try=0; 

    //_________________________________________________________________________________________
    //take a vector in the rest-frame of the A', and boost it to the lab frame
    auto boost_and_rotate = [&pt, this](const TVector3 &p, double E)
    {
        double gamma = fBeam_energy / pt.inv_mass; 
        double beta  = std::sqrt( 1. - 1./(gamma*gamma) ); 

        TVector3 p_boost(
            p[0], 
            p[1], 
            gamma * (p[2] + beta*E)
        ); 
        
        return p_boost; 
    };
    //_________________________________________________________________________________________

    //_________________________________________________________________________________________
    auto check_acceptance = [&](const TVector3& p, bool is_RHRS)
    {   
        Trajectory_t traj_HCS{
            pt.vertex.x() + (p.x()/p.z())*(0. - pt.vertex.z()), 
            pt.vertex.y() + (p.y()/p.z())*(0. - pt.vertex.z()), 
            p.x()/p.z(), 
            p.y()/p.z()
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
    //_________________________________________________________________________________________
    
    do {                
        pt.P_electron = Random_rotation(pt.P_electron, rotation_RMS); 
        pt.P_electron = pt.P_electron.Unit() * (pt.inv_mass/2.);
        
        //rotate the direction of the beam electron
        pt.Pe = boost_and_rotate( pt.P_electron, pt.inv_mass/2.); 
        pt.Pp = boost_and_rotate(-pt.P_electron, pt.inv_mass/2.);

        if (check_acceptance(pt.Pp, true) && check_acceptance(pt.Pe, false)) { 
            pt.in_acceptance =true; 
            break; 
        }

    } while (++i_try < max_tries); 
    
#if !VERBOSE 
    if ((pt.in_acceptance==false) && (fState.in_acceptance==true)) { return; }
#endif 
    
    pt.amplitude = BetheHeitler_pairprod_amplitude( pt.Pp, pt.Pe, fBeam_energy );

    bool accepted=false; 
    
    double old_amp = fState.amplitude; 
    if (pt.in_acceptance) {
        if ( pt.amplitude > fState.amplitude || pt.amplitude/fState.amplitude > Rand() ) {

            fState = pt; 
            accepted=true; 
        }
    }
#if VERBOSE

    printf(
        "~~~~~~~~~~~~~~~~~~~~~ phase space ~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "Positron:\n"
        "    P           (%+.2f, %+.2f, %.2f) MeV/c\n"
        "   |P|          %.2f MeV/c\n"
        "\n"
        "Electron:\n"
        "    P           (%+.2f, %+.2f, %.2f) MeV/c\n"
        "   |P|          %.2f MeV/c\n"
        "\n"
        "amplitude      %.4e\n"
        "in acceptance? %s\n"
        "P-accept:      %.5f\n"
        "accepted?      %s\n"
        "\n",
        pt.Pp[0], pt.Pp[1], pt.Pp[2],
        pt.Pp.Mag(),
        pt.Pe[0], pt.Pe[1], pt.Pe[2],
        pt.Pe.Mag(),
        pt.amplitude,
        pt.in_acceptance ? "yes" : "no", 
        pt.in_acceptance ? min( pt.amplitude / old_amp, 1. ) : 0., 
        (accepted && pt.in_acceptance) ? "yes" : "no"
    ); 

    if ((pt.in_acceptance==false) && (fState.in_acceptance==true)) { return; }
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
double BetheHeitlerGenerator::BetheHeitler_pairprod_amplitude(const TVector3& Pp, const TVector3& Pe, const double k) const
{
    //compute some quantities we will need.
    // we assume that the incoming electron is relativistic, while the outgoing electron is not. 
    
    //photon 3-momentum
    const TVector3 K( 0., 0., k );

    //magnitude of vectors
    double pp = Pp.Mag(); 
    double pe = Pe.Mag(); 
    
    // energy of particles 
    double ep = std::sqrt( fMass_e*fMass_e + pp*pp ); 
    double ee = std::sqrt( fMass_e*fMass_e + pe*pe ); 
    
    //four-vector dot product of K and P1,P0 
    double K_dot_Pp = k*ep - (K*Pp); 
    double K_dot_Pe = k*ee - (K*Pe); 

    //compute cos_phi (angle between (K,P0) and (K,P1) plane)
    TVector3 Pp_x_K = Pp.Cross(K); 
    TVector3 Pe_x_K = Pe.Cross(K); 
    double cos_phi  = (Pp_x_K * Pe_x_K) / ( Pp_x_K.Mag() * Pe_x_K.Mag() ); 

    //momentum transfer 
    double q2 = (K - Pp - Pe).Mag2(); 

    double cos_theta_p = K * Pp / (k*pp);
    double cos_theta_e = K * Pe / (k*pe);
    
    double sin_theta_p = std::sqrt( 1. - cos_theta_p*cos_theta_p ); 
    double sin_theta_e = std::sqrt( 1. - cos_theta_e*cos_theta_e ); 
    
    double amplitude=0.;  
        
    amplitude += (pp*pp*sin_theta_p*sin_theta_p)*( 4.*ee*ee - q2 )/(K_dot_Pp*K_dot_Pp); 

    amplitude += (pe*pe*sin_theta_e*sin_theta_e)*( 4.*ep*ep - q2 )/(K_dot_Pe*K_dot_Pe); 

    amplitude += 2.*(pe*sin_theta_e)*(pp*sin_theta_p)*cos_phi*( 4.*ee*ep + q2 )/(K_dot_Pe*K_dot_Pp); 

    amplitude += -2.*k*k*( pp*pp*sin_theta_p*sin_theta_p + pe*pe*sin_theta_e*sin_theta_e  
                            + 2.*(pe*sin_theta_e)*(pp*sin_theta_p)*cos_phi )/(K_dot_Pe*K_dot_Pp); 
    
    amplitude *= - sin_theta_e*sin_theta_p * pp*pe/(k * q2*q2); 

    return amplitude; 
}
//_________________________________________________________________________________________________
//_________________________________________________________________________________________________


