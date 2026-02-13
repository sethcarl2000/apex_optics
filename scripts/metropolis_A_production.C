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

class AprimeGenerator {
private: 
    
    double fBeam_energy = 2200;                 /// beam energy, MeV
    double fMass_Aprime;                        /// Rest mass of A', MeV/c^2
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
    double fSpectrometer_p0 = 1108;         /// Central momentum of spectrometer, MeV/c 

    const double x_sweep_range     = 2.0e-5;    //maximum change in x   
    const double theta_sweep_range = 0.300e-3;  //maximum change in A' theta
    const double phi_sweep_range   = pi/18.;    //maximum change in A' phi
    const double foil_jump_prob    = 0.20;      //probablility to 'jump' from the current foil to another foil
    const double rotation_RMS      = 0.020/std::sqrt(3.);  /// proportional to electron's in angle, in radians (in A' rest frame). 
    
    
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
        double x=1.;                  /// E_A / E_0 
        double theta=0.;              /// Angle between A' and beam, lab frame 
        double phi=0.;                /// azimuthal angle of A' around beam
        
        TVector3 P_electron;          /// momentum of electron in A' rest frame
        TVector3 Pe_boosted;          /// momentum of electorn in the lab-frame
        TVector3 Pp_boosted;          /// momentum of positron in the lab-frame
        int foil_num=4;               /// the number of the production foil we're generating from 

        double amplitude=0.;          /// 'Amplitude' of the PDF(x,theta) at current values of x,theta
        bool in_acceptance=false;     /// is this phase-point in the acceptance?  
    }; 

    State_t fState; /// current state of the system

    /// @param  i_foil the index of production foil to use
    /// @return z-position of production foil
    double Get_foil_z(int i_foil) { return -238.8e-3  + ((double)i_foil)*55e-3; }

    /// perform random rotation of given TVector3, with given RMS rotation magnitude. preserve magnitude of vector.  
    TVector3 Random_rotation(const TVector3& v, double RMS_rotation_mag); 

    //given the ratio x := Ea / E0 (A' energy to intial beam energy, in Lab), and 'theta_A' (angle between beam and A' in Lab), 
    //this returns the rel. differential cross section. 
    double A_production_amplitude(double x, double theta_A); 

public: 

    AprimeGenerator(double _m_A, double _beam_E=2200., double _spec_p0=1108.); 

    ~AprimeGenerator() {};

    void SetRange_x_sv(bool is_RHRS, double min, double max);
    void SetRange_y_sv(bool is_RHRS, double min, double max);
    
    void SetRange_p0(double _p0) { fRange_dp=_p0; }

    /// perform metropolis update 
    void Update(); 

    /// @return true if current state is inside acceptance
    inline bool Inside_Acceptance() const { return fState.in_acceptance; }

    /// @return event vertex, in HCS (m)
    inline TVector3 GetVertex() const { return fState.vertex; }; 

    /// @return 3-momentum of electron, in HCS (MeV/c)
    inline TVector3 GetPe() const { return fState.Pe_boosted; }; 

    /// @return 3-momentum of positron, in HCS (MeV/c)
    inline TVector3 GetPp() const { return fState.Pp_boosted; }; 
};


namespace {
    constexpr double beam_energy = 2.200;    /// beam energy, GeV
    constexpr double m_A         = 0.320;    /// Rest mass of A', GeV/c2
    constexpr double m_e         = 0.501e-3; /// e+/e- rest mass, GeV/c2 
    constexpr double E_electron = m_A/2.;    /// electron/positrion energy in A' rest frame, GeV
        
    constexpr double beam_E2 = beam_energy*beam_energy; 
    constexpr double m_A2 = m_A*m_A; 
    constexpr double m_e2 = m_e*m_e; 

    constexpr double pi = 3.14159265359; 
    constexpr double sqrt_2 = 1.41421356237; 

    constexpr int N_foils = 10;     /// number of production foils

    constexpr double spectrometer_p0 = 1.108; //GeV/c 
    
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

//given the ratio x := Ea / E0 (A' energy to intial beam energy, in Lab), and 'theta_A' (angle between beam and A' in Lab), 
//this returns the rel. differential cross section. 
double A_production_amplitude(double x, double theta_A) 
{
    if (x < 0. || x > 1.) return 0.;

    double th_A2 = theta_A*theta_A; 

    double U = beam_E2 * x * th_A2   +   m_A2*(1. - x)/x   +   m_e2*x; 

    return beam_E2*x*(  1. - x + (x*x/2.) - x*x*(1-x)*m_A2*beam_E2*th_A2/(U*U)  )/(U*U); 
}

int metropolis_A_production(const bool is_RHRS)
{
    printf("peak amplitude: %.3e\n", A_production_amplitude(1. - 1.e-5, 0.)); 
    
    const double p_electron_mag = sqrt(E_electron*E_electron - m_e2); 
    
    struct MetPhasePoint_t {

        double x_hcs = 0.;            /// react vertex x (HCS), m       
        double y_hcs = 0.;            /// react vertex y (HCS), m
        double x=1.;                  /// E_A / E_0 
        double theta=0.;              /// Angle between A' and beam, lab frame 
        double phi=0.;                /// azimuthal angle of A' around beam
        
        TVector3 P_electron;          /// momentum of electron in A' rest frame
        TVector3 Pe_boosted;          /// momentum of electorn in the lab-frame
        TVector3 Pp_boosted;          /// momentum of positron in the lab-frame
        int    foil_num=4;            /// the number of the production foil we're generating from 

        double amplitude=-1.;          /// 'Amplitude' of the PDF(x,theta) at current values of x,theta
        bool in_acceptance=false;
    }; 


    double x_hcs = 0.;            /// react vertex x (HCS), m       
    double y_hcs = 0.;            /// react vertex y (HCS), m
    double x=1.;                        /// E_A / E_0 
    double theta=0.;                    /// Angle between A' and beam, lab frame 
    double phi=0.;                      /// azimuthal angle of A' around beam
    double amplitude=0.;                /// 'Amplitude' of the PDF(x,theta) at current values of x,theta
    int    foil_num=4;                  /// the number of the production foil we're generating from 
    
    /// orientation of the electron in its own rest-frame 
    
    const double x_sweep_range     = 2.0e-5;    //maximum change in x   
    const double theta_sweep_range = 0.300e-3;  //maximum change in A' theta
    const double phi_sweep_range   = pi/4.;     //maximum change in A' phi
    const double foil_jump_prob    = 0.20;      //probablility to 'jump' from the current foil to another foil

    /// ---- geometrical / momentum acceptance of both detectors
    const double range_L_x_sv[] = {-85.e-3, +85.e-3}; 
    const double range_L_y_sv[] = {-70.e-3, +80.e-3}; 

    const double range_R_x_sv[] = {-85.e-3, +85.e-3}; 
    const double range_R_y_sv[] = {-80.e-3, +70.e-3}; 

    const double range_dp       = +0.06; 

    //we divide by sqrt(3) because we compose this random vector with 3 gaussian random vars.
    const double rotation_RMS   = 0.020/sqrt(3.);  /// proportional to electron's in angle, in radians (in A' rest frame). 
                             
    //initalize the twister
    random_device rd; 
    twister = mt19937(rd()); 
    
    /*/perform a metropolis update
    auto metropolis_update = [  is_RHRS,

                                x_sweep_range,      //metropolis random-walk parameters 
                                theta_sweep_range, 
                                phi_sweep_range, 
                                foil_jump_prob, 
                                rotation_RMS, 

                                range_R_x_sv, range_R_y_sv, //detector acceptances
                                range_L_x_sv, range_L_y_sv, 
                                range_dp, 

                                p_electron_mag  //electron momentum magnitude in A' rest frame
                            ](MetPhasePoint_t& point) 
    {

        auto new_point = point; 

        new_point.x     += (1. - 2.*Rand())*x_sweep_range; 
        new_point.theta += (1. - 2.*Rand())*theta_sweep_range; 
        new_point.phi   += (1. - 2.*Rand())*phi_sweep_range; 
        
        new_point.amplitude = A_production_amplitude(new_point.x, new_point.theta);  

        double jump_r = Rand(); 
        if (jump_r < foil_jump_prob) {
            //decide whether to jump to the next or previous foil
            if (jump_r / foil_jump_prob < 0.5) { new_point.foil_num += +1; } else { new_point.foil_num += -1; } 
            //fix the new foil num if its out-of-range
            if (new_point.foil_num < 0) new_point.foil_num = 0; 
            if (new_point.foil_num >= N_foils) new_point.foil_num = N_foils-1; 
        }
        
        //get a new electorn direction (in A' rest frame)
        new_point.P_electron = Random_rotation(new_point.P_electron, rotation_RMS); 

        new_point.P_electron = new_point.P_electron.Unit() * p_electron_mag; 

        double E_a = new_point.x * beam_energy; 
        double p_A_mag = sqrt( (E_a*E_a) - m_A2 ); 

        TVector3 p_A(
            p_A_mag * cos(new_point.theta) * cos(new_point.phi), 
            p_A_mag * sin(new_point.theta) * sin(new_point.phi), 
            p_A_mag * cos(new_point.theta)
        ); 

        //now, in the rest-frame of the A', we generate the directions of the positron / electron randomly
        double gamma_A = E_a / m_A; 
        double beta_A  = sqrt( 1. - 1./(gamma_A*gamma_A) ); 

        //_________________________________________________________________________________________
        //take a vector in the rest-frame of the A', and boost it to the lab frame
        auto boost_and_rotate = [gamma_A, beta_A, &new_point](const TVector3 &p, double E)
        {
            TVector3 p_boost(
                p[0], 
                p[1], 
                gamma_A * (p[2] + beta_A*E)
            ); 
            
            //so, now that we've boosed the momentum back to the lab, 
            //let's rotate this vector, so that it's boost is in line with the momentum of the A' 
            p_boost.RotateX( new_point.theta );
            p_boost.RotateZ( new_point.phi );
            
            return p_boost; 
        };
        //_________________________________________________________________________________________
        
        //randomly generate the vertex, in line with the production target geometry        
        new_point.x_hcs = 2.e-3*( 1. - 2.*Rand() ); 
        new_point.y_hcs = 2.e-3*( 1. - 2.*Rand() ); 

        TVector3 vertex(
            new_point.x_hcs, 
            new_point.y_hcs,
            Get_foil_z(new_point.foil_num)
        ); 

        //check the acceptances of the electron and positron
        new_point.Pe_boosted = boost_and_rotate(new_point.P_electron, E_electron); 
        new_point.Pp_boosted = boost_and_rotate(-1.*new_point.P_electron, E_electron); 

        new_point.in_acceptance = true; 
        
        //reject update if the new electron is outside the momentum range
        if (fabs(new_point.Pe_boosted.Mag() - spectrometer_p0)/spectrometer_p0 > range_dp) { new_point.in_acceptance=false; } 
        if (fabs(new_point.Pp_boosted.Mag() - spectrometer_p0)/spectrometer_p0 > range_dp) { new_point.in_acceptance=false; } 

        //___________________________________________________________________________________________________
        auto check_acceptance = [&](const TVector3& p, bool is_RHRS)
        {   
            Trajectory_t traj_HCS{
                vertex.x() + (p.x()/p.z())*(0. - vertex.z()), 
                vertex.y() + (p.y()/p.z())*(0. - vertex.z()), 
                p.x()/p.z(), 
                p.y()/p.z(), 
                p.Mag()
            }; 

            auto traj_SCS = ApexOptics::HCS_to_SCS(is_RHRS, traj_HCS);
            
            double xmin = is_RHRS ? range_R_x_sv[0] : range_L_x_sv[0]; 
            double xmax = is_RHRS ? range_R_x_sv[1] : range_L_x_sv[1]; 

            if (traj_SCS.x > xmax || traj_SCS.x < xmin) return false;   
         
            double ymin = is_RHRS ? range_R_y_sv[0] : range_L_y_sv[0]; 
            double ymax = is_RHRS ? range_R_y_sv[1] : range_L_y_sv[1];     

            if (traj_SCS.y > ymax || traj_SCS.y < ymin) return false; 

            return true; 
        };
        //___________________________________________________________________________________________________

        if (check_acceptance(new_point.Pp_boosted, true)  == false) { new_point.in_acceptance=false; } 
        if (check_acceptance(new_point.Pe_boosted, false) == false) { new_point.in_acceptance=false; } 

        
#if VERBOSE
        Trajectory_t tj_e_HCS{
            vertex.x() + (new_point.Pe_boosted.x()/new_point.Pe_boosted.z())*(0. - vertex.z()), 
            vertex.y() + (new_point.Pe_boosted.y()/new_point.Pe_boosted.z())*(0. - vertex.z()), 
            new_point.Pe_boosted.x()/new_point.Pe_boosted.z(), 
            new_point.Pe_boosted.y()/new_point.Pe_boosted.z(), 
            new_point.Pe_boosted.Mag()
        }; 

        auto tj_e_SCS = ApexOptics::HCS_to_SCS(false, tj_e_HCS);

        Trajectory_t tj_p_HCS{
            vertex.x() + (new_point.Pp_boosted.x()/new_point.Pp_boosted.z())*(0. - vertex.z()), 
            vertex.y() + (new_point.Pp_boosted.y()/new_point.Pp_boosted.z())*(0. - vertex.z()), 
            new_point.Pp_boosted.x()/new_point.Pp_boosted.z(), 
            new_point.Pp_boosted.y()/new_point.Pp_boosted.z(), 
            new_point.Pp_boosted.Mag()
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
            new_point.x_hcs*1e3,
            new_point.y_hcs*1e3, 
            new_point.foil_num,
            new_point.P_electron.x()*1e3, new_point.P_electron.y()*1e3, new_point.P_electron.z()*1e3, 
            1. - new_point.x, 
            new_point.theta, 
            new_point.phi * (180./pi), 
            tj_e_SCS.x*1e3, 
            tj_e_SCS.y*1e3, 
            tj_p_SCS.x*1e3, 
            tj_p_SCS.y*1e3, 
            new_point.in_acceptance ? "yes" : "no"
        );
#endif 

        //reject update
        if ((!new_point.in_acceptance) && (point.in_acceptance)) { return; }

        //if we've gotten here, we can consider an update. 
        if ( new_point.amplitude > point.amplitude || new_point.amplitude/point.amplitude > Rand() ) {

            point = new_point; 
        }
        return; 
    }; */ 

    AprimeGenerator Aprime_generator(320., 2200., 1108); 

    TH2D* x_and_theta  = new TH2D("h", "A' produced;x;#theta_{A}", 200, 1. - 16e-5, 1., 200, -1.5, +1.5); 

    TH1D* h_e_momentum = new TH1D("h_momentum", "electron / positron momentum;GeV/c;", 200, -1., 1.+beam_energy); 

    TH2D* h_dx_dy_hall = new TH2D("h_dx_dy", "electron/positron angle in lab;dx/dz (mrad);dy/dz (mrad)", 200, -200,+200, 200, -200,+200); 

    TH2D* h_x_y = new TH2D(
        "h_x_y_sv", 
        Form("electrons on LHRS sieve face - m_{A} = %.1f MeV/c^{2};x_{sv} mm (lab-vertical);y_{sv} mm (lab-horizontal)",m_A*1.e3), 
        200,-120.,+120., 
        200,-120.,+120.
    ); 

    MetPhasePoint_t phase_point{
        .x_hcs      = 0., 
        .y_hcs      = 0., 
        .x          = 1.,
        .theta      = 0.,
        .phi        = 0., 
        .P_electron = TVector3(+p_electron_mag, 0., 0.),
        .foil_num   = 4, 
        .amplitude  = A_production_amplitude(1., 0.) /// a negative amplitude indicates that this state is invalid  
    };

    const long int n_stats = 1e6;
    long int n_generated=0; 

    const long int max_starting_attempts = 1e6; 
    long int i_start=0; 

    printf("generating %li events...", n_stats); cout << flush; 
    while (n_generated < n_stats) {

        Aprime_generator.Update(); 
        
        if (!Aprime_generator.Inside_Acceptance()) { 
            if (++i_start > max_starting_attempts) break; 
            continue; 
        }
        //metropolis_update(phase_point); 

        TVector3 vertex = Aprime_generator.GetVertex(); 

        auto p = Aprime_generator.GetPe(); 

        Trajectory_t traj_HCS{
            vertex.x() + (p.x()/p.z())*(0. - vertex.z()), 
            vertex.y() + (p.y()/p.z())*(0. - vertex.z()), 
            p.x()/p.z(), 
            p.y()/p.z(), 
            p.Mag()
        }; 

        auto traj_SCS = ApexOptics::HCS_to_SCS(is_RHRS, traj_HCS);

        h_x_y->Fill( traj_SCS.x*1e3, traj_SCS.y*1e3 ); 

        n_generated++; 
    }
    cout << "done" << endl; 

    new TCanvas; 
    h_dx_dy_hall->Draw("col");

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
AprimeGenerator::AprimeGenerator(double _mass_A, double _beam_E, double _spec_p0)
    : fMass_Aprime{_mass_A},
    fBeam_energy{_beam_E},
    fSpectrometer_p0{_spec_p0},
    fE_electron{fMass_Aprime/2.}
{
    std::random_device rd;

    fTwister = std::mt19937(rd()); 

    fState = State_t{
        .vertex     = TVector3(0., 0., Get_foil_z(4)), 
        .x          = 1.,
        .theta      = 0.,
        .phi        = 0., 
        .P_electron = TVector3(+sqrt(fE_electron*fE_electron - fMass_e*fMass_e), 0., 0.),
        .foil_num   = 4, 
        .amplitude  = A_production_amplitude(1., 0.), /// a negative amplitude indicates that this state is invalid  
        .in_acceptance=false
    };
}
//_________________________________________________________________________________________________
TVector3 AprimeGenerator::Random_rotation(const TVector3& v, double RMS_rotation_mag) 
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
double AprimeGenerator::A_production_amplitude(double x, double theta_A)
{
    if (x < 0. || x > 1.) return 0.;

    double th_A2   = theta_A*theta_A; 
    double beam_E2 = fBeam_energy*fBeam_energy; 
    double m_A2    = fMass_Aprime*fMass_Aprime; 

    double U = beam_E2 * x * th_A2   +   m_A2*(1. - x)/x   +  (fMass_e*fMass_e)*x; 

    return beam_E2*x*(  1. - x + (x*x/2.) - x*x*(1-x)*m_A2*beam_E2*th_A2/(U*U)  )/(U*U); 
}
//_________________________________________________________________________________________________
void AprimeGenerator::SetRange_x_sv(bool is_RHRS, double min, double max)
{
    if (is_RHRS) { fRange_R_x_sv[0] =min; fRange_R_x_sv[1] =max; }
    else         { fRange_L_x_sv[0] =min; fRange_L_x_sv[1] =max; }
}
//_________________________________________________________________________________________________
void AprimeGenerator::SetRange_y_sv(bool is_RHRS, double min, double max)
{
    if (is_RHRS) { fRange_R_y_sv[0] =min; fRange_R_y_sv[1] =max; }
    else         { fRange_L_y_sv[0] =min; fRange_L_y_sv[1] =max; }
}
//_________________________________________________________________________________________________
void AprimeGenerator::Update()
{
    auto new_point = fState; 

    new_point.x     += (1. - 2.*Rand())*x_sweep_range; 
    new_point.theta += (1. - 2.*Rand())*theta_sweep_range; 
    new_point.phi   += (1. - 2.*Rand())*phi_sweep_range; 

    double jump_r = Rand(); 
    if (jump_r < foil_jump_prob) {
        //decide whether to jump to the next or previous foil
        if (jump_r / foil_jump_prob < 0.5) { new_point.foil_num += +1; } else { new_point.foil_num += -1; } 
        //fix the new foil num if its out-of-range
        if (new_point.foil_num < 0) new_point.foil_num = 0; 
        if (new_point.foil_num >= N_foils) new_point.foil_num = N_foils-1; 
    }
    
    //get a new electorn direction (in A' rest frame)
    new_point.P_electron = Random_rotation(new_point.P_electron, rotation_RMS);  

    double p_electron_mag = std::sqrt( fE_electron*fE_electron - fMass_e*fMass_e ); 

    new_point.P_electron = new_point.P_electron.Unit() * p_electron_mag; 

    double E_a = new_point.x * fBeam_energy; 
    double p_A_mag = sqrt( (E_a*E_a) - fMass_Aprime*fMass_Aprime ); 

    TVector3 p_A(
        p_A_mag * cos(new_point.theta) * cos(new_point.phi), 
        p_A_mag * sin(new_point.theta) * sin(new_point.phi), 
        p_A_mag * cos(new_point.theta)
    ); 

    //now, in the rest-frame of the A', we generate the directions of the positron / electron randomly
    double gamma_A = E_a / fMass_Aprime; 
    double beta_A  = sqrt( 1. - 1./(gamma_A*gamma_A) ); 

    //_________________________________________________________________________________________
    //take a vector in the rest-frame of the A', and boost it to the lab frame
    auto boost_and_rotate = [gamma_A, beta_A, &new_point](const TVector3 &p, double E)
    {
        TVector3 p_boost(
            p[0], 
            p[1], 
            gamma_A * (p[2] + beta_A*E)
        ); 
        
        //so, now that we've boosed the momentum back to the lab, 
        //let's rotate this vector, so that it's boost is in line with the momentum of the A' 
        p_boost.RotateX( new_point.theta );
        p_boost.RotateZ( new_point.phi );
        
        return p_boost; 
    };
    //_________________________________________________________________________________________
    
    //randomly generate the vertex, in line with the production target geometry        
    new_point.vertex = TVector3( 
        2.e-3*( 1. - 2.*Rand() ), 
        2.e-3*( 1. - 2.*Rand() ), 
        Get_foil_z(new_point.foil_num)
    );
    
    //check the acceptances of the electron and positron
    new_point.Pe_boosted = boost_and_rotate(new_point.P_electron, fE_electron); 
    new_point.Pp_boosted = boost_and_rotate(-1.*new_point.P_electron, fE_electron); 

    new_point.in_acceptance = true; 
    
    //reject update if the new electron is outside the momentum range
    if (fabs(new_point.Pe_boosted.Mag() - fSpectrometer_p0)/fSpectrometer_p0 > fRange_dp) { new_point.in_acceptance=false; } 
    if (fabs(new_point.Pp_boosted.Mag() - fSpectrometer_p0)/fSpectrometer_p0 > fRange_dp) { new_point.in_acceptance=false; } 

    //___________________________________________________________________________________________________
    auto check_acceptance = [&](const TVector3& p, bool is_RHRS)
    {   
        Trajectory_t traj_HCS{
            new_point.vertex.x() + (p.x()/p.z())*(0. - new_point.vertex.z()), 
            new_point.vertex.y() + (p.y()/p.z())*(0. - new_point.vertex.z()), 
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

    if (check_acceptance(new_point.Pp_boosted, true)  == false) { new_point.in_acceptance=false; } 
    if (check_acceptance(new_point.Pe_boosted, false) == false) { new_point.in_acceptance=false; } 
    
#if VERBOSE
    Trajectory_t tj_e_HCS{
        new_point.vertex.x() + (new_point.Pe_boosted.x()/new_point.Pe_boosted.z())*(0. - new_point.vertex.z()), 
        new_point.vertex.y() + (new_point.Pe_boosted.y()/new_point.Pe_boosted.z())*(0. - new_point.vertex.z()), 
        new_point.Pe_boosted.x()/new_point.Pe_boosted.z(), 
        new_point.Pe_boosted.y()/new_point.Pe_boosted.z(), 
        new_point.Pe_boosted.Mag()
    }; 

    auto tj_e_SCS = ApexOptics::HCS_to_SCS(false, tj_e_HCS);

    Trajectory_t tj_p_HCS{
        new_point.vertex.x() + (new_point.Pp_boosted.x()/new_point.Pp_boosted.z())*(0. - new_point.vertex.z()), 
        new_point.vertex.y() + (new_point.Pp_boosted.y()/new_point.Pp_boosted.z())*(0. - new_point.vertex.z()), 
        new_point.Pp_boosted.x()/new_point.Pp_boosted.z(), 
        new_point.Pp_boosted.y()/new_point.Pp_boosted.z(), 
        new_point.Pp_boosted.Mag()
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
        new_point.vertex.x()*1e3,
        new_point.vertex.y()*1e3, 
        new_point.foil_num,
        new_point.P_electron.x(), new_point.P_electron.y(), new_point.P_electron.z(), 
        1. - new_point.x, 
        new_point.theta * 1.e3, 
        new_point.phi * (180./pi), 
        tj_e_SCS.x*1e3, 
        tj_e_SCS.y*1e3, 
        tj_p_SCS.x*1e3, 
        tj_p_SCS.y*1e3, 
        new_point.in_acceptance ? "yes" : "no"
    );
#endif 

    //reject update
    if ((!new_point.in_acceptance) && (fState.in_acceptance)) { return; }
    
    new_point.amplitude = A_production_amplitude(new_point.x, new_point.theta);  

    //if we've gotten here, we can consider an update. 
    if ( new_point.amplitude > fState.amplitude || new_point.amplitude/fState.amplitude > Rand() ) {

        fState = new_point; 
    }
    return; 
}; 
//_________________________________________________________________________________________________
//_________________________________________________________________________________________________
//_________________________________________________________________________________________________
//_________________________________________________________________________________________________


