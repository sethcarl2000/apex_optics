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

namespace {
    constexpr double beam_energy = 2.200;    //GeV/c2
    constexpr double m_A         = 0.320;    //GeV/c2
    constexpr double m_e         = 0.501e-3; //GeV/c2 

    constexpr double beam_E2 = beam_energy*beam_energy; 
    constexpr double m_A2 = m_A*m_A; 
    constexpr double m_e2 = m_e*m_e; 

    constexpr double pi = 3.14159265359; 
    constexpr double sqrt_2 = 1.41421356237; 

    constexpr double spectrometer_p0 = 1.108; //GeV/c 
}; 

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

    double x=1.; 
    double theta=0.; 
    double amplitude=0.; 

    const double x_sweep_range = 2.0e-5;        //maximum change in x
    const double theta_sweep_range = 0.300e-3;  //maximum change in theta

    random_device rd; 
    mt19937 twister(rd()); 
    uniform_real_distribution<double> rand_generator(0., 1.); 

    //picks a random foil
    uniform_int_distribution<int> random_foil(0, 9); 

    //return a random foil z-position 
    auto random_prod_foil_z = [&random_foil, &twister]()
    {       
        return -238.8e-3  + ((double)random_foil(twister))*55e-3; 
    }; 

    normal_distribution gaus_generator(0., 1.); 
    
    auto get_uniform = std::bind(rand_generator, twister); 
    auto get_gaus    = std::bind(gaus_generator, twister); 

    auto metropolis_update = [&x, &theta, &amplitude, &get_uniform, x_sweep_range, theta_sweep_range]() 
    {
        double x_new     = x     + (1. - 2.*get_uniform())*x_sweep_range; 
        double theta_new = theta + (1. - 2.*get_uniform())*theta_sweep_range; 

        double new_amplitude = A_production_amplitude(x_new, theta_new); 

        if (new_amplitude > amplitude || get_uniform() < new_amplitude/amplitude ) {
            amplitude = new_amplitude; 
            x         = x_new; 
            theta     = theta_new; 
        }
    }; 

    TH2D* x_and_theta = new TH2D("h", "A' produced;x;#theta_{A}", 200, 1. - 16e-5, 1., 200, -1.5, +1.5); 

    TH1D* h_e_momentum = new TH1D("h_momentum", "electron / positron momentum;GeV/c;", 200, -1., 1.+beam_energy); 

    TH2D* h_dx_dy_hall = new TH2D("h_dx_dy", "electron/positron angle in lab;dx/dz (mrad);dy/dz (mrad)", 200, -200,+200, 200, -200,+200); 

    TH2D* h_x_y = new TH2D(
        "h_x_y_sv", 
        Form("electrons on LHRS sieve face - m_{A} = %.1f MeV/c^{2};x_{sv} mm (lab-vertical);y_{sv} mm (lab-horizontal)",m_A*1.e3), 
        200,-70,+70, 200,-70,+70); 

    const long int n_stats = 3e5;
    long int n_generated=0; 

    printf("generating %li events...", n_stats); cout << flush; 
    while (n_generated < n_stats) {

        metropolis_update(); 

        x_and_theta->Fill(x, theta*1e3);

        //now, let's create a positron / electron pair, and boost them to the lab frame. 

        //frist, create the direction of the A' (randomly)
        double phi = 2.*pi*get_uniform(); 

        double E_a = x * beam_energy; 
        double p_A_mag = sqrt( (E_a*E_a) - m_A2 ); 

        TVector3 p_A(
            p_A_mag * cos(theta) * cos(phi), 
            p_A_mag * sin(theta) * sin(phi), 
            p_A_mag * cos(theta)
        ); 

        //now, in the rest-frame of the A', we generate the directions of the positron / electron randomly
        double gamma_A = E_a / m_A; 
        double beta_A = sqrt( 1. - 1./(gamma_A*gamma_A) ); 

        //energy of positron & electron in A' rest frame. this is very close to the energy of the electron 
        double E_electron = m_A/2.; 
        double p_electron_mag = sqrt( (E_electron*E_electron) - m_e2 ); 

        TVector3 p_electron( get_gaus(), get_gaus(), get_gaus() );
        p_electron = p_electron_mag * p_electron.Unit();  

        TVector3 p_positron = -1. * p_electron; 
        
        //_________________________________________________________________________________________
        //take a vector in the rest-frame of the A', and boost it to the lab frame
        auto boost_and_rotate = [gamma_A, beta_A, theta, phi, E_electron](TVector3 &p, double E)
        {
            //E    = gamma_A * (E    + beta_A*p[2]); 
            p[2] = gamma_A * (p[2] + beta_A*E   ); 
            
            //so, now that we've boosed the momentum back to the lab, 
            //let's rotate this vector, so that it's boost is in line with the momentum of the A' 
            p.RotateX( -theta ); 
            p.RotateZ( phi ); 
        };
        //_________________________________________________________________________________________

        boost_and_rotate( p_electron, E_electron ); 
        boost_and_rotate( p_positron, E_electron ); 

        //check momentum of positron and electron
        if (fabs(p_electron.Mag() - spectrometer_p0)/spectrometer_p0 > 0.05) continue; 
        if (fabs(p_positron.Mag() - spectrometer_p0)/spectrometer_p0 > 0.05) continue; 

        n_generated++; 

        h_e_momentum->Fill( sqrt( p_electron.Mag2() + m_e2 ) );
        h_dx_dy_hall->Fill( 1e3*p_electron.x()/p_electron.z(), 1e3*p_electron.y()/p_electron.z() ); 
        
        h_e_momentum->Fill( sqrt( p_positron.Mag2() + m_e2 ) ); 
        h_dx_dy_hall->Fill( 1e3*p_positron.x()/p_positron.z(), 1e3*p_positron.y()/p_positron.z() ); 
    
        //randomly generate the vertex, in line with the production target geometry
        TVector3 vertex(
            2.e-3*( 1. - 2.*get_uniform()), 
            2.e-3*( 1. - 2.*get_uniform()), 
            random_prod_foil_z()
        ); 

        //electron's trajectory in Hall Coordinate System
        Trajectory_t traj_electron_HCS{
            vertex.x() + (p_electron.x()/p_electron.z())*(0. - vertex.z()), 
            vertex.y() + (p_electron.y()/p_electron.z())*(0. - vertex.z()), 
            p_electron.x()/p_electron.z(), 
            p_electron.y()/p_electron.z(), 
            p_electron.Mag() 
        }; 

        auto traj_electron_SCS = ApexOptics::HCS_to_SCS(is_RHRS, traj_electron_HCS); 

        h_x_y->Fill( traj_electron_SCS.x*1e3, traj_electron_SCS.y*1e3 ); 
    }
    cout << "done" << endl; 

    new TCanvas; 
    x_and_theta->Draw("col");

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

