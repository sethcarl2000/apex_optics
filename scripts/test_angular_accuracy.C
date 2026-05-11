#include <ROOT/RDataFrame.hxx> 
#include <vector> 
#include <cmath> 
#include <utility> 
#include <ApexOptics.h> 
#include <TH1D.h> 
#include <TCanvas.h>
#include <TPad.h> 

using namespace std; 

const double sieve_angle_L = ApexOptics::Get_sieve_angle(false); 
const double sieve_angle_R = ApexOptics::Get_sieve_angle(true); 

//coefficient of the uncertainty in theta
double delta_theta_L(double theta_R, double theta_L)
{
    double val = cos(theta_R) * sin(theta_L); 
    return val * val; 
}

//coefficient of the uncertainty of phi
double delta_phi_L(double phi_R, double phi_L)
{
    phi_R += sieve_angle_R; 
    phi_L += sieve_angle_L; 
    double val = cos(phi_R) * sin(phi_L); 
    return val * val; 
}

//coefficient of the uncertainty of phi
double delta_phi(double phi_R, double phi_L, double theta_R, double theta_L)
{
    phi_R += sieve_angle_R; 
    phi_L += sieve_angle_L; 
    double val = 0.; 

    val = pow( cos(phi_R)*sin(phi_L), 2 ) + pow( cos(phi_L)*sin(phi_R), 2 ); 
    val *= 1./pow( sin(theta_R)*sin(theta_L) + sin(phi_L)*sin(phi_R), 2 ); 

    return val; 
}

//coefficient of the uncertainty of theta
double delta_theta(double phi_R, double phi_L, double theta_R, double theta_L)
{
    phi_R += sieve_angle_R; 
    phi_L += sieve_angle_L; 
    double val = 0.; 

    val = pow( cos(theta_R)*sin(theta_L), 2 ) + pow( cos(theta_L)*sin(theta_R), 2 ); 
    val *= 1./pow( sin(theta_R)*sin(theta_L) + sin(phi_L)*sin(phi_R), 2 ); 

    return val; 
}

int test_angular_accuracy()
{
    const char* path_infile = "data/mc/mc_L_production.root";
    const char* tree_name = "tracks_fp"; 
    
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(tree_name, path_infile);

    struct PonintPair_t { double theta,phi; };

    auto events = *df 
        .Define("theta_phi", [](double theta_sv, double phi_sv)
        {
            return PonintPair_t{theta_sv,phi_sv};
        }, {"dxdz_sv", "dydz_sv"})
        
        .Take<PonintPair_t>("theta_phi"); 

    auto hist_delta_theta = new TH1D("h_theta", "#Delta_{#theta}", 200, 0., 1500); 
    auto hist_delta_phi   = new TH1D("h_phi",   "#Delta_{#phi}",   200, 0., 1500); 

    for (int i=0; i<events.size()-1; i++) {

        auto ptL = events[i]; 
        auto ptR = events[i]; 

        //symmetry flip between the arms
        ptR.phi *= -1.; 

        hist_delta_theta->Fill( delta_theta(ptR.phi, ptL.phi, ptR.theta, ptL.theta) );

        hist_delta_phi  ->Fill( delta_phi(ptR.phi, ptL.phi, ptR.theta, ptL.theta) );
    }

    new TCanvas; 
    gPad->SetLogy(1); 
    

    hist_delta_theta->SetLineColor(kRed); 
    hist_delta_theta->SetFillColor(kRed); 
    hist_delta_theta->SetFillStyle(3004); 

    hist_delta_theta->Draw();

    hist_delta_phi->SetLineColor(kBlack); 
    hist_delta_phi->SetFillColor(kBlack); 
    hist_delta_phi->SetFillStyle(3005); 
    
    hist_delta_phi->Draw("SAME");

    return 0; 
}