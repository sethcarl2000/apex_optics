#include <TF2.h> 
#include <TH2D.h> 
#include <TAxis.h> 
#include <TFitResult.h>
#include <TFitResultPtr.h> 
#include <TCanvas.h> 
#include <TPad.h>
#include <TRandom3.h>  
#include <TGraph.h> 
#include <cmath> 
#include <vector> 
#include <iostream> 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <functional> 
#include <limits> 
#include <optional> 
#include "include/Fit_fcn_to_TH2D.h"

using namespace std; 

//a function describing a 2d-gaussian function
double fcn_2d_gauss(const double *X, const double* par) 
{
    const double amplitude = fabs(par[0]); 
    const double& x0   = par[1]; 
    const double& y0   = par[2];
    const double a_xx  = fabs(par[3]); 
    const double a_xy  = par[4];
    const double a_yy  = fabs(par[5]);
    const double& offset = par[6];
    
    const double x = X[0] - x0;
    const double y = X[1] - y0; 

    return amplitude * exp( -0.5*((x*x*a_xx) + (x*y*a_xy) + (y*y*a_yy)) ) + offset; 
}


int gaussfit_2d()
{
    const char* const here = "gaussfit_2d";

    auto tf2 = new TF2(
        "gauss_fcn",        //name of fcn 
        fcn_2d_gauss,       //ptr to fcn that will do calculation
        -5.,5.,             // x min and max
        -5.,5.,             // y min and max
        6                   // number of parameters
    ); 

    
    auto hist = new TH2D("h", "test", 50, -5, 5, 50, -5, 5); 
    TRandom3 rand; 
    for (int i=0; i<1e5; i++) {
        double x = rand.Gaus(); 
        double y = 2. + (rand.Gaus()*0.8 + x*0.8); 
        hist->Fill( x, y ); 
    }

    for (int i=0; i<1e6; i++) { hist->Fill( -5. +  rand.Rndm()*10., -5. + rand.Rndm()*10. ); }


    //look at fcn above for detailed def.
    double params[] = {
        hist->GetMaximum(),     //amplitude
        0.,     //x-mean
        0.,     //y-mean
        1.,     // a_xx
        0.,     // a_xy
        1.,     // a_yy
    }; 
    tf2->SetParameters(params); 
    
    const vector<FitVariable_t> fit_vars{
        {"amplitude",   hist->GetMaximum()},
        {"x0",          0.},
        {"y0",          1.5},
        {"a_xx",        1.},
        {"a_xy",        0.},
        {"a_yy",        1.},
        {"offset",      400.}
    };

    auto fit_params = Fit_fcn_to_TH2D(hist, fcn_2d_gauss, fit_vars); 

    //hist->GetZaxis()->SetRangeUser(0., hist->GetMaximum()*1.25); 
    hist->Draw("colz"); 

    if (fit_params) {

        int i_var=0; 
        for (auto var : fit_params.value()) tf2->SetParameter(i_var++, var.val); 
        
        tf2->Draw("SAME"); 
    }

    return 0; 
}