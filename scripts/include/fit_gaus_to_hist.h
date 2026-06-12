#ifndef fit_gaus_to_hist_H
#define fit_gaus_to_hist_H 

#include <TH1.h> 
#include <TF1.h>
#include <TAxis.h>  
#include <stdexcept>
#include <TFitResultPtr.h> 
#include <TFitResult.h>
#include <TString.h> 
#include <cmath> 

/// @brief Attempts to fit a gaussian to a histogram, with a constant background
/// @param hist histogram to fit
/// @param radius radius of fit
/// @param center center of gaus, to be optimized by fit (value of arg. overwritten)
/// @param sigma sigma of gaus, to be optimized by fit (value of arg. overwritten) 
/// @param amplitude amplitude of gaussian peak
/// @param base constant vertical offset of gaussian profile 
/// @param do_draw if true, then the fit will be drawn 
void fit_gaus_to_hist( 
    TH1 *hist, 
    double radius, 
    double &center, 
    double &sigma, 
    double &amplitude, 
    double &base,
    bool do_draw=false
) 
{
    //first guess for the center of the fit. we assume here that the gaus we wanna fit is the tallest peak in the hist. 
    center = hist->GetBinCenter( hist->GetMaximumBin() ); 
    
    //first guess for sigma 
    sigma = radius/4; 

    auto x_ax = hist->GetXaxis(); 
    
    double x_min = std::max( center-radius, x_ax->GetXmin() );    
    double x_max = std::min( center+radius, x_ax->GetXmax() );

    //the guess for the 'base' will be based on the bins on either extreme end of the fit
    base = 0.5*(
        hist->GetBinContent( x_ax->GetBinCenter(x_min) ) +
        hist->GetBinContent( x_ax->GetBinCenter(x_max) )
    ); 

    //amplitude over background
    amplitude = hist->GetMaximum() - base; 
    
    //our function to fit the histogram with   
    auto gauss_with_offset = [](double *x, double *par) {
        double arg = (x[0] - par[1])/par[2];
        return par[0]*std::exp( -arg*arg/2. ) + par[3]; 
    };

    auto gaus_fit = new TF1("gaus_fit", gauss_with_offset, x_min, x_max, 4); 
  
    //set the parameters 
    gaus_fit->SetParameter( 0, amplitude ); 
    gaus_fit->SetParameter( 1, center ); 
    gaus_fit->SetParameter( 2, sigma ); 
    gaus_fit->SetParameter( 3, base ); 
    
    auto fit_ptr = hist->Fit("gaus_fit", (do_draw ? "L R S Q" : "L N R S Q")); 
    
    if (fit_ptr->IsValid()==false) {
        throw std::logic_error("<fit_gaus_to_hist>: fit failed."); 
        return; 
    }

    center = fit_ptr->Parameter(1); 
    sigma  = std::fabs(fit_ptr->Parameter(2)); 
    
    return; 
}

#endif 