#ifndef fit_gaus_to_hist_H
#define fit_gaus_to_hist_H 

#include <TH1.h> 
#include <TF1.h>
#include <TAxis.h>  
#include <stdexcept>
#include <TFitResultPtr.h> 
#include <Math/ProbFuncMathCore.h>
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

    double dx = (x_max - x_min)/((double)x_ax->GetNbins());

    //the guess for the 'base' will be based on the bins on either extreme end of the fit
    base = 0.5*(
        hist->GetBinContent( x_ax->FindBin(center + radius) ) +
        hist->GetBinContent( x_ax->FindBin(center - radius) )
    ) / dx; 

    //amplitude over background
    amplitude = (hist->GetBinContent(x_ax->FindBin(center)))/dx - base; 
    
    //our function to fit the histogram with   
    auto gauss_with_offset = [dx](double *x, double *par) {

        double sigma = std::fabs(par[2]);
        double x0 = par[1];
        double arg = (x[0] - x0)/sigma; 
        
        double expect = ROOT::Math::normal_cdf(x[0]+dx/2., sigma, x0) - ROOT::Math::normal_cdf(x[0]-dx/2., sigma, x0);
        expect *= sigma*2.50662827463; 
        return par[0]*expect + par[3]*dx; 
    };

    auto gaus_fit = new TF1("gaus_fit", gauss_with_offset, x_min, x_max, 4); 
  
    //set the parameters 
    gaus_fit->SetParameter( 0, amplitude ); 
    gaus_fit->SetParameter( 1, center ); 
    gaus_fit->SetParameter( 2, sigma ); 
    gaus_fit->SetParameter( 3, base ); 
    
    auto fit_ptr = hist->Fit("gaus_fit", (do_draw ? "L R S" : "L N R S Q")); 
    
    if (fit_ptr->IsValid()==false) {
        throw std::logic_error("<fit_gaus_to_hist>: fit failed."); 
        return; 
    }

    amplitude   = fit_ptr->Parameter(0); 
    center      = fit_ptr->Parameter(1); 
    sigma       = std::fabs(fit_ptr->Parameter(2)); 
    base        = fit_ptr->Parameter(3); 
    
    return; 
}

#endif 