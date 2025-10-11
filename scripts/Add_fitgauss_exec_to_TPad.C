#include "TROOT.h"
#include <cmath> 
#include <TPad.h> 
#include <iostream> 
#include <cstdio> 
#include <TH1D.h>
#include <TFitResultPtr.h> 
#include <TFitResult.h> 
#include <TF1.h>

using namespace std; 



void fit_gaussian_to_hist(const double radius); 


int Add_fitgauss_exec_to_TPad(const double radius, TPad* pad=nullptr)
{
    //if no specific TPad is passed, then just use the 'gPad' global pad. 
    const char* const here = "Add_fitgauss_exec_to_TPad"; 

    if (pad==nullptr) {

        //if no specific pad is passed, then just use the global 'gPad', which should exist 
        // any time there's an active TCanvas. 
        if (gPad==nullptr) {
            
            Error(here, "TPad* passed is null (2nd arg), and gPad is also null... we have to TPads to fit!"); 
            return -1; 
        } else {

            pad = (TPad*)gPad; 
        }
    }
    
    pad->AddExec("exec_fitgauss", Form("fit_gaussian_to_hist(%.8e)", radius)); 
    
    return 1; 
}

void fit_gaussian_to_hist(const double radius) 
{
    const char* const here = "exec_fitgauss"; 

    if (!gPad) return;
    if (!gPad->GetSelected()) return;  
    //if (gPad->GetSelected()->IsA() != TClass::GetClass<TH1D>()) return; 

    if (gPad->GetEvent() != 11) return; 

    //find hits
    TH1D* hist=nullptr; 
    auto plist = gPad->GetListOfPrimitives(); 
    for (TObject* obj : *plist) { if (obj->InheritsFrom("TH1")) { hist = (TH1D*)obj; break; } }
    if (!hist) return; 

    double x = gPad->AbsPixeltoX( gPad->GetEventX() ); 
    
    auto gaus_fcn = [](double *x, double *par) {

        double arg = (x[0] - par[1])/par[2]; 
        return par[0] * TMath::Exp( - 0.5 * arg * arg ) + par[3]; 
    }; 


    auto fit = new TF1("gaussfit", "gaus(0) + [3]", x - radius, x + radius ); 

    auto x_ax = hist->GetXaxis(); 

    int bin_low  = x_ax->FindBin( x - radius ); 
    int bin_high = x_ax->FindBin( x + radius ); 

    double max=0.; int max_bin=0;  
    for (int i=bin_low; i<=bin_high; i++) {
        if (max < hist->GetBinContent(i)) {
            max = hist->GetBinContent(i);
            max_bin = i; 
        }
    } 

    double x0_guess = x_ax->GetBinCenter(max_bin); 

    double x_low_val  = hist->GetBinContent( bin_low  );
    double x_high_val = hist->GetBinContent( bin_high ); 

    double const_offset_guess = 0.5 * ( x_low_val + x_high_val ); 
    
    fit->SetParameter(0, max - const_offset_guess); 
    fit->SetParameter(1, x0_guess); 
    fit->SetParameter(2, radius/2.);
    fit->SetParameter(3, const_offset_guess); 

    // options (second character argument):
    //  S - store the hit in the returned FitResult object 
    //  R - only fit to the function rainge 
    //  L - use 'log liklihood' method to fit, which is appropriate for histograms where bins represent fit-counts
    //  Q - don't print information about the fit (we'll do that oursevles.)
    //  
    auto fitresult = hist->Fit("gaussfit", "S R L Q"); 
    
    if (!fitresult->IsValid()) {
        Warning(here, "Fit failed - TFitResult returned is invalid."); 
        return; 
    }

    const double *fit_params = fitresult->GetParams(); 
    const double *fit_errors = fitresult->GetErrors(); 

    if (fit_params==nullptr || fit_errors==nullptr) {
        Warning(here, "Fit failed - array of parameters returned is null."); 
        return; 
    }

    printf(
        "\n"
        "~~~ fit successful. parameters of fit: ~~~~~~~~~~~~~~~~~~~~~~\n"
        " - mean        = %+.6e +/- %.6e\n"
        " - sigma       = %.6e +/- %.6e\n"
        " - offset      = %.6e +/- %.6e\n"
        " - amplitude   = %.6e +/- %.6e\n"
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n", 
        fit_params[1], fit_errors[1],
        fabs(fit_params[2]), fit_errors[2], //we call 'fabs' here because sometimes the fit finds a 'best-fit' value for this with sigma < 0. 
        fit_params[3], fit_errors[3],
        fit_params[0], fit_errors[0]
    ); cout << flush; 
    
    fit->Draw("SAME");
}