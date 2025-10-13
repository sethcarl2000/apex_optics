#ifndef Fit_gauss_to_TH2D_H
#define Fit_gauss_to_TH2D_H

#include <cmath> 
#include "Fit_fcn_to_TH2D.h"
#include <TH2D.h>
#include <TAxis.h> 
#include <functional>
#include <optional> 
#include <vector> 
#include <cstdio> 
#include <TF2.h> 
#include <TBox.h> 
#include <TPad.h> 

namespace Gauss2D {

    //_______________________________________________________________________________________________________________________
    struct FitResult_t {
        //these are the domains of the fit
        double 
            x0, x1, 
            y0, y1; 

        double 
            amplitude,
            cent_x, cent_y, 
            a_xx, a_xy, a_yy,
            offset; 

        double sigma_x()    const { return 1./std::sqrt(a_xx); }
        double sigma_y()    const { return 1./std::sqrt(a_yy); }
        double sigma_xy()   const { return 1./std::sqrt(a_xy); }

        //angle this ellipse makes with the horizontal axis
        double theta()      const { return 0.5*std::atan( a_xy/(a_xx - a_yy)); }

        //sigma of major axis of ellipse
    };
    //_______________________________________________________________________________________________________________________
    
    
    //a function describing a 2d-gaussian function
    //_______________________________________________________________________________________________________________________
    double fcn(const double *X, const double* par) 
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
    //_______________________________________________________________________________________________________________________
    
    
    //generate a 2d fcn (TF2) based off of a (valid) fit result
    //_______________________________________________________________________________________________________________________
    TF2 *Make_TF2_from_fit(const Gauss2D::FitResult_t& fit, const char* name = "gaus2d") {

        //create a new TF2
        auto tf2 = new TF2(name, Gauss2D::fcn, fit.x0,fit.x1, fit.y0,fit.y1, 7);
        
        //set the parameters
        double params[] = {
            fit.amplitude, 
            fit.cent_x, 
            fit.cent_y,
            fit.a_xx, 
            fit.a_xy,
            fit.a_yy, 
            fit.offset
        }; 
        tf2->SetParameters(params); 

        //return the new TF2    
        return tf2; 
    }
    //_______________________________________________________________________________________________________________________
    

    //_______________________________________________________________________________________________________________________
    std::optional<Gauss2D::FitResult_t> Fit_TH2D(
        TH2D* hist, 
        double fit_domain_x0 = -1.e100,
        double fit_domain_y0 = -1.e100,
        double fit_domain_x1 = +1.e100,
        double fit_domain_y1 = +1.e100
    )
    {
        //here, we're going to load in some defualt values
        if (hist==nullptr) return std::nullopt; 

        //find the minimum and maximum values for the hist IN THIS SUBRANGE 
        auto x_ax = hist->GetXaxis(); 
        auto y_ax = hist->GetYaxis(); 

        //if no fit domain is given, take the range limits to be that of the histogram 
        fit_domain_x0 = max<double>( fit_domain_x0, x_ax->GetXmin() ); 
        fit_domain_x1 = min<double>( fit_domain_x1, x_ax->GetXmax() ); 
        
        fit_domain_y0 = max<double>( fit_domain_y0, y_ax->GetXmin() ); 
        fit_domain_y1 = min<double>( fit_domain_y1, y_ax->GetXmax() ); 

        double subrange_max=-1.e30; 
        double subrange_min=+1.e30; 

        for (int bx=x_ax->FindBin(fit_domain_x0); bx<=x_ax->FindBin(fit_domain_x1); bx++) {
            for (int by=y_ax->FindBin(fit_domain_y0); by<=y_ax->FindBin(fit_domain_y1); by++) {

               subrange_max = max<double>( subrange_max, hist->GetBinContent(bx, by) ); 
               subrange_min = min<double>( subrange_min, hist->GetBinContent(bx, by) ); 
            }
        }

        const double offset     = subrange_min; 
        const double amplitude  = subrange_max - offset; 

        const double cent_x     = 0.5 * (fit_domain_x1 + fit_domain_x0); 
        const double cent_y     = 0.5 * (fit_domain_y1 + fit_domain_y0); 

        //this corresponds to a starting sigma of 1/2 the fit range of each axis
        const double a_xx   = 1. / pow(0.5 * fabs(fit_domain_x1 - fit_domain_x0), 2); 
        const double a_yy   = 1. / pow(0.5 * fabs(fit_domain_y1 - fit_domain_y0), 2); 

        auto first_guess = Gauss2D::FitResult_t{
            .x0 = fit_domain_x0, 
            .x1 = fit_domain_x1,
            .y0 = fit_domain_y0, 
            .y1 = fit_domain_y1,

            .amplitude  = amplitude, 
            .cent_x     = cent_x,
            .cent_y     = cent_y,
            .a_xx       = a_xx,
            .a_xy       = 0.,
            .a_yy       = a_yy,
            .offset     = offset
        };
        
        /*auto tf2_first_guess = Gauss2D::Make_TF2_from_fit(first_guess); 
        //tf2_first_guess->SetLineColor(kBlack); 
        tf2_first_guess->Draw("SAME"); */ 

        std::vector<FitVariable_t> vars{ 
            {"amplitude", amplitude},
            {"center_x",  cent_x},
            {"center_y",  cent_y},
            {"a_xx",      a_xx},
            {"a_xy",      0.},
            {"a_yy",      a_yy},
            {"offset",    offset}
        };

        auto fit_result_opt = Fit_fcn_to_TH2D(hist, Gauss2D::fcn, vars, 
            fit_domain_x0, fit_domain_y0,
            fit_domain_x1, fit_domain_y1
        ); 
        
        /*
        printf("fit range: x = [%+.3f,%+3.f], y = [%+.3f,%+.3f]\n", 
            fit_domain_x0, fit_domain_x1, 
            fit_domain_y0, fit_domain_y1
        ); 
        if (gPad) {
            auto box = new TBox(fit_domain_x0,fit_domain_x1, fit_domain_y0,fit_domain_y1); 
            box->SetLineColor(kRed); 
            
            if (fit_result_opt.has_value()) { box->SetLineStyle(kSolid); } else { box->SetLineStyle(kDotted); }
            box->SetFillStyle(0); 
            box->Draw("SAME");
        }*/ 

        if (!fit_result_opt) return std::nullopt; 

        auto fit_result = fit_result_opt.value(); 

        return Gauss2D::FitResult_t{
            .x0 = fit_domain_x0, 
            .x1 = fit_domain_x1,
            .y0 = fit_domain_y0, 
            .y1 = fit_domain_y1,

            .amplitude  = fit_result[0].val, 
            .cent_x     = fit_result[1].val,
            .cent_y     = fit_result[2].val,
            .a_xx       = fit_result[3].val,
            .a_xy       = fit_result[4].val,
            .a_yy       = fit_result[5].val,
            .offset     = fit_result[6].val
        }; 

    }
    //_______________________________________________________________________________________________________________________
    

};



#endif 