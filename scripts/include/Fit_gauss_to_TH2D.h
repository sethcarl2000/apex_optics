#ifndef Fit_gauss_to_TH2D_H
#define Fit_gauss_to_TH2D_H

#include <cmath> 
#include "Fit_fcn_to_TH2D.h"
#include <TH2D.h>
#include <functional>
#include <optional> 
#include <vector> 
#include <TF2.h> 

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

        const double offset     = hist->GetMinimum(); 
        const double amplitude  = hist->GetMaximum() - offset; 

        const double cent_x     = 0.5 * (fit_domain_x1 + fit_domain_x0); 
        const double cent_y     = 0.5 * (fit_domain_y1 + fit_domain_y0); 

        const double a_xx   = 0.5 * (fit_domain_x1 - fit_domain_x0); 
        const double a_yy   = 0.5 * (fit_domain_y1 - fit_domain_y0); 

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
    

};



#endif 