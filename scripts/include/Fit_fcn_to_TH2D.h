#ifndef Fit_TH2D_H
#define Fit_TH2D_H

#include <TF2.h> 
#include <TH2D.h> 
#include <TAxis.h> 
#include <cmath> 
#include <vector> 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <functional> 
#include <limits> 
#include <optional> 

//________________________________________________________________________________________________________
//
//  This takes a std::function as input (with signature double(const double* X, const double* par)), 
//  and fits a (sub-range) of a TH2D.  
//  args: 
//  
//  Fit_fcn_to_TH2D( 
//      TH2D *hist,         
//          - histogram
//      std::function<double(const double*, const double*)>  fcn,   
//          - the function to fit. the form should be double* X, and double* par. 
//      double fit_domain_... 
//              ...x0   - min x-range. optional: if none is provided, the histo min x-val is used. 
//              ...x1   - max x-range. 
//              ...y0   - min y-range.
//              ...y1   - max y-range
//  )
//
//________________________________________________________________________________________________________

//template<typename T> struct MinMaxVal_t { T min, max; }; 
static constexpr double NaN_double = std::numeric_limits<double>::quiet_NaN();  

struct FitVariable_t { 
    std::string name; 
    double val; 
    double error=NaN_double; 
}; 


std::optional<std::vector<FitVariable_t>> Fit_fcn_to_TH2D(
    TH2D* hist, 
    std::function<double(const double*,const double*)> fcn, 
    const std::vector<FitVariable_t>& varlist,
    double fit_domain_x0 = -1.e100,
    double fit_domain_y0 = -1.e100,
    double fit_domain_x1 = +1.e100,
    double fit_domain_y1 = +1.e100
)
{   
    //first, collect all the bin contents we need
    auto x_ax = hist->GetXaxis();
    auto y_ax = hist->GetYaxis();

    //make a copy of the input fcn we can mess with
    const int n_fcn_params = varlist.size();    

    fit_domain_x0 = std::max<double>( fit_domain_x0, x_ax->GetXmin() ); 
    fit_domain_y0 = std::max<double>( fit_domain_y0, y_ax->GetXmin() ); 
    
    fit_domain_x1 = std::min<double>( fit_domain_x1, x_ax->GetXmax() ); 
    fit_domain_y1 = std::min<double>( fit_domain_y1, y_ax->GetXmax() ); 
    
    const int 
        bin_x_min{ x_ax->FindBin(fit_domain_x0) },
        bin_x_max{ x_ax->FindBin(fit_domain_x1) }; 

    const int 
        bin_y_min{ y_ax->FindBin(fit_domain_y0) },
        bin_y_max{ y_ax->FindBin(fit_domain_y1) }; 
    
    //get the bin contents from the histogram
    const int n_bins_x = bin_x_max - bin_x_min; 
    const int n_bins_y = bin_y_max - bin_y_min; 

    struct BinContent_t { double x, y, n; }; 
    std::vector<BinContent_t> bin_contents; bin_contents.reserve(n_bins_x * n_bins_y); 

    for (int ix=bin_x_min; ix<=bin_x_max; ix++) {
        for (int iy=bin_y_min; iy<=bin_y_max; iy++) {

            bin_contents.push_back({
                x_ax->GetBinCenter(ix),
                y_ax->GetBinCenter(iy), 
                hist->GetBinContent(ix, iy) 
            }); 
        }
    }

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    
    if (!minimizer) return std::nullopt; 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(2);    

    //simple function to compute log liklihood
    //____________________________________________________________________________________________________
    auto LogLiklihood = [&bin_contents, &fcn](const double* par) {

        double log_liklihood = 0.;

        //array of x-vals which we will pass as the fcn argument
        double X[] = {0., 0.}; 

        //loop thru all bins
        for (const auto& bin : bin_contents) {
            
            X[0] = bin.x; 
            X[1] = bin.y; 

            double fcn_val  = fcn(X, par);// //fcn->Eval(bin.x, bin.y);  //fcn_2d_gauss(X, par);

            //histogram bin counting starts from '1' and goes to 'ninbs', inclusive. 

            //use stirling's approx 
            double log_nfact = bin.n > 0. ? log_nfact = bin.n * log(bin.n) - bin.n : 0.; 

            //this is just the log of the poisson dist, taken with n=hist_val, and <n>=fcn_val
            log_liklihood += (bin.n * log(fcn_val)) - fcn_val - log_nfact; 
        }
        return -1.*log_liklihood; 
    };
    //___________________________________________________________________________________________________________

    auto f_minimizer = ROOT::Math::Functor(LogLiklihood, n_fcn_params);
    
    minimizer->SetFunction(f_minimizer);

    //set the list of variables
    int i_var=0; 
    for (const auto& var : varlist) {
        minimizer->SetVariable(i_var++, var.name.c_str(), var.val, 1e-4);
    }

    bool fit_status = minimizer->Minimize();

    if (!fit_status) return std::nullopt; 

    const double* par_result = minimizer->X(); 
    const double* par_errors = minimizer->Errors(); 

    if (!par_result) return std::nullopt; 

    std::vector<FitVariable_t> vals{}; 
    for (int i=0; i<n_fcn_params; i++) {
        
        vals.push_back({
            varlist[i].name,
            par_result[i],
            par_errors[i]
        }); 
    }
    
    return vals; 
} 

#endif 