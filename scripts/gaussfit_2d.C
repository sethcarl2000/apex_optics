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
    
    const double x = X[0] - x0;
    const double y = X[1] - y0; 

    return amplitude * exp( -0.5*((x*x*a_xx) + (x*y*a_xy) + (y*y*a_yy)) );
}

//template<typename T> struct MinMaxVal_t { T min, max; }; 

vector<double> Fit_TF2_to_TH2D(TH2D* hist, TF2* fcn_inp, const std::vector<std::pair<std::string,double>>& varlist)
{   
    //first, collect all the bin contents we need
    auto x_ax = hist->GetXaxis();
    auto y_ax = hist->GetYaxis();

    //make a copy of the input fcn we can mess with
    auto fcn = (TF2*)fcn_inp->Clone(); 

    const int n_fcn_params = fcn->GetNpar(); 

    double fit_domain_x0, fit_domain_x1; 
    double fit_domain_y0, fit_domain_y1; 

    fcn->GetRange( 
        fit_domain_x0, fit_domain_x1, 
        fit_domain_y0, fit_domain_y1 
    ); 

    fit_domain_x0 = min<double>( fit_domain_x0, x_ax->GetXmin() ); 
    fit_domain_y0 = min<double>( fit_domain_y0, y_ax->GetXmin() ); 
    
    fit_domain_x1 = max<double>( fit_domain_x1, x_ax->GetXmax() ); 
    fit_domain_y1 = max<double>( fit_domain_y1, y_ax->GetXmax() ); 
    
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
    vector<BinContent_t> bin_contents; bin_contents.reserve(n_bins_x * n_bins_y); 

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
    
    if (!minimizer) return {}; 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(0);    

    //simple function to compute log liklihood
    //____________________________________________________________________________________________________
    auto LogLiklihood = [&bin_contents, &fcn](const double* par) {

        fcn->SetParameters(par);

        double log_liklihood = 0.;

        //array of x-vals which we will pass as the fcn argument
        double X[] = {0., 0.}; 

        //loop thru all bins
        for (const auto& bin : bin_contents) {
            
            X[0] = bin.x; 
            X[1] = bin.y; 

            double fcn_val  = fcn_2d_gauss(X, par);// //fcn->Eval(bin.x, bin.y);  //fcn_2d_gauss(X, par);

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
        minimizer->SetVariable(i_var++, var.first.c_str(), var.second, 1e-4);
    }

    auto status = minimizer->Minimize();

    cout << "minimizer status: " << std::boolalpha << status << endl; 

    const double* par_result = minimizer->X(); 

    if (!par_result) return {}; 

    vector<double> vals{}; 
    for (int i=0; i<n_fcn_params; i++) vals.push_back(par_result[i]); 
    
    return vals; 
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

    
    auto hist = new TH2D("h", "test", 60, -5, 5, 60, -5, 5); 
    TRandom3 rand; 
    for (int i=0; i<1e5; i++) {
        double x = rand.Gaus(); 
        double y = rand.Gaus() - 0.4*x; 
        hist->Fill( x, y ); 
    }

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
    
    auto fit_params = Fit_TF2_to_TH2D(hist, tf2, {
        {"amplitude",   hist->GetMaximum()},
        {"x0",          1.},
        {"y0",          0.},
        {"a_xx",        1.},
        {"a_xy",        0.},
        {"a_yy",        1.}
    });
    
    tf2->SetParameters(fit_params.data()); 

    hist->Draw("col2"); 
    //tf2->SetParameter(0, hist->GetMaximum()); 
    tf2->Draw("SAME"); 

    return 0; 
}