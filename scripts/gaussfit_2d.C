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
#include <iostream> 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>

using namespace std; 

//a function describing a 2d-gaussian function
double fcn_2d_gauss(const double *X, const double* par) 
{
    const double amplitude = fabs(par[0]); 
    const double& x0   = par[1]; 
    const double& y0   = par[2];
    const double a_xx  = fabs(par[3]); 
    const double a_xy  = fabs(par[4]);
    const double a_yy  = fabs(par[5]);
    
    const double x = X[0] - x0;
    const double y = X[1] - y0; 

    return amplitude * exp( -0.5*((x*x*a_xx) + (x*y*a_xy) + (y*y*a_yy)) );
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

    //look at fcn above for detailed def.
    double params[] = {
        1.,     //amplitude
        0.,     //x-mean
        0.,     //y-mean
        1.,     // a_xx
        0.,     // a_xy
        1.,     // a_yy
    }; 
    tf2->SetParameters(params); 

    auto hist = new TH2D("h", "test", 60, -5, 5, 60, -5, 5); 
    TRandom3 rand; 
    for (int i=0; i<1e5; i++) {
        hist->Fill( rand.Gaus(), rand.Gaus() ); 
    }

    //compute the log-liklihood of a hist with a given TF2
    auto LogLiklihood = [](TH2D* hist, TF2* fcn) 
    {
        double log_liklihood = 0.; 

        const int nbins_x = hist->GetXaxis()->GetNbins(); 
        const int nbins_y = hist->GetYaxis()->GetNbins(); 

        double* bin_x_ptr = new double[nbins_x]; hist->GetXaxis()->GetCenter(bin_x_ptr); 
        double* bin_y_ptr = new double[nbins_y]; hist->GetYaxis()->GetCenter(bin_y_ptr); 

        for (int ix=0; ix<nbins_x; ix++) {

            const double x = bin_x_ptr[ix];

            for (int iy=0; iy<nbins_y; iy++) {

                const double y = bin_y_ptr[iy];

                double fcn_val  = fcn->Eval(x, y); 

                //histogram bin counting starts from '1' and goes to 'ninbs', inclusive. 
                double hist_val = hist->GetBinContent(ix+1, iy+1); 

                //use stirling's approx 
                double log_nfact = hist_val > 0. ? log_nfact = hist_val * log(hist_val) - hist_val : 0.; 

                //this is just the log of the poisson dist, taken with n=hist_val, and <n>=fcn_val
                log_liklihood += (hist_val * log(fcn_val)) - fcn_val - log_nfact; 
            }
        }
        return log_liklihood; 
    };

    hist->Draw("col2"); 

    vector<double> pts_x, pts_LL; 

    double x0 = -2.; 
    for (int i=0; i<201; i++) {

        tf2->SetParameter(1, x0);

        pts_x.push_back(x0); 

        double LL = LogLiklihood(hist, tf2);

        pts_LL.push_back(LL);

        x0 += 0.02; 

        printf("x %+.3f, LL: %.3e", x0, LL); cout << endl;      
    }

    tf2->SetParameter(0, hist->GetMaximum()); 
    tf2->Draw("SAME"); 

    new TCanvas; 
    auto graph = new TGraph(pts_x.size(), pts_x.data(), pts_LL.data()); 

    graph->SetTitle("Log Liklihood vs x_{0};x_{0};L.L.");
    graph->Draw(); 
    
    return 0; 
}