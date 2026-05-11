#include <TCanvas.h>
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <functional> 
#include <TRandom3.h> 
#include <TH1D.h> 
#include <cmath> 
#include <TAxis.h> 
#include <Math/ProbFunc.h>
#include <iostream> 
#include <TGraph.h> 
#include <TLine.h> 
#include <TStyle.h> 
#include <THStack.h>
#include <TScatter.h> 
#include <TH2D.h> 
#include <ROOT/RVec.hxx> 
#include <RMatrix.h> 
#include <TF1.h> 

using namespace std; 


const double width_background = 1.0; 

double gaus_pdf(double x, double sigma=1., double mean=0.)
{  
    double arg = (x - mean)/sigma; 
    return 0.398942280401 * exp( -0.5 * arg * arg ) / sigma; 
}

struct HistPoint_t{ double x,N; }; 

vector<HistPoint_t> compute_q0(TH1D* hist, double signal_sigma)
{
    using namespace ROOT::Math; 

    //make a mass hypothesis for each bin, starting from bin 10 and going to bin N-10.  

    auto ax = hist->GetXaxis(); 
    
    const double total_stats = hist->Integral(); 

    vector<HistPoint_t> points; points.reserve(ax->GetNbins()); 

    for (int b=1; b<=ax->GetNbins(); b++) points.push_back({ ax->GetBinCenter(b), hist->GetBinContent(b) }); 

    vector<HistPoint_t> points_background;
    points_background.reserve(points.size()); 
    
    const double dx = points[1].x - points[0].x; 
    
    double normalization_mult = 
        normal_cdf(points.back().x  + dx/2., width_background, 0.) - 
        normal_cdf(points.front().x - dx/2., width_background, 0.); 

    //cout << "Normalization mult: " << normalization_mult << endl; 

    for (const auto& point : points) {

        double x = point.x;

        double expect = normal_cdf(x + dx/2., width_background) - normal_cdf(x - dx/2., width_background); 
        expect *= total_stats / normalization_mult; 

        points_background.push_back({ x, expect });
    }

    vector<HistPoint_t> points_lambda;  

    for (int i=10; i<ax->GetNbins()-10; i++) {

        const double signal_mean = ax->GetBinCenter(i); 
            
        vector<HistPoint_t> points_signal;
        points_signal.reserve(points.size()); 

        
        for (const auto& point : points) {

            double x = point.x;
            double expect = normal_cdf(x + dx/2., signal_sigma, signal_mean) - normal_cdf(x - dx/2., signal_sigma, signal_mean); 
            expect *= total_stats; 
            points_signal.push_back({ x, expect });
        }
        
        //computes negative log liklihood for signal/noise hypothesis 'mu' 
        auto NLL = [&points_background, &points_signal, &points, signal_mean](const double *par) 
        {
            const auto& mu = par[0];    
            double nll=0.; 

            for (int b=0; b<points_background.size(); b++) {
                double expect = points_background[b].N + (mu * points_signal[b].N); 
                double actual = points[b].N; 

                nll += (actual * log(expect)) - expect - ( actual < 1. ? 0. : (actual*log(actual)) - actual );
            }
            return -nll; 
        };
            
        ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
        
        minimizer->SetMaxFunctionCalls(1e7); 
        minimizer->SetMaxIterations(1e6); 
        minimizer->SetTolerance(1e-3);
        minimizer->SetPrintLevel(0);    
            
        auto f_minimizer = ROOT::Math::Functor(NLL, 1);
        
        minimizer->SetFunction(f_minimizer);

        //set the list of variables
        minimizer->SetVariable(0, "mu", 0., 1e-4);

        bool fit_status = minimizer->Minimize();
        if (!fit_status) continue; 

        const double* par_result = minimizer->X(); 
        const double* par_errors = minimizer->Errors(); 

        double mu_null[] = {0.}; 
        double q0 = -2.*log(NLL(mu_null)/NLL(par_result)); 

        if (par_result[0] < 0.) q0 *= -1.; 

        points_lambda.push_back({ signal_mean, q0 }); 
    }   

    return points_lambda; 
}

const double Z2 = 0.02275; 
const double Z1 = 0.15866;

int pretend_bump()
{
    const long int default_stats = 1e5; 

    const double mean_signal  = +1.0; 
    const double width_signal = 0.05; 

    const double signal_background_ratio = 1e-2; 

    auto hist_q0 = new TH1D("h_q0", "Distribution of q_{0} statistic", 200, -1.,+1.); 
    
    TRandom3 randGen; 

    int nbins = 150; 

    double signal_background = 0.01; //
    double signal_background_dx = 0.0001; //

    
    double min_mult = -0.50;
    double max_mult = +0.50;

    auto hist_mult = new TH2D("h_mult", "Modifier: #sigma #times M,    stats #times 1;Multiplier (M);q_{0}", 50, min_mult, max_mult, 50, -3., +3.); 

    cout << "starting experiment..." << endl; 

    const double signal_mean_range[] = {-2.5, +2.5};

    using namespace ROOT::VecOps; 
    double A{0.}, B{0.};

    auto h_sigtest = TH1D("h_st", "", nbins, -3.,+3);
    for (int i=0; i<3000; i++) {

        const double signal_mean = randGen.Uniform(signal_mean_range[0], signal_mean_range[1]);
        double multiplier = randGen.Uniform(min_mult, max_mult);

        //run an experiment with the given parameters, and reutrn the min_q0 found
        auto collect_data = [&h_sigtest,signal_mean,width_signal,&randGen](const long int N_total, const double SNR, const double signal_sigma)
        {
            h_sigtest.Reset(); 
            for (long int i=0; i<N_total; i++) {
                h_sigtest.Fill( randGen.Gaus() * width_background ); 
            }
            for (long int i=0; i<N_total*SNR; i++) {
                h_sigtest.Fill( randGen.Gaus() * signal_sigma + signal_mean ); 
            }
            auto Q0 = compute_q0(&h_sigtest, signal_sigma); 
            double min_q0 = +1e30; 
            for (auto q : Q0) min_q0 = min<double>(q.N, min_q0); 
            return min_q0; 
        };   

        //first, conduct the 'original' experiment
        double minq0_original = collect_data(default_stats, signal_background, width_signal); 
        
        double minq0_modified = collect_data(default_stats, signal_background, width_signal * exp(multiplier)); 

        double y = minq0_modified - minq0_original; 
        A += multiplier * multiplier; 
        B += y * multiplier; 

        hist_mult->Fill( multiplier, minq0_modified - minq0_original); 
    }

    double coeff = B / A; 
    printf("k_{sigma} = %.8f\n", coeff); 
    auto tf1 = new TF1("coeff", "x * [0]", min_mult, max_mult);
    tf1->SetParameter(0, coeff); 

    new TCanvas; 
    
    hist_mult->SetTitle("min q_{0}(#sigma_{s} * e^{A}, N) - min q_{0}(#sigma_{s}, N);A;#Delta min q_{0}"); 
    hist_mult->Draw("col");     

    auto line = new TLine(min_mult,0., max_mult,0.);
    line->Draw(); 
    line = new TLine(0.,-3., 0.,+3.); 
    line->Draw(); 

    tf1->Draw("SAME"); 
    return 0; 

#if 0 
    for (int t=0; t<10e3; t++) {
        
        auto h = TH1D("hh", "test", 100, -3,+3); 

        for (long int i=0; i<stats; i++) h.Fill( randGen.Gaus() * width_background ); 
        //for (long int i=0; i<stats*signal_background_ratio; i++) h.Fill( randGen.Gaus() * width_signal + mean_signal );

        auto q0 = compute_q0(&h, width_signal); 

        for (auto test : q0) hist_q0->Fill( test.N ); 
    }

    //hist_q0->Draw(); 

    //return the x-value 
    auto Get_inverse_cdf = [&hist_q0](double P)
    {
        auto ax = hist_q0->GetXaxis(); 

        P *= hist_q0->Integral(); 
        double val=0.; 
        int b=1; const int b_max = ax->GetNbins(); 
        while (val < P) {
            if (b >= b_max) break; 
            val += hist_q0->GetBinContent(++b); 
        }   
        return ax->GetBinCenter(b); 
    };


    auto h_signal_background = new TH1D("h_sb", "signal + background", 100, -3,+3); 
    auto h_signal            = new TH1D("h_s",  "signal",              100, -3,+3); 
    auto h_background        = new TH1D("h_b",  "background",          100, -3,+3); 

    for (long int i=0; i<stats; i++) {
        double bg = randGen.Gaus() * width_background; 
        h_signal_background->Fill( bg );
        h_background       ->Fill( bg );
    }
    for (long int i=0; i<stats*signal_background_ratio; i++) {
        double sig = randGen.Gaus() * width_signal + mean_signal; 
        h_signal_background->Fill( sig );
        h_signal->Fill( sig );
    }

    auto q0 = compute_q0(h_signal_background, width_signal); 

    TCanvas *c; 
    c = new TCanvas("c", "Comparison", 1400, 450); 
    c->Divide(2,1);
    gStyle->SetOptStat(0);
    

    h_signal_background->SetTitle("'Real' data;"); 
    c->cd(1); h_signal_background->Draw("E");

    auto stack = new THStack("stack", "Background + Signal Comparison");

    h_background->SetFillColor(kWhite);
    //h_background->Draw("SAME"); 
    stack->Add(h_background);

    h_signal->SetFillColor(kGreen);
    //h_signal_background->Draw(); 
    stack->Add((TH1*)h_signal->Clone("h_sig_copy"));

    c->cd(2); stack->Draw(); 
    //h_background->Draw();
    h_signal->SetFillStyle(3004);
    //h_signal->SetLineColor(kGreen);
    h_signal->DrawCopy("SAME"); 
    //*/
    
    new TCanvas; 

    TGraph *graph; 

    vector<double> q0_x, q0_val; 
    for (auto point : q0) { q0_x.push_back(point.x); q0_val.push_back(point.N); }
    graph = new TGraph(q0_x.size(), q0_x.data(), q0_val.data());
    graph->SetMarkerStyle(kOpenCircle); 
    graph->SetMarkerSize(0.70);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->Draw("ALP"); 

    //how many sigmas to draw? 
    TLine* line; 

    const int n_stddev = 3; 
    for (int i=0; i<n_stddev; i++) {

        double x_sig = Get_inverse_cdf(ROOT::Math::normal_cdf((double)i+1,1.,0.)); 

        line = new TLine( q0_x.front(),-x_sig, q0_x.back(),-x_sig ); line->SetLineStyle(kDashed); line->Draw();
        line = new TLine( q0_x.front(),+x_sig, q0_x.back(),+x_sig ); line->SetLineStyle(kDashed); line->Draw();
    }
    line = new TLine( q0_x.front(),0., q0_x.back(),0. ); line->SetLineStyle(kDotted); line->Draw();
    //return 0; 

    return 0; 
#endif 
}