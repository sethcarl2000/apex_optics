#include "include/TestAngleReco.h"
#include "include/RDFNodeAccumulator.h"
#include "include/Get_TParameter_from_TFile.h"
#include <NPolyArrayChain.h>
#include <NPolyArray.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <ROOT/RDataFrame.hxx>
#include <ApexOptics.h>
#include <vector>
#include <ROOT/RVec.hxx>
#include <ApexOptics.h> 
#include <cstdio> 
#include <TGraph.h>
#include <TCanvas.h> 
#include <TF1.h> 
#include <TVirtualPad.h> 

using namespace std; 
using ApexOptics::Trajectory_t; 

//what order should our Xsv => Xsv 'spoiling' polynomial be? 
const int spoiling_poly_order = 1; 

using namespace ROOT::VecOps; 
namespace {
    //random number generator to be used by different fcns 
    TRandom3 *randGen=nullptr; 
};

struct EvaluatedAngleResult_t { 
    TestAngleReco::AngleFitResult_t dxdz, dydz; 
    double dxdz_actual, dydz_actual; 
};

using AngleEvaluator = std::function<EvaluatedAngleResult_t(ROOT::RDF::RNode)>; 

EvaluatedAngleResult_t Spoil_data_and_measure(const double gaus_noise, ROOT::RDF::RNode df, const AngleEvaluator& evaluator);


int test_TestAngleReco_with_montecarlo(const char* path_infile, const char* target_name)
{
    const char* const here = "test_TestAngleReco_with_montecarlo"; 

    randGen = new TRandom3; 
    
    auto param_bool = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS"); 
    if (!param_bool.has_value()) {
        Error(here, "Unable to extract 'is_RHRS' TParameter, from file '%s'.", path_infile); 
        return -1; 
    }
    const bool is_RHRS = param_bool.value(); 

    auto target = ApexOptics::GetTarget(target_name);

    //first, we must make a random matrix

    ROOT::EnableImplicitMT();
    RDFNodeAccumulator rna(ROOT::RDataFrame("tracks_fp", path_infile)); 

    rna.Define("Xsv_rvec", [](double x, double y, double dxdz, double dydz, double dpp)
    {
        return RVec<double>{x,y,dxdz,dydz,dpp}; 
    }, {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"});

    //turn off drawing & printing 
    //TestAngleReco::kDrawing = false; 
    //TestAngleReco::kQuiet   = true; 

    auto evaluator = [is_RHRS, &target](ROOT::RDF::RNode df) {
        
        EvaluatedAngleResult_t result; 
        
        auto result_dydz = TestAngleReco::Evaluate(
            is_RHRS, TestAngleReco::kDydz, df, target,
            4, 13, 
            0, 10, 
            1.25, 0.50, 0.20,
            "dxdz_spoil", "dydz_spoil",
            -0.045, +0.055,
            -0.035, +0.020
        ); 

        auto result_dxdz = TestAngleReco::Evaluate(
            is_RHRS, TestAngleReco::kDxdz, df, target,
            4, 13, 
            0, 10, 
            1.25, 0.50, 1.00,
            "dxdz_spoil", "dydz_spoil",
            -0.045, +0.055,
            -0.035, +0.020
        ); 
        if (!(result_dxdz.has_value() && result_dydz.has_value())) return result; 

        result.dxdz = result_dxdz.value();
        result.dydz = result_dydz.value();

        return result; 
    };

#if 0 
    const double noise_range[] = {0., 3e-4}; 

    vector<double> 
        actual_dxdz, measure_dxdz, 
        actual_dydz, measure_dydz;  
    
    //evaluating... 
    const int n_tests = 5; 

    auto vector_min = [](const vector<double>& v) {
        double min = +1e30; 
        for (double x : v) if (x<min) min=x; 
        return min; 
    };
    auto vector_max = [](const vector<double>& v) {
        double max = -1e30; 
        for (double x : v) if (x>max) max=x; 
        return max; 
    };

    for (int i=0; i<n_tests; i++) {

        cout << "Evaluating..." << flush; 
        auto result = Spoil_data_and_measure(randGen->Uniform(noise_range[0], noise_range[1]), rna.Get(), evaluator); 

        measure_dxdz.push_back(result.dxdz.RMS_position); 
        actual_dxdz.push_back(result.dxdz_actual); 

        measure_dydz.push_back(result.dydz.RMS_position); 
        actual_dydz.push_back(result.dydz_actual); 

        cout << "done." << endl; 
    }

    auto c = new TCanvas("C", "canv", 1400,650); 
    c->Divide(2,1);


    auto g_dxdz = new TGraph(n_tests, actual_dxdz.data(), measure_dxdz.data());
    auto g_dydz = new TGraph(n_tests, actual_dydz.data(), measure_dydz.data());

    g_dxdz->SetTitle("Measured vs actual error of dx/dz;dx/dz RMS (actual);dx/dz RMS (measured)");
    g_dydz->SetTitle("Measured vs actual error of dy/dz;dy/dz RMS (actual);dy/dz RMS (measured)");

    g_dxdz->SetMarkerSize(0.75);
    g_dxdz->SetMarkerStyle(kOpenCircle); 

    g_dydz->SetMarkerSize(0.75);
    g_dydz->SetMarkerStyle(kOpenCircle); 

    TVirtualPad *pad; 

    pad = c->cd(1);
    pad->SetLeftMargin(0.12);
    g_dxdz->Draw("AP");
    g_dxdz->GetYaxis()->SetRangeUser(0., vector_max(measure_dxdz)*1.1);
    auto tf1_1 = new TF1("line_x", "x", vector_min(actual_dxdz), vector_max(actual_dxdz)); 
    tf1_1->SetLineStyle(kDashed); 
    tf1_1->Draw("SAME");
    
    
    pad = c->cd(2);
    pad->SetLeftMargin(0.12); 
    g_dydz->Draw("AP"); 
    g_dxdz->GetYaxis()->SetRangeUser(0., vector_max(measure_dydz)*1.1);
    auto tf1_2 = new TF1("line_y", "x", vector_min(actual_dydz), vector_max(actual_dydz)); 
    tf1_2->SetLineStyle(kDashed); 
    tf1_2->Draw("SAME");


    return 0; 
#endif
    
    const double gaus_noise = 0.;//1e-3; 

    const RVec<double> normalizer = {0.03,0.03,0.03,0.03, 1.}; 
            
    auto spoiling_array = NPolyArrayChain::CreateIdentityNPolyArray(5, spoiling_poly_order); 
    for (int i=0; i<5; i++) { 
        auto poly = spoiling_array.Get_poly(i); 
        
        auto normailizations = poly->Eval_noCoeff(normalizer);
        for (int e=0; e<poly->Get_nElems(); e++) 
            poly->Get_elem(e)->coeff += randGen->Gaus()*gaus_noise / normailizations[e]; 
    }
    //spoiling_array.Print(); 

    
    auto df_spoil = rna.Get()

        .Define("Xsv_spoil", [&spoiling_array](const RVec<double>& Xsv)
        {
            return spoiling_array.Eval(Xsv); 
        }, {"Xsv_rvec"})

        .Define("dxdz_spoil", [](const RVec<double>& v){ return v[2]; }, {"Xsv_spoil"})    
        .Define("dydz_spoil", [](const RVec<double>& v){ return v[3]; }, {"Xsv_spoil"}); 


    auto result_dxdz = TestAngleReco::Evaluate(
        is_RHRS, TestAngleReco::kDxdz, df_spoil, target,
        4, 13, 
        1, 10, 
        1.25, 0.50, 1.00,
        "dxdz_spoil", "dydz_spoil",
        -0.045, +0.055,
        -0.035, +0.020
    ); 

    auto result_dydz = TestAngleReco::Evaluate(
        is_RHRS, TestAngleReco::kDydz, df_spoil, target,
        4, 13, 
        1, 10, 
        1.25, 0.50, 0.25,
        "dxdz_spoil", "dydz_spoil",
        -0.045, +0.055,
        -0.035, +0.020
    ); 

    auto hist_spoil = df_spoil
        .Histo2D({"h_angles", "dx/dz_{sv} vs dy/dz_{sv};dx/dx_{sv};dy/dz_{sv}", 200, -0.05, 0.06, 200, -0.04, 0.03}, "dxdz_spoil", "dydz_spoil"); 
    
    //hist_spoil->DrawCopy("col2"); 

    return 0; 
}

 
EvaluatedAngleResult_t Spoil_data_and_measure(const double gaus_noise, ROOT::RDF::RNode df, const AngleEvaluator& evaluator)
{   
    using namespace ROOT::VecOps; 
    //first, generate a random polynomial array with noise added
    auto spoiling_array = NPolyArrayChain::CreateIdentityNPolyArray(5, spoiling_poly_order); 
    for (int i=0; i<5; i++) { auto poly = spoiling_array.Get_poly(i); 
        for (int e=0; e<poly->Get_nElems(); e++) poly->Get_elem(e)->coeff += randGen->Gaus()*gaus_noise; 
    }

    auto df_spoil = df

        .Define("Xsv_spoil", [&spoiling_array](const RVec<double>& Xsv)
        {
            return spoiling_array.Eval(Xsv); 
        }, {"Xsv_rvec"})

        .Define("dxdz_spoil", [](const RVec<double>& v){ return v[2]; }, {"Xsv_spoil"})    
        .Define("dydz_spoil", [](const RVec<double>& v){ return v[3]; }, {"Xsv_spoil"});
    
    auto result = evaluator(df_spoil); 

    //2-element r-vec of quadrature-summed errors
    auto df_err = df_spoil
        .Define("err_dxdz", [](double dxdz, double dxdz_spoil)
        {
            return pow(dxdz - dxdz_spoil, 2); 
        }, {"dxdz_sv", "dxdz_spoil"})
        .Define("err_dydz", [](double dydz, double dydz_spoil)
        {
            return pow(dydz - dydz_spoil, 2); 
        }, {"dydz_sv", "dydz_spoil"});

    const double count = (double)*df_err.Count(); 

    double err_dxdz = *df_err.Sum("err_dxdz") / count;
    double err_dydz = *df_err.Sum("err_dydz") / count;
    

    result.dxdz_actual = sqrt(err_dxdz); 
    result.dydz_actual = sqrt(err_dydz); 

    printf(
        "Results: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        " Dx/dz (lab-vertical):\n"
        " measured: %.4e\n"
        " actual:   %.4e\n"
        "\n"
        " Dy/dz (lab-horizontal):\n"
        " measured: %.4e\n"
        " actual:   %.4e\n",
        result.dxdz.RMS_position, 
        result.dxdz_actual, 
        result.dydz.RMS_position, 
        result.dydz_actual
    ); cout << flush; 

    return result; 
}; 