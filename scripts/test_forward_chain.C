
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector> 
#include <map>
#include <utility> 
#include "include/RDFNodeAccumulator.h"
#include "include/TestAngleReco.h"
#include "include/ChainedOpticsModel.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include <TParameter.h>
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TGraph.h> 
#include <TCanvas.h> 
#include <NPolyArrayChain.h> 
#include <limits> 
#include <Interactive3dHist.hxx>
#include <TGClient.h> 
#include <TStopwatch.h> 

using namespace std; 
using namespace ROOT::VecOps; 
using ApexOptics::OpticsTarget_t; 
using ApexOptics::Trajectory_t; 

#define USE_MULTITHREADDING false

namespace {
        
    constexpr double range_dxdz[] = {-0.05, +0.05}; 
    constexpr double range_dydz[] = {-0.04, +0.04}; 

    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

    const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
    const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

    const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
    const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
}; 

//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_forward_chain( const char* path_infile ="data/replay/real_L_V2.root",
                        const char* target_name ="V2",
                        const char* path_dbfile ="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
                        const char* tree_name   ="tracks_fp" ) 
{

    const char* const here = "test_forward_chain"; 

    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }
    
    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(here, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the proper TTree
    TTree* tree = (TTree*)infile->Get(tree_name); 
    if (!tree) {
        Error(here, "could not find TTree '%s'", tree_name); 
        return 1; 
    }

    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 
       
    delete tree; 
    infile->Close(); 
    delete infile; 

    //try to get the target we need. 
    OpticsTarget_t target; 
    try { 
        target = ApexOptics::GetTarget(string(target_name)); 
    
    } catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }

    //in the syntax below, an '<=' arrow represents an input or output of a polynomial. 
    // if a polynomial is written '[Poly]' then it is trained soley on monte-carlo data. 
    // if a polynomial is written '_Poly_', then it is trained on real data, which is fed into monte-carlo polynomials 
    //      (in order to reconstruct XQ1).
    // uncomment whichever configuration you want to use.
    // 

    struct NPolyArrayConstructor_t { string path{}; vector<string> coords{}; int input_DoF{0}; }; 

    
    NPolyArray parr_correction = ApexOptics::Parse_NPolyArray_from_file("data/poly/fits_21Dec/V123_sv-yhcs_sv_L_4ord.dat", {
        "x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"
    }, 6);

    ChainedOpticsModel* model = new ChainedOpticsModel(is_RHRS); 
    model->CreateChainRev({

        // sv <= [Poly] <= fp
        //{"data/poly/fits_5Dec/V123_fp_sv_4ord.dat", branches_sv, 4} //*/ 
        {"data/poly/fits_21Dec/V123_fp_sv_L_4ord.dat", branches_sv, 4} //*/ 

        /*/ sv <= [Poly] <= fp-fwd <= _Poly7_ <= fp
        {"data/csv/poly_fits_fp_fp-fwd_L_4ord.dat", branches_fwd_fp, 4}, 
        {"data/csv/poly_prod_fp_sv_L_4ord.dat",     branches_sv,     4} //*/ 

        /*/ sv <= [Poly] q1 <= [Poly] <= fp-fwd <= _Poly_ <= fp 
        {"data/csv/poly_fits_fp_fp-fwd_L_4ord.dat", branches_fwd_fp, 4}, 
        {"data/csv/poly_prod_fp_q1_L_4ord.dat",     branches_q1,     4},
        {"data/csv/poly_prod_q1_sv_L_4ord.dat",     branches_sv,     5} //*/ 

        /*/ sv <= [Poly] q1-fwd <= _Poly_ <= fp 
        {"data/poly/fits_6Nov/V123_fp_q1-fwd_4ord.dat", branches_fwd_q1, 4}, 
        {"data/csv/poly_prod_q1_sv_L_4ord.dat",     branches_sv,     5} //*/ 

        /*/ sv <= _Poly_ q1-rev <= [Poly] <= fp 
        {"data/csv/poly_prod_fp_q1_L_4ord.dat",     branches_q1,     4}, 
        {"data/csv/poly_fits_q1-rev_sv_L_4ord.dat", branches_sv,     5} //*/ 
    }); 
    model->CreateChainFwd({

        /*/ sv => [Poly] => fp
        {"data/csv/poly_prod_sv_fp_L_3ord.dat", branches_fp, 5} //*/ 

        /*/ sv => _Poly_ => fp
        {"data/poly/fits_21Dec/V123_sv_fp_L_4ord.dat", branches_fp, 5} //*/ 

        /*/ sv => [Poly] q1-fwd => _Poly_ => fp 
        {"data/poly/fits_6Nov/V123_fp_q1-fwd_4ord.dat", branches_fwd_q1, 4}, 
        {"data/csv/poly_prod_q1_sv_L_4ord.dat",         branches_sv,     5} //*/ 

        /*/ sv => _Poly_ q1-rev => [Poly] => fp 
        {"data/csv/poly_prod_fp_q1_L_4ord.dat",     branches_q1,     4}, 
        {"data/csv/poly_fits_q1-rev_sv_L_4ord.dat", branches_sv,     5} //*/ 
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_L_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_21Dec/V123_fp-fwd_fp_L_3ord.dat", branches_fp, 4} 
        //*/ 
    
    }); 

    /*NPolyArrayChain chain_rev, chain_fwd; 
    
    //try to parse all polys from files
    try {
        
        for (const auto& path_and_coord : path_and_coords_rev) {

            const char* path  = path_and_coord.path.c_str(); 
            const auto  coord = path_and_coord.coords;
            const int   DoF   = path_and_coord.input_DoF;  

            chain_rev.AppendArray( ApexOptics::Parse_NPolyArray_from_file(path, coord, DoF) ); 
        }

        for (const auto& path_and_coord : path_and_coords_fwd) {

            const char* path  = path_and_coord.path.c_str(); 
            const auto  coord = path_and_coord.coords;
            const int   DoF   = path_and_coord.input_DoF;  

            chain_fwd.AppendArray( ApexOptics::Parse_NPolyArray_from_file(path, coord, DoF) ); 
        }

    } catch (const std::exception& e) {

        Error(here, "Something went wrong parsing one of the NPolyArrays.\n what(): %s", e.what()); 
        return -1; 
    }*/


    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 

#if USE_MULTITHREADDING
    ROOT::EnableImplicitMT(); 
    Info(__func__, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 
#else 
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
    Info(__func__, "Multi-threadding is disabled; running in single-thread mode."); 
#endif

    RDFNodeAccumulator rna(tree_name, path_infile); 

    //total events in this RDataFrame
    const long int n_events = *rna.Get().Count(); 

    printf("Processing %li events from file '%s'...", n_events, path_infile); cout << endl; 

    //probably not the most elegant way to do this, but here we are. 
    const long int max_benchmark_events = 7e5; //n_events; 

    rna = rna.Get().Range(max_benchmark_events); 
    
    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ x, y, dxdz, dydz };
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Overwrite("Xsv_first_guess", [&model](Trajectory_t Xfp)
        {
            return model->Compute_Xsv_first_guess(Xfp); 
        }, {"Xfp"});

    rna.Overwrite("Xsv_corrected", [&parr_correction](Trajectory_t Xsv_fg, TVector3 vtx)
        {
            RVec<double> Xsv_fg_rvec = ApexOptics::Trajectory_t_to_RVec(Xsv_fg); 
            Xsv_fg_rvec.push_back(vtx.y()); 
            
            RVec<double>&& Xsv_corrected = parr_correction.Eval(Xsv_fg_rvec); 
            return ApexOptics::RVec_to_Trajectory_t(Xsv_corrected); 

        }, {"Xsv_first_guess", "position_vtx"}); 

    rna.Overwrite("Xsv_reco", [&model](Trajectory_t Xfp, TVector3 vtx_hcs)
        {
            return model->Compute_Xsv(Xfp, vtx_hcs);
        }, {"Xfp", "position_vtx"});

    rna.Define("Xhcs_reco", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
        
        }, {"Xsv_reco"});

    rna.Overwrite("z_reco_horizontal", [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {   
            return - ( Xhcs.y - vtx.y() ) / Xhcs.dydz; 
        }, {"Xhcs_reco", "position_vtx"});

    rna.Overwrite("z_reco_vertical",   [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 
        }, {"Xhcs_reco", "position_vtx"});

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }
    
    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xsv_first_guess", {
        {"fg_x_sv",    &Trajectory_t::x},
        {"fg_y_sv",    &Trajectory_t::y},
        {"fg_dxdz_sv", &Trajectory_t::dxdz},
        {"fg_dydz_sv", &Trajectory_t::dydz}
    }); 

    rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        "reco_x_sv",
        "reco_y_sv",
        "reco_dxdz_sv",
        "reco_dydz_sv"
    }); 

    rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_corrected", {
        "corr_x_sv",
        "corr_y_sv",
        "corr_dxdz_sv",
        "corr_dydz_sv"
    }); 

    rna.Define("y_pos", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}); 

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    const double min_cerenkov_sum = 1.5e3; 
    const double min_Esh_Eps_sum = 0.5; 

    
    auto hist_z_dy_nocut = rna.Get()
        .Histo2D<double>({"h_z_dy_nocut", "No PID cut;z_{tg};dy/dz_{sv}", 200, -0.4, +0.4, 200, -0.04, 0.03}, "z_reco_vertical", "reco_dydz_sv"); 

   

    auto hist_z = rna.Get()
        .Histo1D<double>({"h_z", "Z_{hcs} reconstruction;z_{tg};", 200, -0.4, +0.4}, "z_reco_vertical"); 

    
    //create both histograms
    auto hist_xy = rna.Get()
        .Histo2D({"h_xy",     "x_{sv} vs y_{sv} (linear dp/p correction);x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "reco_x_sv", "reco_y_sv");
    
    auto hist_angles = rna.Get()
        .Histo2D({"h_angles", "dx/dz_{sv} vs dy/dz_{sv} (linear dp/p correction);dx/dx_{sv};dy/dz_{sv}", 200, range_dxdz[0],range_dxdz[1], 200, range_dydz[0],range_dydz[1]}, "reco_dxdz_sv", "reco_dydz_sv"); 
    

    //create both histograms
    auto hist_fg_xy = rna.Get()
        .Histo2D({"h_fg_xy",     "x_{sv} vs y_{sv};x_{sv};y_{sv}", 200, -0.040, 0.045, 200, -0.045, 0.010}, "fg_x_sv", "fg_y_sv");

    auto hist_fg_angles = rna.Get()
        .Histo2D({"h_fg_angles", "dx/dz_{sv} vs dy/dz_{sv};dx/dx_{sv};dy/dz_{sv}", 200, range_dxdz[0],range_dxdz[1], 200, range_dydz[0],range_dydz[1]}, "fg_dxdz_sv", "fg_dydz_sv");

    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    /*/Measure lab-vertical angle (dxdz)
    TestAngleReco::Evaluate(is_RHRS, TestAngleReco::kDxdz, rna.Get(), target, 
        4,13, 
        1,8, 
        1.35, 0.50, 1.00,
        "reco_dxdz_sv", "reco_dydz_sv",
        false, 
        range_dxdz[0], range_dxdz[1],
        range_dydz[0], range_dydz[1]
    ); 

    //Measure lab-horizontal angle (dydz)
    TestAngleReco::Evaluate(is_RHRS, TestAngleReco::kDydz, rna.Get(), target, 
        4,13, 
        2,10, 
        1.25, 0.65, 0.20,
        "reco_dxdz_sv", "reco_dydz_sv",
        true, 
        range_dxdz[0], range_dxdz[1],
        range_dydz[0], range_dydz[1]
    ); 

    return 0;  //*/ 

    TCanvas *c; 
    c = new TCanvas("c_z_reco", c_title, 700, 700); 

    //this is the timing benchmark.
    const int n_benchmarks = 1; 
    //const long int max_benchmark_events = 1e5; //n_events; 

    vector<double> vec_times, vec_n_events; 

    cout << "Starting benchmark.........." << endl; 

    for (int i=0; i<n_benchmarks; i++) {

        const long int n_benchmark_events = (i+1)*max_benchmark_events/n_benchmarks; 

        printf("Benchmarking %li events...", n_benchmark_events); cout << flush;  

        TStopwatch timer; 

        double dummy_sum = *rna.Get()
            .Range(n_benchmark_events)
            .Define("dummy_sum", [](const Trajectory_t& Xsv){ return Xsv.x; }, {"Xsv_reco"}) 
            .Sum("dummy_sum");
            
        double elapsed_time = timer.RealTime(); 

        printf("done. time: %.3f s", elapsed_time); cout << endl; 

        vec_times.push_back(elapsed_time); 
        vec_n_events.push_back((double)n_benchmark_events); 
    }

    double time_per_event;
    if (n_benchmarks > 1) {
        //now, find the slope
        double sum_xx{0.}, sum_xy{0.}, sum_x{0.}, sum_y{0.}; 
        for (int i=0; i<vec_times.size(); i++) {
            double x = vec_n_events[i];
            double y = vec_times[i]; 
            sum_xx += x*x; 
            sum_xy += x*y; 
            sum_x  += x; 
            sum_y  += y; 
        } 
        time_per_event = ( sum_xy - sum_x*sum_y )/( sum_xx - sum_x*sum_x );
        
        auto g_benchmark = new TGraph(vec_times.size(), vec_n_events.data(), vec_times.data()); 

        g_benchmark->SetTitle("Timing benchmark;N. benchmark events.;time (s)"); 
        g_benchmark->SetMarkerStyle(kOpenCircle);
        g_benchmark->Draw("PLC"); 

    } else {

        time_per_event = vec_times.back()/vec_n_events.back(); 
    }
    
    printf("extrapolating from slope, the time per event is:    %.4f ms/event \n", time_per_event*1.e3);

    return 0; 
    

    hist_z->DrawCopy(); 


    c = new TCanvas("c1", c_title, 1200, 600); 
    c->SetLeftMargin(0.12); c->SetRightMargin(0.05); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1); hist_fg_xy->DrawCopy("col2"); 
    c->cd(2); hist_xy->DrawCopy("col2");
    
    c = new TCanvas("c2", c_title, 1200, 600); 
    c->Divide(2,1, 0.01,0.01); 

    c->cd(1); hist_fg_angles->DrawCopy("col2"); 
    c->cd(2); hist_angles->DrawCopy("col2"); 

    return 0; 
}