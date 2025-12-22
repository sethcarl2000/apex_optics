//std lib headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include <stdexcept> 

//ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TVector3.h> 
#include <TSystem.h>
#include <TParameter.h> 
#include <TStyle.h>
#include <TCanvas.h> 

//APEX optics headers
#include <ApexOptics.h> 
#include <PolynomialCut.h> 
#include "include/RDFNodeAccumulator.h"
#include "include/ChainedOpticsModel.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include "include/Add_TParameter_to_TFile.h"
#include <NPolyArrayChain.h> 

using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 
using ApexOptics::OpticsTarget_t; 

namespace {
    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

    const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
    const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

    const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
    const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};

    inline bool is_nan(double _val) { return _val != _val; }
 
};



//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int cleanup_optics_replay(  const bool is_RHRS          = false,
                            const char* path_infile     = "data/replay/replay.4768.root",
                            const char* target_name     = "V2",
                            const char* path_outfile    = "data/replay/real_L_V2.root",
                            const char* path_polycut    = "",//"data/replay/fits_6Nov/polycut-V3.dat",
                            const char* tree_name       = "track_data" ) 
{
    const char* const here = "cleanup_optics_replay";
    
    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df(tree_name, path_infile); 

    //minimum sum of all cerenkov ADCs. 
    const double min_cerenkov_sum = 1.5e3;  
    const double min_Esh_Eps_sum  = 0.45; 

    const double central_momentum = 1104.;  //in MeV

    //this accounts for the fact that the y-raseter was not calibrated correctly for these runs. 
    const double y_correction = +1.72835e-3;

    //select the target chosen
    OpticsTarget_t target; 
    
    if (string(target_name)=="") { 

        target = ApexOptics::GetTarget("none"); 
    
    } else {

        try { target = ApexOptics::GetTarget(string(target_name)); } 

        catch (const std::exception& e) {

            Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
            return -1; 
        }
    }
    
    printf("Optics target chosen is '%s'.\n", target.name.c_str()); 


    //try to create the polynomial cu
    PolynomialCut* polycut = nullptr; 

    const bool use_polycut = (string(path_polycut) != ""); 

    
    if (use_polycut) {
        printf("We are applying a polynomial cut to the data: '%s'\n", path_polycut); 
    } else {
        printf("We are NOT going to apply a polynomial cut to filter events - none was specified.\n"); 
    }

    if (use_polycut) {
        try { 
            cout << "using polycut file: " << path_polycut << endl; 

            polycut = new PolynomialCut; 
            polycut->Parse_dbfile(path_polycut); 
        } 
        catch (const PolynomialCut::DBFileException& e) {
            Error(here, "Problem parsing cutfile. Exception:\n %s", e.what());
            return -1;  
        }
    }
    const char* cutbranch_x = "z_reco_vertical"; 
    const char* cutbranch_y = "dydz_sv"; 

    //try to initialize the optics model
    ChainedOpticsModel* model = new ChainedOpticsModel(is_RHRS); 

    try {
        model->CreateChainRev({ // sv <= [Poly] q1-fwd <= _Poly_ <= fp 
            {"data/poly/fits_6Nov/V123_fp_sv_4ord.dat", branches_sv, 4}
        }); 

        model->CreateChainFwd({ // sv => [Poly] => fp
            {"data/csv/poly_prod_sv_fp_L_4ord.dat",         branches_fp, 5},
            {"data/poly/fits_6Nov/V123_fp-fwd_fp_1ord.dat", branches_fp, 4} 
        }); 
    
    } catch (const std::exception& e) {
        Error(here, "Something went wrong trying to build the polynomial chains.\n what(): %s", e.what()); 
        return -1; 
    }
    
    //we need to reverse-engineer some informaiton to 'fix' the y-raster. 
    //  this constant here (0.070302882883/2) is the measured ratio: (raster-timing-delay / 1-rast.-period).
    //  computed by looking at the graphs of h-wire runs, we can look at the 'peaks' when we plot y-rast vs. y-BPM  
    //                  /\<=|
    //      ___________/  \__________________________ <= bpm
    //                      |=>/\
    //      __________________/  \___________________ => bpm 
    //
    //                   |------|  <= 2*raster-timing-delay
    //
    //      |----------------------------------------| <= total raster amplitude
    //
    const double y_raster_amplitude_offset = 0.070302882883 / 2.; 
    
    const double y_min = *df.Define("y_vtx", [](TVector3 vtx){ return vtx.y(); }, {"react_vertex"}).Min("y_vtx"); 
    const double y_max = *df.Define("y_vtx", [](TVector3 vtx){ return vtx.y(); }, {"react_vertex"}).Max("y_vtx"); 

    //the actual amount by which the y-vtx value will be changed. 
    const double y_vtx_offset = y_raster_amplitude_offset * ( y_max - y_min ); 
    
    auto df_rast_info = df

        .Define("y_vtx", [](TVector3 vtx){ return vtx.y(); }, {"react_vertex"})

        //average y-value of both BPMs. 
        .Define("y_BPM", [](TVector2 BPMA, TVector2 BPMB)
        {
            return ( BPMA.Y() + BPMB.Y() )/2.; 
        }, {"r_BPMA", "r_BPMB"})

        //this 'parameter' is a way to normalize the y-raster value on the interval [0,1), 
        .Define("raster_parameter", [=](double y_vtx)
        {
            return (y_vtx - y_min) / (y_max - y_min); 
        }, {"y_vtx"});

    const double min_rast_y_bpm = *df_rast_info 
        //the events with the smallest raster; what is their average y-bpm value? 
        .Filter([](double rast_param){ return rast_param < 0.01; }, {"raster_parameter"})
        .Mean("y_BPM"); 

    const double max_rast_y_bpm = *df_rast_info 
        //the events with the smallest raster; what is their average y-bpm value? 
        .Filter([](double rast_param){ return rast_param > 0.99; }, {"raster_parameter"})
        .Mean("y_BPM"); 
    
    //now, we're ready to create a 'function' which can tell us if an event is 'above' or 'below' this line, which
    // will tell us whether we should add or subtract the phase. 
    auto Phase_correction = [&](double y_vtx, double y_bpm) 
    {
        double raster_param = (y_vtx - y_min) / (y_max - y_min); 

        double phase_line_value = min_rast_y_bpm + (max_rast_y_bpm - min_rast_y_bpm) * raster_param; 

        return (y_bpm > phase_line_value) ? -y_vtx_offset : y_vtx_offset;  
    }; 

    cout << "Creating output file..." << flush; 

    const char* br_n_tracks = (is_RHRS ? "R_tr_n" : "L_tr_n"); 

    //handle the fact that the branches have different names between the right / left arms
    const char* branch_cerenkov_adc     = (is_RHRS ? "R_cer_a_c" : "L_cer_a_c" ); 
    const char* branch_preshower_adc    = (is_RHRS ? "R_ps_a_c"  : "L_prl1_a_c"); 
    const char* branch_shower_adc       = (is_RHRS ? "R_sh_a_c"  : "L_prl2_a_c"); 

    //select events with 1 (and only 1) track reconstructed. 
    printf(
        "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "input file:            '%s'\n"
        "output file:           '%s'\n"
        "n. events (pre-cut):   %llu\n",
        path_infile, 
        path_outfile,
        *df.Filter([](int n){ return n==1; }, {br_n_tracks}).Count()
    );
    cout << flush; 
    
    RDFNodeAccumulator rna( df.Filter([](int n){ return n==1; }, {br_n_tracks}) ); 
    
    //rename the focal-plane coords
    vector<pair<string, string>> fp_coords = {
        {"x_fp",    "_tr_tra_x" },
        {"y_fp",    "_tr_tra_y" },
        {"dxdz_fp", "_tr_tra_th"},
        {"dydz_fp", "_tr_tra_ph"},
    }; 

    for (auto& branch : fp_coords) { 

        //
        const char* br_new = branch.first.c_str(); 

        //add an 'R' or 'L' to the start of the branch's name, to indicate which arm its refering to. 
        branch.second = (is_RHRS ? "R" : "L") + branch.second;
        const char* br_old = branch.second.c_str(); 

        rna.Define(br_new, [](const RVec<double>& v){ return v.at(0); }, {br_old});
    }

    //now, actually do some calculations.

    //perform the conversion from 'transport' coordinates (tra) to 'focal-plane' coordinates (fp)
    rna.Overwrite("dxdz_fp", [](double x_tra, double dxdz_tra)
        {
            return dxdz_tra - x_tra/6.; 
        }, {"x_fp", "dxdz_fp"});

    //apply the phase correction to the vertex
    rna.Define("y_BPM", [](TVector2 BPMA, TVector2 BPMB){ return (BPMA.Y() + BPMB.Y())/2.; }, {"r_BPMA", "r_BPMB"}); 

    rna.Define("y_hcs_corrected", [y_correction, &Phase_correction](TVector3 vtx, double y_BPM)
        {
            return vtx.y() + Phase_correction(vtx.y(), y_BPM) + y_correction; 
        }, {"react_vertex", "y_BPM"});

    rna.Define("position_vtx", [&](TVector3 vtx, double y_BPM)
        {
            //fix the y-value, unless this is an H-wire run, in which case we should NOT fix it. 
            if (!is_nan(target.x_hcs)) { vtx[0] = target.x_hcs; }  
            if (!is_nan(target.y_hcs)) { vtx[1] = target.y_hcs; } else { vtx[1] += Phase_correction(vtx.y(), y_BPM) + y_correction; }  
            if (!is_nan(target.z_hcs)) { vtx[2] = target.z_hcs; }  

            return vtx; 

        }, {"react_vertex", "y_BPM"}); 

    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{x, y, dxdz, dydz};    
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_reco", [&model](Trajectory_t Xfp, TVector3 vtx_hcs)
        {
            return model->Compute_Xsv(Xfp, vtx_hcs);
        }, {"Xfp", "position_vtx"}); 
        

    //now, add the constant y-offset
    rna.Overwrite("position_vtx", [y_correction, &target](TVector3 vtx_hcs, double y_hcs_corrected)
        { 
            //only do this if this ISN'T a h-wire run
            if (target.y_hcs != target.y_hcs) vtx_hcs[1] = y_hcs_corrected; 
            return vtx_hcs; 
        }, {"position_vtx", "y_hcs_corrected"}); 

    rna.Define("position_vtx_scs", [is_RHRS](TVector3 vtx_hcs)
        {
            return ApexOptics::HCS_to_SCS(is_RHRS, vtx_hcs); 
        }, {"position_vtx"}); 


    //perform a very basic cut on sum of cerenkov ADCs 
    rna.Define("cer_sum", [](const ROOT::RVec<double>& cer_a_c)
        {
            double val=0.; 
            for (double adc : cer_a_c) if (fabs(adc) < 5e4) val += adc; 
            return val; 
        }, {branch_cerenkov_adc}); 

        //for now, this is actually only the MAXIMUM block. 
    rna.Define("E_ps", [](const ROOT::RVec<double>& blocks_adc)
        {
            double val=0.; 
            for (double adc : blocks_adc) val = max<double>(val, adc);
            return val; 
        }, {branch_preshower_adc});

        //for now, this is actually only the MAXIMUM block. 
    rna.Define("E_sh", [](const ROOT::RVec<double>& blocks_adc)
        {
            double val=0.; 
            for (double adc : blocks_adc) val = max<double>(val, adc); //val += adc; 
            return val; 
        }, {branch_shower_adc}); 

       
        //ratio of shower energy over momentum 
    rna.Define("E_sh_p_ratio", [central_momentum](double E, Trajectory_t Xsv)
        {
            return E / (central_momentum * ( 1. + Xsv.dpp ));  
        }, {"E_sh", "Xsv_reco"}); 

        //ratio of pre-shower energy over momentum 
    rna.Define("E_ps_p_ratio", [central_momentum](double E, Trajectory_t Xsv)
        {
            return E / (central_momentum * ( 1. + Xsv.dpp ));  
        }, {"E_ps", "Xsv_reco"}); 
    
    //add the sieve branches
    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"x_sv",    &Trajectory_t::x},
        {"y_sv",    &Trajectory_t::y},
        {"dxdz_sv", &Trajectory_t::dxdz},
        {"dydz_sv", &Trajectory_t::dydz},
        {"dpp_sv",  &Trajectory_t::dpp}
    }); 

    //make a cut on z_vertical & dydz, if we're using a polycut
    if (use_polycut) {

        rna.Define("z_reco_vertical", [is_RHRS](TVector3 vtx, Trajectory_t Xsv)
        {
            auto Xhcs = ApexOptics::SCS_to_HCS(is_RHRS, Xsv);
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz; 
        }, {"position_vtx", "Xsv_reco"}); 

        //apply the polynomial cut
        rna = rna.Get().Filter([&polycut](double dydz_sv, double z_reco_vertical)
            {
                return polycut->IsInside(z_reco_vertical, dydz_sv); 
            }, {"dydz_sv", "z_reco_vertical"});
    }
    
    //do a cut on cerenkov values
    auto df_out = rna.Get()

        //do a cut on cerenkov sum
        .Filter([min_cerenkov_sum](double cer_sum){ return cer_sum > min_cerenkov_sum; }, {"cer_sum"})

        .Filter([min_Esh_Eps_sum](double Eps_p, double Esh_p){ return Eps_p + Esh_p > min_Esh_Eps_sum; }, {"E_ps_p_ratio", "E_sh_p_ratio"}); 
    
    printf(
        "n. events (post-cut):  %llu\n"
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
        *df_out.Count()
    );  
    

    df_out.Snapshot("tracks_fp", path_outfile, {

            //focal-plane coordinates
            "x_fp",         
            "y_fp",
            "dxdz_fp",
            "dydz_fp",

            //sieve coordinates
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",

            //vertex/target information
            "position_vtx",         //the position vtx, in APEX HCS
            "position_vtx_scs",     //the position vtx, in SCS
            "y_hcs_corrected",      //the y-position from raster info.
                                    // this will be no different from position_vtx.y() for non-H-wire runs, but for H-wire runs 
                                    // (where position_vtx.y() is overwritten with the real H-wire pos from survey data, this is 
                                    // a way to keep the raster info.) 
            "y_BPM",                // average of both BPMs for this event.
            "r_BPMA",               // event-wise BPMA info
            "r_BPMB",               // event-wise BPMB info
            "Raster2_current_x",    // raw raster current - x
            "Raster2_current_y",    // raw raster current - y  

            //pid information
            "cer_sum",      
            "E_sh",
            "E_ps",
            "E_ps_p_ratio",
            "E_sh_p_ratio"
        });

    cout << "done." << endl; 

    cout << "Adding parameters..." << flush; 
    
    auto file = new TFile(path_outfile, "UPDATE"); 

    Add_TParameter_to_TFile("is_RHRS",          is_RHRS); 
    Add_TParameter_to_TFile("min_cerenkov_sum", min_cerenkov_sum); 
    Add_TParameter_to_TFile("min_Esh_Eps_sum",  min_Esh_Eps_sum); 

    file->Close(); 
    
    cout << "done." << endl; 

    return 0; 
}