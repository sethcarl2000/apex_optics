#include <NPolyArray.h> 
#include <stdexcept> 
#include <ROOT/RDataFrame.hxx> 
#include <ApexOptics.h> 
#include <ROOT/RVec.hxx>
#include "include/RDFNodeAccumulator.h"
#include <iostream> 
#include <TSystem.h> 
#include <TFile.h> 
#include <TParameter.h> 

using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 


//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Trajectory_t struct. 
ROOT::RDF::RNode add_branch_from_Trajectory_t(ROOT::RDF::RNode df, const char* branch_in, map<string, double Trajectory_t::*> branches)
{
    const int n_nodes = branches.size() + 1; 
    RVec<ROOT::RDF::RNode> nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Trajectory_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Trajectory_t::coord' 
        nodes.push_back( nodes.back()
            
            .Define(branch_name, [coord](const Trajectory_t& track) { return track.*coord; }, {branch_in})
        ); 
    }
    
    return nodes.back(); 
}


int make_fitdata_with_MC_model(const char* path_infile, const char* path_outfile) 
{
    const char* const here = "make_fitdata_with_MC_model"; 
    
    const char* tree_name = "tracks_fp"; 

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }
    
    //check if the infile could be opened
    auto infile = new TFile(path_infile, "READ"); 
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

    //this polynomial will only be used to compute dp/p, as a funciton of sieve coordinates. 
    const char* path_NPolyArray_dpp   = "data/csv/poly_prod_fp_sv_L_5ord.dat";  // dp/p <======= fp

    const char* path_NPolyArray_sv_q1 = "data/csv/poly_prod_sv_q1_L_5ord.dat";  //   sv => q1 
    const char* path_NPolyArray_q1_fp = "data/csv/poly_prod_q1_fp_L_5ord.dat";  //         q1 => fp 
    const char* path_NPolyArray_sv_fp = "data/csv/poly_prod_sv_fp_L_5ord.dat";  //   sv =======> fp 

    const char* path_NPolyArray_fp_sv = "data/csv/poly_prod_fp_sv_L_5ord.dat";  //   sv <======= fp 
    const char* path_NPolyArray_fp_q1 = "data/csv/poly_prod_fp_q1_L_5ord.dat";  //         q1 <= fp 
    const char* path_NPolyArray_q1_sv = "data/csv/poly_prod_q1_sv_L_5ord.dat";  //   sv <= q1 


    //try to parse the polynomials from files. catch any exceptions, if they occur. 
    NPolyArray parr_fp_dpp; 
    
    NPolyArray parr_sv_q1, parr_q1_fp, parr_sv_fp; 
    NPolyArray parr_fp_sv, parr_fp_q1, parr_q1_sv;  

    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}; 
    const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}; 
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"}; 

    
    try {
        parr_fp_dpp = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_dpp,   {"dpp_sv"}, 4); 

        parr_sv_q1  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_sv_q1, branches_q1, 5); 
        parr_q1_fp  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_q1_fp, branches_fp, 5);
        parr_sv_fp  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_sv_fp, branches_fp, 5);

        parr_fp_sv  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_fp_sv, branches_sv, 4); 
        parr_fp_q1  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_fp_q1, branches_q1, 4);
        parr_q1_sv  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_q1_sv, branches_sv, 5);
    
    } catch (const std::exception& e) {
        
        Error(here, "Something went wrong trying to parse NPolyArray files.\n what(): %s", e.what()); 
        return -1; 
    }

    
    //now, try to construct the RDataFrame. 
    ROOT::RDataFrame *df = nullptr; 

    try {        
        ROOT::EnableImplicitMT();         
        df = new ROOT::RDataFrame("tracks_fp", path_infile); 
    
    } catch (const std::exception& e) {

        Error(here, "Someting went wrong trying to construct the RDataFrame.\n what(): %s", e.what()); 
        return -1; 
    } 

    RDFNodeAccumulator rna(*df); 

    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
    {
        return Trajectory_t{ x, y, dxdz, dydz }; 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 

    rna.Define("Xsv", [&parr_fp_dpp](double x, double y, double dxdz, double dydz, Trajectory_t Xfp)
    {
        return Trajectory_t{
            .x    = x, 
            .y    = y,
            .dxdz = dxdz, 
            .dydz = dydz, 
            .dpp  = parr_fp_dpp.Eval(ApexOptics::Trajectory_t_to_RVec(Xfp)).at(0)
        }; 

    }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv", "Xfp"}); 

    //propagate the 'real' Xsv to XQ1 and Xfp
    rna.Define("Xq1_fwd", [&parr_sv_q1](Trajectory_t Xsv)
    {
        return ApexOptics::RVec_to_Trajectory_t(
            parr_sv_q1.Eval(ApexOptics::Trajectory_t_to_RVec(Xsv))
        ); 
    }, {"Xsv"});

    rna.Define("Xfp_fwd", [&parr_q1_fp](Trajectory_t Xq1)
    {
        return ApexOptics::RVec_to_Trajectory_t(
            parr_q1_fp.Eval(ApexOptics::Trajectory_t_to_RVec(Xq1))
        ); 
    }, {"Xq1_fwd"});

    //propagate the 'real' Xfp to XQ1 and Xsv
    rna.Define("Xq1_rev", [&parr_fp_q1](Trajectory_t Xfp)
    {
        return ApexOptics::RVec_to_Trajectory_t(
            parr_fp_q1.Eval(ApexOptics::Trajectory_t_to_RVec(Xfp))
        ); 
    }, {"Xfp"});

    rna.Define("Xsv_rev", [&parr_q1_sv](Trajectory_t Xq1)
    {
        return ApexOptics::RVec_to_Trajectory_t(
            parr_q1_sv.Eval(ApexOptics::Trajectory_t_to_RVec(Xq1))
        ); 
    }, {"Xq1_rev"});

    //add the forward reconstructions 
    rna = add_branch_from_Trajectory_t(rna.Get(), "Xq1_fwd", {
        {"fwd_x_q1",    &Trajectory_t::x}, 
        {"fwd_y_q1",    &Trajectory_t::y},
        {"fwd_dxdz_q1", &Trajectory_t::dxdz},
        {"fwd_dydz_q1", &Trajectory_t::dydz},
        {"fwd_dpp_q1",  &Trajectory_t::dpp} 
    });

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xfp_fwd", {
        {"fwd_x_fp",    &Trajectory_t::x}, 
        {"fwd_y_fp",    &Trajectory_t::y},
        {"fwd_dxdz_fp", &Trajectory_t::dxdz},
        {"fwd_dydz_fp", &Trajectory_t::dydz}
    }); 


    //add the backward reconstructions
    rna = add_branch_from_Trajectory_t(rna.Get(), "Xq1_rev", {
        {"rev_x_q1",    &Trajectory_t::x}, 
        {"rev_y_q1",    &Trajectory_t::y},
        {"rev_dxdz_q1", &Trajectory_t::dxdz},
        {"rev_dydz_q1", &Trajectory_t::dydz},
        {"rev_dpp_q1",  &Trajectory_t::dpp} 
    });

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_rev", {
        {"rev_x_sv",    &Trajectory_t::x}, 
        {"rev_y_sv",    &Trajectory_t::y},
        {"rev_dxdz_sv", &Trajectory_t::dxdz},
        {"rev_dydz_sv", &Trajectory_t::dydz},
        {"rev_dpp_sv",  &Trajectory_t::dpp} 
    });

    cout << "Making the snapshot..." << flush; 

    //try to make the snapshot. 
    try {
        
        rna.Get().Snapshot("tracks_fp", path_outfile, {
            //real data branches
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",
    
            "x_fp", 
            "y_fp",
            "dxdz_fp",
            "dydz_fp",

            //forward reconstructions
            "fwd_x_q1",
            "fwd_y_q1",
            "fwd_dxdz_q1",
            "fwd_dydz_q1",
            "fwd_dpp_q1",

            "fwd_x_fp",
            "fwd_y_fp",
            "fwd_dxdz_fp",
            "fwd_dydz_fp",

            //reverse reconstructions
            "rev_x_q1",
            "rev_y_q1",
            "rev_dxdz_q1",
            "rev_dydz_q1",
            "rev_dpp_q1",

            "rev_x_sv",
            "rev_y_sv",
            "rev_dxdz_sv",
            "rev_dydz_sv",
            "rev_dpp_sv"
        }); 

    } catch (const std::exception& e) {

        cout << endl; 
        Error(here, "Something went wrong making the output file.\n what(): %s", e.what()); 
        return -1; 
    }
    cout << "done.\nOuptut file created: '" << path_outfile << "'" << endl; 

    //add the 'is_RHRS' parameter to the file
    infile = new TFile(path_outfile, "UPDATE"); 
    
    param_is_RHRS = new TParameter<bool>("is_RHRS", is_RHRS); 
    param_is_RHRS->Write(); 
    
    infile->Close(); 
    delete infile; 

    return 0; 
}