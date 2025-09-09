#include <NPolyArray.h> 
#include <stdexcept> 
#include <ROOT/RDataFrame.hxx> 
#include <ApexOptics.h> 
#include <ROOT/RVec.hxx>
#include "include/RDFNodeAccumulator.h"
#include <iostream> 

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

    //this polynomial will only be used to compute dp/p, as a funciton of sieve coordinates. 
    const char* path_NPolyArray_dpp   = "data/csv/poly_prod_fp_sv_L_6ord.dat"; 
    const char* path_NPolyArray_sv_q1 = "data/csv/poly_prod_sv_q1_L_6ord.dat"; 
    const char* path_NPolyArray_q1_fp = "data/csv/poly_prod_q1_fp_L_6ord.dat"; 
    
    //try to parse the polynomials from files. catch any exceptions, if they occur. 
    NPolyArray parr_fp_dpp, parr_sv_q1, parr_q1_fp; 
    
    try {
        parr_fp_dpp = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_dpp,   {"dpp_sv"}, 4); 
        parr_sv_q1  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_sv_q1, {"x_q1", "y_q1", "dxdz_q1", "dydz_q1", "dpp_q1"}, 5); 
        parr_q1_fp  = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray_q1_fp, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}, 5); 
    
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

    rna.Define("Xq1_reco", [&parr_sv_q1](Trajectory_t Xsv)
    {   
        RVec<double>&& Xsv_rvec = ApexOptics::Trajectory_t_to_RVec(Xsv); 

        Trajectory_t&& Xq1 = ApexOptics::RVec_to_Trajectory_t(parr_sv_q1.Eval(Xsv_rvec)); 

        return Xq1; 

    }, {"Xsv"}); 

    rna.Define("Xfp_reco", [&parr_q1_fp](Trajectory_t Xq1)
    {   
        RVec<double>&& Xq1_rvec = ApexOptics::Trajectory_t_to_RVec(Xq1); 

        Trajectory_t&& Xfp = ApexOptics::RVec_to_Trajectory_t(parr_q1_fp.Eval(Xq1_rvec)); 

        return Xfp; 

    }, {"Xq1_reco"}); 

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xq1_reco", {
        {"reco_x_q1",    &Trajectory_t::x}, 
        {"reco_y_q1",    &Trajectory_t::y},
        {"reco_dxdz_q1", &Trajectory_t::dxdz},
        {"reco_dydz_q1", &Trajectory_t::dydz},
        {"reco_dpp_q1",  &Trajectory_t::dpp} 
    }); 

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xfp_reco", {
        {"reco_x_fp",    &Trajectory_t::x}, 
        {"reco_y_fp",    &Trajectory_t::y},
        {"reco_dxdz_fp", &Trajectory_t::dxdz},
        {"reco_dydz_fp", &Trajectory_t::dydz}
    }); 



    cout << "Making the snapshot..." << flush; 

    //try to make the snapshot. 
    try {
        
        rna.Get().Snapshot("tracks_fp", path_outfile, {
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",

            "reco_x_q1",
            "reco_y_q1",
            "reco_dxdz_q1",
            "reco_dydz_q1",
            "reco_dpp_q1",

            "reco_x_fp",
            "reco_y_fp",
            "reco_dxdz_fp",
            "reco_dydz_fp",

            "x_fp", 
            "y_fp",
            "dxdz_fp",
            "dydz_fp"
        }); 

    } catch (const std::exception& e) {

        cout << endl; 
        Error(here, "Something went wrong making the output file.\n what(): %s", e.what()); 
        return -1; 
    }
    cout << "done.\nOuptut file created: '" << path_outfile << "'" << endl; 

    return 0; 
}