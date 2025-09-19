
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include "include/RDFNodeAccumulator.h"
#include <TParameter.h>
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h> 
#include <TRandom3.h> 
//this last one should eventually be moved to the ApexOptics namespace, as it has general use outside of the 'isolate_sieveholes' app. 
#include <SieveHole.h>
#include <TLine.h> 
#include <TAxis.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
#include <memory>
#include <TRandom3.h> 
#include <TBox.h> 
#include "include/TestAngleReco.h"

using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 
using ApexOptics::OpticsTarget_t; 

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

//the first in the pair is the x-value,
//and the second in the pair is the y-value. 
//this is evenly-spaced x-points, with the y-points being that of a normalizedd gaussian with sigma=1. 
//the first point is x=0, the last point goes out to x > 0.; 
const size_t gauss_points = 200; 
const double gauss_max_x  = 10.; 

inline double fcn_gauss(double sigma, double x) {
    return 0.398942280401 * exp( -0.5 * pow( x/sigma, 2 ) ) / sigma; 
}

//semicirle with height 1 and radius 'rad'
inline double fcn_semicircle(double rad, double x) {
    return fabs(x) > rad ? 0. : sqrt( 1. - pow( x/rad, 2 ) );  
} 

//a gaussian fit-function which is an approximate convolution between a gaussian and semicircle
double gaussian_semicircle(double *x, double *par) {
    //parameters works this way: 
    // par[0] = x0 - center of semicircle (and of whole fit)
    // par[1] = radius of semicircle
    // par[2] = sigma of gaussian convolution
    // par[3] = overall magnitude ( =1 when sigma=0 )
    // par[4] = const. offset
    double val=0.; 
    
    double dx = x[0] - par[0]; 

    const double& rad   = fabs(par[1]); 
    const double& sigma = fabs(par[2]); 

    double gauss_dx = gauss_max_x / ((double)gauss_points - 1); 

    double t=0.; 

    //perform the numerical convolution
    for (size_t i=0; i<gauss_points; i++) {

        val += gauss_dx * fcn_gauss(1., t ) * ( fcn_semicircle(rad, dx - t*sigma ) + fcn_semicircle(rad, dx + t*sigma ) ); 
        t   += gauss_dx; 
    }

    return val * par[3] + par[4]; 
}

using ApexOptics::OpticsTarget_t; 



//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
// if the argument path_dXsv == "", then the linear extrapolation technique will NOT be used. 
int test_angular_accuracy(  const char* path_infile ="data/replay/real_L_V2.root",
                            std::string target_name ="V2",
                            const char* path_dXsv   ="data/csv/poly_dXsv_L_2ord.dat",
                            const char* path_dbfile ="data/csv/poly_WireAndFoil_fp_sv_L_4ord.dat",  
                            const char* tree_name   ="tracks_fp" ) 
{
    const char* const here = "test_angular_accuracy"; 

    //try to get the target we need. 
    OpticsTarget_t target; 
    try { 
        target = ApexOptics::GetTarget(string(target_name)); 
    
    } catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }

    //get the arm we need. 
    TFile *infile; 
    TParameter<bool> *param_is_RHRS; 

    try {
        infile = new TFile(path_infile, "READ");

        //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
        param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
        if (!param_is_RHRS) {
            Error(here, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
            return 1; 
        }
    
    } catch (const std::exception& e) {

        Error(here, "Execption caught trying to open file.\n what(): %s", e.what()); 
        return -1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 
       
    infile->Close(); 
    delete infile; 

    //are we testing the vertical or horizontal angle? 
    // if x_hcs is NOT nan, that means that this is a vertical wire run. 
    // otherwise, its a horizontal-wire run. 
    const bool test_horizontal_angle = (target.x_hcs == target.x_hcs); 


    const vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    }; 

    const vector<string> branches_fp{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const vector<string> branches_dXsv{
        "d_x_sv",
        "d_y_sv",
        "d_dxdz_sv",
        "d_dydz_sv"
    }; 

    const bool use_linear_extrapolation = (string(path_dXsv) != ""); 

    NPolyArray parr_fp_sv, parr_dXsv;  

    try { 
        parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file(path_dbfile, branches_sv, 4); 
        if (use_linear_extrapolation)
            parr_dXsv  = ApexOptics::Parse_NPolyArray_from_file(path_dXsv, branches_dXsv, 5); 
    
    } catch (const std::exception& e) {

        Error(here, "Error caught trying to parse NPolyArray from file.\n what(): %s", e.what()); 
        return -1; 
    }

    //choose which row/column to start at, and how many rows/columns to do in either direction. 

    const int center_row = 8; 
    const int n_side_rows = 4; 

    const int center_col = 5; 
    const int n_side_cols = 4; 



    cout << "parsing done." << endl; 
    //now, we're ready to deal with the data. 

    
    //now, initialize the semicircle-gaussian convolution function



    const double y_rast_rms = 0.113e-3; 

    const double d_dpp = 1e-3; 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //get the average react-vertex. 
    TVector3 react_vtx_hcs(
        0.,
        0.,
        target.z_hcs
    ); 

    //figure out which coordinates are well-defined for this target (they are 'nan' if not). 
    if (target.x_hcs == target.x_hcs) {
        react_vtx_hcs[0] = target.x_hcs; 
    } else {
        react_vtx_hcs[0] = *df.Define("x", [](TVector3 vtx){ return vtx.x(); }, {"position_vtx"}).Mean("x"); 
    }

    if (target.y_hcs == target.y_hcs) {
        react_vtx_hcs[1] = target.y_hcs; 
    } else {
        react_vtx_hcs[1] = *df.Define("y", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}).Mean("y"); 
    }
    
    //get the react vertex in Sieve Coordinates
    const TVector3 react_vtx_scs = ApexOptics::HCS_to_SCS(is_RHRS, react_vtx_hcs); 

    printf("React vertex (SCS): {%+3.1f, %+3.1f, %+3.1f}\n", 
        1e3 * react_vtx_hcs.x(), 
        1e3 * react_vtx_hcs.y(), 
        1e3 * react_vtx_hcs.z()
    ); 

    //this is a little helper class which is meant to avoid some of the awkward syntax typically associated with RDataFrames creation. 
    auto rna = RDFNodeAccumulator(df); 

    //this branch would only be defined if this is monte-carlo data. 
    const bool is_MonteCarlo = false; //rna.IsBranchDefined("x_sv"); 

    cout << (is_MonteCarlo ? "Monte-carlo mode" : "Real-data mode") << endl; 

    TRandom3 rand; 

    //probably not the most elegant way to do this, but here we are. 
    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz }; 

        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_first_guess",  [&parr_fp_sv](const Trajectory_t& Xfp)
        {
            return ApexOptics::RVec_to_Trajectory_t( parr_fp_sv.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz}) ); 

        }, {"Xfp"});

    //define a y-react pos. if it's monte-carlo data, then we want to add a realistic RMS. 
    rna.Define("y_vtx", [y_rast_rms, &rand, is_MonteCarlo](TVector3 vtx)
        {
            return vtx.y() + ( is_MonteCarlo ? rand.Gaus() * y_rast_rms : 0. ); 

        }, {"position_vtx"}); 
    
    auto Compute_dp = [d_dpp, is_RHRS](Trajectory_t Xsv, Trajectory_t dXsv, double y_vtx, double z_vtx)
    {
        dXsv.dpp = 1.; 

        Trajectory_t Xsv_dp = Xsv + (dXsv * d_dpp); 

        Trajectory_t Xhcs    = ApexOptics::SCS_to_HCS(is_RHRS, Xsv); 
        Trajectory_t Xhcs_dp = ApexOptics::SCS_to_HCS(is_RHRS, Xsv_dp); 

        // now, project the track onto the target's z-plane.
        double y1 = Xhcs_dp.y  +  Xhcs_dp.dydz * z_vtx; 
        double y0 = Xhcs.y     +  Xhcs.dydz    * z_vtx;

        return ( y_vtx - y0 ) * d_dpp / ( y1 - y0 ); 
    };

    if (use_linear_extrapolation) {
        //if we're using the linear extrapolation model

        rna.Define("dXsv",     [&parr_dXsv](Trajectory_t Xsv)
            {
                return ApexOptics::RVec_to_Trajectory_t( parr_dXsv.Eval({Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp}) ); 
            }, {"Xsv_first_guess"}); 

        rna.Define("dp", [&Compute_dp](Trajectory_t Xsv, Trajectory_t dXsv, double y_vtx, TVector3 vtx)
            {   
                return Compute_dp(Xsv, dXsv, y_vtx, vtx.z()); 

            }, {"Xsv_first_guess", "dXsv", "y_vtx", "position_vtx"}); 

        rna.Define("Xsv_reco", [](Trajectory_t Xsv, Trajectory_t dXsv, double dp)
            {
                return Xsv + ( dXsv * dp ); 

            }, {"Xsv_first_guess", "dXsv", "dp"});
    
    } else {

        rna.Define("Xsv_reco", [](Trajectory_t Xsv)
            {
                return Xsv; 

            }, {"Xsv_first_guess"});
    }

    rna.Define("Xhcs_reco", [is_RHRS](const Trajectory_t& Xsv)
        {
            return ApexOptics::SCS_to_HCS(is_RHRS, Xsv);  

        }, {"Xsv_reco"});

    rna.Overwrite("z_reco_horizontal", [](const Trajectory_t& Xhcs, double y_vtx)
        {   
            return - ( Xhcs.y - y_vtx ) / Xhcs.dydz; 

        }, {"Xhcs_reco", "y_vtx"});

    rna.Overwrite("z_reco_vertical",   [](const Trajectory_t& Xhcs, const TVector3& vtx)
        {
            return - ( Xhcs.x - vtx.x() ) / Xhcs.dxdz;

        }, {"Xhcs_reco", "position_vtx"});

    //check the status of the RDFNodeAccumulator obejct before proceeding
    if (rna.GetStatus() != RDFNodeAccumulator::kGood) {
        Error(here, "RDFNodeAccumulator reached error status when defining branches.\n Message: %s", rna.GetErrorMsg().c_str()); 
        return -1; 
    }

    const double hole_smearing = 0.400e-3; 

    //if it's montecarlo, then act like the 'real' answers are the 'reconstructed' answers. 
    if (is_MonteCarlo) {

        rna.Overwrite("Xsv_reco", [&rand, hole_smearing](double x, double y, double dxdz, double dydz)
            { 
                return Trajectory_t{ 
                    x, 
                    y, 
                    dxdz + hole_smearing * rand.Gaus(), 
                    dydz + hole_smearing * rand.Gaus()}; 
            }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv"});
    }

    rna = add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"reco_x_sv",    &Trajectory_t::x},
        {"reco_y_sv",    &Trajectory_t::y},
        {"reco_dxdz_sv", &Trajectory_t::dxdz},
        {"reco_dydz_sv", &Trajectory_t::dydz}
    });  




    auto fit_result_ret = TestAngleReco(is_RHRS, rna.Get(), target, center_row, n_side_rows, center_col, n_side_cols); 

    AngleFitResult_t fit; 

    if (fit_result_ret) {
        fit = fit_result_ret.value(); 
    } else {
        Error(here, "Something went wrong with the fitresult."); 
        return -1; 
    }
    cout << "fit_result (sigma overall):" << fit.sigma_dydz_overall << endl; 

    return 0; 
}
