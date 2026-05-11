#include <memory>
#include <string> 
#include <map>
#include <fstream>
#include <iomanip> 
#include <TFile.h>
#include <TVector3.h>
#include <TParameter.h>
#include <TSystem.h> 
#include <TDatime.h> 
#include <memory> 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx> 
#include "include/RDFNodeAccumulator.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include <ApexOptics.h> 
#include <sstream> 

using namespace std;
using namespace ROOT::VecOps;
using ApexOptics::Trajectory_t; 

//sieve coords  {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}
//q1 coords     {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}
//q1-fwd coords {"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"} 
//fp coords     {"x_fp","y_fp","dxdz_fp","dydz_fp"}

//  this creates a 'raster' polynomial, which takes the 'first pass' sieve-coordinates and outputs a 
//  'corrected' set of sieve coordinates, which uses information from the vertical raster to 'correct' its first guess. 
int make_raster_polynomial( const int poly_order =2,
                            const char* path_infile  ="data/sieve_holes/fits_21Dec/V123.root",
                            const char* stem_outfile ="data/poly/fits_21Dec/V123_sv-yhcs_sv_L_2ord.dat",  
                            const char* tree_name    ="tracks_fp") 
{
    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(__func__, "libApexOptics could not be loaded."); 
        return 1; 
    }

    auto infile = new TFile(path_infile, "READ");

    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(__func__, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(__func__, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 

    infile->Close(); 
    delete infile; 

    ROOT::EnableImplicitMT(); 
    Info(__func__, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame. 
    //first, let's produce the data. 

    //this is the polynomial we will use to make the 'first pass' sieve coordiantes 
    const char* path_first_pass_poly = "data/poly/fits_21Dec/V123_fp_sv_L_4ord.dat";
    
    //parse the first-pass polynomial array 
    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}; 
    const vector<string> branches_fg_sv{"fg_x_sv","fg_y_sv","fg_dxdz_sv","fg_dydz_sv","fg_dpp_sv", "y_hcs"}; 
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"}; 
    

    NPolyArray parr_forward = ApexOptics::Parse_NPolyArray_from_file(path_first_pass_poly, branches_sv, 4);                                 

    auto df_ptr = unique_ptr<ROOT::RDataFrame>(nullptr); 
    try { df_ptr = unique_ptr<ROOT::RDataFrame>(new ROOT::RDataFrame(tree_name, path_infile)); }
    catch (const std::invalid_argument& e) {
        Error(__func__, "Trying to create RDataFrame threw std::invalid_argument exception.\n what(): %s", e.what()); 
        return -1; 
    } 
    ROOT::RDataFrame& df = *df_ptr; 

    auto rna = RDFNodeAccumulator(df); 

    rna.Define("Xsv_first_pass", [&parr_forward](double x, double y, double dxdz, double dydz)
        {
            RVec<double> Xsv = parr_forward.Eval({x,y,dxdz,dydz}); 
            return Trajectory_t{Xsv[0],Xsv[1],Xsv[2],Xsv[3],Xsv[4]}; 
        }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}); 

    rna.Define("y_hcs", [](TVector3 vtx_hcs){ return vtx_hcs.y(); }, {"position_vtx"}); 

    rna = Add_branch_from_Trajectory_t(rna.Get(), "Xsv_first_pass", {
        {"fg_x_sv", &Trajectory_t::x},
        {"fg_y_sv", &Trajectory_t::y},
        {"fg_dxdz_sv", &Trajectory_t::dxdz},
        {"fg_dydz_sv", &Trajectory_t::dydz},
        {"fg_dpp_sv", &Trajectory_t::dpp}
    }); 

    cout << "Creating polynomials for inputs => outputs..." << flush; 
    //now, let's create the first-pass reconstruction of the sieve-coordinates

    // we're going to map back onto Xsv, using the raster info 
    map<string,NPoly*> polymap = ApexOptics::Create_NPoly_fit(rna.Get(), poly_order, branches_fg_sv, branches_sv);
    
    cout << "done." << endl;

    string path_outfile(stem_outfile);

    //create the output file ___________________________________________
    if (path_outfile != "") {
        cout << "Creating dbfiles for polynomials..." << flush; 
        
        // first, we're going to write some of our own data to the file.
        // the 'trunc' flag means that we're going to wipe the contents of the file and start over.  
        fstream outfile(path_outfile, ios::out | ios::trunc); 

        if (!outfile.is_open()) {
            Error(__func__, "Unable to open file '%s'.", path_outfile.c_str()); 
            return -1; 
        }

        //now that we can open the file, let's add some metadata to it. 
        ostringstream metadata_format; 
        metadata_format << 
            "# Polynomial metadata: \n"
            "# -- Order: "<< poly_order <<" \n"
            "# -- Trained on data: '"<< path_infile <<"'\n"
            "# -- first-guess (fg) polynomial: '"<< path_first_pass_poly <<"'\n"
            "# -- Input branches: { "; 

        for (const auto& br : branches_fg_sv) metadata_format << (br + " "); 
        
        TDatime dtime; 
        metadata_format << 
            "}\n"
            "# -- Time created: "<< dtime.AsString() <<"\n"; 

        //const int md_len = metadata_format.size(); 
        //char metadata[md_len];
        //sprintf(metadata, metadata_format.c_str(), poly_order, path_infile, path_first_pass_poly, dtime.AsString()); 
        cout << metadata_format.str() << endl; 

        outfile << metadata_format.str() << endl; 
        outfile.close(); 

        ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap, true); 

        cout << "done." << endl; 

    } else { 

        cout << "Skipping db-file creation." << endl; 
    }
    //____________________________________________________________________________

    //delete our poly models
    //for (auto it = polymap.begin(); it != polymap.end(); it++ ) delete it->second;
    
    return 0;
}