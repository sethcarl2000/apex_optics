#include "include/NPolyArray_fit.h" 
#include "include/Get_TParameter_from_TFile.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include "include/RDFNodeAccumulator.h"
#include <NPolyArray.h>
#include <ApexOptics.h> 
#include <TDatime.h> 
#include <fstream> 
#include <sstream> 
#include <vector> 
#include <string> 
#include <cstdio>
#include <map> 

//sieve coords  {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}
//q1 coords     {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}
//q1-fwd coords {"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"} 
//fp coords     {"x_fp","y_fp","dxdz_fp","dydz_fp"}

using namespace std;

/***
 * @brief creates monte-carlo based polynomial optics models
 * 
 */
int create_fp_correction_poly(  const int poly_order     =5,
                                const char* path_infile  ="data/sieve_holes/fits_21Dec/V123.root",
                                const char* path_outfile ="data/poly/mc_fp_fp-fwd_L_5ord.dat", 
                                const char* path_sv_fp_poly="data/poly/mc_sv_fp_L_5ord.dat") 
{
    const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS").value_or(false); 

    RDFNodeAccumulator rna("tracks_fp", path_infile); 

    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"}; 
    const vector<string> branches_fp_fwd{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"}; 

    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    
    NPolyArray parr_sv_fp = ApexOptics::Parse_NPolyArray_from_file(path_sv_fp_poly, branches_fp, 5);     

    //this minor shift in dp/p makes the monte-carlo and real data line up a bit better (for the LHRS)
    const double dpp_correction = -0.0005; 


    rna.Define("Xfp_reco", [&parr_sv_fp, dpp_correction](double x, double y, double dxdz, double dydz, double dpp)
        {
            return ApexOptics::RVec_to_Trajectory_t(
                parr_sv_fp.Eval({x,y,dxdz,dydz,dpp + dpp_correction})            
            ); 
        }, branches_sv);

    rna = Add_branches_from_Trajectory_t(rna.Get(), "Xfp_reco", 
        {"fwd_x_fp", "fwd_y_fp", "fwd_dxdz_fp", "fwd_dydz_fp"}
    ); 
    
    // first, we're going to write some of our own data to the file.
    // the 'trunc' flag means that we're going to wipe the contents of the file and start over.  
    fstream outfile(path_outfile, ios::out | ios::trunc); 

    if (!outfile.is_open()) {
        Error(__func__, "Unable to open file '%s'.", path_outfile); 
        return -1; 
    }

    //now that we can open the file, let's add some metadata to it. 
    ostringstream metadata; 
    
    // Muon's comment (21 Jan 26): )_------------------------------- 

    //add basic metadata 
    metadata <<   
        "# -- Polynomial metadata: \n"
        "# -- Order: " << poly_order << "\n"
        "# -- Trained on data: " << path_infile << "\n"
        "# -- Monte-carlo sv=>fp polynomial: " << path_sv_fp_poly << "\n"  
        "# -- dp/p_sv correction from training data: " << dpp_correction << "\n"; 
        
    //add the branches; 
    metadata << "# -- Input branches: { "; 
    for (const auto& br : branches_fp) metadata << br << " "; 
    metadata << "}\n"; 

    //add a timestamp
    TDatime dtime; 
    metadata << 
        "# -- Time created: " << dtime.AsString() << "\n"; 

    cout << metadata.str() << endl; 

    outfile << metadata.str() << endl; 
    outfile.close(); 

    //now, actually perform the fit
    map<string,NPoly*> polymap = ApexOptics::Create_NPoly_fit(rna.Get(), poly_order,
        branches_fp, branches_fp_fwd 
    );

    //save the new polynomial 
    ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap, true); 

    return 0; 
}