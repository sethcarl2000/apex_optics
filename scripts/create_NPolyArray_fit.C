#include "include/NPolyArray_fit.h"

//sieve coords  {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}
//q1 coords     {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}
//q1-fwd coords {"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"} 
//fp coords     {"x_fp","y_fp","dxdz_fp","dydz_fp"}

using namespace std; 

int create_NPolyArray_fit(  const int poly_order     =5,
                            const char* path_infile  ="data/mc/mc_R_production.root",
                            const char* stem_outfile ="data/poly/mc_q1_fp_R_5ord.dat",  
                            const vector<string> inputs  ={"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"},
                            const vector<string> outputs ={"x_fp","y_fp","dxdz_fp","dydz_fp"},
                            const double pruning_constant=-1., 
                            const char* tree_name    ="tracks_fp") 
{
    const int poly_order     =5;
    const char* path_infile  ="data/mc/mc_R_production.root";
    const char* stem_outfile ="data/poly/mc_q1_fp_R_5ord.dat";  
    const vector<string> inputs  {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const vector<string> outputs {"x_fp","y_fp","dxdz_fp","dydz_fp"};
    const double pruning_constant   =-1.;
    const char* tree_name    ="tracks_fp";

    return NPolyArray_fit(
        poly_order,
        path_infile, 
        stem_outfile, 
        inputs, outputs, 
        pruning_constant, 
        tree_name
    ); 
}