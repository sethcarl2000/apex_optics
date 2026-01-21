#include "include/NPolyArray_fit.h" 
#include "include/Get_TParameter_from_TFile.h"
#include <sstream> 
#include <cstdio>

//sieve coords  {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}
//q1 coords     {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}
//q1-fwd coords {"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"} 
//fp coords     {"x_fp","y_fp","dxdz_fp","dydz_fp"}

using namespace std;

/***
 * @brief creates monte-carlo based polynomial optics models
 * 
 */
int create_MC_polynomial_models(const int poly_order     =5,
                                const char* path_infile  ="data/mc/mc_L_production.root",
                                const char* stem_outfile ="data/poly/mc") 
{
    vector<string> inputs, outputs; 
    char path_outfile[1000]; 
    
    const bool is_RHRS = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS").value_or(false); 

    //   sv => q1 
    inputs  = vector<string>{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    outputs = vector<string>{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}; 
    sprintf(path_outfile, "%s_sv_q1_%c_%iord.dat", stem_outfile, is_RHRS?'R':'L', poly_order); 
    if (NPolyArray_fit(poly_order, path_infile, path_outfile, inputs, outputs) != 0) {
        Error(__func__, "Something went wrong trying to create 'sv => q1' polynomial."); 
        return -1; 
    }
    //         q1 => fp 
    inputs  = vector<string>{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    outputs = vector<string>{"x_fp","y_fp","dxdz_fp","dydz_fp"}; 
    sprintf(path_outfile, "%s_q1_fp_%c_%iord.dat", stem_outfile, is_RHRS?'R':'L', poly_order); 
    if (NPolyArray_fit(poly_order, path_infile, path_outfile, inputs, outputs) != 0) {
        Error(__func__, "Something went wrong trying to create 'sv => q1' polynomial."); 
        return -1; 
    }
    //   sv =======> fp 
    inputs  = vector<string>{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    outputs = vector<string>{"x_fp","y_fp","dxdz_fp","dydz_fp"}; 
    sprintf(path_outfile, "%s_sv_fp_%c_%iord.dat", stem_outfile, is_RHRS?'R':'L', poly_order); 
    if (NPolyArray_fit(poly_order, path_infile, path_outfile, inputs, outputs) != 0) {
        Error(__func__, "Something went wrong trying to create 'sv => q1' polynomial."); 
        return -1; 
    }

    //   sv <======= fp 
    inputs  = vector<string>{"x_fp","y_fp","dxdz_fp","dydz_fp"};
    outputs = vector<string>{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}; 
    sprintf(path_outfile, "%s_fp_sv_%c_%iord.dat", stem_outfile, is_RHRS?'R':'L', poly_order); 
    if (NPolyArray_fit(poly_order, path_infile, path_outfile, inputs, outputs) != 0) {
        Error(__func__, "Something went wrong trying to create 'sv => q1' polynomial."); 
        return -1; 
    }
    //         q1 <= fp 
    inputs  = vector<string>{"x_fp","y_fp","dxdz_fp","dydz_fp"};
    outputs = vector<string>{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}; 
    sprintf(path_outfile, "%s_fp_q1_%c_%iord.dat", stem_outfile, is_RHRS?'R':'L', poly_order); 
    if (NPolyArray_fit(poly_order, path_infile, path_outfile, inputs, outputs) != 0) {
        Error(__func__, "Something went wrong trying to create 'sv => q1' polynomial."); 
        return -1; 
    }
    //   sv <= q1 
    inputs  = vector<string>{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    outputs = vector<string>{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}; 
    sprintf(path_outfile, "%s_q1_sv_%c_%iord.dat", stem_outfile, is_RHRS?'R':'L', poly_order); 
    if (NPolyArray_fit(poly_order, path_infile, path_outfile, inputs, outputs) != 0) {
        Error(__func__, "Something went wrong trying to create 'sv => q1' polynomial."); 
        return -1; 
    }

    return 0; 
}