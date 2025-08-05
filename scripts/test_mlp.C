#include "TROOT.h"


using namespace std; 
using namespace ROOT::VecOps; 

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int test_mlp(  const char* path_infile="",
                const char* path_dbfile="data/csv/mlp_prod_sv_fp_L.dat",  
                const char* tree_name="tracks_fp" )
{
    const char* const here = "train_mlp"; 
    
    vector<string> branches_input   = {
        "x_sv", 
        "y_sv", 
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };
    const size_t n_inputs = branches_input.size(); 

    vector<string> branches_output  = {
        "x_q1",
        "y_q1",
        "dxdz_q1",
        "dydz_q1"
    }; 
    const size_t n_outputs = branches_output.size(); 

    MultiLayerPerceptron *mlp = ApexOptics::Parse_mlp_from_file(path_dbfile); 
    
    if (!mlp) {
        Error(here, "Unable to parse MLP from file '%s'", path_dbfile); 
        return -1; 
    }

    const double hrs_momentum = 1104.; 

    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df(tree_name, path_infile); 

    auto df_input = df
        .Define("x_sv",     [](TVector3 v){ return v.x(); }, {"position_sieve"})
        .Define("y_sv",     [](TVector3 v){ return v.y(); }, {"position_sieve"})
        .Define("dxdz_sv",  [](TVector3 v){ return v.x()/v.z(); }, {"momentum_sieve"})
        .Define("dydz_sv",  [](TVector3 v){ return v.y()/v.z(); }, {"momentum_sieve"})
        .Define("dpp_sv",   [hrs_momentum](TVector3 v){ return (v.Mag()-hrs_momentum)/hrs_momentum; }, {"momentum_sieve"}) 

        .Define("x_q1",     [](TVector3 v){ return v.x(); }, {"position_Q1"})
        .Define("y_q1",     [](TVector3 v){ return v.y(); }, {"position_Q1"})
        .Define("dxdz_q1",  [](TVector3 v){ return v.x()/v.z(); }, {"momentum_Q1"})
        .Define("dydz_q1",  [](TVector3 v){ return v.y()/v.z(); }, {"momentum_Q1"})
        .Define("dpp_q1",   [hrs_momentum](TVector3 v){ return (v.Mag()-hrs_momentum)/hrs_momentum; }, {"momentum_Q1"})

    
        .Define("X_inp", [n_inputs ](){ RVec<double> v{}; v.reserve(n_inputs);  return v; }, {})
        .Define("Z_out", [n_outputs](){ RVec<double> v{}; v.reserve(n_outputs); return v; }, {}); 

    vector<ROOT::RDF::RNode> nodes{ df_input }; 

    //for each event, create a RVec<double> for the inputs and outputs.
    //We do this by adding each input to the vector, one-at-a-time, in order. 
    for (const auto& str : branches_input) {
        nodes.push_back( 
            nodes.back().Redefine("X_inp", [](RVec<double>& v, double x){ v.push_back(x); return v; }, {"X_inp", str.data()})
        ); 
    }
    //now, do the same for the outputs. 
    for (const auto& str : branches_output) {
        nodes.push_back( 
            nodes.back().Redefine("Z_out", [](RVec<double>& v, double x){ v.push_back(x); return v; }, {"Z_out", str.data()})
        ); 
    }
    
    nodes.push_back( 
        nodes.back()
        .Define("Z_model",  [mlp](const RVec<double>& X)
        {
            return mlp->Eval(X); 
        }, {"X_inp"}) 
        
        .Define("Z_error", [](const RVec<double>& Z_out, const RVec<double>& Z_model)
        {
            return Z_out - Z_model; 
        }, {"Z_out", "Z_model"})
    ); 

    
    //now, acually measure the error of each output. 
    int i_output=0;    
    for (const auto& str : branches_output) {
        
        char buff[100]; 
        sprintf(buff, "err_%s", str.data()); 

        nodes.push_back(
            nodes.back().Define(buff, [i_output](const RVec<double>& Z_err)
            {
                return Z_err.at(i_output); 
            }, {"Z_error"})
        ); 
        ++i_output; 
    }

    auto df_output = nodes.back(); 


    auto h_x_q1     = df_output.Histo1D<double>({"h_x",     ";x_{Q1};",     200, -1, -1}, "err_x_q1"); 
    auto h_y_q1     = df_output.Histo1D<double>({"h_y",     ";y_{Q1};",     200, -1, -1}, "err_y_q1"); 
    auto h_dxdz_q1  = df_output.Histo1D<double>({"h_dxdz",  ";dx/dz_{Q1};", 200, -1, -1}, "err_dxdz_q1"); 
    auto h_dydz_q1  = df_output.Histo1D<double>({"h_dydz",  ";dy/dz_{Q1};", 200, -1, -1}, "err_dydz_q1"); 

    auto c = new TCanvas("c", "Error of mlp", 1200, 800); 

    c->Divide(2,2, 0.01,0.01); 

    gStyle->SetOptStat(0);
    gStyle->SetPalette(0); 

    c->cd(1); h_x_q1   ->DrawCopy("col2"); 
    c->cd(2); h_y_q1   ->DrawCopy("col2"); 
    c->cd(3); h_dxdz_q1->DrawCopy("col2"); 
    c->cd(4); h_dydz_q1->DrawCopy("col2"); 

    return 0; 
}