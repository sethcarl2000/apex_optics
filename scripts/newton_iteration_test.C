#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

struct Track_t {
    double x,y,dxdz,dydz,dpp; 
};

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int newton_iteration_test(  bool is_RHRS=false,
                            const char* path_infile="",
                            const char* path_dbfile="data/csv/db_prod_mc_sv_fp_L_3ord.dat",  
                            const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fitpoints_mc_sv_fp"; 

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

    delete tree; 
    infile->Close(); 
    delete infile; 

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 


    const double hrs_momentum = 1104.0; 

    vector<string> branches_fp = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const int DoF_sv = 5; 

    //parse all sv=>fp polynomials from file
    vector<NPoly> poly_vec; 

    for (const string& str : branches_fp) {

        NPoly poly(DoF_sv); 

        ApexOptics::Parse_NPoly_from_file(path_dbfile, str.data(), &poly); 

        poly_vec.push_back(poly); 
    }

    //now that we have created the polys, we can create the NPolyModel object
    NPolyModel *pmod = new NPolyModel(poly_vec); 


    auto Iterate_to_Xfp = [pmod](const Track_t& Xfp, const Track_t& Xsv) 
    {
        return Track_t{}; //noop
    };


    //Define all of the branches we want to create models to map between
    auto df_output = df

        //this is the only difference between VDC TRANSPORT COORDINATES (tra) and FOCAL PLANE COORDINATES (fp)
        .Redefine("dxdz_fp", [](double x_tra, double dxdz_tra){return dxdz_tra - x_tra/6.;}, {"x_fp", "dxdz_fp"})

        .Define("x_sv",      [](TVector3 v){ return v.x(); },        {"position_sieve"})
        .Define("y_sv",      [](TVector3 v){ return v.y(); },        {"position_sieve"})
        .Define("dxdz_sv",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_sieve"})
        .Define("dydz_sv",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_sieve"})
        .Define("dpp_sv",    [hrs_momentum](TVector3 v){ return (v.z()-hrs_momentum)/hrs_momentum; }, {"momentum_sieve"})

        //define the Track_t structs that we will need to use for the newton iteration
        .Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Track_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz, .dpp=0. };  
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})

        .Define("Xsv", [](double x, double y, double dxdz, double dydz, double dpp)
        {
            return Track_t{ .x=x, .y=y, .dxdz=dxdz, .dydz=dydz, .dpp=dpp };  
        }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"})

        .Define("X_fp_model",  [pmod](double x, double y, double dxdz, double dydz, double dpp)
        {
            return pmod->Eval({x,y,dxdz,dydz,dpp});         
        }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"})   

        .Define("err_x_fp",     [](RVec<double>& Xfp, double x){ return Xfp[0]-x; }, {"X_fp_model", "x_fp"})
        .Define("err_y_fp",     [](RVec<double>& Xfp, double x){ return Xfp[1]-x; }, {"X_fp_model", "y_fp"})
        .Define("err_dxdz_fp",  [](RVec<double>& Xfp, double x){ return Xfp[2]-x; }, {"X_fp_model", "dxdz_fp"})
        .Define("err_dydz_fp",  [](RVec<double>& Xfp, double x){ return Xfp[3]-x; }, {"X_fp_model", "dydz_fp"});

    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords: %s", path_dbfile); 
    auto c = new TCanvas("c", b_c_title, 1200, 800); 

    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); 
    auto h_x    = df_output.Histo1D({"h_x",    "Error of x_fp;mm", 200, -10e-3, 10e-3},    "err_x_fp"); 
    h_x->DrawCopy(); 
    c->cd(2); 
    auto h_y    = df_output.Histo1D({"h_y",    "Error of y_fp;mm", 200, -10e-3, 10e-3},    "err_y_fp"); 
    h_y->DrawCopy(); 
    c->cd(3); 
    auto h_dxdz = df_output.Histo1D({"h_dxdz", "Error of dxdz_fp;mrad", 200, -2e-3, 2e-3}, "err_dxdz_fp"); 
    h_dxdz->DrawCopy(); 
    c->cd(4); 
    auto h_dydz = df_output.Histo1D({"h_dydz", "Error of dydz_fp;mrad", 200, -2e-3, 2e-3}, "err_dydz_fp"); 
    h_dydz->DrawCopy(); 


    delete pmod;    

    return 0;
}