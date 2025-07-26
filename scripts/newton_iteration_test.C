#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

struct Track_t {
    double x,y,dxdz,dydz,dpp; 

    //implicit conversion from this type to RVec<double> 
    operator RVec<double>() const { return RVec<double>{x,y,dxdz,dydz,dpp}; }    
};

//#define EVENT_RANGE 200

//given a DB file to look in, and a list of output polynomial names, will return an NPolyArray object with all the relevant polys filled in
NPolyArray Parse_NPolyArray_from_file(const char* path_dbfile, vector<string> output_names) 
{
    const char* const here = "Parse_NPolyArray_from_file";

    //parse all relevant polys from file
    vector<NPoly> poly_vec; 

    for (const string& str : output_names) {

        NPoly poly(DoF_sv); 

        ApexOptics::Parse_NPoly_from_file(path_dbfile, str.data(), &poly); 

        if (poly.Get_nElems()==0) {
            Warning(here, "Poly '%s' did not find have elements in file '%s'.", str.data(), path_dbfile); 
        }

        poly_vec.push_back(poly); 
    }

    return NPolyArray(poly_vec); 
}   


//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int newton_iteration_test(  const char* path_infile="",
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

#ifdef EVENT_RANGE
    if (ROOT::IsImplicitMTEnabled()) {
        ROOT::DisableImplicitMT(); 
    }
    Info(here, "Multi-threadding is disabled. Events to run : %i", EVENT_RANGE); 
#else 
    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size : %i", ROOT::GetThreadPoolSize()); 
#endif 

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
    const int DoF_fp = 4; 

    //parse all sv=>fp polynomials from file
    vector<NPoly> poly_vec; 

    for (const string& str : branches_fp) {

        NPoly poly(DoF_sv); 

        ApexOptics::Parse_NPoly_from_file(path_dbfile, str.data(), &poly); 

        poly_vec.push_back(poly); 
    }

    //now that we have created the polys, we can create the NPolyArray object
    NPolyArray *parr = new NPolyArray(poly_vec); 

    //check that each poly has at least some elements (otherwise, there has been some sort of file-open error)
    for (int i=0; i<parr->Get_DoF_out(); i++) {
        if (parr->Get_poly(i)->Get_nElems() <= 0) {
            Error(here, "NPoly found without elements. Something has gone wrong with the db file..."); 
            return 1; 
        }
    }

    const int n_iterations = 7; 

    auto rv_dot = [](const RVec<double>& u, const RVec<double>& v) {
        double ret(0.); for (const double& x : u * v) ret += x; return ret; 
    };
    
    auto rv_mag = [](const RVec<double>& u) { 
        double ret(0.); for (const double& x : u * u) ret += x; return sqrt(ret); 
    };

    //number of microns to compute vdc-smearing by
    double vdc_smearing_um = 0.; 

    //what we're effectivley doing here is using newton's method for iterating towrads the root of a nonlinear system.
    // the system we're trying to solve is the following least-square problem: 
    //  - the 'chi-square' in this case is the square error between the model's value for Xfp, and the actual value.
    //    this is given by the d_Xfp vector above. 
    //  - Therefore, the 'F' vector is our evaluation of the gradient of this function, which will be zero at the 
    //    minimum error value (if it isn't a local, false minima.)
    //  - to use newton's method, we need to compute the Jacobian of our 'F' funciton. this is what 'J' will be. 
    //
    auto find_next_Xsv = [parr, rv_dot, rv_mag, DoF_sv, DoF_fp](RVec<double>& Xfp, RVec<double>& Xsv) {
        
        //Get the difference between the model's evaluation of Xfp, and the actual value. 
        RVec<double> d_Xfp{ parr->Eval(Xsv) - Xfp }; 

        RMatrix dGi_dxj = parr->Jacobian(Xsv); 
        RVec<RMatrix> dGi_dxj_dxk = parr->HessianTensor(Xsv); 
 
        //Compute the 'F' vector and the 'J' matrix
        RMatrix J(DoF_sv, DoF_sv, 0.); 
        RVec<double> F(DoF_sv, 0.); 
        
        for (int i=0; i<DoF_fp; i++) {

            for (int j=0; j<DoF_sv; j++) {
            
                F.at(j) += d_Xfp.at(i) * dGi_dxj.at(i,j);
                
                for (int k=0; k<DoF_sv; k++) {
                    J.at(j,k) += 
                        (dGi_dxj.at(i,j) * dGi_dxj.at(i,k))   +   (d_Xfp.at(i) * dGi_dxj_dxk[i].at(j,k)); 
                }
            }
        }
        
        double det = J.Determinant(); 
        
        //printf(" det: %+.3e\n", J.Determinant()); 
        if ((det != det) || (fabs(det) < 1e-30)) return RVec<double>{}; 

        return  -1. * J.Solve( F ); 
    };

    auto Iterate_to_Xfp = [find_next_Xsv, n_iterations, rv_mag, parr](const Track_t& Xfp, const Track_t& Xsv) 
    {   
        RVec<double> Xfp_vec{ 
            Xfp.x,
            Xfp.y,
            Xfp.dxdz,
            Xfp.dydz 
        };

        RVec<double> Xsv_vec{
            Xsv.x,
            Xsv.y,
            Xsv.dxdz,
            Xsv.dydz,
            Xsv.dpp
        }; 
        
        for (int i=0; i<n_iterations; i++) {
            
            //printf("it %3i error: % .9f\n", i, rv_mag(Xfp_vec - parr->Eval(Xsv_vec)));

            RVec<double> d_Xsv = find_next_Xsv(Xfp_vec, Xsv_vec);
            
            //cout << "d_Xsv: " << d_Xsv << endl;
            //if the 'J' matrix computed above is singular, the 'find_next_Xsv' fcn return an empty RVec
            if (d_Xsv.size() == 0) break; 
            
            Xsv_vec += d_Xsv; 
        }

        //printf("final error: % .9f\n", rv_mag(Xfp_vec - parr->Eval(Xsv_vec)));
        return Xsv_vec; 
    };


    //Define all of the branches we want to create models to map between
    auto df_output = df

#ifdef EVENT_RANGE
        .Range(EVENT_RANGE)
#endif 
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

        .Define("Xsv_model", [Iterate_to_Xfp](Track_t& Xfp, Track_t& Xsv)
        {
            return Iterate_to_Xfp(Xfp, Xsv); 
//            return Track_t{ .x=ret[0], .y=ret[1], .dxdz=ret[2], .dydz=ret[3], .dpp=ret[4]};         
        }, {"Xfp", "Xsv"})

        .Define("reco_x_sv",        [](const RVec<double>& X){ return X[0]; },    {"Xsv_model"})
        .Define("reco_y_sv",        [](const RVec<double>& X){ return X[1]; },    {"Xsv_model"})
        .Define("reco_dxdz_sv",     [](const RVec<double>& X){ return X[2]; },    {"Xsv_model"})
        .Define("reco_dydz_sv",     [](const RVec<double>& X){ return X[3]; },    {"Xsv_model"})
        .Define("reco_dpp_sv",      [](const RVec<double>& X){ return X[4]; },    {"Xsv_model"})
        
        .Define("Xfp_model",  [parr](const RVec<double> &Xsv)
        {
            return parr->Eval(Xsv); 
        }, {"Xsv_model"})   

        .Define("err_x_fp",     [](RVec<double>& Xfp, double x){ return (Xfp[0]-x)*1e3; }, {"Xfp_model", "x_fp"})
        .Define("err_y_fp",     [](RVec<double>& Xfp, double x){ return (Xfp[1]-x)*1e3; }, {"Xfp_model", "y_fp"})
        .Define("err_dxdz_fp",  [](RVec<double>& Xfp, double x){ return (Xfp[2]-x)*1e3; }, {"Xfp_model", "dxdz_fp"})
        .Define("err_dydz_fp",  [](RVec<double>& Xfp, double x){ return (Xfp[3]-x)*1e3; }, {"Xfp_model", "dydz_fp"})

        .Define("err_x_sv",     [](RVec<double>& Xsv, double x){ return (Xsv[0]-x)*1e3; }, {"Xsv_model", "x_sv"})
        .Define("err_y_sv",     [](RVec<double>& Xsv, double x){ return (Xsv[1]-x)*1e3; }, {"Xsv_model", "y_sv"})
        .Define("err_dxdz_sv",  [](RVec<double>& Xsv, double x){ return (Xsv[2]-x)*1e3; }, {"Xsv_model", "dxdz_sv"})
        .Define("err_dydz_sv",  [](RVec<double>& Xsv, double x){ return (Xsv[3]-x)*1e3; }, {"Xsv_model", "dydz_sv"})
        .Define("err_dpp_sv",   [](RVec<double>& Xsv, double x){ return (Xsv[4]-x)*1e3; }, {"Xsv_model", "dpp_sv"});
    
    //this histogram will be of the actual sieve-coords
    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords. db:'%s', data:'%s'", path_dbfile, path_infile); 
    
    char buff_hxy_title[200];  
    sprintf(buff_hxy_title, "Reconstructed sieve coordinates. VDC smearing: %.1f um;x_sv;y_sv", vdc_smearing_um); 
    
    gStyle->SetPalette(kSunset);
    gStyle->SetOptStat(0); 
    auto c2         = new TCanvas("c2", b_c_title); 
    auto h_xy_sieve = df_output.Histo2D({"h_xy_sieve", buff_hxy_title, 250, -45e-3, 45e-3, 250, -45e-3, 45e-3}, "reco_x_sv", "reco_y_sv"); 
    h_xy_sieve->DrawCopy("col2"); 
    return 0; 


    auto c = new TCanvas("c1", b_c_title, 1200, 800); 

    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); 
    auto h_x    = df_output.Histo1D({"h_x",    "Error of x_fp;mm", 200, -5, 5},    "err_x_sv"); 
    h_x->DrawCopy(); 
    c->cd(2); 
    auto h_y    = df_output.Histo1D({"h_y",    "Error of y_fp;mm", 200, -5, 5},    "err_y_sv"); 
    h_y->DrawCopy(); 
    c->cd(3); 
    auto h_dxdz = df_output.Histo1D({"h_dxdz", "Error of dxdz_fp;mrad", 200, -5, 5}, "err_dxdz_sv"); 
    h_dxdz->DrawCopy(); 
    c->cd(4); 
    auto h_dydz = df_output.Histo1D({"h_dydz", "Error of dydz_fp;mrad", 200, -5, 5}, "err_dydz_sv"); 
    h_dydz->DrawCopy(); 
    
    
    

    delete parr;    

    return 0;
}