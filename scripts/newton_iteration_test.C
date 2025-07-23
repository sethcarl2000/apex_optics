#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

struct Track_t {
    double x,y,dxdz,dydz,dpp; 

    //implicit conversion from this type to RVec<double> 
    operator RVec<double>() const { return RVec<double>{x,y,dxdz,dydz,dpp}; }    
};

#define EVENT_RANGE 1

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

    //now that we have created the polys, we can create the NPolyModel object
    NPolyModel *pmod = new NPolyModel(poly_vec); 

    const int n_iterations = 10; 

    auto rv_dot = [](const RVec<double>& u, const RVec<double>& v) {
        double ret(0.); for (const double& x : u * v) ret += x; return ret; 
    };
    
    auto rv_mag = [](const RVec<double>& u) { 
        double ret(0.); for (const double& x : u * u) ret += x; return sqrt(ret); 
    };

    auto find_next_Xsv = [pmod, rv_dot, rv_mag, DoF_sv, DoF_fp](RVec<double>& Xfp, RVec<double>& Xsv) {
        
        //Get the difference between the model's evaluation of Xfp, and the actual value. 
        RVec<double> d_Xfp{ pmod->Eval(Xsv) - Xfp }; 

        const NPoly* polys[DoF_fp]; 

        //there will be 1 Hessian matrix per polynomial
        RMatrix hess_i[DoF_fp];

        //there will also be 1 gradient per polynomial
        RVec<double> grad_i[DoF_fp];

        for (int i=0; i<DoF_fp; i++) {
            
            //get a polynomial ptr, and const-cast it to make sure we can't modify it (thread safety!)
            polys[i] = const_cast<const NPoly*>(pmod->Get_poly(i)); 
            
            //now, compute the hessian matrix & gradients
            hess_i[i] = polys[i]->Hessian(Xsv); 
            grad_i[i] = polys[i]->Gradient(Xsv); 
        }

        //what we're effectivley doing here is using newton's method for iterating towrads the root of a nonlinear system.
        // the system we're trying to solve is the following least-square problem: 
        //  - the 'chi-square' in this case is the square error between the model's value for Xfp, and the actual value.
        //    this is given by the d_Xfp vector above. 
        //  - Therefore, the 'F' vector is our evaluation of the gradient of this function, which will be zero at the 
        //    minimum error value (if it isn't a local, false minima.)
        //  - to use newton's method, we need to compute the Jacobian of our 'F' funciton. this is what 'J' will be. 

        RMatrix J(DoF_sv, DoF_sv, 0.); 
        RVec<double> F(DoF_sv, 0.); 
        
        //fill the 'F' vector and the 'J' matrix
        for (int i=0; i<DoF_fp; i++) {

            for (int j=0; j<DoF_sv; j++) {
            
                F.at(j) += d_Xfp.at(i) * grad_i[i].at(j);  
                
                for (int k=0; k<DoF_sv; k++) {
                    J.at(j,k) += (grad_i[i].at(j) * grad_i[i].at(k)) + (d_Xfp.at(i) * hess_i[i].at(j,k)); 
                }
            }
        }
        
        double det = J.Determinant(); 
        
        printf(" det: %+.3e\n", J.Determinant()); 
        if ((det != det) || (fabs(det) < 1e-30)) return RVec<double>{}; 

        return  -1.* J.Solve( F ); 
    };

    auto Iterate_to_Xfp = [find_next_Xsv, n_iterations, rv_mag, pmod](const Track_t& Xfp, const Track_t& Xsv) 
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

        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        for (int i=0; i<n_iterations; i++) {
            
            printf("it %3i error: % .9f\n", i, rv_mag(Xfp_vec - pmod->Eval(Xsv_vec)));

            auto d_Xsv = find_next_Xsv(Xfp_vec, Xsv_vec);
            cout << "d_Xsv: " << d_Xsv << endl;
            if (d_Xsv.size() == 0) break; 
            
            Xsv_vec += d_Xsv; 
        }

        printf("final error: % .9f\n", rv_mag(Xfp_vec - pmod->Eval(Xsv_vec)));
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
            auto ret = Iterate_to_Xfp(Xfp, Xsv); 
            return Track_t{ .x=ret[0], .y=ret[1], .dxdz=ret[2], .dydz=ret[3], .dpp=ret[4]};         
        }, {"Xfp", "Xsv"})

        .Define("Xfp_model",  [pmod](Track_t &Xsv)
        {
            return pmod->Eval({
                Xsv.x,
                Xsv.y,
                Xsv.dxdz,
                Xsv.dydz,
                Xsv.dpp
            });         
        }, {"Xsv_model"})   

        .Define("err_x_fp",     [](RVec<double>& Xfp, double x){ return Xfp[0]-x; }, {"Xfp_model", "x_fp"})
        .Define("err_y_fp",     [](RVec<double>& Xfp, double x){ return Xfp[1]-x; }, {"Xfp_model", "y_fp"})
        .Define("err_dxdz_fp",  [](RVec<double>& Xfp, double x){ return Xfp[2]-x; }, {"Xfp_model", "dxdz_fp"})
        .Define("err_dydz_fp",  [](RVec<double>& Xfp, double x){ return Xfp[3]-x; }, {"Xfp_model", "dydz_fp"});

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