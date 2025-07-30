#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

struct Track_t {
    double x,y,dxdz,dydz,dpp; 

    //implicit conversion from this type to RVec<double> 
    operator RVec<double>() const { return RVec<double>{x,y,dxdz,dydz,dpp}; }  
    
};

//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Track_t struct. 
ROOT::RDF::RNode add_branch_from_Track_t(ROOT::RDF::RNode df, const char* branch_in, map<string, double Track_t::*> branches)
{
    const int n_nodes = branches.size() + 1; 
    RVec<ROOT::RDF::RNode> df_nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Track_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Track_t::coord' 
        auto new_node = df_nodes.back()

            .Define(branch_name, [coord](const Track_t& track) { return track.*coord; }, {branch_in}); 

        df_nodes.push_back(new_node); 
    }
    
    return df_nodes.back(); 
}

//#define EVENT_RANGE 2000

//given a DB file to look in, and a list of output polynomial names, will return an NPolyArray object with all the relevant polys filled in
NPolyArray Parse_NPolyArray_from_file(const char* path_dbfile, vector<string> output_names, const int DoF) 
{
    const char* const here = "Parse_NPolyArray_from_file";

    //parse all relevant polys from file
    vector<NPoly> poly_vec; 

    for (const string& str : output_names) {

        NPoly poly(DoF); 

        ApexOptics::Parse_NPoly_from_file(path_dbfile, str.data(), &poly); 

        if (poly.Get_nElems()==0) {
            Warning(here, "Poly '%s' did not find have elements in file '%s'.", str.data(), path_dbfile); 
        }

        poly_vec.push_back(poly); 
    }

    return NPolyArray(poly_vec); 
}   


//transform coordinates from Sieve Coordinate System (SCS) to Hall coordinate system (HCS)
void SCS_to_HCS(Track_t& track, const bool is_RHRS) 
{
    //direction (SCS)
    auto dir = TVector3( track.dxdz, track.dydz, 1. );

    auto pos = TVector3( track.x, track.y, 0. ) + ApexOptics::Get_sieve_pos(is_RHRS); 

    //rotate both the position and the direction
    dir.RotateZ( -TMath::Pi()/2. ); 
    dir.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

    pos.RotateZ( -TMath::Pi()/2. ); 
    pos.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

    //compute the new slopes
    track.dxdz = dir.x() / dir.z(); 
    track.dydz = dir.y() / dir.z(); 

    //use these new slopes to project the track onto the z=0 plane in HCS 
    track.x = pos.x() - track.dxdz * pos.z(); 
    track.y = pos.y() - track.dydz * pos.z(); 
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

    vector<string> branches_q1 = {
        "x_q1",
        "y_q1",
        "dxdz_q1",
        "dydz_q1",
        "dpp_q1"
    };

    vector<string> branches_sv = {
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };

    const int DoF_sv = 5; 
    const int DoF_fp = 4; 

    //now that we have created the polys, we can create the NPolyArray object
    
    //this file has elements for all the polynomials which map from the SIEVE to the Q1
    const char* path_db_sv_q1 = "data/csv/db_prod_sv_q1_fp_L_2-2-ord.dat"; 

    //this file has elements for all the polynomials which map from the Q1 to the FP 
    const char* path_db_q1_fp = "data/csv/db_prod_sv_q1_fp_L_2-3-ord.dat"; 
    
    //This is the poly-array which gives us a 'starting point' to use
    const char* path_db_fp_sv = "data/csv/db_center_fp_sv_L_3ord.dat"; 
    NPolyArray poly_array_fp_sv = Parse_NPolyArray_from_file(path_db_fp_sv, branches_sv, DoF_fp); 
    
    NPolyArray* parr_forward = &poly_array_fp_sv; 

    //this is the poly-array which maps directly from the sieve to the focal plane. 
    const char* path_db_sv_fp = "data/csv/db_prod_sv_fp_L_3ord.dat"; 
    NPolyArray poly_array_sv_fp = Parse_NPolyArray_from_file(path_db_sv_fp, branches_fp, DoF_sv);
    

    NPolyArray poly_array_sv_q1 = Parse_NPolyArray_from_file(path_db_sv_q1, branches_q1, DoF_sv);
    NPolyArray poly_array_q1_fp = Parse_NPolyArray_from_file(path_db_q1_fp, branches_fp, DoF_sv);

    
    //if we 'nest' the models, we use the model which is the sv=>q1 model fed into the q1=>fp model. like: fp(q1(sv)). 
    cout << " -- Nesting arrays..." << flush;  

    //NPolyArray arr_model = NPolyArray::Nest( poly_array_q1_fp, poly_array_sv_q1 );
    NPolyArray arr_model = poly_array_sv_fp; 

    //count how many elements there are in this monstrosity
    cout << "done.\n"; 
    int n_elems_total(0); 
    for (int i=0; i<arr_model.Get_DoF_out(); i++) {
        int n_elems = arr_model.Get_poly(i)->Get_nElems(); 
        printf(" > poly %i - n_elems, maxPower: %5i, %i\n", i, n_elems, arr_model.Get_poly(i)->Get_maxPower());
        n_elems_total += n_elems;  
    }

    printf(" -- total number of elements: %i\n", n_elems_total); cout << flush; 


    //if we use this model, we just use the model that maps directly from the SIEVE to the FOCAL PLANE
    NPolyArray *parr = &arr_model; 


    //check that each poly has at least some elements (otherwise, there has been some sort of file-open error)
    for (int i=0; i<parr->Get_DoF_out(); i++) {
        if (parr->Get_poly(i)->Get_nElems() <= 0) {
            Error(here, "NPoly found without elements. Something has gone wrong with the db file..."); 
            return 1; 
        }
    }

    cout << " -- NPolyArray size: " << parr->Get_DoF_in() << " x " << parr->Get_DoF_out() << endl; 


    const int n_iterations = 7; 

    auto rv_dot = [](const RVec<double>& u, const RVec<double>& v) {
        double ret(0.); for (const double& x : u * v) ret += x; return ret; 
    };
    
    auto rv_mag = [](const RVec<double>& u) { 
        double ret(0.); for (const double& x : u * u) ret += x; return sqrt(ret); 
    };

    //number of microns to compute vdc-smearing by
    double vdc_smearing_um = 0.; 

#if 0 
    //'fan out' from the central found trajcetory, to adjacent trajectories. see which will be best. 
    const int n_trajectories=50; 
    const double trajectory_spacing=0.50e-3; 

    auto Find_trajectories = [n_trajectories, trajectory_spacing, parr](const Track_t& Xfp, const Track_t& xsv) 
    {
        //search tracjectories 'forward'
        RVec<Track_t> traj; traj.reserve(n_trajectories); 

        for (int i=0; i<n_trajectories; i++) { 

            RMatirx J( parr->Jacobian(Xsv) ); 
        }
    };
#endif 


    //Define all of the branches we want to create models to map between
    auto df_output = df

#ifdef EVENT_RANGE
        .Range(EVENT_RANGE)
#endif 
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

        .Define("Xsv_first_guess", [parr_forward](Track_t& Xfp)
        {
            auto v = parr_forward->Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz});  
            
            return Track_t{ 
                .x      = v[0],
                .y      = v[1],
                .dxdz   = v[2],
                .dydz   = v[3],
                .dpp    = v[4]
            };
        }, {"Xfp"})

        .Define("Xsv_model", [parr](Track_t& Xfp, Track_t& Xsv)
        {
            RVec<double> Xsv_rvec{
                Xsv.x,
                Xsv.y,
                Xsv.dxdz,
                Xsv.dydz,
                Xsv.dpp    
            }; 

            parr->Iterate_to_root(Xsv_rvec, {Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz}, 8);
            
            return Track_t{
                .x      = Xsv_rvec[0],
                .y      = Xsv_rvec[1],
                .dxdz   = Xsv_rvec[2],
                .dydz   = Xsv_rvec[3],
                .dpp    = Xsv_rvec[4]
            }; 

        }, {"Xfp", "Xsv_first_guess"})

        .Define("Xsv_trajectories",     [](const Track_t& Xfp, const Track_t& Xsv)
        {

        })

        .Define("reco_x_sv",        [](const Track_t& X){ return X.x; },    {"Xsv_model"})
        .Define("reco_y_sv",        [](const Track_t& X){ return X.y; },    {"Xsv_model"})
        .Define("reco_dxdz_sv",     [](const Track_t& X){ return X.dxdz; },    {"Xsv_model"})
        .Define("reco_dydz_sv",     [](const Track_t& X){ return X.dydz; },    {"Xsv_model"})
        .Define("reco_dpp_sv",      [](const Track_t& X){ return X.dpp; },    {"Xsv_model"})
        
        .Define("Xfp_model",  [parr](const Track_t&Xsv)
        {
            auto Xfp = parr->Eval({Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp});
            return Track_t{.x=Xfp[0], .y=Xfp[1], .dxdz=Xfp[2], .dydz=Xfp[3], .dpp=Xfp[4]}; 

        }, {"Xsv_model"})   

        .Define("err_x_fp",     [](Track_t& Xfp, double x){ return (Xfp.x-x)*1e3; }, {"Xfp_model", "x_fp"})
        .Define("err_y_fp",     [](Track_t& Xfp, double x){ return (Xfp.y-x)*1e3; }, {"Xfp_model", "y_fp"})
        .Define("err_dxdz_fp",  [](Track_t& Xfp, double x){ return (Xfp.dxdz-x)*1e3; }, {"Xfp_model", "dxdz_fp"})
        .Define("err_dydz_fp",  [](Track_t& Xfp, double x){ return (Xfp.dydz-x)*1e3; }, {"Xfp_model", "dydz_fp"})

        .Define("err_x_sv",     [](Track_t& Xsv, double x){ return (Xsv.x-x)*1e3; }, {"Xsv_model", "x_sv"})
        .Define("err_y_sv",     [](Track_t& Xsv, double x){ return (Xsv.y-x)*1e3; }, {"Xsv_model", "y_sv"})
        .Define("err_dxdz_sv",  [](Track_t& Xsv, double x){ return (Xsv.dxdz-x)*1e3; }, {"Xsv_model", "dxdz_sv"})
        .Define("err_dydz_sv",  [](Track_t& Xsv, double x){ return (Xsv.dydz-x)*1e3; }, {"Xsv_model", "dydz_sv"})
        .Define("err_dpp_sv",   [](Track_t& Xsv, double x){ return (Xsv.dpp-x)*1e3; }, {"Xsv_model", "dpp_sv"})

        .Define("Xhcs", [is_RHRS](const Track_t& Xsv)
        {
            Track_t Xhcs{Xsv}; 
            SCS_to_HCS(Xhcs, is_RHRS); 
            return Xhcs;    
        }, {"Xsv_model"}) 

        .Define("Xhcs_first_guess", [is_RHRS](const Track_t& Xsv)
        {
            Track_t Xhcs{Xsv}; 
            SCS_to_HCS(Xhcs, is_RHRS); 
            return Xhcs;    
        }, {"Xsv_first_guess"});


    auto df_hcs = add_branch_from_Track_t(df_output, "Xhcs", {
        {"x_hcs",       &Track_t::x},
        {"y_hcs",       &Track_t::y},
        {"dxdz_hcs",    &Track_t::dxdz},
        {"dydz_hcs",    &Track_t::dydz},
        {"dpp_hcs",     &Track_t::dpp}
    }); 

    auto df_hcs_guess = add_branch_from_Track_t(df_output, "Xhcs_first_guess", {
        {"x_hcs_fg",       &Track_t::x},
        {"y_hcs_fg",       &Track_t::y},
        {"dxdz_hcs_fg",    &Track_t::dxdz},
        {"dydz_hcs_fg",    &Track_t::dydz},
        {"dpp_hcs_fg",     &Track_t::dpp}
    });

    //book the histograms we need. 
    char buff_hxy_title[200];  
    sprintf(buff_hxy_title, "Reconstructed sieve coordinates. VDC smearing: %.1f um;x_sv;y_sv", vdc_smearing_um); 
    
    auto h_xy_sieve = df_output.Histo2D<double>({"h_xy_sieve", buff_hxy_title, 250, -45e-3, 45e-3, 250, -45e-3, 45e-3}, "reco_x_sv", "reco_y_sv"); 
    

    auto h_xy_hcs = df_hcs
        .Histo2D<double>({"h_xy_hcs", "Projection of sieve-coords onto z_HCS=0;x_hcs;y_hcs", 250, -25e-3,25e-3, 250, -25e-3,25e-3}, "x_hcs", "y_hcs"); 


    auto h_x    = df_output.Histo1D<double>({"h_x",    "Error of x_fp;mm", 200, -5, 5},      "err_x_sv"); 
    auto h_y    = df_output.Histo1D<double>({"h_y",    "Error of y_fp;mm", 200, -5, 5},      "err_y_sv"); 
    auto h_dxdz = df_output.Histo1D<double>({"h_dxdz", "Error of dxdz_fp;mrad", 200, -5, 5}, "err_dxdz_sv"); 
    auto h_dydz = df_output.Histo1D<double>({"h_dydz", "Error of dydz_fp;mrad", 200, -5, 5}, "err_dydz_sv"); 
    

    //this histogram will be of the actual sieve-coords
    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords. db:'%s', data:'%s'", path_dbfile, path_infile); 
    
    
    gStyle->SetPalette(kSunset);
    //gStyle->SetOptStat(0); 
    new TCanvas("c2", b_c_title); 
    h_xy_sieve->DrawCopy("col2"); 


    new TCanvas("c3", b_c_title); 
    h_xy_hcs->DrawCopy("col2"); 

    auto c = new TCanvas("c1", b_c_title, 1200, 800); 
    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); 
    h_x->DrawCopy(); 
    c->cd(2); 
    h_y->DrawCopy(); 
    c->cd(3); 
    h_dxdz->DrawCopy(); 
    c->cd(4); 
    h_dydz->DrawCopy(); 
    
      

    return 0;
}