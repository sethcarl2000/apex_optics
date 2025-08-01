#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

struct Track_t {
    double x,y,dxdz,dydz,dpp; 

    //implicit conversion from this type to RVec<double> 
    operator RVec<double>() const { return RVec<double>{x,y,dxdz,dydz,dpp}; }  
    
    Track_t operator-(const Track_t& rhs) const {   
        return Track_t{
            .x      = x     - rhs.x,
            .y      = y     - rhs.y, 
            .dxdz   = dxdz  - rhs.dxdz,
            .dydz   = dydz  - rhs.dydz,
            .dpp    = dpp   - rhs.dpp
        }; 
    }; 

    Track_t operator+(const Track_t& rhs) const {   
        return Track_t{
            .x      = x     + rhs.x,
            .y      = y     + rhs.y, 
            .dxdz   = dxdz  + rhs.dxdz,
            .dydz   = dydz  + rhs.dydz,
            .dpp    = dpp   + rhs.dpp
        }; 
    }; 


    Track_t operator*(double mult) const {   
        return Track_t{
            .x      = x     * mult,
            .y      = y     * mult, 
            .dxdz   = dxdz  * mult,
            .dydz   = dydz  * mult,
            .dpp    = dpp   * mult
        }; 
    };

};

class TrajectoryCollection : public TObject {
public: 
    TrajectoryCollection(const ROOT::RVec<Track_t>& _trajectories) 
        : fN_elems{(int)_trajectories.size()}, fTrajectories{_trajectories} {

        dx = ((double)fN_elems-1); 
    };
    
    ~TrajectoryCollection() {}; 

    Track_t Get(double idx) const {
        if (idx >= 1.00 || idx < 0.) { Error("Get", "Invalid index given: %f, must be [0,1).", idx);
            return Track_t{};
        } 

        double remainder = dx * idx; 
        int index = (int)remainder; 
        remainder += -1.*(double)index; 

        const Track_t & t1 = fTrajectories[index+1]; 
        const Track_t & t0 = fTrajectories[index]; 

        return t0 + ((t1 - t0) * remainder);  
    }; 

private:    
    int fN_elems; 
    RVec<Track_t> fTrajectories; 

    double dx; 
};

//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Track_t struct. 
void add_branch_from_Track_t(   std::vector<ROOT::RDF::RNode>& df_nodes, 
                                const char* branch_in, 
                                map<string, double Track_t::*> branches )
{
    const int n_nodes = branches.size() + 1; 
    
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Track_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Track_t::coord' 
        auto new_node = df_nodes.back()

            .Define(branch_name, [coord](const Track_t& track) { return track.*coord; }, {branch_in}); 

        df_nodes.push_back(new_node); 
    }
    
    return; 
}

#define EVENT_RANGE 30000

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
                            const char* path_dbfile="data/csv/db_prod_sv_fp_L_3ord.dat",  
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
    Info(here, "Multi-threadding is disabled. Max events to run: %i", EVENT_RANGE); 
#else 
    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size : %i", ROOT::GetThreadPoolSize()); 
#endif 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

#ifdef EVENT_RANGE
    const int n_events_max = *df.Count(); 
    const int n_events = min<int>(EVENT_RANGE, n_events_max); 
#else 
    const int n_events = *df.Count(); 
#endif

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

    auto rv_unit = [&rv_mag](const RVec<double>& u) {
        return u / rv_mag(u); 
    }; 

    //number of microns to compute vdc-smearing by
    double vdc_smearing_um = 0.; 

    const double vertex_uncertainty_x = 1.00e-3;
    const double vertex_uncertainty_y = 1.00e-3;  

 
    //'fan out' from the central found trajcetory, to adjacent trajectories. see which will be best. 
    const int n_trajectories=50; 
    const double trajectory_spacing=0.250e-3; 

    auto Find_trajectories = [  n_trajectories,     
                                trajectory_spacing, 
                                parr, 
                                DoF_fp, 
                                DoF_sv, 
                                &rv_mag ](const Track_t& Xfp, const Track_t& Xsv) 
    {
        
        RVec<double> Xfp_rv{ Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz }; 

        auto Get_next_trajectory = [parr, trajectory_spacing, &rv_mag, &Xfp_rv](RVec<Track_t>* t_vec, RVec<double>& Xsv, double oreintation)
        {
            RVec<double> J_arr{ *(parr->Jacobian(Xsv).Data()) }; 
            int i_elem=0; 

            //this jacobian is 4x5, as we're mapping from R^5 => R^4. we're gonna 'peel off' the last column, and
            //then solve the resultant 4x4 matrix. 
            RVec<double> Ji_arr;    Ji_arr.reserve(DoF_fp*DoF_fp);  
            RVec<double> J0;        J0.reserve(DoF_fp); 

            for (int i=0; i<DoF_fp; i++) {
                for (int j=0; j<DoF_fp; j++) Ji_arr.push_back( J_arr[i_elem++] ); 
            
                J0.push_back( J_arr[i_elem++] ); 
            }

            RMatrix Ji(DoF_fp, DoF_fp, Ji_arr); Ji.ReportSingular()=false; 

            auto dX = Ji.Solve( J0 ); 

            //chekc for null vector
            if (dX.size() != DoF_fp) return 0; 
            //check for NaN vals
            for (double& x : dX) if ( x != x ) return 0; 
            
            dX.push_back( -1. ); 
            dX *= oreintation * trajectory_spacing/rv_mag(dX);  

            //now that we have a new trajectory 'suggestion', lets use it: 
            Xsv += dX;

            //now, we perform a few iterations to 'fix' it
            parr->Iterate_to_root( Xsv, Xfp_rv, 3 ); 

            t_vec->push_back({
                .x      = Xsv[0], 
                .y      = Xsv[1],
                .dxdz   = Xsv[2],
                .dydz   = Xsv[3],
                .dpp    = Xsv[4]
            }); 
            return 1; 
        };

        //0bvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv q  AA1
        RVec<double> Xsv_rv{ Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp };

        ///get 'forward' trajectories
        RVec<Track_t> traj; traj.reserve(n_trajectories);  
        int i_traj=0; 
        while ( ++i_traj < n_trajectories ) {
            if (Get_next_trajectory(&traj, Xsv_rv, 1. ) != 1) break; 
        } 

        ///get 'backward' trajectories
        Xsv_rv = Track_t{ Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp }; 
        RVec<Track_t> back_traj; back_traj.reserve(n_trajectories); 
        
        i_traj=0;
        while ( ++i_traj < n_trajectories ) {
            if (Get_next_trajectory(&back_traj, Xsv_rv, -1.) != 1) break; 
        }

        RVec<Track_t> trajectories; 
        trajectories.reserve(traj.size() + back_traj.size() + 1); 

        //add the 'back aspect' tracks
        for (int i=back_traj.size()-1; i>=0; i--) trajectories.push_back( back_traj[i] ); 

        //add the 'central' (starting) tracks
        trajectories.push_back( Xsv ); 

        //add the 'front aspect' tracks
        copy( traj.begin(), traj.end(), back_inserter(trajectories) ); 

        return trajectories; 
    };


#ifdef EVENT_RANGE
    printf("Processing %i events...", min<int>(EVENT_RANGE, n_events)); 
#else 
    printf("Processing %i events...", n_events); 
#endif
    cout << flush; 

    //Define all of the branches we want to create models to map between
    auto df_proc = df

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

        .Define("Xsv_first_guess", [parr_forward, parr](Track_t& Xfp)
        {
            auto Xsv = parr_forward->Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz});  
            
            parr->Iterate_to_root(Xsv, {Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz}, 8);
            
            return Track_t{
                .x      = Xsv[0],
                .y      = Xsv[1],
                .dxdz   = Xsv[2],
                .dydz   = Xsv[3],
                .dpp    = Xsv[4]
            }; 

        }, {"Xfp"})

        .Define("Xsv_trajectories",     [&Find_trajectories](const Track_t& Xfp, const Track_t& Xsv)
        {
            return Find_trajectories(Xfp, Xsv); 

        }, {"Xfp", "Xsv_first_guess"})  

        .Define("Xsv_best",     [ is_RHRS,
                                  vertex_uncertainty_x,
                                  vertex_uncertainty_y ](const RVec<Track_t>& traj, TVector3 vtx)
        {
            const Track_t *Xsv_best = nullptr;
            double err2_best = 1e30; 

            double vx = vtx.x() + gRandom->Gaus() * vertex_uncertainty_x; 
            double vy = vtx.y() + gRandom->Gaus() * vertex_uncertainty_y; 

            //look for the best track
            for (const Track_t& track_scs : traj) {
                
                Track_t track_hcs{track_scs}; 

                SCS_to_HCS(track_hcs, is_RHRS); 
                
                //project our 'track' onto the z-plane of the react vertex
                double err2(0.); 
                err2 += pow( ((track_hcs.x  +  track_hcs.dxdz * vtx.z())    -   vx)/vertex_uncertainty_x, 2 );
                err2 += pow( ((track_hcs.y  +  track_hcs.dydz * vtx.z())    -   vy)/vertex_uncertainty_y, 2 ); 

                if (err2 < err2_best) {
                    Xsv_best = &track_scs;
                    err2_best = err2; 
                }
            }

            if (!Xsv_best) {
                Error("Define(Xsv_best)", "No good sieve coord found!"); return Track_t{}; 
            }

            return *Xsv_best; 

        }, {"Xsv_trajectories", "position_vtx"})

        .Define("n_trajectories", [](const RVec<Track_t>& traj){ return (int)traj.size(); }, {"Xsv_trajectories"})

        .Define("Xfp_model",  [parr](const Track_t&Xsv)
        {
            auto Xfp = parr->Eval({Xsv.x, Xsv.y, Xsv.dxdz, Xsv.dydz, Xsv.dpp});
            return Track_t{.x=Xfp[0], .y=Xfp[1], .dxdz=Xfp[2], .dydz=Xfp[3], .dpp=Xfp[4]}; 

        }, {"Xsv_best"})   

        .Define("err_x_fp",     [](Track_t& Xfp, double x){ return (Xfp.x-x)*1e3; }, {"Xfp_model", "x_fp"})
        .Define("err_y_fp",     [](Track_t& Xfp, double x){ return (Xfp.y-x)*1e3; }, {"Xfp_model", "y_fp"})
        .Define("err_dxdz_fp",  [](Track_t& Xfp, double x){ return (Xfp.dxdz-x)*1e3; }, {"Xfp_model", "dxdz_fp"})
        .Define("err_dydz_fp",  [](Track_t& Xfp, double x){ return (Xfp.dydz-x)*1e3; }, {"Xfp_model", "dydz_fp"})

        .Define("Xsv_error", [](Track_t Xsv_mod, Track_t Xsv)
        {   
            return (Xsv_mod - Xsv) * 1e3; 
        }, {"Xsv_best", "Xsv"})

        .Define("Xhcs", [is_RHRS](const Track_t& Xsv)
        {
            Track_t Xhcs{Xsv};
            SCS_to_HCS(Xhcs, is_RHRS);
            return Xhcs;
        }, {"Xsv_best"}) 

        .Define("Xhcs_first_guess", [is_RHRS](const Track_t& Xsv)
        {
            Track_t Xhcs{Xsv};
            SCS_to_HCS(Xhcs, is_RHRS);
            return Xhcs;
        }, {"Xsv_first_guess"});


    vector<ROOT::RDF::RNode> output_nodes{ df_proc };  

    add_branch_from_Track_t( output_nodes, "Xhcs", {
        {"x_hcs",       &Track_t::x},
        {"y_hcs",       &Track_t::y},
        {"dxdz_hcs",    &Track_t::dxdz},
        {"dydz_hcs",    &Track_t::dydz},
        {"dpp_hcs",     &Track_t::dpp}
    }); 

    add_branch_from_Track_t( output_nodes, "Xsv_first_guess", {
        {"x_sv_fg",       &Track_t::x},
        {"y_sv_fg",       &Track_t::y},
        {"dxdz_sv_fg",    &Track_t::dxdz},
        {"dydz_sv_fg",    &Track_t::dydz},
        {"dpp_sv_fg",     &Track_t::dpp}
    }); 

    add_branch_from_Track_t( output_nodes, "Xsv_best", {
        {"reco_x_sv",       &Track_t::x},
        {"reco_y_sv",       &Track_t::y},
        {"reco_dxdz_sv",    &Track_t::dxdz},
        {"reco_dydz_sv",    &Track_t::dydz},
        {"reco_dpp_sv",     &Track_t::dpp}
    }); 

    add_branch_from_Track_t( output_nodes, "Xsv_error", {
        {"err_x_sv",       &Track_t::x},
        {"err_y_sv",       &Track_t::y},
        {"err_dxdz_sv",    &Track_t::dxdz},
        {"err_dydz_sv",    &Track_t::dydz},
        {"err_dpp_sv",     &Track_t::dpp}
    }); 

    auto df_output = output_nodes.back(); 

    //book the histograms we need. 
    char buff_hxy_title[200];  
    sprintf(buff_hxy_title, "Reconstructed sieve coordinates. VDC smearing: %.1f um;x_sv;y_sv", vdc_smearing_um); 
    
    auto h_xy_sieve = df_output
        .Histo2D<double>({"h_xy_sieve", "Sieve coordinates (best fit);x_sv;y_sv", 250, -45e-3, 45e-3, 250, -45e-3, 45e-3}, "reco_x_sv", "reco_y_sv"); 

    auto h_xy_sieve_fg = df_output
        .Histo2D<double>({"h_xy_sieve", "Sieve coordinates (first guess);x_sv;y_sv", 250, -45e-3, 45e-3, 250, -45e-3, 45e-3}, "x_sv_fg", "y_sv_fg"); 

    auto h_xy_hcs = df_output
        .Histo2D<double>({"h_xy_hcs", "Projection of sieve-coords onto z_HCS=0;x_hcs;y_hcs", 250, -25e-3,25e-3, 250, -25e-3,25e-3}, "x_hcs", "y_hcs"); 

    auto h_n_trajectories = df_output
        .Histo1D<int>({"h_n_traj", "Number of trajectories generated", 121, -0.5, 120.5}, "n_trajectories"); 


    auto h_x    = df_output.Histo1D<double>({"h_x",    "Error of x_fp;mm", 200, -5, 5},      "err_x_sv"); 
    auto h_y    = df_output.Histo1D<double>({"h_y",    "Error of y_fp;mm", 200, -5, 5},      "err_y_sv"); 
    auto h_dxdz = df_output.Histo1D<double>({"h_dxdz", "Error of dxdz_fp;mrad", 200, -5, 5}, "err_dxdz_sv"); 
    auto h_dydz = df_output.Histo1D<double>({"h_dydz", "Error of dydz_fp;mrad", 200, -5, 5}, "err_dydz_sv"); 
    

    //this histogram will be of the actual sieve-coords
    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords. db:'%s', data:'%s'", path_dbfile, path_infile); 


    
    gStyle->SetPalette(kSunset);
    //gStyle->SetOptStat(0); 
    new TCanvas("c4", b_c_title); 
    
    TStopwatch timer; 
    
    h_xy_sieve->DrawCopy("col2"); 
    
    double realtime( timer.RealTime() ); 
    printf("done.\nprocessed %i events in %f seconds (%f ms/event) - %i thread(s)\n", 
        n_events, 
        realtime, 
        1e3 * realtime / ((double)n_events), 
        (ROOT::GetThreadPoolSize()==0 ? 1 : ROOT::GetThreadPoolSize()) 
    );
    cout << endl;  
    
    
    new TCanvas("c2", b_c_title); 
    h_xy_sieve_fg->DrawCopy("col2"); 


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