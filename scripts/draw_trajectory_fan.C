#include <TGraph.h> 
#include <TLine.h> 
#include <TCanvas.h>
#include <TStyle.h> 
#include <TH2D.h> 
#include <NPolyArray.h> 
#include <ApexOptics.h> 
#include <vector>
#include <string> 
#include <TVector3.h> 
#include <stdexcept> 
#include <ROOT/RVec.hxx>
#include <cmath> 
#include <stdexcept> 
#include <sstream> 
#include <limits> 
#include <utility> 

using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 
using ApexOptics::OpticsTarget_t; 

using ApexOptics::RVec_to_Trajectory_t; 

//compute the magnitude of a ROOT::RVec<double> vector. 
double rvec_mag(const RVec<double>& v) {
    double err=0.; 
    for (double x : v) err += x*x; 
    return sqrt(err); 
}

#define VWIRE 1

int draw_trajectory_fan(const string target_name = "V2", const bool is_RHRS=false)
{
    const char* const here = "draw_trajectory_fan"; 

    OpticsTarget_t target; 
    try { 
    
        target = ApexOptics::GetTarget(target_name); 
    
    } catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }

    const vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    }; 

    const vector<string> branches_fp{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const char* path_poly_fp_sv = "data/csv/poly_prod_fp_sv_L_6ord.dat"; 
    const char* path_poly_sv_fp = "data/csv/poly_prod_sv_fp_L_6ord.dat"; 

    //try to parse the NPolyArray-s from the given db-files
    NPolyArray parr_fp_sv, parr_sv_fp; 
    try {
        parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file(path_poly_fp_sv, branches_sv, 4); 
        parr_sv_fp = ApexOptics::Parse_NPolyArray_from_file(path_poly_sv_fp, branches_fp, 5); 
    
    } catch (const std::exception& e) {

        Error(here, "Exception caught trying to parse NPolyArray db-files.\n what(): %s", e.what()); 
        return -1; 
    }


    //this function takes a starting vertex (in HCS), and a target position on the sieve-face, and given dp/p momentum parameter, 
    //and generates a 'Trajectory_t' in SCS 
    auto Generate_test_trajectory = [is_RHRS](TVector3 react_vertex, double x_sv, double y_sv, double dpp)
    {
        //first, convert the TVector3 to SCS. 
        react_vertex.RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
        react_vertex.RotateZ( TMath::Pi()/2. ); 

        react_vertex += -ApexOptics::Get_sieve_pos(is_RHRS);

        //now compute the 'angles' in scs, using the intercept with the sieve-face
        double dxdz = ( x_sv - react_vertex.x() )/( 0. - react_vertex.z() ); 
        double dydz = ( y_sv - react_vertex.y() )/( 0. - react_vertex.z() ); 

        //return a Trajectory_t struct with all this information 
        return Trajectory_t{ 
            .x    = x_sv, 
            .y    = y_sv, 
            .dxdz = dxdz, 
            .dydz = dydz, 
            .dpp  = dpp
        }; 
    }; 

    const double center_dpp = +0.000; 

    if (target_name == "none") target = ApexOptics::GetTarget("V2"); 

    vector<Trajectory_t> test_trajectories{};

    TVector3 react_vtx_center(
        target.x_hcs,
        0.,
        target.z_hcs
    ); 

    if (target.name == "V1") {

        //V1-wire coords 
        test_trajectories = vector<Trajectory_t>{
            Generate_test_trajectory(react_vtx_center, -0.025, -0.0200, center_dpp),     //lower row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.0200, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.025, -0.0200, center_dpp),

            Generate_test_trajectory(react_vtx_center, -0.025, -0.0075, center_dpp),     //middle row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.0075, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.025, -0.0075, center_dpp),

            Generate_test_trajectory(react_vtx_center, -0.025, +0.0050, center_dpp),     //upper row
            Generate_test_trajectory(react_vtx_center,  0.000, +0.0050, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.025, +0.0050, center_dpp)
        };
    }

    if (target.name == "V2") {

        //V2-wire coords 
        test_trajectories = vector<Trajectory_t>{
            Generate_test_trajectory(react_vtx_center, -0.025, -0.025, center_dpp),     //lower row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.025, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.025, -0.025, center_dpp),

            Generate_test_trajectory(react_vtx_center, -0.025, -0.015, center_dpp),     //middle row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.015, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.025, -0.015, center_dpp),

            Generate_test_trajectory(react_vtx_center, -0.025, -0.005, center_dpp),     //upper row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.005, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.025, -0.005, center_dpp)
        }; 
    }

    if (target.name == "V3") {

        //V3-wire coords
        test_trajectories = vector<Trajectory_t>{
            Generate_test_trajectory(react_vtx_center, -0.020, -0.035, center_dpp),     //lower row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.035, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.020, -0.035, center_dpp),

            Generate_test_trajectory(react_vtx_center, -0.020, -0.025, center_dpp),     //middle row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.025, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.020, -0.025, center_dpp),

            Generate_test_trajectory(react_vtx_center, -0.020, -0.015, center_dpp),     //upper row
            Generate_test_trajectory(react_vtx_center,  0.000, -0.015, center_dpp),
            Generate_test_trajectory(react_vtx_center, +0.020, -0.015, center_dpp)
        };
    } 

    //the max +/- level of dp/p to search. 
    const double dpp_search_range = 0.50e-3; 

    //double this number +1 is the nubmer of 'trajectory points' to draw
    const size_t half_trajectory_points = 20; 

    //set up & draw the histogram
    auto hist_x_y   = new TH2D("h_xy",   "x_{sv} vs y_{sv}",         300, -0.04, 0.04, 300, -0.04, 0.03); 
    auto hist_dx_dy = new TH2D("h_dxdy", "dx/dz_{sv} vs dy/dz_{sv}", 300, -0.04, 0.04, 300, -0.04, 0.03); 

    auto canvas = new TCanvas("c", Form("trajectories; dpp = % .4f +/- %.4f", center_dpp, dpp_search_range), 1400, 600); 
    canvas->Divide(2,1, 0.01,0.01); 

    //draw the 'center' of each trajectory
    for (auto traj : test_trajectories) {
        hist_x_y  ->Fill( traj.x, traj.y ); 
        hist_dx_dy->Fill( traj.dxdz, traj.dydz ); 
    }

    hist_x_y  ->SetMarkerStyle(kOpenCircle);
    hist_dx_dy->SetMarkerStyle(kOpenCircle);

    hist_x_y  ->SetMarkerSize(1.); 
    hist_dx_dy->SetMarkerSize(1.); 

    canvas->cd(1); hist_x_y->Draw("SCAT P"); 
    canvas->cd(2); hist_dx_dy->Draw("SCAT P"); 

    gStyle->SetOptStat(0); 


    //now, lets draw our trajectories. 

    
    //if the difference between an 'actual' and 'desired' Xsv is smaller than this euclidean dist, then quit iterations. 
    const double error_threshold = 1e-10; 

    //let's loop through all of them one-at-a-time.
    for (Trajectory_t traj : test_trajectories) {

        //first, try to find the 'focal plane coordinate' cooresponding to this trajectory. 

        auto Xfp_rvec = parr_sv_fp.Eval({traj.x, traj.y, traj.dxdz, traj.dydz, traj.dpp});  

        auto Xsv_fwd_rvec = parr_fp_sv.Eval(Xfp_rvec); 

        const RVec<double> Xsv_rvec_actual{ traj.x, traj.y, traj.dxdz, traj.dydz, traj.dpp }; 

        Trajectory_t traj_fwd_model{
            .x      = Xsv_fwd_rvec[0],
            .y      = Xsv_fwd_rvec[1],
            .dxdz   = Xsv_fwd_rvec[2],
            .dydz   = Xsv_fwd_rvec[3],
            .dpp    = Xsv_fwd_rvec[4]
        }; 

        cout << "new traj.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
        printf("traj (real vs fwd-model:)  {x,y,dxdz,dydz,dpp}:\n"); 
        printf(" -- {%+.4f,%+.4f,%+.4f,%+.4f,%+.4f}\n", traj.x, traj.y, traj.dxdz, traj.dydz, traj.dpp); 
        printf(" -- {%+.4f,%+.4f,%+.4f,%+.4f,%+.4f}\n", traj_fwd_model.x, traj_fwd_model.y, traj_fwd_model.dxdz, traj_fwd_model.dydz, traj_fwd_model.dpp); 
        
        //now, iterate to a Xsv coordinate which 'matches' the Xfp coordinate. 
        auto Xsv_rvec = Xsv_fwd_rvec;
        
        printf(" after %i iterations, final error: % .3e\n", 
            parr_sv_fp.Iterate_to_root(Xsv_rvec, Xfp_rvec, 15, error_threshold), 
            rvec_mag(Xfp_rvec - parr_sv_fp.Eval(Xsv_rvec))
        ); 

        //now, let's do our little dance of finding 'adjacent' trajectories to draw. 
        auto Find_next_traj = [&](RVec<double> Xsv, const double d_dpp) 
        {
            auto J = parr_sv_fp.Jacobian(Xsv_rvec); 

            //the last column of this jacobian corresponds to dp/p. since this is the variable we will be varying, we will 
            RVec<double> J4{ 
                -J.at(0,4),
                -J.at(1,4), 
                -J.at(2,4), 
                -J.at(3,4)
            }; 

            RMatrix Ji(4,4, {
                J.at(0,0), J.at(0,1), J.at(0,2), J.at(0,3),
                J.at(1,0), J.at(1,1), J.at(1,2), J.at(1,3),
                J.at(2,0), J.at(2,1), J.at(2,2), J.at(2,3),
                J.at(3,0), J.at(3,1), J.at(3,2), J.at(3,3)
            }); 

            auto dX = Ji.Solve( J4 ); 

            if (dX.empty()) return RVec<double>{}; 
            for (double x : dX) { if (x != x) return RVec<double>{}; }
            
            dX *= d_dpp; 

            dX.push_back( d_dpp ); 

            Xsv += dX; 

            //now, iterate to 'correct' whatever small error may be in our linear extrapolation 'Xsv += dX' 
            parr_sv_fp.Iterate_to_root(Xsv, Xfp_rvec, 15, 1e-10); 

            return Xsv; 
        }; 
        //_________________________________________________________________________________________________
        
        
        const double d_dpp = dpp_search_range/((double)half_trajectory_points-1); 

        //if the next vector is further than this, then quit finding new points (something is wrong!)
        const double max_dist_to_next = d_dpp * 20.; 

        vector<Trajectory_t> trajectories_forward; 

        auto Xsv_next = Xsv_rvec; 
        
        for (size_t i=0; i<half_trajectory_points; i++) {
            
            Xsv_next = Find_next_traj(Xsv_next, d_dpp); 

            //this happens if the 'Find_next_traj' fails for some reason. if so, stop finding new points. 
            if (Xsv_next.size() != 5) break; 

            trajectories_forward.push_back( RVec_to_Trajectory_t(Xsv_next) ); 
        }

        //now, find the 'backward' trajectories
        vector<Trajectory_t> trajectories_backward;
        
        Xsv_next = Xsv_rvec; 

        for (size_t i=0; i<half_trajectory_points; i++) {

            auto Xsv_new = Find_next_traj(Xsv_next, -1.*d_dpp); 

            //this happens if the 'Find_next_traj' fails for some reason. if so, stop finding new points. 
            if (Xsv_next.size() != 5) break; 
            
            //check if the new vector has changed too much. if so, something has gone wrong in the iteration process
            if (rvec_mag(Xsv_next - Xsv_new) > max_dist_to_next) break; 

            //check if the momentum is 'out of range' of the maximum search range
            if (fabs(Xsv_next[4] - traj.dpp) > dpp_search_range) break; 

            Xsv_next = Xsv_new; 

            trajectories_backward.push_back( RVec_to_Trajectory_t(Xsv_next) );
        }

        //now we're done. make the vectors. 
        cout << "number of trajectories (fwd/back): " 
             << trajectories_forward.size() << " / " 
             << trajectories_backward.size() << endl; 

        //now, combine them all in one vector. 
        vector<Trajectory_t> trajectories;
        trajectories.reserve( trajectories_forward.size() + trajectories_backward.size() + 1 );
        
        for (int i=trajectories_backward.size()-1; i>=0; i--) 
            trajectories.push_back( trajectories_backward.at(i) ); 
        
        trajectories.push_back( RVec_to_Trajectory_t(Xsv_rvec) ); 

        for (int i=0; i<trajectories_forward.size(); i++) 
            trajectories.push_back( trajectories_forward.at(i) ); 
        
        //now, we have all the trajectories, both backward and forward. 
        const size_t n_pts = trajectories.size(); 

        double *pts_x    = new double[n_pts];
        double *pts_y    = new double[n_pts];
        double *pts_dxdz = new double[n_pts];
        double *pts_dydz = new double[n_pts];
        
        int i=0; 
        for (Trajectory_t traj : trajectories) {

            pts_x[i] = traj.x;
            pts_y[i] = traj.y;
            pts_dxdz[i] = traj.dxdz; 
            pts_dydz[i] = traj.dydz; 

            i++; 
        }

        auto graph_x_y      = new TGraph(n_pts, pts_x, pts_y);
        auto graph_dx_dy    = new TGraph(n_pts, pts_dxdz, pts_dydz);

        graph_x_y   ->SetLineColor(kBlack); 
        graph_dx_dy ->SetLineColor(kBlack); 

        canvas->cd(1); graph_x_y->Draw("SAME");
        canvas->cd(2); graph_dx_dy->Draw("SAME"); 
    }
    cout << endl; 

    return 0; 
}