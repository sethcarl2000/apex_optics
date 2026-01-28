

//APEX-otpics headers
#include <AngleRecoTester.h> 
#include <ApexOptics.h> 
#include <RMatrix.h>
//ROOT headers
#include <TBox.h> 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <ROOT/RVec.hxx>
#include <TParameter.h>
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h>
#include <TVirtualPad.h> 
#include <SieveHole.h>
#include <TLine.h> 
#include <TAxis.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
//C++ std-lib headers
#include <iostream>
#include <sstream>
#include <string> 
#include <memory>
#include <limits> 
#include <optional> 
#include <functional> 
#include <stdexcept> 
#include <sstream> 
#include <cmath> 
#include <iostream> 

using namespace std; 
using ApexOptics::Trajectory_t; 
using ApexOptics::OpticsTarget_t; 
using namespace ROOT::VecOps; 

namespace {
    //returns true if a floating-point type is NAN
    template<typename T> bool is_nan(T x) { return (x != x); }

    constexpr double NaN_double = std::numeric_limits<double>::quiet_NaN(); 



    //These are for the 'AngleRecoTester::MeasureSlope' alg, we need to invert a 2x2 matrix. this is the cutoff for a 
    //
    //singular matrix. 
    constexpr double singular_determinant_cutoff = 1e-8;  
    //
    //the max number of newton iterations (for slope-fitting) before ceasing
    constexpr int max_newton_iterations = 10; 
    //
    //minimum height ratio between the max for a hole-peak, and the height of any of the veritcal slices 
    constexpr double min_height_ratio = 0.15; 
}


//____________________________________________________________________________________________________________________________
AngleRecoTester::AngleRecoTester(bool _is_RHRS, ApexOptics::OpticsTarget_t _target, ROOT::RDF::RNode _dataframe)
    : 
    fIs_RHRS{_is_RHRS}, 
    fTarget{_target},
    fNode{_dataframe}, 
    //flags which are enabled by default.
    fStateFlag{kDraw_slopes | kDraw_boxes | kDraw_centroid_lines}
{

}

//____________________________________________________________________________________________________________________________
std::optional<AngleFitResult_t> AngleRecoTester::Measure(
    int32_t flag, 
    const int first_row,                                    
    const int last_row,      
    const int first_col,
    const int last_col,
    const double row_cut_width,
    const double col_cut_width,
    const double bg_cut_width
)
{
    //see if we're measuring dydz
    const bool test_dydz = (flag & kDydz);
    
    //see if we're measuring slopes
    const bool measure_slope = (flag & kSlopes); 

    //measuring slopes while also measuring dx/dz is no longer supported
    if ((flag & kDxdz) && (flag & kSlopes)) { 
        Error(__func__, "Measuring slope wile measuring dx/dz error is no longer supported. Switch to measuring dy/dz to also measure slopes."); 
        return std::nullopt; 
    } 

    const char* const test_name = test_dydz ? "dydz" : "dxdz"; 

    char c_title[255]; 
    
    //check that the fTarget given is a V-wire
    if (!(fTarget.name == "V1" || 
        fTarget.name == "V2" || 
        fTarget.name == "V3"))  {
        
        Error(__func__, "Target provided is not a V-wire target: %s", fTarget.name.c_str()); 
        return std::nullopt; 
    }

    //RDFNodeAccumulator rna(df); 

    auto node_vtx = fNode

        .Define("vtx_scs", [this](const TVector3& vtx_hcs)
        {
            return ApexOptics::HCS_to_SCS(fIs_RHRS, vtx_hcs); 
        }, {"position_vtx"})

        .Define("vtx_u", [test_dydz](const TVector3& v){ return test_dydz ? v.y() : v.x(); }, {"position_vtx_scs"});

    const double vtx_z = *node_vtx.Define("vtx_z", [](const TVector3& v){ return v.z(); }, {"position_vtx_scs"}).Mean("vtx_z"); 

    const double vtx_y_hcs = *node_vtx.Define("y_hcs", [](const TVector3& v){ return v.y(); }, {"position_vtx"}).Mean("y_hcs"); 

    const TVector3 react_vtx_scs = ApexOptics::HCS_to_SCS( fIs_RHRS, TVector3(fTarget.x_hcs, vtx_y_hcs, fTarget.z_hcs) ); 

    const double vtx_u = test_dydz ? react_vtx_scs.y() : react_vtx_scs.x(); 
    
    const double min_vtx_u = *node_vtx.Min("vtx_u");
    const double max_vtx_u = *node_vtx.Max("vtx_u");

    //set the color pallete
    if (fDo_drawing) {
        gStyle->SetPalette(kSunset); 
        gStyle->SetOptStat(0); 
    }

    //now, we can compute where each sieve-hole SHOULD be
    std::vector<SieveHole> sieve_holes = ApexOptics::ConstructSieveHoles(fIs_RHRS); 

    auto FindSieveHole = [&sieve_holes](int row, int col) 
    {
        auto it = std::find_if( 
            sieve_holes.begin(),
            sieve_holes.end(), 
            [row,col](const SieveHole& elem){ return row==elem.row && col==elem.col; }
        ); 
        if (it == sieve_holes.end()) {
            return SieveHole{-1,-1}; 
        } 
        return SieveHole{*it}; 
    }; 

    //all units in meters
    const double sieve_thickness_z = 12e-3;         //sieve thickness (m)
    const double sieve_hole_spacing_y = 4.826e-3;   //sieve-hole spacing, m
    const double sieve_hole_spacing_x = 5.842e-3;   //sieve-hole spacing, m
    
    //in radians
    const double dxdz_cut_width = (sieve_hole_spacing_x * row_cut_width) / fabs(vtx_z); 
    const double dydz_cut_width = (sieve_hole_spacing_y * col_cut_width) / fabs(vtx_z); 

    // 'u' is the direction we're evaluating accuracy in, and 'v' is the other direction we're not fitting in. 
    const double du_cut_width = test_dydz ? dydz_cut_width : dxdz_cut_width; 
    const double dv_cut_width = test_dydz ? dxdz_cut_width : dydz_cut_width; 
    
    if (fDo_drawing) {
        gStyle->SetPalette(kBird);
    }

    
    //Stores a single dx/dz and dy/dz data point
    struct AnglePair_t { double dxdz,dydz; }; 

    //make a vector of all events 
    //if (!fQuiet) std::cout << "Making copy of data..." << flush; 

    std::vector<AnglePair_t> data = *fNode 

        .Define("event", [](double dxdz, double dydz)
        {
            return AnglePair_t{dxdz,dydz}; 
        }, {fBranch_dxdz, fBranch_dydz})
        
        .Take<AnglePair_t>("event");  

    //if (!fQuiet) std::cout << "done." << endl; 

        
    TCanvas* c{nullptr};
    
    //come up with an appropriate histogram title
    const char* dxdy_hist_name = test_dydz
        ? "Test of dy/dz (theta-tg, lab-horizontal);dx/dz (rad);dy/dz (rad)" 
        : "Test of dx/dz (phi-tg, lab-vertical);dx/dz (rad);dy/dz (rad)"; 
    

    //create, fill, and draw the 2d-hist
    auto h_dx_dy = new TH2D(Form("h_2d_%s",test_name), dxdy_hist_name, 200, fDxdz_min, fDxdz_max, 200, fDydz_min, fDydz_max); 
    
    for (const auto& ev : data) { h_dx_dy->Fill( ev.dxdz, ev.dydz ); }

    //we use these variables do easily switch between our different sub-pads
    TVirtualPad *c_2d, *c_profiles;  
    if (fDo_drawing) {

        c = new TCanvas(Form("ctest_%s",test_name), "Angle Reco - dy/dz", 1600, 800); 
        c->Divide(2,1);

        c_2d         = c->cd(1);
        c_profiles   = c->cd(2);
        
        c_2d->cd(); 
        c_2d->SetLeftMargin(0.15); 
        c_2d->SetRightMargin(0.02); 
        h_dx_dy->DrawCopy("col"); 
    }

    size_t n_holes_measured =0; 
    
    int i_canv=1;

    //row counter
    
    //order of polynomial to fit the backgroudn with.  
    const int background_polynomial_order = 10; 
    
    AngleFitResult_t fit_result; 
    
    const int u0 = test_dydz ? first_col : first_row; 
    const int u1 = test_dydz ? last_col  : last_row; 

    const int v0 = test_dydz ? first_row : first_col; 
    const int v1 = test_dydz ? last_row  : last_col; 

    const double du_max = test_dydz ? fDydz_max : fDxdz_max; 
    const double du_min = test_dydz ? fDydz_min : fDxdz_min;
    
    const double dv_max = test_dydz ? fDxdz_max : fDydz_max; 
    const double dv_min = test_dydz ? fDxdz_min : fDydz_min;

    //this allows us to generalize some of the code below. 
    double AnglePair_t::*du_ap = test_dydz ? &AnglePair_t::dydz : &AnglePair_t::dxdz; 
    double AnglePair_t::*dv_ap = test_dydz ? &AnglePair_t::dxdz : &AnglePair_t::dydz; 

    double SieveHole::*u_sh = test_dydz ? &SieveHole::y : &SieveHole::x; 
    double SieveHole::*v_sh = test_dydz ? &SieveHole::x : &SieveHole::y;

    double Trajectory_t::*du_trj = test_dydz ? &Trajectory_t::dydz : &Trajectory_t::dxdz; 
    double Trajectory_t::*dv_trj = test_dydz ? &Trajectory_t::dxdz : &Trajectory_t::dydz; 
    
    //each one of these vectors is a different slice to fit
    std::vector<std::vector<SieveHole>> hole_slices; 

    for (int iv = v0; iv <= v1; iv++) {

        //Backstory - I made a terrible decision a few years ago when I decided that sieve-holes on adjacent, staggered rows 
        // ought to have the same column-index, and so we need to enage in this funny business of fitting TWO histogram slices per row, 
        // but ONLY if we're cheking dx/dz (because for dy/dz all holes are in the same column)
        if (test_dydz) {
            //for testing dy/dz, all holes are in-line with each other
            vector<SieveHole> slice{}; 
            for (int iu = u0; iu <= u1; iu++) {
            
                int row = test_dydz ? iv : iu;
                int col = test_dydz ? iu : iv;
                
                auto hole = FindSieveHole(row,col); 
                if (hole.row>=0 && hole.col>=0) slice.push_back(hole);               
            }
            hole_slices.push_back(slice);    
        } else {
            //for fitting dxdz, the holes in the same column are staggered. 
            vector<SieveHole> slice_even, slice_odd{}; 
            for (int iu = u0; iu <= u1; iu++) {
            
                int row = test_dydz ? iv : iu;
                int col = test_dydz ? iu : iv;
                
                auto hole = FindSieveHole(row,col); 
                if (hole.row>=0 && hole.col>=0) {
                    if (hole.row%2==0) { slice_even.push_back(hole); } else { slice_odd.push_back(hole); }
                }               
            }
            hole_slices.push_back(slice_even); 
            hole_slices.push_back(slice_odd); 
        }
    }

    if (fDo_drawing) c_profiles->Divide( 1, hole_slices.size(), 0.01,0.); 

    // Fit all holes - dydz
    int icanv_profile=hole_slices.size(); 
    for (const auto& holes : hole_slices) {

        vector<Trajectory_t> hole_angles;
        for (auto& hole : holes) { 
            hole_angles.emplace_back(Trajectory_t{
                hole.x, 
                hole.y,
                ( hole.x - react_vtx_scs.x() )/( 0. - vtx_z ),
                ( hole.y - react_vtx_scs.y() )/( 0. - vtx_z ) 
            }); 
        }
                
        // get the dudz-value of the row. for fitting dy/dz, all holes have the same 'dv'. 
        const double center_dv = hole_angles.front().*dv_trj; 
        // this is needed for testing dx/dz 

        // create and fill the row histogram
        auto hist_slice = unique_ptr<TH1D>(new TH1D(Form("h_slice_%i",icanv_profile), "", 200, du_min, du_max)); 
        
        auto hist_slice_bg = (TH1D*)hist_slice->Clone("h_row_bg"); 

        for (const auto& ev : data) {
            
            //all holes are in a straight line
            if (fabs(ev.*dv_ap - center_dv) < dv_cut_width/2.)              hist_slice->Fill( ev.*du_ap ); 
            if (fabs(ev.*dv_ap - center_dv) < bg_cut_width*dv_cut_width/2.) hist_slice_bg->Fill( ev.*du_ap ); 
        }

        //if (!fQuiet) printf("Drawing slice %i...",icanv_profile); cout << flush; 
        
        if (fDo_drawing) {
            hist_slice->GetYaxis()->SetNdivisions(0); 
            c_profiles->cd(icanv_profile); hist_slice->DrawCopy(); 
        }
        auto x_axis = hist_slice->GetXaxis(); 
        
        // now, try to fit a polynomial to the background.
        // do do this, we must make a list of 'excluded' regions, where we know the holes ought to be. 
        vector<array<double,2>> exclusions;
        for (int i=0; i<hole_angles.size(); i++) {

            exclusions.push_back({
                hole_angles[i].*du_trj - du_cut_width/2., 
                hole_angles[i].*du_trj + du_cut_width/2.
            });
        }
        //if (!fQuiet) cout << "Fitting background..." << flush; 

        //fit the background with a polynomial, 
        auto background_fcn_unnormalized = FitBackground(hist_slice_bg, exclusions, background_polynomial_order);
        auto background_fcn = [background_fcn_unnormalized, bg_cut_width](double *x, double *par){ return background_fcn_unnormalized(x,par)/bg_cut_width; };
        auto tf1_bg = new TF1(Form("bg_poly_%i",icanv_profile), background_fcn, x_axis->GetXmin(), x_axis->GetXmax(), 0); 
        delete hist_slice_bg; 
        tf1_bg->SetLineStyle(kDashed);
        tf1_bg->SetLineColor(kRed); 
        if (fDo_drawing) { c_profiles->cd(icanv_profile); tf1_bg->DrawCopy("SAME"); } 
        

        //if (!fQuiet) cout << "done." << endl; 

        for (size_t i=0; i<holes.size(); i++) {

            //the 'SieveHole' struct associated with this hole
            const auto& hole  = holes[i]; 

            //the 'Trajectory_t' struct associated with this hole
            const auto& angle = hole_angles[i]; 

            auto fcn_gauss_offset = [background_fcn](double *x, double *par) { 
                double sigma = fabs(par[2]); 
                return par[0] * exp( -0.5 * pow( (x[0] - par[1])/sigma, 2 ) ) + background_fcn(x,nullptr); 
            }; 

            //get the maximum value in this range
            int bin_min = x_axis->FindBin( angle.*du_trj - du_cut_width/2. ); 
            int bin_max = x_axis->FindBin( angle.*du_trj + du_cut_width/2. ); 

            double maxval = -1.; 
            int max_bin = -1; 
            for (int bin=bin_min; bin<=bin_max; bin++) {

                if (hist_slice->GetBinContent(bin) > maxval) {
                    maxval = hist_slice->GetBinContent(bin); 
                    max_bin = bin; 
                }
            }
            //cout << "Max bin val: " << maxval << endl; 
            double x0 = x_axis->GetBinCenter(max_bin); 

            //double offset = (hist_row->GetBinContent(bin_min) + hist_row->GetBinContent(bin_max))/2.; 

            //upper edge of sieve hole, accounting for paralax, assuming that the tungsten material of the sieve is impenitrable
            // (which is known not to be true!!)
            const double R_hole = hole.radius_front; 
            
            double hole_du_min = (hole.*u_sh - R_hole - max_vtx_u) / ( (max_vtx_u < hole.*u_sh - R_hole ? 0. : sieve_thickness_z) - vtx_z ); 
            double hole_du_max = (hole.*u_sh + R_hole - min_vtx_u) / ( (min_vtx_u < hole.*u_sh + R_hole ? sieve_thickness_z : 0.) - vtx_z ); 
            

            //aparent width of the hole (if sieve is impermiable)
            const double hole_angle_width = hole_du_max - hole_du_min;
            
            //aparent centroid of this hole (under the assumption that the sieve is impenitrable)
            const double hole_du = (hole_du_max + hole_du_min) / 2.;
                
            //auto fit = unique_ptr<TF1>(new TF1("holefit", gaussian_semicircle, hole_dydz -2.5e-3, hole_dydz +2.5e-3, 5));              
            
            if (hole_du - du_cut_width/2. > (test_dydz?fDydz_max:fDxdz_max)) continue; 
            if (hole_du + du_cut_width/2. < (test_dydz?fDydz_min:fDxdz_min)) continue; 

            //our TF1 which we will fit our hole-peak with. 
            auto fit = unique_ptr<TF1>(new TF1(
                "holefit", 
                fcn_gauss_offset, 
                //take a close look here at the fit-ranges of our function. if it's a big-hole, we extend the fit range by 25%. 
                max<double>( du_min, hole_du - (hole.is_big ? 1.25 : 1.00)*du_cut_width/2. ), 
                min<double>( du_max, hole_du + (hole.is_big ? 1.25 : 1.00)*du_cut_width/2. ), 
                3
            ));              

            double bin_center[] = {hole_du}; 

            fit->SetParameter(0, maxval - background_fcn(bin_center,nullptr)); 
            fit->SetParameter(1, hole_du);
            fit->SetParameter(2, hole_angle_width / 2.5);

            auto fitresult = hist_slice->Fit("holefit", "S R B Q N L"); 

            //if the fit failed, then skip. 
            if (!fitresult.Get() || !fitresult->IsValid()) continue; 

            //error of the dy/dz position of the hole
            double       hole_du_fit        = fitresult->Parameter(1);
            const double hole_du_staterr    = fitresult->ParError(1);
            const double hole_sigma_fit     = fabs(fitresult->Parameter(2)); 
            const double hole_amplitude_fit = fitresult->Parameter(0); 

            //check to see if this fit is reasonable. 
            //if the height of the gaussian peak is less than 5% of the histogram max, then discard it
            if (hole_amplitude_fit < hist_slice->GetMaximum() * 0.05) continue; 

            //if the relative error of 'sigma' is more than 5%, then discard it. 
            if (fabs(fitresult->ParError(2) / hole_sigma_fit) > 0.05 ) continue; 

            //jnhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhqw ,8
            // -Muon's comment (25 Jan 26)


            const double du_bin = (x_axis->GetXmax() - x_axis->GetXmin()) / ((double)x_axis->GetNbins() - 1); 

            //approximate the number of signal events for this hole, using the parameters of the gaussian fit
            //the const number out front is sqrt(2*pi)
            //
            double hole_n_events = 2.50662827463 * hole_sigma_fit * hole_amplitude_fit * du_bin; 


            //now, let's try to measure the slope 
            SlopeFit_t slope; 
            if (measure_slope && test_dydz) {

                //slope measured for this hole
                slope = MeasureSlope(
                    h_dx_dy, 
                    center_dv, 
                    hole_du, 
                    { center_dv - dv_cut_width/2., center_dv + dv_cut_width/2. }, 
                    { hole_du   - du_cut_width/2., hole_du   + du_cut_width/2. }, 
                    hole_sigma_fit,
                    background_fcn, 
                    c_2d
                );

                //don't record this hole if its slope measurement failed 
                if (is_nan(slope.m)) continue; 
            
                //fix the 'hole_du_fit' variable to use the offset measured in the slope-fitting algorithm 
                hole_du_fit = slope.b + slope.m*center_dv; 

            }
            //cout << "slope: " << slope.m << endl; //" +/- " << slope.m_err << endl; 

            //we're counting this hole as having been 'measured'
            n_holes_measured++;

            fit->SetLineColor(kRed); 
            
            if (fDo_drawing) {
                c_profiles->cd(icanv_profile); fit->DrawCopy("SAME");

                TLine *line; 
                line = new TLine(hole_du_fit, 0., hole_du_fit, hist_slice->GetMaximum()); 
                line->SetLineColor(kBlack); 
                line->Draw(); 

                //draw a line of where the hole SHOULD BE 
                line = new TLine(hole_du, 0., hole_du, hist_slice->GetMaximum()); 
                line->SetLineColor(kRed); 
                line->Draw(); 

                double x0 = test_dydz ? center_dv : hole_du_fit; 
                double y0 = test_dydz ? hole_du_fit : center_dv;  

                //draw the slope-line
                if (measure_slope && (is_nan(slope.m)==false) && (fStateFlag & kDraw_slopes)) {
                    
                    auto slope_line = new TF1("slope_line", 
                        [slope, x0,y0](double *X, double *par)
                        {
                            return slope.m*X[0] + slope.b;
                        }, 
                        angle.dxdz - dxdz_cut_width/2., 
                        angle.dxdz + dxdz_cut_width/2., 
                        0
                    ); 

                    slope_line->SetLineWidth(2); 
                    slope_line->SetLineColor(kBlack);
                        
                    c_2d->cd();
                    slope_line->DrawCopy("SAME"); 
                }
                
                //draw the box that represents the cut used 
                if (fStateFlag & kDraw_boxes) {
                    auto box = new TBox(
                        angle.dxdz - dxdz_cut_width/2., 
                        angle.dydz - (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.,
                        angle.dxdz + dxdz_cut_width/2., 
                        angle.dydz + (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.
                    );
                    //make the box have no fill (transparent), and draw it.

                    c_2d->cd();
                    box->SetFillStyle(0); 
                    box->Draw(); 
                }
                
                //draw the lines that represent the ideal centroid & measured centroid
                if (fStateFlag & kDraw_centroid_lines) {
                    
                    auto line = new TLine; 
                    c_2d->cd();

                    //this is the thin black line that is the centroid's ideal position 
                    if (test_dydz) {
                        line->SetLineColor(kBlack);
                        line->SetLineWidth(1); 
                        line->DrawLine(
                            angle.dxdz - dxdz_cut_width/2., 
                            hole_du,
                            angle.dxdz + dxdz_cut_width/2., 
                            hole_du
                        );

                    } else {
                        //this is the thin black line that is the centroid's ideal position 
                        line->SetLineColor(kBlack);
                        line->SetLineWidth(1); 
                        line->DrawLine(
                            hole_du, 
                            angle.dydz - (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.,
                            hole_du,
                            angle.dydz + (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.
                        );

                        //this is the thick black line that's the measured hole's position 
                        line->SetLineColor(kBlack);
                        line->SetLineWidth(2); 
                        line->DrawLine(
                            hole_du_fit, 
                            angle.dydz - (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.,
                            hole_du_fit,
                            angle.dydz + (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.
                        );
                    }//if (test_dydz)

                }//if (fStateFlag & kDraw_centroid_lines)  -- drawing centroid lines
            }//if (fDo_drawing)

            //add this result to the fitresult 
            fit_result.fits.push_back({
                .hole           = hole, 
                .angle_real     = hole_du, 
                .angle_fit      = hole_du_fit, 
                .angle_sigma    = hole_sigma_fit,
                .angle_slope    = measure_slope ? slope.m : numeric_limits<double>::quiet_NaN(),
                .amplitude      = hole_n_events, 
                .angle_fit_staterr = hole_du_staterr
            });
        }

        //if (!fQuiet) printf("Done with row %2i/%i\n", icanv_profile, v1-v0+1); cout << flush; 

        icanv_profile--; 
    }
    //this will actually store the errors computed

    double pos_offset_RMS             = 0.; 
    double avg_peak_width             = 0.; 
    double avg_pos_offset_uncertainty = 0.; 

    double max_pos_err =0.;
    for (auto holefit : fit_result.fits) {
        double err = fabs(holefit.angle_fit - holefit.angle_real); 
        if (err > max_pos_err) max_pos_err = err; 
    }

    auto h_pos_err_dist = new TH1D(Form("h_pos_dist_%s",test_name), Form("Distribution of hole RMS values (%s);hole RMS (mrad);",test_dydz?"dy/dz":"dx/dz"), 35, 0., 1.1); 

    //spacing bewtween sieve holes in the dx/dz (theta) direction
    const double dTheta = sieve_hole_spacing_x / vtx_z; 

    for (auto holefit : fit_result.fits) {

        if (test_dydz && measure_slope) { 
            //if we're testing phi (dy/dz), we must account for hole-slopes
            double phi0 = holefit.angle_fit - holefit.angle_real; 
            double m    = holefit.angle_slope; 

            //this is the result of integrating the square-error over the range [-dTheta/2, +dTheta/2]. 
            //where the error of any region is [phi0 + m*dTheta]
            double hole_square_error = (pow(phi0 + m*dTheta/2., 3) - pow(phi0 - m*dTheta/2., 3))/(3.*m*dTheta);

            pos_offset_RMS += hole_square_error;

            h_pos_err_dist->Fill(1.e3*sqrt(hole_square_error)); 
        } else {

            //otherwise, just compute the RMS of each hole's centroid
            double hole_square_error = pow(holefit.angle_fit - holefit.angle_real, 2); 
            pos_offset_RMS += hole_square_error; 

            h_pos_err_dist->Fill(1e3*sqrt(hole_square_error)); 
        }
        avg_peak_width  += holefit.angle_sigma; 
        avg_pos_offset_uncertainty += holefit.angle_fit_staterr;   
        
    }

    if (fDo_drawing) {
        new TCanvas(Form("h_pos_err_dist_%s",test_dydz?"dydz":"dxdz"), "dist. of hole-position error"); 
        h_pos_err_dist->SetStats(1);
        h_pos_err_dist->Draw(); 
    }

    const double n_holes_measured_d = (double)n_holes_measured; 

    pos_offset_RMS = sqrt( pos_offset_RMS / n_holes_measured_d );
    avg_peak_width *= 1. / n_holes_measured_d;
    avg_pos_offset_uncertainty *= 1. / n_holes_measured_d; 
    //fit_result.sigma_dydz_position = sqrt( pos_error / amplitude_total ); 
    //fit_result.sigma_dydz_smearing = smear_error / (double)n_holes_measured; 

    if (!fQuiet) {
        printf(
            "%s:\n"
            " ~~ %s %.6f (mrad)\n"
            " ~~ Mean peak width: ................. %.6f (mrad)\n",
            
            (test_dydz?"dy/dz_sv {Phi-tg}":"dx/dz_sv {Theta-tg}"),
            (test_dydz?"centroid offset + slope RMS: .....":"centroid offset RMS: ............."),
            pos_offset_RMS * 1e3,
            avg_peak_width * 1e3
        );
    }
    
    fit_result.RMS_position = pos_offset_RMS; 
    fit_result.RMS_smearing = avg_peak_width; 

    return fit_result; 
}   

//____________________________________________________________________________________________________________________________
function<double(double*,double*)> 
        AngleRecoTester::FitBackground(TH1D* hist, const vector<array<double,2>>& exclusions, const int order)
{
    //lets get a vector of all bins / positions
    struct HistPoint_t { double x,y,err; };

    auto x_ax = hist->GetXaxis(); 
    vector<HistPoint_t> bin_points; 

    HistPoint_t max_point{0., -1.e30, 0.}; 

    //for
    for (int bin=1; bin<=x_ax->GetNbins(); bin++) {

        double bin_x = x_ax->GetBinCenter(bin); 

        //check to see if this bin is inside an 'exclusion zone' 
        bool exclude=false; 
        for (const auto& zone : exclusions) if ( bin_x > zone[0] && bin_x < zone[1] ) { exclude=true; break; } 
        if (exclude) continue; 

        HistPoint_t new_point{bin_x, hist->GetBinContent(bin), hist->GetBinError(bin)};

        //so, this point is not excluded. let's add it to the list of points to fit.    
        bin_points.push_back(new_point); 

        if (new_point.y > max_point.y) max_point = new_point; 
    }

    //now, let's do a chi-square fit first. 
    RVec<double> B(order+1, 0.); 
    RMatrix A(order+1,order+1, 0.); 

    for (const auto& pt : bin_points) {

        RVec<double> Xmu{1.}; 
        for (int i=1; i<=order; i++) Xmu.push_back( Xmu.back() * pt.x ); 

        for (int i=0; i<=order; i++) {

            B[i] += pt.y * Xmu[i] / max<double>(pt.err * pt.err,1.); 
            
            for (int j=0; j<=order; j++) A.get(i,j) += Xmu[i] * Xmu[j] / max<double>(pt.err * pt.err,1.); ; 
        }
    }

    auto&& coeffs = A.Solve(B); 

    auto poly_fcn = [coeffs](double *x, double *par){
        double val=coeffs.back(); 
        for (int i=coeffs.size()-2; i>=0; i--) val = coeffs[i] + (x[0]*val); 
        return val; 
    };
    //auto tf1 = new TF1("poly_test", poly_fcn, x_ax->GetXmin(), x_ax->GetXmax(), 0); 

    return poly_fcn; 

    double normalization = coeffs.back(); 
    for (int i=coeffs.size()-2; i>=0; i--) normalization = coeffs[i] + (max_point.x*normalization); 

    normalization = max_point.y/exp(normalization); 

    //_________________________________________________________________________________________________
    auto result_fcn = [coeffs,normalization](double* xptr, double* par) {
        
        //use homer's method to compute the argument of the exponential 

        double arg = coeffs.back();
        for (int i=coeffs.size()-2; i>=0; i--) arg = coeffs[i] + (xptr[0]*arg); 
        
        return normalization * exp( arg ); 
    };
    //_________________________________________________________________________________________________

    return result_fcn; 

}
//____________________________________________________________________________________________________________________________
AngleRecoTester::SlopeFit_t AngleRecoTester::MeasureSlope(
    TH2D* h_data, 
    const double x0, 
    const double y0, 
    const std::array<double,2>& xlim, 
    const std::array<double,2>& ylim, 
    const double hole_sigma,
    const std::function<double(double*,double*)>& bg_fcn, 
    TVirtualPad *pad_2d
) 
{
    auto xax = h_data->GetXaxis(); 
    auto yax = h_data->GetYaxis(); 

    const int bins_x[] = { xax->FindBin(xlim[0]), xax->FindBin(xlim[1]) }; 
    const int bins_y[] = { yax->FindBin(ylim[0]), yax->FindBin(ylim[1]) }; 

    //we have to re-nomralized the background function, as we're looking at a single bin-width
    const double bg_normalization = ( (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()-1) ) / (xlim[1] - xlim[0]); 

    struct GausFitPoints_t { double x,y; };   

    vector<GausFitPoints_t> slope_fitpoints; 

    //find the maximum point in a histogram's subrange 
    double max_height_diff=-1.; 
    for (int bx=bins_x[0]; bx<=bins_x[1]; bx++) {
        for (int by=bins_y[0]; by<=bins_y[1]; by++) { 

            double y_ptr[] = { yax->GetBinCenter(by) }; 
            double height_diff = h_data->GetBinContent(bx,by) - bg_normalization*bg_fcn(y_ptr, nullptr); 

            if (height_diff > max_height_diff) max_height_diff = height_diff; 
        }
    }

    //now, we go thru each bin, and fit a gaussian to the data. 
    for (int bx=bins_x[0]; bx<=bins_x[1]; bx++) {

        //now, we fit a gaussian to this 'slice' 

        //the y-pos of the bins
        vector<double> Y;
        //the number of hits in the bin
        vector<double> Z;
        //the difference between the stats of this bin, and the background expcetation 
        vector<double> Z_diff; 

        for (int by=bins_y[0]; by<=bins_y[1]; by++) {
            
            Y.push_back( yax->GetBinCenter(by) ); 
            
            double bin_content = h_data->GetBinContent(bx,by); 

            Z.push_back(bin_content); 
            
            //now, we will find the difference between the entries in this bin, and the expectation from the 
            // background. 

            //this is necessary, because our background function only accepts pointers for its argument
            double y_ptr[] = {Y.back()}; 

            //now, this is the difference between the background expectation, and the actual bin content. 
            bin_content += -bg_normalization*bg_fcn(y_ptr, nullptr); 

            Z_diff.push_back(bin_content); 
        }

        //now, we're ready to fit a gaussian to this bin.
        
        double y_ptr[] { yax->GetBinCenter(yax->FindBin(y0)) }; 

        //get the height of the histogram at the 'central bin' (this will be our first guess for the height of the gaussian)
        double A_first_guess = h_data->GetBinContent( bx, yax->FindBin(y0) ) - bg_normalization*bg_fcn(y_ptr, nullptr); 

        //height of gaussian above bg
        double A{A_first_guess};
        //center of gaussian 
        double y_cent{y0}; 

        //now, let's iterate. 
        const int npts = Y.size();
        const double sig2 = hole_sigma*hole_sigma;  
        
        //____________________________________________________________________________________________________________________________________________
        //this helper function tries to optimize the NLL fit w/r/t 'a' and 'ybar', and returns the NLL. 
        //returns NaN if the iteration fails 
        auto Newton_iterate = [h_data, bx, &Y, &Z, &Z_diff, npts, sig2](double& a, double &ybar) 
        {
            double dEps_dA{0.};     // \partial \epsilon^2 / \partial A 
            double d2Eps_dA2{0.};   // \partial^2 \epsilon^2 / \partial A^2
            double dEps_dy{0.};     // \partial \epsilon / \partial ybar
            double d2Eps_dydA{0.};  // \partial^2 \epsilon^2 / \partial ybar \partial A
            double d2Eps_dy2{0.};   // \partial^2 \epsilon^2 / \partial^2 ybar

            for (int i=0; i<=npts; i++) {

                double z_i    = Z[i];
                double zdif_i = Z_diff[i]; 
                double dy_i   = Y[i] - ybar; 
            
                double eta_i = exp( -0.5 * dy_i*dy_i/sig2 ); 
                double eps_i = a*eta_i - zdif_i;

                //first derivs
                dEps_dA += eps_i*eta_i;
     
                dEps_dy += a * eps_i*eta_i * dy_i/sig2; 

                //second derivs
                d2Eps_dA2  += eta_i*eta_i; 

                d2Eps_dydA += ( a*eta_i + eps_i )*eta_i * dy_i/sig2; 

                d2Eps_dy2  += ( a*eta_i + eps_i ) * (a*eta_i) * (dy_i/sig2)*(dy_i/sig2)  -  (a/sig2) * (eps_i*eta_i); 
            }

            //now, invert the following matrix: 
            /*
            *       [ a   b ]^-1  =  [  d  -b ] 
            *       [ b   d ]        [ -b   a ]/(ad - bb)
            *
            * a = d2Eps_dA2; 
            * b = d2Eps_dydA; 
            * d = d2Eps_dy2; 
            */

            double detA = (d2Eps_dA2*d2Eps_dy2) - (d2Eps_dydA*d2Eps_dydA); 

            //check if the matrix is singular
            if (is_nan(detA) || fabs(detA) < singular_determinant_cutoff) return NaN_double;  

            a    += - ( d2Eps_dy2 *dEps_dA  - d2Eps_dydA*dEps_dy)/detA; 
            ybar += - (-d2Eps_dydA*dEps_dA  + d2Eps_dA2 *dEps_dy)/detA; 

            //now, re-measure the chi2. 
            double chi2=0.; 
            for (int i=0; i<=npts; i++) {

                double zdif_i = Z_diff[i]; 
                double dy_i   = Y[i] - ybar;

                chi2 += pow(zdif_i - a*exp( -0.5 * dy_i*dy_i/sig2 ), 2); 
            }

            //printf(" chi2: %.3e     A: %+.1e  y-bar: %+.1e\n", chi2/2., a, ybar); std::cout << std::flush; 
            
            return chi2/2.; 
            
            
        };
        //____________________________________________________________________________________________________________________________________________
        
        //now, we do our iterations. break when certain conditions are met.
        // starting iterations 
        //printf("starting iterations ~~~~~~~~~~~~~~~~~~~~~\n"); std::cout << std::flush;  
        int i_it=0;
        double chi2;  
        do {

            //printf("it %3i ", i_it); std::cout << std::flush; 
            chi2 = Newton_iterate(A, y_cent); 
            
        } while ( (is_nan(chi2)==false) && ++i_it < max_newton_iterations ); 
        //std::cout << "\n"; 

        //now, do some basic checks

        

        //if the center we found is out of range, throw it out.
        if (y_cent < ylim[0] || y_cent > ylim[1]) continue; 

        //printf("ratio of A to max height for this hole: %.4f\n", A / max_height_diff); std::cout << std::flush; 

        //if the max height found is too low, throw it out. 
        if (A / max_height_diff < min_height_ratio) continue; 
         
        double x_bin = xax->GetBinCenter(bx); 

        slope_fitpoints.push_back({ .x = x_bin, .y = y_cent }); 

        //printf("%+.3f, %+.3f\n", (x_bin-x0)*1.e3, (y_cent-y0)*1.e3); std::cout << std::flush; 
    }   

    //printf(" ------------------------- number of fit-point for this hole: %zi\n", slope_fitpoints.size()); 

    //now, find the slope and offset. 
    //draw the boxes, if that's relevant. 
    if (fDo_drawing && (fStateFlag & kDraw_slope_points)) {

        const double bin_halfwidth_x = 0.5*(xax->GetXmax()-xax->GetXmin())/((double)xax->GetNbins()-1.);         
        const double bin_halfwidth_y = 0.5*(yax->GetXmax()-yax->GetXmin())/((double)yax->GetNbins()-1.); 

        pad_2d->cd(); 
        for (const auto& pt : slope_fitpoints) {
            
            auto box = new TBox(
                pt.x - bin_halfwidth_x, 
                pt.y - bin_halfwidth_y, 
    
                pt.x + bin_halfwidth_x, 
                pt.y + bin_halfwidth_y
            ); 
            box->SetFillStyle(0); 
            box->Draw();
        }
    }
    
    //min number of points to attempt fitting the slope 
    if (slope_fitpoints.size() < 3) return { NaN_double }; 
    
    //this just your everyday least-squares for mx + b = y
    double sum_x{0.}, sum_y{0.}, sum_xx{0.}, sum_xy{0.}; 
    for (const auto& pt : slope_fitpoints) {
        sum_x += pt.x; 
        sum_y += pt.y; 
        sum_xx += pt.x*pt.x; 
        sum_xy += pt.x*pt.y; 
    } 
    double N = (double)slope_fitpoints.size(); 
    
    double m = ( N*sum_xy - sum_x*sum_y ) / ( N*sum_xx - sum_x*sum_x ); 
    double b = (sum_y - m*sum_x)/N; 

    return { m, 0., b }; 

};

//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________

ClassImp(AngleRecoTester); 