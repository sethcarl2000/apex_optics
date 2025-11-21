#ifndef TestHorizontalAngle_H
#define TestHorizontalAngle_H

//C++ std-lib headers
#include <iostream>
#include <sstream>
#include <string> 
#include <memory>
#include <limits> 
#include <optional> 
#include <functional> 
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
#include <SieveHole.h>
#include <TLine.h> 
#include <TAxis.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
//APEX-otpics headers
#include <ApexOptics.h> 
#include <RMatrix.h>

//this stores the results of the fit 
struct AngleFit_t {

    //the sieve hole associated with this fit
    SieveHole hole;
    //the position that this angle should be, based on target used and sieve-hole position 
    double angle_real;
    //the position of the peak, from the gaussian fit
    double angle_fit; 
    //the sigma of the gaussian fit
    double angle_sigma; 
    //the amplitude of the gaussian fit (with background subtracted). 
    double amplitude; 
}; 

//
struct AngleFitResult_t {
    
    double RMS_position;
    double RMS_smearing; 
    double RMS_overall; 

    std::vector<AngleFit_t> fits{}; 
}; 

std::function<double(double*,double*)> FitBackground(TH1D* hist, const vector<pair<double,double>>& exclusions, const int order); 

//specify which angle to test
constexpr bool kTestAngle_dxdz = false; 
constexpr bool kTestAngle_dydz = true;

//test the horizontal angle reconstruction efficiency
std::optional<AngleFitResult_t> TestAngleReco(  const bool is_RHRS,         //arm to use
                                                const bool test_dydz,       //which angle to test
                                                ROOT::RDF::RNode df,        //RDataFrame to get our data from 
                                                ApexOptics::OpticsTarget_t target,      //optics target to use  
                                                const int first_row,        //first row to start at
                                                const int last_row,         //Last row to end at       
                                                const int first_col,        //first column to start at
                                                const int last_col,         //Last column to end at
                                                const double row_cut_width=1.2,     //width to make vertical cut, expressed as a fraction of the row spacing 
                                                const double col_cut_width=0.80,    //width to make cut around hole, in fraction of inter-column spacing
                                                const double bg_cut_width=1.00,     //fraction of row_cut_width to use to fit the background with
                                                const char* name_dxdz="reco_dxdz_sv", 
                                                const char* name_dydz="reco_dydz_sv", 
                                                const double dxdz_min = -0.04,
                                                const double dxdz_max = +0.04, 
                                                const double dydz_min = -0.04, 
                                                const double dydz_max = +0.04, 
                                                const char* name_react_vtx="position_vtx" )
{
    using namespace std; 
    using ApexOptics::Trajectory_t; 
    using ApexOptics::OpticsTarget_t; 

    const char* const here = "TestHorizontalAngle"; 

    char c_title[255]; 
    
    //check that the target given is a V-wire
    if (!(target.name == "V1" || 
          target.name == "V2" || 
          target.name == "V3"))  {
        
        Error(here, "Target provided is not a V-wire target: %s", target.name.c_str()); 
        return std::nullopt; 
    }

    //get the react vertex
    const double vtx_y_hcs = *df.Define("y", [](TVector3 v){ return v.y(); }, {name_react_vtx}).Mean("y"); 

    const TVector3 react_vtx_scs = ApexOptics::HCS_to_SCS( is_RHRS, TVector3(target.x_hcs, vtx_y_hcs, target.z_hcs) ); 

    const double vtx_u = test_dydz ? react_vtx_scs.y() : react_vtx_scs.x(); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 
    gStyle->SetOptStat(0); 

    //now, we can compute where each sieve-hole SHOULD be
    vector<SieveHole> sieve_holes = ApexOptics::ConstructSieveHoles(is_RHRS); 

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
    const double sieve_thickness_z = 12e-3; 
    const double sieve_hole_spacing_y = 4.826e-3; 
    const double sieve_hole_spacing_x = 5.842e-3; 
    
    //in radians
    const double dxdz_cut_width = (sieve_hole_spacing_x * row_cut_width) / fabs(react_vtx_scs.z()); 
    const double dydz_cut_width = (sieve_hole_spacing_y * col_cut_width) / fabs(react_vtx_scs.z()); 

    // 'u' is the direction we're evaluating accuracy in, and 'v' is the other direction we're not fitting in. 
    const double du_cut_width = test_dydz ? dydz_cut_width : dxdz_cut_width; 
    const double dv_cut_width = test_dydz ? dxdz_cut_width : dydz_cut_width; 
    
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0); 

    //Stores a single dx/dz and dy/dz data point
    struct AnglePair_t { double dxdz,dydz; }; 

    //make a vector of all events 
    std::cout << "Making copy of data..." << flush; 

    vector<AnglePair_t> data = *df 

        .Define("event", [](double dxdz, double dydz)
        {
            return AnglePair_t{dxdz,dydz}; 
        }, {name_dxdz, name_dydz})
        
        .Take<AnglePair_t>("event");  

    std::cout << "done." << endl; 

    auto c = new TCanvas("ctest_dy", "Angle Reco - dy/dz", 1600, 800); 
    c->Divide(2,1);

    auto c_2d         = c->cd(1);
    auto c_profiles   = c->cd(2);

    //create, fill, and draw the 2d-hist
    auto h_dx_dy = new TH2D("h", "dx/dz_{sv} vs dy/dz_{sv} (rad)", 200, dxdz_min, dxdz_max, 200, dydz_min, dydz_max); 
    
    for (const auto& ev : data) {
        h_dx_dy->Fill( ev.dxdz, ev.dydz ); 
    }
    c_2d->cd(); h_dx_dy->DrawCopy("col"); 
    
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

    const double du_max = test_dydz ? dydz_max : dxdz_max; 
    const double du_min = test_dydz ? dydz_min : dxdz_min;
    
    const double dv_max = test_dydz ? dxdz_max : dydz_max; 
    const double dv_min = test_dydz ? dxdz_min : dydz_min;

    //this allows us to generalize some of the code below. 
    double AnglePair_t::*du_ap = test_dydz ? &AnglePair_t::dydz : &AnglePair_t::dxdz; 
    double AnglePair_t::*dv_ap = test_dydz ? &AnglePair_t::dxdz : &AnglePair_t::dydz; 

    double SieveHole::*u_sh = test_dydz ? &SieveHole::y : &SieveHole::x; 
    double SieveHole::*v_sh = test_dydz ? &SieveHole::x : &SieveHole::y;

    double Trajectory_t::*du_trj = test_dydz ? &Trajectory_t::dydz : &Trajectory_t::dxdz; 
    double Trajectory_t::*dv_trj = test_dydz ? &Trajectory_t::dxdz : &Trajectory_t::dydz; 
    
    //each one of these vectors is a different slice to fit
    vector<vector<SieveHole>> hole_slices; 
    for (int iv = v0; iv <= v1; iv++) {

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
    
    c_profiles->Divide( 1, hole_slices.size(), 0.01,0.); 

    // Fit all holes - dydz
    int icanv_profile=hole_slices.size(); 
    for (const auto& holes : hole_slices) {

        vector<Trajectory_t> hole_angles;
        for (auto& hole : holes) { 
            hole_angles.emplace_back(Trajectory_t{
                hole.x, 
                hole.y,
                ( hole.x - react_vtx_scs.x() )/( 0. - react_vtx_scs.z() ),
                ( hole.y - react_vtx_scs.y() )/( 0. - react_vtx_scs.z() ) 
            }); 
        }
                
        // get the dudz-value of the row. for fitting dy/dz, all holes have the same 'dv'. 
        const double center_dv = hole_angles.front().*dv_trj; 
        // this is needed for testing dx/dz 

        // create and fill the row histogram
        auto hist_slice = new TH1D(Form("h_slice_%i",icanv_profile), "", 200, du_min, du_max); 
        
        auto hist_slice_bg = (TH1D*)hist_slice->Clone("h_row_bg"); 

        //Backstory - I made a terrible decision a few years ago when I decided that sieve-holes on adjacent, staggered rows 
        // ought to have the same column-index, and so we need to enage in this funny business of fitting TWO histogram slices per row, 
        // but ONLY if we're cheking dx/dz (because all holes in the same column )
        for (const auto& ev : data) {
            
            //all holes are in a straight line
            if (fabs(ev.*dv_ap - center_dv) < dv_cut_width/2.)              hist_slice->Fill( ev.*du_ap ); 
            if (fabs(ev.*dv_ap - center_dv) < bg_cut_width*dv_cut_width/2.) hist_slice_bg->Fill( ev.*du_ap ); 
        }

        printf("Drawing slice %i...",icanv_profile); cout << flush; 
        
        
        hist_slice->GetYaxis()->SetNdivisions(0); 
        c_profiles->cd(icanv_profile); hist_slice->DrawCopy(); 
        auto x_axis = hist_slice->GetXaxis(); 

        
        // now, try to fit a polynomial to the background.
        // do do this, we must make a list of 'excluded' regions, where we know the holes ought to be. 
        vector<pair<double,double>> exclusions;
        for (int i=0; i<hole_angles.size(); i++) {

            exclusions.push_back({
                hole_angles[i].*du_trj - du_cut_width/2., 
                hole_angles[i].*du_trj + du_cut_width/2.
            });
        }
        cout << "Fitting background..." << flush; 

        //fit the background with a polynomial, 
        auto background_fcn_unnormalized = FitBackground(hist_slice_bg, exclusions, background_polynomial_order);
        auto background_fcn = [background_fcn_unnormalized, bg_cut_width](double *x, double *par){ return background_fcn_unnormalized(x,par)/bg_cut_width; };
        auto tf1_bg = new TF1(Form("bg_poly_%i",icanv_profile), background_fcn, x_axis->GetXmin(), x_axis->GetXmax(), 0); 
        delete hist_slice_bg; 
        tf1_bg->SetLineStyle(kDashed);
        tf1_bg->SetLineColor(kRed); 
        c_profiles->cd(icanv_profile); tf1_bg->DrawCopy("SAME"); 
        

        cout << "done." << endl; 

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
            const double du_hi = ( hole.*u_sh - vtx_u ) / ( sieve_thickness_z - react_vtx_scs.z() ); 

            const double du_lo = ( hole.*u_sh - vtx_u ) / ( 0.                - react_vtx_scs.z() ); 

            const double hole_radius = hole.radius_front / fabs(react_vtx_scs.z());
            const double hole_du     = (du_hi + du_lo) / 2.;
                
            //auto fit = unique_ptr<TF1>(new TF1("holefit", gaussian_semicircle, hole_dydz -2.5e-3, hole_dydz +2.5e-3, 5));              
            
            if (hole_du - du_cut_width/2. > (test_dydz?dydz_max:dxdz_max)) continue; 
            if (hole_du + du_cut_width/2. < (test_dydz?dydz_min:dxdz_min)) continue; 

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
            fit->SetParameter(2, hole_radius / 2.5);

            auto fitresult = hist_slice->Fit("holefit", "S R B Q N L"); 

            //if the fit failed, then skip. 
            if (!fitresult.Get() || !fitresult->IsValid()) continue; 

            //error of the dy/dz position of the hole
            const double hole_du_fit        = fitresult->Parameter(1);
            const double hole_sigma_fit     = fabs(fitresult->Parameter(2)); 
            const double hole_amplitude_fit = fitresult->Parameter(0); 

            //check to see if this fit is reasonable. 
            //if the height of the gaussian peak is less than 5% of the histogram max, then discard it
            if (hole_amplitude_fit < hist_slice->GetMaximum() * 0.05) continue; 

            //if the relative error of 'sigma' is more than 5%, then discard it. 
            if (fabs(fitresult->ParError(2) / hole_sigma_fit) > 0.05 ) continue; 


            const double du_bin = (x_axis->GetXmax() - x_axis->GetXmin()) / ((double)x_axis->GetNbins() - 1); 

            //approximate the number of signal events for this hole, using the parameters of the gaussian fit
            //the const number out front is sqrt(2*pi)
            //
            double hole_n_events = 2.50662827463 * hole_sigma_fit * hole_amplitude_fit * du_bin; 

            //we're counting this hole as having been 'measured'
            n_holes_measured++;


            fit->SetLineColor(kRed); 
            c_profiles->cd(icanv_profile); fit->DrawCopy("SAME");  

            //draw a line of where the hole IS 
            TLine *line; 
            line = new TLine(hole_du_fit, 0., hole_du_fit, hist_slice->GetMaximum()); 
            line->SetLineColor(kBlack); 
            line->Draw(); 

            //draw a line of where the hole SHOULD BE 
            line = new TLine(hole_du, 0., hole_du, hist_slice->GetMaximum()); 
            line->SetLineColor(kRed); 
            line->Draw(); 

            //draw the box that represents the cut used 
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

            //add this result to the fitresult 
            fit_result.fits.push_back({
                .hole           = hole, 
                .angle_real     = angle.*du_trj, 
                .angle_fit      = hole_du_fit, 
                .angle_sigma    = hole_sigma_fit,
                .amplitude      = hole_n_events
            });
        }

        printf("Done with row %2i/%i\n", icanv_profile, v1-v0+1); cout << flush; 

        icanv_profile--; 
    }
    //this will actually store the errors computed

    double amplitude_total  = 0.; 
    double pos_error        = 0.; 
    double smear_error      = 0.; 

    for (auto holefit : fit_result.fits) {

        pos_error       += pow(holefit.angle_fit - holefit.angle_real, 2) * holefit.amplitude; 
        amplitude_total += holefit.amplitude; 

        smear_error     += holefit.angle_sigma; 
    }

    pos_error = sqrt( pos_error / amplitude_total );
    smear_error = smear_error / (double)n_holes_measured; 
    //fit_result.sigma_dydz_position = sqrt( pos_error / amplitude_total ); 
    //fit_result.sigma_dydz_smearing = smear_error / (double)n_holes_measured; 

    printf(
        "%s:\n"
        " ~~ position error: %.4e\n"
        " ~~ smearing error: %.4e\n",

        (test_dydz?"dy/dz_sv {Phi-tg}":"dx/dz_sv {Theta-tg}"),
        pos_error,
        smear_error
    );

    fit_result.RMS_position = pos_error; 
    fit_result.RMS_smearing = smear_error; 

    return fit_result; 
}

std::function<double(double*,double*)> FitBackground(
    TH1D* hist, 
    const vector<pair<double,double>>& exclusions, 
    const int order) 
{
    using namespace std; 
    using namespace ROOT::VecOps; 

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
        for (const auto& zone : exclusions) if ( bin_x > zone.first && bin_x < zone.second ) { exclude=true; break; } 
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

#endif 
