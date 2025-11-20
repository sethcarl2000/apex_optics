#ifndef TestHorizontalAngle_H
#define TestHorizontalAngle_H

#include <iostream>
#include <sstream>
#include <string> 
#include <TParameter.h>
#include <ApexOptics.h> 
#include <TVector3.h> 
#include <TSystem.h> 
#include <TStyle.h> 
#include <TColor.h> 
#include <TCanvas.h>
//this last one should eventually be moved to the ApexOptics namespace, as it has general use outside of the 'isolate_sieveholes' app. 
#include <SieveHole.h>
#include <TLine.h> 
#include <TAxis.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
#include <functional> 
#include <memory>
#include <TBox.h> 
#include <optional> 
#include <Math/Factory.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <RMatrix.h> 
#include <ROOT/RVec.hxx>
#include <limits> 

using ApexOptics::OpticsTarget_t; 
using ApexOptics::Trajectory_t; 

//this stores the results of the fit 
struct AngleFit_t {

    SieveHole hole; 
    double angle_real;
    double angle_fit; 
    double angle_sigma; 
}; 

//
struct AngleFitResult_t {
    
    double sigma_dxdz_position;
    double sigma_dxdz_smearing; 
    double sigma_dxdz_overall; 

    double sigma_dydz_position;
    double sigma_dydz_smearing; 
    double sigma_dydz_overall; 

    std::vector<AngleFit_t> fits_dxdz{}; 
    std::vector<AngleFit_t> fits_dydz{}; 
}; 

std::function<double(double*,double*)> FitBackground(TH1D* hist, const vector<pair<double,double>>& exclusions, const int order); 

//test the horizontal angle reconstruction efficiency
std::optional<AngleFitResult_t> TestAngleReco(  const bool is_RHRS,         //arm to use
                                                ROOT::RDF::RNode df,        //RDataFrame to get our data from 
                                                OpticsTarget_t target,      //optics target to use  
                                                const int first_row,     //first row to start at
                                                const int last_row,      //Last row to end at       
                                                const int first_col,     //first column to start at
                                                const int last_col,      //Last column to end at
                                                const double row_cut_width=1.2,     //width to make vertical cut, expressed as a fraction of the row spacing 
                                                const double col_cut_width=0.80,    //width to make cut around hole, in fraction of inter-column spacing
                                                const double col_bg_cut_width=1.00, //fraction of col_cut_width to use to fit the background with
                                                const char* name_dxdz="reco_dxdz_sv", 
                                                const char* name_dydz="reco_dydz_sv", 
                                                const double dxdz_min = -0.04,
                                                const double dxdz_max = +0.04, 
                                                const double dydz_min = -0.04, 
                                                const double dydz_max = +0.04, 
                                                const char* name_react_vtx="position_vtx" )
{
    using namespace std; 

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

    auto c_dyTest = new TCanvas("c_holes", "Angle Reco", 1600, 800); 
    c_dyTest->Divide(2,1);

    auto cdyTest_2d         = c_dyTest->cd(1);
    auto cdyTest_profiles   = c_dyTest->cd(2);
    
    //let's give names to each subcanvas: 
    const int canvId_dx_dy = 1; 

    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0); 

    //Stores a single dx/dz and dy/dz data point
    struct AnglePair_t { double dxdz,dydz; }; 

    //make a vector of all events 
    cout << "Making copy of data..." << flush; 

    vector<AnglePair_t> data = *df 

        .Define("event", [](double dxdz, double dydz)
        {
            return AnglePair_t{dxdz,dydz}; 
        }, {name_dxdz, name_dydz})
        
        .Take<AnglePair_t>("event");  

    cout << "done." << endl; 

    
    //create, fill, and draw the 2d-hist
    auto h_dx_dy = new TH2D("h", "dx/dz_{sv} vs dy/dz_{sv} (rad)", 200, dxdz_min, dxdz_max, 200, dydz_min, dydz_max); 
    
    for (const auto& ev : data) {
        h_dx_dy->Fill( ev.dxdz, ev.dydz ); 
    }
    cdyTest_2d->cd(); 
    h_dx_dy->Draw("col"); 


    size_t n_holes_measured =0; 
    
    //here, we're testing vertical wires.     
    cdyTest_profiles->Divide( 1, last_row - first_row + 1, 0.01,0.); 

    int i_canv=1;

    //row counter
    int i_row=1; 

    //order of polynomial to fit the backgroudn with.  
    const int background_polynomial_order = 10; 
    
    AngleFitResult_t fit_result; 

    //the elements of this pair are 'dxdz' and 'dydz'; they represent which holes were fit successfully (and therefore should be drawn). 
    vector<AnglePair_t> holes_to_draw; 

    //loop thru each row 
    for (int row = first_row; row <= last_row; row++) {

        vector<SieveHole> row_holes; 
        for (int col=first_col; col<=last_col; col++) {
            
            auto hole = FindSieveHole(row, col); 
            if (hole.row>=0 && hole.col>=0) row_holes.push_back(hole);                
        }

        vector<Trajectory_t> hole_angles;
        for (auto& hole : row_holes) { 
            hole_angles.emplace_back(Trajectory_t{
                hole.x, 
                hole.y,
                ( hole.x - react_vtx_scs.x() )/( 0. - react_vtx_scs.z() ),
                ( hole.y - react_vtx_scs.y() )/( 0. - react_vtx_scs.z() ) 
            }); 
        }
                
        // get the dxdz-value of the row
        const double row_dxdz = hole_angles.front().dxdz; 
        
        // create and fill the row histogram
        auto hist_row = new TH1D(Form("h_row_%i",row-first_row), "dy/dz_{sv}", 200, dydz_min, dydz_max); 
        auto hist_row_bg = (TH1D*)hist_row->Clone("h_row_bg"); 
        
        for (const auto& ev : data) {
            
            if (fabs(ev.dxdz - row_dxdz) < dxdz_cut_width/2.) hist_row->Fill( ev.dydz ); 
            if (fabs(ev.dxdz - row_dxdz) < col_bg_cut_width*dxdz_cut_width/2.) hist_row_bg->Fill( ev.dydz ); 
        }

        printf("Drawing row %i...",row); cout << flush; 
        
        cdyTest_profiles->cd(i_row); 
        hist_row->GetYaxis()->SetNdivisions(0); 
        hist_row->DrawCopy(); 
        auto x_axis = hist_row->GetXaxis(); 

        
        // now, try to fit a polynomial to the background.
        // do do this, we must make a list of 'excluded' regions, where we know the holes ought to be. 
        vector<pair<double,double>> exclusions; 
        for (const auto& angle : hole_angles) {
            exclusions.push_back({
                angle.dydz - dydz_cut_width/2., 
                angle.dydz + dydz_cut_width/2.
            }); 
        }
        cout << "Fitting background..." << flush; 
        auto background_fcn_unnormalized = FitBackground(hist_row_bg, exclusions, background_polynomial_order);
        delete hist_row_bg; 
        auto background_fcn = [background_fcn_unnormalized, col_bg_cut_width](double *x, double *par){ return background_fcn_unnormalized(x,par)/col_bg_cut_width; };
        auto tf1_bg = new TF1(Form("bg_poly_%i",row-first_row), background_fcn, x_axis->GetXmin(), x_axis->GetXmax(), 0); 
        tf1_bg->SetLineStyle(kDashed); 
        tf1_bg->Draw("SAME"); 

        cout << "done." << endl; 

        int i_col =0; 
        for (size_t i=0; i<row_holes.size(); i++) {

            //the 'SieveHole' struct associated with this hole
            const auto& hole  = row_holes[i]; 

            //the 'Trajectory_t' struct associated with this hole
            const auto& angle = hole_angles[i]; 

            auto fcn_gauss_offset = [background_fcn](double *x, double *par) { 
                double sigma = fabs(par[2]); 
                return par[0] * exp( -0.5 * pow( (x[0] - par[1])/sigma, 2 ) ) + background_fcn(x,nullptr); 
            }; 

            //get the maximum value in this range
            int bin_min = x_axis->FindBin( angle.dydz - dydz_cut_width/2. ); 
            int bin_max = x_axis->FindBin( angle.dydz + dydz_cut_width/2. ); 

            double maxval = -1.; 
            int max_bin = -1; 
            for (int bin=bin_min; bin<=bin_max; bin++) {

                if (hist_row->GetBinContent(bin) > maxval) {
                    maxval = hist_row->GetBinContent(bin); 
                    max_bin = bin; 
                }
            }
            //cout << "Max bin val: " << maxval << endl; 
            double x0 = x_axis->GetBinCenter(max_bin); 

            //double offset = (hist_row->GetBinContent(bin_min) + hist_row->GetBinContent(bin_max))/2.; 

            //upper edge of sieve hole, accounting for paralax, assuming that the tungsten material of the sieve is impenitrable
            // (which is known not to be true!!)
            const double dydz_hi = ( (hole.y + hole.radius_front) - react_vtx_scs.y() ) / ( sieve_thickness_z - react_vtx_scs.z() ); 

            const double dydz_lo = ( (hole.y - hole.radius_front) - react_vtx_scs.y() ) / ( 0. - react_vtx_scs.z() ); 

            const double hole_radius = (dydz_hi - dydz_lo) / 2.;
            const double hole_dydz   = (dydz_hi + dydz_lo) / 2.;
                
            //auto fit = unique_ptr<TF1>(new TF1("holefit", gaussian_semicircle, hole_dydz -2.5e-3, hole_dydz +2.5e-3, 5));              
            
            if (hole_dydz - dydz_cut_width/2. > dydz_max) continue; 
            if (hole_dydz + dydz_cut_width/2. < dydz_min) continue; 

            //our TF1 which we will fit our hole-peak with. 
            auto fit = unique_ptr<TF1>(new TF1(
                "holefit", 
                fcn_gauss_offset, 
                //take a close look here at the fit-ranges of our function. if it's a big-hole, we extend the fit range by 25%. 
                max<double>( dydz_min, hole_dydz - (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2. ), 
                min<double>( dydz_max, hole_dydz + (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2. ), 
                3
            ));              

            double bin_center[] = {hole_dydz}; 

            fit->SetParameter(0, maxval - background_fcn(bin_center,nullptr)); 
            fit->SetParameter(1, hole_dydz);
            fit->SetParameter(2, hole_radius / 2.5);

            auto fitresult = hist_row->Fit("holefit", "S R B Q N L"); 

            //if the fit failed, then skip. 
            if (!fitresult.Get() || !fitresult->IsValid()) continue; 

            //error of the dy/dz position of the hole
            const double hole_dydz_fit  = fitresult->Parameter(1);
            const double hole_sigma_fit = fabs(fitresult->Parameter(2)); 
            const double hole_amplitude_fit = fitresult->Parameter(0); 

            //check to see if this fit is reasonable. 
            //if the height of the gaussian peak is less than 5% of the histogram max, then discard it
            if (hole_amplitude_fit < hist_row->GetMaximum() * 0.05) continue; 

            //if the relative error of 'sigma' is more than 5%, then discard it. 
            if (fabs(fitresult->ParError(2) / hole_sigma_fit) > 0.05 ) continue; 


            const double dx_hist = (x_axis->GetXmax() - x_axis->GetXmin()) / ((double)x_axis->GetNbins() - 1); 

            //approximate the number of signal events for this hole, using the parameters of the gaussian fit
            //the const number out front is sqrt(2*pi)
            //
            double hole_n_events = 2.50662827463 * hole_sigma_fit * hole_amplitude_fit * dx_hist; 

            //we're counting this hole as having been 'measured'
            n_holes_measured++;

            cdyTest_profiles->cd(i_row);

            fit->SetLineColor(kRed); 
            fit->DrawCopy("SAME");  

            //draw a line of where the hole SHOULD BE 
            auto line = new TLine(hole_dydz, 0., hole_dydz, hist_row->GetMaximum()); 
            line->SetLineColor(kBlack); 
            line->Draw(); 

            //draw a line of where the hole SHOULD BE 
            line = new TLine(hole_dydz_fit, 0., hole_dydz_fit, hist_row->GetMaximum()); 
            line->SetLineColor(kRed); 
            line->Draw(); 

            holes_to_draw.push_back({angle.dxdz, angle.dydz});

            //draw the box that represents the cut used 
            auto box = new TBox(
                angle.dxdz - dxdz_cut_width/2., 
                angle.dydz - (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.,
                angle.dxdz + dxdz_cut_width/2., 
                angle.dydz + (hole.is_big ? 1.25 : 1.00)*dydz_cut_width/2.
            );
            //make the box have no fill (transparent), and draw it.
             
            cdyTest_2d->cd();

            box->SetFillStyle(0); 
            box->Draw(); 

            //add this result to the fitresult 
            fit_result.fits_dydz.push_back({
                .hole           = hole, 
                .angle_real     = angle.dydz, 
                .angle_fit      = hole_dydz_fit, 
                .angle_sigma    = hole_sigma_fit
            });
        }

        /*line = new TLine( row_dxdz + dxdz_cut_width/2., -0.04,  row_dxdz + dxdz_cut_width/2., +0.04 ); 
        line->Draw("SAME");

        line = new TLine( row_dxdz - dxdz_cut_width/2., -0.04,  row_dxdz - dxdz_cut_width/2., +0.04 ); 
        line->Draw("SAME");*/ 

        printf("Done with row %2i/%i\n", i_row, last_row - first_row + 1); cout << flush; 

        i_row++; 
    }
    //this will actually store the errors computed
    double dydz_pos_error = 0.; 
    double dydz_smear_error = 0.; 

    for (auto holefit : fit_result.fits_dydz) {

        dydz_pos_error   += pow(holefit.angle_fit - holefit.angle_real, 2); 
        dydz_smear_error += holefit.angle_sigma; 
    }

    fit_result.sigma_dydz_position = sqrt( dydz_pos_error / ((double)n_holes_measured)); 
    fit_result.sigma_dydz_smearing = dydz_smear_error / (double)n_holes_measured; 

    printf("position error: %.4e\n", fit_result.sigma_dydz_position );
    printf("smearing error: %.4e\n", fit_result.sigma_dydz_smearing );  

    fit_result.sigma_dydz_overall = sqrt( 
        pow(fit_result.sigma_dydz_position, 2) + 
        pow(fit_result.sigma_dydz_smearing, 2) 
    ); 

    printf("Total: %.4e\n", fit_result.sigma_dydz_overall); 

    return fit_result; 
}

std::function<double(double*,double*)> FitBackground(TH1D* hist, const vector<pair<double,double>>& exclusions, const int order) 
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
