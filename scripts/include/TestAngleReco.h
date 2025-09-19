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
#include <memory>
#include <TBox.h> 
#include <optional> 

using namespace std; 
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


//test the horizontal angle reconstruction efficiency
std::optional<AngleFitResult_t> TestAngleReco(  const bool is_RHRS,         //arm to use
                                                ROOT::RDF::RNode df,        //RDataFrame to get our data from 
                                                OpticsTarget_t target,      //optics target to use  
                                                const int center_row,       //central row to start at
                                                const int n_side_rows,      //number of rows IN ADDITION TO THE CENTRAL row on either side to fit       
                                                const int center_col,       //central column to fit
                                                const int n_side_cols,      //number of columns on eiter side to fit
                                                const double row_cut_width=1.2,     //width to make vertical cut, expressed as a fraction of the row spacing 
                                                const double col_cut_width=0.80,    //width to make cut around hole, in fraction of inter-column spacing
                                                const char* name_dxdz="reco_dxdz_sv", 
                                                const char* name_dydz="reco_dydz_sv", 
                                                const char* name_react_vtx="position_vtx" )
{
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
            ostringstream oss; 
            oss << "in <FindSieveHole>: invalid sieve hole requested. row/col: " << row << "/" << col; 
            throw invalid_argument(oss.str()); 
            return SieveHole{}; 
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

    auto c_holes = new TCanvas("c_holes", "Angle Reco", 800, 600); 
    
    const double dydz_max = +0.04; 
    const double dydz_min = -0.04;

    const int canvId_dx_dy = 1; 

    gStyle->SetPalette(kBird);

    auto h_dx_dy = df
        .Histo2D<double>({"h", "dx/dz_{sv} vs dy/dz_{sv} (rad)", 200, -0.04, 0.04, 200, -0.04, 0.04}, name_dxdz, name_dydz); 

    h_dx_dy->DrawCopy("col2"); 


    size_t n_holes_measured =0; 
    
    //here, we're testing vertical wires.     
    auto c_cuts = new TCanvas("c_cuts", "", 1600, 800); 
    c_cuts->Divide( 1, 1 + n_side_rows*2, 0.01,0.); 

    int i_canv=1;

    //row counter
    int i_row=1; 

    
    AngleFitResult_t fit_result; 


    //the elements of this pair are 'dxdz' and 'dydz'; they represent which holes were fit successfully (and therefore should be drawn). 
    vector<pair<double,double>> holes_to_draw; 

    //loop thru eaach column 
    for (int row = center_row-n_side_rows; row <= center_row+n_side_rows; row++) {

        vector<SieveHole> row_holes; 
        for (int col=center_col-n_side_cols; col<center_col+n_side_cols; col++) {
            try {
                row_holes.push_back( FindSieveHole(row, col) ); 
            } catch (const std::exception& e) {
                Error(here, "Error trying to access sieve hole.\n what(): %s", e.what()); 
                return std::nullopt; 
            }
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
                
        //..
        const double row_dxdz = hole_angles.front().dxdz; 
        
        c_cuts->cd(i_row); 

        auto hist_row = df
            
            .Filter([dxdz_cut_width, row_dxdz](double dxdz){ return fabs(dxdz - row_dxdz) < dxdz_cut_width/2.; }, {name_dxdz})

            .Histo1D<double>({Form("h_dydz_%i",i_row), "dy/dz_{sv}", 200, dydz_min, dydz_max}, name_dydz);

        //make it so that the y-axes don't have tick-marks
        hist_row->GetYaxis()->SetNdivisions(0); 
        hist_row->DrawCopy(); 

        auto x_axis = hist_row->GetXaxis(); 

        int i_col =0; 
        for (size_t i=0; i<row_holes.size(); i++) {

            const auto& hole  = row_holes[i]; 
            const auto& angle = hole_angles[i]; 

            auto fcn_gauss_offset = [](double *x, double *par) { 
                double sigma = fabs(par[2]); 
                return par[0] * exp( -0.5 * pow( (x[0] - par[1])/sigma, 2 ) ) + fabs(par[3]); 
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

            double offset = (hist_row->GetBinContent(bin_min) + hist_row->GetBinContent(bin_max))/2.; 

            //upper edge of sieve hole, accounting for paralax, assuming that the tungsten material of the sieve is impenitrable
            // (which is known not to be true!!)
            const double dydz_hi = ( (hole.y + hole.radius_front) - react_vtx_scs.y() ) / ( sieve_thickness_z - react_vtx_scs.z() ); 

            const double dydz_lo = ( (hole.y - hole.radius_front) - react_vtx_scs.y() ) / ( 0. - react_vtx_scs.z() ); 

            const double hole_radius = (dydz_hi - dydz_lo) / 2.;
            const double hole_dydz   = (dydz_hi + dydz_lo) / 2.;
                
            //auto fit = unique_ptr<TF1>(new TF1("holefit", gaussian_semicircle, hole_dydz -2.5e-3, hole_dydz +2.5e-3, 5));              

            auto fit = unique_ptr<TF1>(new TF1(
                "holefit", 
                fcn_gauss_offset, 
                hole_dydz - dydz_cut_width/2., 
                hole_dydz + dydz_cut_width/2., 
                4
            ));              

            fit->SetParameter(0, maxval - offset); 
            fit->SetParameter(1, hole_dydz);
            fit->SetParameter(2, hole_radius / 2.5);

            fit->SetParameter(3, offset); 
            fit->SetParLimits(3, 0., maxval); 

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

            c_cuts->cd(i_row);

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
                angle.dydz - dydz_cut_width/2.,
                angle.dxdz + dxdz_cut_width/2., 
                angle.dydz + dydz_cut_width/2.
            );
            //make the box have no fill (transparent), and draw it.
             
            c_holes->cd();

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

        printf("Done with row %2i/%i\n", i_row, 2*n_side_rows + 1); cout << flush; 

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

#endif 
