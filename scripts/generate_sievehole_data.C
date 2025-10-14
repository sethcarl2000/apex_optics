//apex src headers
#include <SieveHole.h> 
#include <ApexOptics.h> 
//local 'scripts/include' folder headers
#include "include/ChainedOpticsModel.h"
#include "include/RDFNodeAccumulator.h"
#include "include/Add_TParameter_to_TFile.h"
#include "include/Add_branch_from_Trajectory_t.h"
//ROOT headers
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h> 
#include <TStyle.h> 
#include <TBox.h> 
#include <TLine.h> 
#include <TF1.h> 
#include <TFitResult.h> 
#include <TFitResultPtr.h> 
#include <TH1D.h> 
#include <TH2D.h> 
//std lib headers
#include <vector>
#include <string> 
#include <optional> 
#include <stdexcept>
#include <cmath>  
#include <utility> 

//simple 1d gauss fcn with const offset
double fcn_gauss_with_offset(const double* X, const double* par); 

TFitResult* Fit_gauss_to_TH1D(TH1D* hist); 


using namespace std; 
using ApexOptics::Trajectory_t; 
using ApexOptics::OpticsTarget_t; 

const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};

//Given a TH2*, find the maximum bin-value in the subrange given by [x0,x1] & [y0,y1]. 
double GetMaximum_in_subrange(TH2* hist, double x0,double y0, double x1,double y1)
{   
    //check if the hist is null. 
    if (!hist) {
        throw std::invalid_argument("in <Max_bin_in_subrange>: null TH2* object was passed."); 
        return 0.; 
    }

    auto x_ax = hist->GetXaxis(); 
    auto y_ax = hist->GetYaxis(); 

    double max_content=-1.e30; 

    for (int bx=x_ax->FindBin(x0); bx<=x_ax->FindBin(x1); bx++) {
        for (int by=y_ax->FindBin(y0); by<=y_ax->FindBin(y1); by++) {

            max_content = max<double>( max_content, hist->GetBinContent(bx, by) ); 
        }
    }
    return max_content; 
}


//_______________________________________________________________________________________________________________________________________________
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int generate_sievehole_data(const char* path_infile     = "data/replay/real_L_H3-dp.root",
                            const char* target_name     = "H3",
                            const char* path_outfile    = "data/sieve_holes/holes_L_H3-nofit.root",
                            const char* path_graphic    = "histos/holes_L_H3",
                            const char* tree_name       = "tracks_fp" )
{
    const char* const here = "cleanup_optics_replay";
    
    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df(tree_name, path_infile); 

    //select the target chosen
    OpticsTarget_t target; 
    try { target = ApexOptics::GetTarget(string(target_name)); } 

    catch (const std::exception& e) {

        Error(here, "Something went wrong trying to get the target info.\n what(): %s", e.what()); 
        return -1; 
    }
    
    printf("Optics target chosen is '%s'.\n", target.name.c_str()); 

    //try to figure out which arm we're using. quit if something goes wrong. 
    auto is_RHRS_opt = Get_TParameter_from_TFile<bool>(path_infile, "is_RHRS"); 
    if (!is_RHRS_opt.has_value()) {
        Error(here, "Something went wrong trying to get the 'is_RHRS' parameter from the file: '%s'", path_infile); 
        return -1; 
    }    
    const bool is_RHRS = is_RHRS_opt.value(); 

    //try to initialize the optics model
    ChainedOpticsModel* model = new ChainedOpticsModel(is_RHRS); 

    try {
        model->CreateChainRev({ // sv <= [Poly] q1-fwd <= _Poly_ <= fp 
            {"data/csv/poly_fits_fp_q1-fwd_L_4ord.dat", branches_fwd_q1, 4}, 
            {"data/csv/poly_prod_q1_sv_L_4ord.dat",     branches_sv,     5} 
        }); 

        model->CreateChainFwd({ // sv => [Poly] => fp
            {"data/csv/poly_prod_sv_fp_L_4ord.dat", branches_fp, 5} 
        }); 
    
    } catch (const std::exception& e) {
        Error(here, "Something went wrong trying to build the polynomial chains.\n what(): %s", e.what()); 
        return -1; 
    }

    //figure out if we're going to make a graphic 
    const bool make_graphic = (string(path_graphic)!=""); 

    //figure out if we're making an output file 
    const bool make_outfile = (string(path_outfile)!="");


    //try and make the RDataFrame
    RDFNodeAccumulator rna(tree_name, path_infile); 

    rna.Define("Xsv", [](double x, double y, double dxdz, double dydz, double dpp)
        {
            return Trajectory_t{x, y, dxdz, dydz, dpp}; 
        }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv", "dpp_sv"});

    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
        {
            return Trajectory_t{x, y, dxdz, dydz}; 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 
    
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
            return std::optional<SieveHole>(std::nullopt); 
        } else {
            return std::optional<SieveHole>(SieveHole{*it});
        }   
    }; 

    
    //get the react vertex
    //we're going to perfom these calculations to get the (average) dx/dz and dy/dz values for each hole. 


    //first, we must compute the react vertex's position in the Sieve Coordinate System (SCS), if it has not already been done. 
    rna.DefineIfMissing("position_vtx_scs", [is_RHRS](TVector3 vtx_hcs)
        {
            return ApexOptics::HCS_to_SCS(is_RHRS, vtx_hcs); 
        }, {"position_vtx"});

    //now, we compute its average value of the react vertex. 
    const TVector3 position_vtx_scs(
        *(rna.Get().Define("val", [](TVector3 vtx){ return vtx.x(); }, {"position_vtx_scs"}).Mean("val")), 
        *(rna.Get().Define("val", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx_scs"}).Mean("val")), 
        *(rna.Get().Define("val", [](TVector3 vtx){ return vtx.z(); }, {"position_vtx_scs"}).Mean("val"))
    );
    
    //once we have this, we can compute the average slope to exect for each hole. 
    const double sieve_dz = 12.5e-3; //in meters
    const double dx = 5.842e-3; //all units in meters
    const double dy = 4.826e-3;

    struct HolePosWithAngles_t { 
        int row, col; 
        double dxdz_sv, dydz_sv;
    };
    std::vector<HolePosWithAngles_t> sieve_hole_angles; 
    for (const auto hole : sieve_holes) {

        //first, we 'aim' at the front of the sieve 
        double dxdz_front = (hole.x - position_vtx_scs.x())/(0. - position_vtx_scs.z()); 
        double dydz_front = (hole.y - position_vtx_scs.y())/(0. - position_vtx_scs.z()); 

        //then, we 'aim' at the back of the sieve hole's exit 
        double dxdz_back  = (hole.x - position_vtx_scs.x())/(sieve_dz - position_vtx_scs.z()); 
        double dydz_back  = (hole.y - position_vtx_scs.y())/(sieve_dz - position_vtx_scs.z()); 

        //we take the average of these to to determine the best spot to 'guess' the center of the sieve hole's true position 
        sieve_hole_angles.push_back({
            .row = hole.row, 
            .col = hole.col, 
            .dxdz_sv = 0.5*(dxdz_front + dxdz_back),
            .dydz_sv = 0.5*(dydz_front + dydz_back)
        }); 
    }

    

    
    auto hist_dx_dy_ptr = rna.Get().Histo2D<double>({"h", "hist xy", 100, -0.05,0.05, 200, -0.04,0.03}, "dxdz_sv","dydz_sv"); 

    auto hist_dx_dy = (TH2D*)hist_dx_dy_ptr->Clone("hist_dx_dy"); 

    auto c = new TCanvas; 
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(kBird); hist_dx_dy->GetZaxis()->SetNdivisions(200.); 
    hist_dx_dy->Draw("col"); 
    
    c->Modified(); 
    c->Update(); 
    
    auto c_fp = new TCanvas("c_fp", "focal-plane coords", 1500, 700); 
    c_fp->Divide(2,2); 
    c_fp->cd(1); hist_dx_dy->Draw("col"); 
    c_fp->Modified(); 
    c_fp->Update(); 

    //range of columns and rows to try to fit 
    const int row_0 = 1;//3; 
    const int row_1 = 12;//12; 
    
    const int col_0 = 0; 
    const int col_1 = 10; 

    //statistical constraints (to avoid fitting cells which have no stats)
    
    //for each fit subrange, the maximum bin-value in that subrange must be at least this fraction of the 
    // histogram's global maximum. otherwise, a fit is not attempted. 
    const double subrange_max_ratio_min = 0.10;    

    //number of bins to use in the 1d fit histogram 
    const int nbins_fit = 25; 

    //fraction of the spacing to the next hole to make our cut. 
    // '=1.' corresponds to a cut which spans to the center of the next hole on either side. 
    const double cut_width_row = 0.60 * dx / fabs(position_vtx_scs.z()); 
    const double cut_width_col = 0.50 * dy / fabs(position_vtx_scs.z());  

    //number of sigmas in x & y to make final cut on events
    const double cut_sigma_mult = 2.; 

    //max acceptable sigma for fit
    const double sigma_x_max = 0.55 * cut_width_row; 
    const double sigma_y_max = 0.75 * cut_width_col; 

    //min ratio between amplitude of gaussian and offset background
    const double min_amplitude_offset_ratio = 0.0;     

    //how far we're going to let the peak center be from the edge
    const double max_peak_dist_from_center_x = 0.25; 
    const double max_peak_dist_from_center_y = 0.45; 


    //cache the sieve-coords and fp-coords in memory, to make it faster. 
    auto df_cache = rna.Get().Cache({
        "x_sv", 
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv",

        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp", 

        "position_vtx_scs"
    });     

    //this will be where we take our data. 
    vector<pair<Trajectory_t,Trajectory_t>> coordinate_data; 

    TBox* cursor_box = nullptr; 

    int i_frame=0; 
    //loop thru all the sieve holes, and try to fit the ones (which exist). 
    for (int row=row_0; row<=row_1; row++) {
        for (int col=col_0; col<=col_1; col++) {

            //find the corresponding sieve hole. if it doesn't exist, then skip it. 
            auto hole_find = FindSieveHole(row, col);
            if (!hole_find.has_value()) continue; 

            const auto hole = hole_find.value(); 

            //find the sieve hole angles we computed before
            HolePosWithAngles_t hole_angle{.row=-1}; 
            for (const auto& ha : sieve_hole_angles) { 
                if (ha.row == hole.row && ha.col == hole.col) { hole_angle = ha; break; }
            }
            if (hole_angle.row==-1) {
                Error(here, "Something went wrong trying to find the correct 'HolePosWithAngles_t'. row=%i, col=%i", 
                    hole.row, hole.col    
                ); 
                return -1;
            }

            const double cut_halfwidth_x = cut_width_row; 
            const double cut_halfwidth_y = cut_width_col; 

            //get the fit subrange
            const double 
                dxdz0{hole_angle.dxdz_sv - cut_halfwidth_x}, 
                dydz0{hole_angle.dydz_sv - cut_halfwidth_y}, 

                dxdz1{hole_angle.dxdz_sv + cut_halfwidth_x}, 
                dydz1{hole_angle.dydz_sv + cut_halfwidth_y}; 

            
            //draw the box of the hole we're trying to fit, currently 
            if (cursor_box != nullptr) delete cursor_box; 
            cursor_box = new TBox(dxdz0,dydz0, dxdz1,dydz1); 
            cursor_box->SetLineColor(kRed); 
            cursor_box->SetFillStyle(0); 
            c->cd(); cursor_box->Draw(); 
            c->Modified(); c->Update(); 
            

            double subrange_max = GetMaximum_in_subrange(hist_dx_dy, dxdz0,dydz0, dxdz1,dydz1); 
            double global_max   = hist_dx_dy->GetMaximum(); 
            //check the ratio of the subrange max to the global max. if it's too small, then don't attemt a fit. 
            if (subrange_max / global_max < subrange_max_ratio_min) continue; 
            
            //create the 2d hist of the events in this box. 
            auto hist_box = df_cache  
                //make a cut on events in our rectangle
                .Filter([&](double dxdz){ return fabs(dxdz - hole_angle.dxdz_sv) < cut_halfwidth_x; }, {"dxdz_sv"})
                .Filter([&](double dydz){ return fabs(dydz - hole_angle.dydz_sv) < cut_halfwidth_y; }, {"dydz_sv"})
                
                .Histo2D<double>({
                    "h_box", "", 
                    nbins_fit, hole_angle.dxdz_sv - cut_halfwidth_x, hole_angle.dxdz_sv + cut_halfwidth_x,
                    nbins_fit, hole_angle.dydz_sv - cut_halfwidth_y, hole_angle.dydz_sv + cut_halfwidth_y    
                }, "dxdz_sv", "dydz_sv"); 

            //create each 1d-hist that we will attempt to fit
            auto fit_x = Fit_gauss_to_TH1D(hist_box->ProjectionX("h_box_dxdz")); 
            if (!fit_x) continue; 
            const double dxdz_cent  = fit_x->Parameter(1); 
            const double dxdz_sigma = fabs(fit_x->Parameter(2)); 

            //check for reasonability 
            // - center of fit must be inside cut range
            if (fabs(dxdz_cent - hole_angle.dxdz_sv)/cut_halfwidth_x > max_peak_dist_from_center_x) continue; 
            // - amplitude must be at least as large as background * 0.20
            if (fabs(fit_x->Parameter(0)/fit_x->Parameter(3)) < min_amplitude_offset_ratio) continue; 
            // - max acceptable sigma
            if (dxdz_sigma > sigma_x_max) continue; 

            
            auto fit_y = Fit_gauss_to_TH1D(hist_box->ProjectionY("h_box_dydz")); 
            if (!fit_y) continue; 
            const double dydz_cent  = fit_y->Parameter(1); 
            const double dydz_sigma = fabs(fit_y->Parameter(2));
            
            //check for reasonability
            // - center of fit must be inside cut range
            if (fabs(dydz_cent - hole_angle.dydz_sv)/cut_halfwidth_y > max_peak_dist_from_center_y) continue; 
            // - amplitude must be at least as large as background * 0.20 
            if (fabs(fit_y->Parameter(0)/fit_y->Parameter(3)) < min_amplitude_offset_ratio) continue; 
            // - max acceptable sigma
            if (dydz_sigma > sigma_y_max) continue; 

            auto box = new TBox(
                dxdz_cent + dxdz_sigma*cut_sigma_mult, dydz_cent + dydz_sigma*cut_sigma_mult, 
                dxdz_cent - dxdz_sigma*cut_sigma_mult, dydz_cent - dydz_sigma*cut_sigma_mult
            );

            c->cd(); 

            box->SetFillStyle(0);
            box->SetLineColor(kBlack); 
            box->SetLineStyle(kDotted); 
            box->Draw(); 

            TLine *line=nullptr; 
            line = new TLine(dxdz_cent,dydz0, dxdz_cent,dydz1); 
            line->SetLineColor(kBlack); 
            line->Draw(); 

            line = new TLine(dxdz0,dydz_cent, dxdz1,dydz_cent); 
            line->SetLineColor(kBlack); 
            line->Draw(); 

            c->Modified(); 
            c->Update(); 

            //now, lets plot the data in the focal coordinate plane
            auto df_cut = df_cache
                //make a cut on events in our rectangle
                .Filter([&](double dxdz){ return fabs(dxdz - dxdz_cent) < dxdz_sigma*cut_sigma_mult; }, {"dxdz_sv"})
                .Filter([&](double dydz){ return fabs(dydz - dydz_cent) < dydz_sigma*cut_sigma_mult; }, {"dydz_sv"}); 

            c_fp->cd(1);
            box->Draw(); 

            auto h_x_y    = df_cut.Histo2D<double>({"h_x_y",  "x_{fp} vs y_{fp}",     50,-0.6,0.6, 60,-0.05,0.05}, "x_fp", "y_fp"); 
            auto h_x_dxdz = df_cut.Histo2D<double>({"h_x_dy", "x_{fp} vs dx/dz_{fp}", 50,-0.6,0.6, 60,-0.05,0.05}, "x_fp", "dxdz_fp"); 
            auto h_x_dydz = df_cut.Histo2D<double>({"h_x_dx", "x_{fp} vs dy/dz_{fp}", 50,-0.6,0.6, 60,-0.05,0.05}, "x_fp", "dydz_fp"); 
            
            c_fp->cd(2); h_x_y->Draw("col"); 
            c_fp->cd(3); h_x_dxdz->Draw("col"); 
            c_fp->cd(4); h_x_dydz->Draw("col"); 

            c_fp->Modified(); 
            c_fp->Update(); 

            if (make_outfile) {
                // 'take' this data from the RDataFrame (get a vector of all events)
                vector<pair<Trajectory_t,Trajectory_t>> hole_coord_data = *df_cut

                    .Define("Xsv", [](double x, double y, double dxdz, double dydz, double dpp){ return Trajectory_t{x,y,dxdz,dydz,dpp}; }, 
                            {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"})

                    .Define("Xfp", [](double x, double y, double dxdz, double dydz){ return Trajectory_t{x,y,dxdz,dydz}; }, 
                            {"x_fp","y_fp","dxdz_fp","dydz_fp"})
                        
                    .Define("coord_data", [&hole, sieve_dz](Trajectory_t Xsv, Trajectory_t Xfp, TVector3 vtx_scs)
                    {
                        //take an average of the 'back' and 'front' hole opening, 
                        // which is a rough way to approximate the effects of 
                        // parallax (for react vertices which are askew from the axis of this sieve-hole)
                        Trajectory_t Xsv_front{
                            hole.x, 
                            hole.y, 
                            (hole.x - vtx_scs.x())/(0. - vtx_scs.z()), 
                            (hole.y - vtx_scs.y())/(0. - vtx_scs.z()),
                            Xsv.dpp
                        };

                        Trajectory_t Xsv_back{
                            hole.x, 
                            hole.y, 
                            (hole.x - vtx_scs.x())/(sieve_dz - vtx_scs.z()), 
                            (hole.y - vtx_scs.y())/(sieve_dz - vtx_scs.z()),
                            Xsv.dpp
                        };

                        //take the average of the two 
                        auto Xsv_corrected = (Xsv_front + Xsv_back) * 0.5; 

                        //compute dx/dz & dy/dz
                        return pair<Trajectory_t, Trajectory_t>({Xsv_corrected, Xfp}); 
                    }, {"Xsv", "Xfp", "position_vtx_scs"})

                    .Take<pair<Trajectory_t,Trajectory_t>>("coord_data"); 

                //now, we actually 'take' the data from the rdataframe. 
                
                //stick this data into the vector we need. 
                coordinate_data.insert( 
                    coordinate_data.end(), 
                    std::make_move_iterator(hole_coord_data.begin()), 
                    std::make_move_iterator(hole_coord_data.end())                      
                );
            }

            if (make_graphic) {
                if (i_frame==0) { c_fp->SaveAs(Form("%s.pdf(",path_graphic)); }
                else            { c_fp->SaveAs(Form("%s.pdf",path_graphic));  } 
                i_frame++; 
            }
        }
    }

    if (make_outfile) {
        //now, we can stick this data in a file. 
        //we need to run in single-thread mode
        if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 

        size_t i_event=0;
        const size_t n_events = coordinate_data.size();  
        ROOT::RDataFrame df_out(n_events);
        RDFNodeAccumulator rna_out(df_out); 
        
        rna_out.Define("coord_data", [&coordinate_data, &i_event, n_events]()
            {
                if (i_event >= n_events) {
                    throw std::logic_error("Error! we have tried to access events which don't exist!"); 
                    return std::pair<Trajectory_t,Trajectory_t>({}); 
                } 
                return coordinate_data.at(i_event++); 
            }, {});
        rna_out.Define("Xsv", [](std::pair<Trajectory_t,Trajectory_t> pair){ return pair.first;  }, {"coord_data"}); 
        rna_out.Define("Xfp", [](std::pair<Trajectory_t,Trajectory_t> pair){ return pair.second; }, {"coord_data"}); 

        rna_out = Add_branch_from_Trajectory_t(rna_out.Get(), "Xsv", {
            {"x_sv",    &Trajectory_t::x}, 
            {"y_sv",    &Trajectory_t::y}, 
            {"dxdz_sv", &Trajectory_t::dxdz}, 
            {"dydz_sv", &Trajectory_t::dydz}, 
            {"dpp_sv",  &Trajectory_t::dpp} 
        }); 

        rna_out = Add_branch_from_Trajectory_t(rna_out.Get(), "Xfp", {
            {"x_fp",    &Trajectory_t::x}, 
            {"y_fp",    &Trajectory_t::y}, 
            {"dxdz_fp", &Trajectory_t::dxdz}, 
            {"dydz_fp", &Trajectory_t::dydz}, 
        }); 

        cout << "making snapshot. (path='" << path_outfile << ")..." << flush; 

        rna_out.Get().Snapshot("tracks_fp", path_outfile, {
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",

            "x_fp",
            "y_fp",
            "dxdz_fp",
            "dydz_fp"
        });

        cout << "done" << endl; 

        //add the parameter to the file
        auto file = new TFile(path_outfile, "UPDATE");
        Add_TParameter_to_TFile<bool>("is_RHRS", is_RHRS);
        file->Close(); 
        delete file;  
    }

    if (make_graphic) {
        c->cd(); 
        c->SaveAs(Form("%s.pdf(",path_graphic)); 
    }
    
    return 0; 
}

//______________________________________________________________________________________________________________________________
double fcn_gauss_with_offset(const double* X, const double* par)
{
    //params: 
    // par[0]   amplitude
    // par[1]   mean
    // par[2]   sigma
    // par[3]   offset

    double x = (X[0] - par[1])/fabs(par[2]); 
    
    return (par[0] * exp(-0.5 * x*x)) + fabs(par[3]); 
}
//______________________________________________________________________________________________________________________________
TFitResult* Fit_gauss_to_TH1D(TH1D* hist)
{
    if (hist==nullptr) return nullptr; 

    //attempt to fit the histogram. lets' make some reasonable guesses for parameters: 
    double offset       = hist->GetMinimum(); 
    double amplitude    = hist->GetMaximum() - offset; 
    double mean         = hist->GetXaxis()->GetBinCenter( hist->GetMaximumBin() ); 
    double sigma        = 0.5 * (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()); 

    //printf("first guess: %f, %f, %f, %f\n", amplitude, mean, sigma, offset);

    double params[] = { amplitude, mean, sigma, offset }; 

    auto tf1 = new TF1(
        "gauss_fit_with_offset", 
        fcn_gauss_with_offset, 
        hist->GetXaxis()->GetXmin(), 
        hist->GetXaxis()->GetXmax(), 
        4
    );

    tf1->SetParameters(params); 

    // options (second character argument):
    //  S - store the hit in the returned FitResult object 
    //  R - only fit to the function range 
    //  L - use 'log liklihood' method to fit, which is appropriate for histograms where bins represent fit-counts
    //  Q - don't print information about the fit (we'll do that oursevles.)
    //  N - do not draw the result
    //  
    TFitResult* fit_result = (hist->Fit("gauss_fit_with_offset", "R L Q S N")).Get(); 

    /*printf("result : %f, %f, %f, %f\n", 
        fit_result->Parameter(0), 
        fit_result->Parameter(1), 
        fit_result->Parameter(2), 
        fit_result->Parameter(3)
    );*/ 

    if (!fit_result) return nullptr; 

    if (!fit_result->IsValid()) return nullptr; 

    return fit_result; 

}