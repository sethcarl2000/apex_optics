#ifndef measure_model_accuracy_H
#define measure_model_accuracy_H

//APEX headers
#include <ApexOptics.h>
#include <AngleRecoTester.h>
#include "RDFNodeAccumulator.h"
#include "BothArmValue_t.h"
#include <SieveHole.h> 
#include <ArmMode.h>
//ROOT headers
#include <ROOT/RDataFrame.hxx>  
#include <TRandom3.h> 
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TVirtualPad.h> 
#include <TBox.h> 
#include <TLine.h> 
#include <TEllipse.h> 
//stdlib headers
#include <vector>
#include <string> 
#include <limits> 
#include <functional> 
#include <stdexcept> 
#include <algorithm> 
#include <map> 

struct MassResolutionFit_t {
    double aprime_mass, R_RMS_dxdz, R_RMS_dydz, L_RMS_dxdz, L_RMS_dydz; 
}; 

struct HoleAndAngles_t {
    SieveHole hole; 
    double dxdz; 
    double dydz; 
}; 


//give this a general function, and it will perform some operation on both arms 


/**
 * @brief measures & reports reconstruction accuracy for the left or right arm
 */
class MeasureModelAccuracy {
private: 

    bool is_NaN(double x) const { return x != x; }
public: 
    
    ///reduction in angular resolution resulting from multiple-scattering in the target
    const double kRMS_multiple_scattering = 0.36e-3; 

    /// path to all A' montecarlo files
    std::vector<std::string> fPath_Aprime_montecarlo{
        "data/mc/mc_Aprime_140-MeV.root",
        "data/mc/mc_Aprime_150-MeV.root",
        "data/mc/mc_Aprime_160-MeV.root",
        "data/mc/mc_Aprime_170-MeV.root",
        "data/mc/mc_Aprime_180-MeV.root",
        "data/mc/mc_Aprime_190-MeV.root", 
        "data/mc/mc_Aprime_200-MeV.root",
        "data/mc/mc_Aprime_210-MeV.root",
        "data/mc/mc_Aprime_220-MeV.root",
        "data/mc/mc_Aprime_230-MeV.root",
        "data/mc/mc_Aprime_240-MeV.root",
        "data/mc/mc_Aprime_250-MeV.root" 
    };

    std::string fPath_R_V1{"data/replay/real_R_V1-v01.root"}; 
    std::string fPath_R_V2{"data/replay/real_R_V2-v01.root"}; 
    std::string fPath_R_V3{"data/replay/real_R_V3-v01.root"}; 

    std::string fPath_L_V1{"data/replay/real_L_V1.root"}; 
    std::string fPath_L_V2{"data/replay/real_L_V2.root"}; 
    std::string fPath_L_V3{"data/replay/real_L_V3.root"}; 

    std::string fBranch_R_dxdz{"R_dxdz_sv"}; 
    std::string fBranch_R_dydz{"R_dydz_sv"}; 
    
    std::string fBranch_L_dxdz{"R_dxdz_sv"}; 
    std::string fBranch_L_dydz{"R_dydz_sv"}; 
    
    ArmMode::Bit fArm_mode = ArmMode::kBoth; 

    /// drawing palette to use for 2D histograms 
    int fPalette = kBird; 

    /// if true, then J.Williamson's results for the accuracy will be drawn. 
    bool fDraw_prev_results=false; 

    /// drawing parameters for dxdz & dydz fits
    int fLineColor_dxdz=kBlack;  
    int fLineColor_dydz=kBlack; 

    //spacing for sieve-holes
    const double sieve_hole_spacing_x = 5.842e-3;   //sieve-hole spacing, m
    const double sieve_hole_spacing_y = 4.826e-3;   //sieve-hole spacing, m

    MeasureModelAccuracy() {}; 

    std::function<ROOT::RDF::RNode(ArmMode::Bit, ROOT::RDF::RNode)> fModel; 

    int EvaluateModel ( 
        std::function<ROOT::RDF::RNode(ArmMode::Bit,ROOT::RDF::RNode)> optics_model
    ); 
}; 


int MeasureModelAccuracy::EvaluateModel(std::function<ROOT::RDF::RNode(ArmMode::Bit, ROOT::RDF::RNode)> model)    
{        
    using namespace ApexOptics; 
    using namespace std;

    struct TargetAndPath_t { ApexOptics::OpticsTarget_t target; std::string path; };

    const char* tree_name="tracks_fp";

    //check to make sure we have the correct branches

    auto is_nan = [](double x) { return x==x ? false : true; };

    //vertical wiress
    const OpticsTarget_t targets[] = {
        GetTarget("V1"),
        GetTarget("V2"),
        GetTarget("V3")
    };

    BothArmValue_t<ROOT::RDataFrame*> 
        df_V1_in{0, 0}, 
        df_V2_in{0, 0}, 
        df_V3_in{0, 0}; 

    //try to construct the RDataFrame's 
    try { 
        df_V1_in.R = new ROOT::RDataFrame(tree_name, fPath_R_V1.c_str()); 
        df_V2_in.R = new ROOT::RDataFrame(tree_name, fPath_R_V2.c_str()); 
        df_V3_in.R = new ROOT::RDataFrame(tree_name, fPath_R_V3.c_str()); 
        df_V1_in.L = new ROOT::RDataFrame(tree_name, fPath_L_V1.c_str()); 
        df_V2_in.L = new ROOT::RDataFrame(tree_name, fPath_L_V2.c_str()); 
        df_V3_in.L = new ROOT::RDataFrame(tree_name, fPath_L_V3.c_str()); 
    }
    catch (const std::invalid_argument &e) { 
        Error(__func__, "Error trying to open one of the V-wire files.\n what(): %s", e.what()); 
        return -1; 
    }

    //get a list of sieve_holes
    auto c_wires = new TCanvas("c_reco", "Optics model applied to all vertical wires", 1400, 700); 
    c_wires->Divide(3,2, 0.00,0.01); 

    BothArmValue_t<TVirtualPad*> c_wire[3] = {
        { c_wires->cd(1), c_wires->cd(4) }, 
        { c_wires->cd(2), c_wires->cd(5) }, 
        { c_wires->cd(3), c_wires->cd(6) }
    };  

    
    c_wire[0](ArmMode::kRHRS)->cd()->SetLeftMargin(0.12);
    c_wire[0](ArmMode::kLHRS)->cd()->SetLeftMargin(0.12);

    constexpr double NaN_double = std::numeric_limits<double>::quiet_NaN(); 
    
    struct SieveHoleAndWire_t {
        SieveHole hole; 
        Trajectory_t traj; 
        double stats=0.; 
        double rms_dxdz=NaN_double;
        double rms_dydz=NaN_double; 
    }; 
    BothArmValue_t<std::vector<SieveHole>> sieve_holes{
        ApexOptics::ConstructSieveHoles(true), 
        ApexOptics::ConstructSieveHoles(false)
    }; 
    //std::vector<SieveHole> L_sieve_holes = ApexOptics::ConstructSieveHoles(false); 

    //const double average_y_hcs = *mc_Vwire_nodes[0].Define("y_hcs", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}).Mean("y_hcs");

    //a list of sieve-holes for each wire
    BothArmValue_t<vector<SieveHoleAndWire_t>> wire_holes[3];
    //now, we compute the get ready to actually perform the measurements 

    //these are flags for use with the angle tester
    const auto kDxdz   = AngleRecoTester::kDxdz; 
    const auto kDydz   = AngleRecoTester::kDydz; 
    const auto kSlopes = AngleRecoTester::kSlopes;
    
    constexpr bool kBool_RHRS = true;
    constexpr bool kBool_LHRS = false;
    
    //get dataframes with the info we need
    BothArmValue_t<ROOT::RDF::RNode> df_Vwire[3] = {
        { model(ArmMode::kRHRS, *df_V1_in.R), model(ArmMode::kLHRS, *df_V1_in.L) },
        { model(ArmMode::kRHRS, *df_V2_in.R), model(ArmMode::kLHRS, *df_V2_in.L) },
        { model(ArmMode::kRHRS, *df_V3_in.R), model(ArmMode::kLHRS, *df_V3_in.L) }
    }; 

    BothArmValue_t<const double[2]> range_dxdz = { {-0.050, +0.050}, {-0.050, +0.050} }; 
    BothArmValue_t<const double[2]> range_dydz = { {-0.025, +0.040}, {-0.035, +0.020} }; 
    
    AngleFitResult_t angles_dxdz[3]; 
    AngleFitResult_t angles_dydz[3]; 

    //each of these values are for wires 1, 2, 3.
    //                                           RHRS --------------------        LHRS ---------------------
    //                                           {   V1      V2      V3   }       {   V1      V2      V3   }  
    BothArmValue_t<int[3]> row_minx           = {{    4,      4,      5   },      {    4,      4,      5   }};  //min row to fit
    BothArmValue_t<int[3]> row_miny           = {{    4,      4,      5   },      {    4,      4,      5   }};
    BothArmValue_t<int[3]> row_maxx           = {{   12,     12,     11   },      {   12,     12,     11   }};  //max row to fit
    BothArmValue_t<int[3]> row_maxy           = {{   12,     12,     11   },      {   12,     12,     11   }};
    BothArmValue_t<int[3]> col_minx           = {{    3,      1,      0   },      {    3,      1,      0   }};
    BothArmValue_t<int[3]> col_miny           = {{    3,      1,      0   },      {    3,      1,      0   }};
    BothArmValue_t<int[3]> col_maxx           = {{    7,      6,      3   },      {    7,      6,      3   }};
    BothArmValue_t<int[3]> col_maxy           = {{    7,      6,      3   },      {    7,      6,      3   }};
    BothArmValue_t<double[3]> row_cut_widthx  = {{ 1.35,   1.35,   1.35   },      { 1.35,   1.35,   1.35   }};
    BothArmValue_t<double[3]> row_cut_widthy  = {{ 1.35,   1.35,   1.35   },      { 1.35,   1.35,   1.35   }};
    BothArmValue_t<double[3]> col_cut_widthx  = {{ 0.50,   0.50,   0.65   },      { 0.50,   0.50,   0.65   }};
    BothArmValue_t<double[3]> col_cut_widthy  = {{ 0.50,   0.50,   0.65   },      { 0.50,   0.50,   0.65   }};
    BothArmValue_t<double[3]> bg_cut_widthx   = {{ 1.00,   1.00,   1.00   },      { 1.00,   1.00,   1.00   }};
    BothArmValue_t<double[3]> bg_cut_widthy   = {{ 1.00,   1.00,   0.70   },      { 1.00,   1.00,   0.70   }};

    /// RMS for sieve-plane variables from drift-plane resolution
    constexpr double vdc_smearing_dxdz_sv = 0.426e-3; //rad
    constexpr double vdc_smearing_dydz_sv = 0.138e-3; //rad
    
    /// RMS from multiple-scattering in target
    constexpr double RMS_multiple_scattering = 0.360e-3; 

    
    //this counter will keep track of the overall stat-weighting used for each hole
    double dxdz_total_statweight = 0.; 
    double dydz_total_statweight = 0.; 

    double dxdz_total_square_error = 0.; 
    double dydz_total_square_error = 0.; 

    const string wire_names[] = {"V1", "V2", "V3"}; 
    BothArmValue_t<string> arm_names{ "RHRS", "LHRS" }; 

    BothArmValue_t<AngleFitResult_t>  
        angleFits_dxdz[3],
        angleFits_dydz[3]; 

    //dummy name for histogram 
    int i_hist=0; 
    
    //for each sieve-hole which has a dx/dz AND dy/dz RMS measured, keep the sieve hole, along with the RMS for both angles. 
    BothArmValue_t<vector<HoleAndAngles_t>> hole_and_RMS[3]; 

    for (int iwire=0; iwire<3; iwire++) {
        
        cout << "------------ "<< wire_names[iwire] <<"   --------------------------------------------------\n" << flush; 

        string wire_name = wire_names[iwire]; 

        double dxdz_wire_statweight = 0.;
        double dydz_wire_statweight = 0.; 

        double dxdz_wire_square_error = 0.; 
        double dydz_wire_square_error = 0.; 

        //target for this wire
        const auto& target = targets[iwire]; 

        //real-data dataframe for this wire
        auto df_wire = df_Vwire[iwire]; 
        //auto& df_wire = df_Vwire[iwire].L; 

        //setup the angle tester for both arms 
        for (const auto arm : {ArmMode::kRHRS, ArmMode::kLHRS}) {

            const string arm_name = arm_names(arm); 

            bool is_RHRS = (ArmMode::kRHRS & arm); 
            AngleRecoTester angle_tester(is_RHRS, target, df_wire(arm)); 

            //get the z-vertex location
            double sieve_dist  = fabs(ApexOptics::Get_sieve_pos(arm==ArmMode::kRHRS).z()); 

            double sieve_angle = ApexOptics::Get_sieve_angle(arm==ArmMode::kRHRS); 
            
            const double z_vertex = sieve_dist - cos(sieve_angle)*target.z_hcs; 

            angle_tester.SetDrawing(false); 

            angle_tester.SetRange_dxdz(range_dxdz(arm)[0], range_dxdz(arm)[1]); 
            angle_tester.SetRange_dydz(range_dydz(arm)[0], range_dydz(arm)[1]); 

            angle_tester.SetName_dxdz("reco_dxdz_sv");
            angle_tester.SetName_dydz("reco_dydz_sv"); 
            
            cout << "------------ "<<(arm & ArmMode::kRHRS ? "RHRS":"LHRS")<<" --------------------------------------------------\n" << flush; 

            angleFits_dxdz[iwire](arm) = angle_tester.Measure(kDxdz,     //measure dxdz
                row_minx(arm)[iwire],row_maxx(arm)[iwire],               //min row, max row 
                col_minx(arm)[iwire],col_maxx(arm)[iwire],               //min col, max col 
                row_cut_widthx(arm)[iwire], col_cut_widthx(arm)[iwire], bg_cut_widthx(arm)[iwire] ).value(); //row cut width, col cut width, background cut width 
        
            angleFits_dydz[iwire](arm) = angle_tester.Measure(kDydz | kSlopes,     //measure dxdz
                row_miny(arm)[iwire],row_maxy(arm)[iwire],               //min row, max row 
                col_miny(arm)[iwire],col_maxy(arm)[iwire],               //min col, max col 
                row_cut_widthy(arm)[iwire], col_cut_widthy(arm)[iwire], bg_cut_widthy(arm)[iwire] ).value(); //row cut width, col cut width, background cut width 
            
            //draw plots of all the fits which were successful 
            
            //draw a histogram of all events
            auto hist_holes = df_wire(arm) 
                .Histo2D<TH2D>({
                    Form("h_%i",++i_hist), 
                    Form("%s, %s;dx/dz (rad);dy/dz (rad)",wire_name.c_str(), arm_name.c_str()), 
                    200, range_dxdz(arm)[0], range_dxdz(arm)[1], 
                    200, range_dydz(arm)[0], range_dydz(arm)[1]
                }, "reco_dxdz_sv", "reco_dydz_sv");

            c_wire[iwire](arm)->cd();
            gStyle->SetOptStat(0);  
            gStyle->SetPalette(fPalette); 
            hist_holes->DrawCopy("col");

            //now, collect a list of holes for which we have measured **both** dxdz and dydz. 
            hole_and_RMS[iwire](arm).clear(); 

            for (const auto& fit_dxdz : angleFits_dxdz[iwire](arm).fits) {
                for (const auto& fit_dydz : angleFits_dydz[iwire](arm).fits) {
                    if (fit_dxdz.hole == fit_dydz.hole) {

                        //don't keep this hole if the dx/dz or dy/dz measurements are NaN 
                        if (is_NaN(fit_dxdz.overall_RMS) || is_NaN(fit_dydz.overall_RMS)) continue; 

                        double dxdz_actual = fit_dxdz.angle_real; 
                        double dydz_actual = fit_dydz.angle_real; 
                        //this line indicates the bounds of the fit region
                        auto line_bound = new TLine; 
                        line_bound->SetLineStyle(kDashed); 
                        line_bound->SetLineColor(fLineColor_dxdz); 

                        //half-width of hole cuts
                        const double dxdz_cut_halfwidth = 0.5*row_cut_widthx(arm)[iwire]*sieve_hole_spacing_x/z_vertex;     

                        //half-width of hole cuts
                        const double dydz_cut_halfwidth = 0.5*col_cut_widthy(arm)[iwire]*sieve_hole_spacing_y/z_vertex; 

                        //box showing bound of fit
                        const bool is_big_hole = fit_dxdz.hole.is_big; 
                        /*auto box = new TBox(
                            dxdz_actual - dxdz_cut_halfwidth, 
                            dydz_actual - (is_big_hole ? 1.25 : 1.00)*dydz_cut_halfwidth,
                            dxdz_actual + dxdz_cut_halfwidth, 
                            dydz_actual + (is_big_hole ? 1.25 : 1.00)*dydz_cut_halfwidth
                        );
                        box->SetLineStyle(kDashed);
                        box->SetLineWidth(1);  
                        box->SetLineColor(fLineColor_dxdz); 
                        box->SetFillStyle(0); 
                        box->Draw(); */ 

                        auto line_actual = new TLine; 
                        line_actual->SetLineStyle(kDotted); 
                        line_actual->SetLineWidth(1); 
                        line_actual->SetLineColor(fLineColor_dxdz); 

                        //draw the dx/dz centerline
                        line_actual->DrawLine(
                            dxdz_actual, 
                            dydz_actual - (is_big_hole ? 1.25 : 1.00)*dydz_cut_halfwidth,
                            dxdz_actual, 
                            dydz_actual + (is_big_hole ? 1.25 : 1.00)*dydz_cut_halfwidth
                        ); 

                        //draw the dy/dz centerline
                        line_actual->DrawLine(
                            dxdz_actual - dxdz_cut_halfwidth, 
                            dydz_actual,
                            dxdz_actual + dxdz_cut_halfwidth, 
                            dydz_actual
                        ); 

                        //now, draw the actual line
                        auto line_measure = new TLine; 
                        line_actual->SetLineStyle(kSolid); 
                        line_actual->SetLineWidth(2); 
                        line_actual->SetLineColor(fLineColor_dxdz); 

                        double dxdz_fit = fit_dxdz.angle_fit; 
                        double dydz_fit = fit_dydz.angle_fit; 
                        double slope    = fit_dydz.angle_slope; 

                        //draw the dx/dz centerline
                        line_measure->DrawLine(
                            dxdz_fit, 
                            dydz_actual - (is_big_hole ? 1.25 : 1.00)*dydz_cut_halfwidth,
                            dxdz_fit, 
                            dydz_actual + (is_big_hole ? 1.25 : 1.00)*dydz_cut_halfwidth
                        ); 

                        //draw the dy/dz centerline
                        line_measure->DrawLine(
                            dxdz_actual - dxdz_cut_halfwidth, 
                            dydz_fit    - dxdz_cut_halfwidth*slope,
                            dxdz_actual + dxdz_cut_halfwidth, 
                            dydz_fit    + dxdz_cut_halfwidth*slope
                        ); 

                        hole_and_RMS[iwire](arm).push_back({
                            .hole = fit_dxdz.hole, 
                            .dxdz = fit_dxdz.overall_RMS, 
                            .dydz = fit_dydz.overall_RMS
                        }); 
                    } 
                }
            }               
            c_wires->Modified(); c_wires->Update(); 
 
                        
            /*for (auto& hRMS : hole_and_RMS[iwire](arm)) {
                hRMS.dxdz = sqrt( hRMS.dxdz*hRMS.dxdz + vdc_smearing_dxdz_sv*vdc_smearing_dxdz_sv ); 
                hRMS.dydz = sqrt( hRMS.dydz*hRMS.dydz + vdc_smearing_dydz_sv*vdc_smearing_dydz_sv ); 
            }*/

        }; 

        printf(
            " ~~ Number of holes fit:\n"
            "       R - (dxdz/dydz/both): %zi / %zi / %zi\n"
            "       L - (dxdz/dydz/both): %zi / %zi / %zi\n", 
            angleFits_dxdz[iwire](ArmMode::kRHRS).fits.size(), 
            angleFits_dydz[iwire](ArmMode::kRHRS).fits.size(), 
            hole_and_RMS[iwire](ArmMode::kRHRS).size(), 

            angleFits_dxdz[iwire](ArmMode::kLHRS).fits.size(), 
            angleFits_dydz[iwire](ArmMode::kLHRS).fits.size(), 
            hole_and_RMS[iwire](ArmMode::kLHRS).size()
        );

    }

    vector<double> pts_mass, pts_error_total, pts_error_optics, pts_error_MS, pts_error_vdc; 

    //now that we have an estimate for the error oat each hole position, we need to use the monte-carlo data 
    // to estimate its error. 
    for (const auto& path_mass : fPath_Aprime_montecarlo) {
        
        double mA; 

        //We have to go thru this rigamarole of using a pointer, because RDataFrame's don't have a default 
        // constructer we can use, and this one is out of scope. 
        ROOT::RDataFrame *df_mc_in{nullptr}; 
        try {
            df_mc_in = new ROOT::RDataFrame(tree_name, path_mass.c_str());
        } 
        catch (const std::invalid_argument &e) {
            Error(__func__, "Error trying to open the monte-carlo file.\n what(): %s", e.what()); 
            return -1; 
        }
        auto df_mc = *df_mc_in; delete df_mc_in; 

        constexpr double V1_zcut =  -96.2e-3; 
        constexpr double V3_zcut = +103.8e-3; 

        //______________________________________________________________________________
        auto find_bestfit_hole = [](
            bool is_RHRS,
            const OpticsTarget_t& target, 
            const vector<HoleAndAngles_t>& holes_and_RMS, 
            double dxdz_sv, double dydz_sv, const TVector3& vtx_hcs)
        {
            int i_hole_best = 0; 
            double min_dist = 1.e30; 

            auto vtx = ApexOptics::HCS_to_SCS(is_RHRS, vtx_hcs); 

            int i_hole=0; 
            for (const auto& hrms : holes_and_RMS) {
                double dxdz_hole = (hrms.hole.x - vtx.x())/(0. - vtx.z()); 
                double dydz_hole = (hrms.hole.y - vtx.y())/(0. - vtx.z()); 

                double dist 
                    = (dxdz_sv - dxdz_hole)*(dxdz_sv - dxdz_hole) + 
                      (dydz_sv - dydz_hole)*(dydz_sv - dydz_hole); 
                
                if (dist < min_dist) {
                    i_hole_best = i_hole; 
                    min_dist = dist; 
                } 
                ++i_hole; 
            }
            return i_hole_best; 
        }; 
        //______________________________________________________________________________
        TRandom3 rand; 

        const double hrs_momentum = 1104.; 

        ROOT::RDF::RResultPtr<TH1D> hist_error;  

        try {
            RDFNodeAccumulator rna(df_mc); 
            
            //define some arm-specific quantities
            for (const auto arm : {ArmMode::kRHRS, ArmMode::kLHRS}) {
                
                const string arm_prefix = (arm & ArmMode::kRHRS) ? "R_" : "L_"; 
                const bool arm_bool = (arm & ArmMode::kRHRS); 

                rna.Define(arm_prefix+"traj_hcs_smeared_full", [&rand, &hole_and_RMS, &targets, arm, arm_bool, vdc_smearing_dxdz_sv, vdc_smearing_dydz_sv, &find_bestfit_hole, V1_zcut, V3_zcut, this](
                    double x, 
                    double y, 
                    double dxdz, 
                    double dydz, 
                    double dpp, 
                    const TVector3& vtx_hcs
                ) {
                    //find which wire this event is closest to (based on z-vetex position)
                    int i_wire=1;  
                    if (vtx_hcs.z() > V3_zcut) i_wire=2; //V3 
                    if (vtx_hcs.z() < V1_zcut) i_wire=0; //V1
                    
                    //find which measured sieve-hole this event is closest to (based on dx/dz and dy/dz in sieve coordinates)
                    int best_hole_id = find_bestfit_hole(true, targets[i_wire], hole_and_RMS[i_wire](arm), dxdz, dydz, vtx_hcs);

                    const auto& hRMS = hole_and_RMS[i_wire](arm).at(best_hole_id); 

                    double rms_dx = sqrt( 
                        hRMS.dxdz*hRMS.dxdz + 
                        kRMS_multiple_scattering*kRMS_multiple_scattering + 
                        vdc_smearing_dxdz_sv*vdc_smearing_dxdz_sv 
                    );  
                    double rms_dy = sqrt( 
                        hRMS.dydz*hRMS.dydz + 
                        kRMS_multiple_scattering*kRMS_multiple_scattering + 
                        vdc_smearing_dydz_sv*vdc_smearing_dydz_sv 
                    );

                    Trajectory_t traj_scs{
                        x, 
                        y, 
                        dxdz + rand.Gaus()*rms_dx,
                        dydz + rand.Gaus()*rms_dy,
                        dpp 
                    };

                    return ApexOptics::SCS_to_HCS(arm_bool, traj_scs); 
                }, {
                    arm_prefix+"x_sv", 
                    arm_prefix+"y_sv", 
                    arm_prefix+"dxdz_sv", 
                    arm_prefix+"dydz_sv", 
                    arm_prefix+"dpp_sv", 
                    arm_prefix+"position_vtx"
                }); 

                rna.Define(arm_prefix+"traj_hcs_smeared_optics", [&rand, &hole_and_RMS, &targets, arm, arm_bool, &find_bestfit_hole, V1_zcut, V3_zcut, this](
                    double x, 
                    double y, 
                    double dxdz, 
                    double dydz, 
                    double dpp, 
                    const TVector3& vtx_hcs
                ) {
                    //find which wire this event is closest to (based on z-vetex position)
                    int i_wire=1;  
                    if (vtx_hcs.z() > V3_zcut) i_wire=2; //V3 
                    if (vtx_hcs.z() < V1_zcut) i_wire=0; //V1
                    
                    //find which measured sieve-hole this event is closest to (based on dx/dz and dy/dz in sieve coordinates)
                    int best_hole_id = find_bestfit_hole(true, targets[i_wire], hole_and_RMS[i_wire](arm), dxdz, dydz, vtx_hcs);

                    const auto& hRMS = hole_and_RMS[i_wire](arm).at(best_hole_id); 

                    Trajectory_t traj_scs{
                        x, 
                        y, 
                        dxdz + rand.Gaus()*hRMS.dxdz,
                        dydz + rand.Gaus()*hRMS.dydz,
                        dpp 
                    };

                    return ApexOptics::SCS_to_HCS(arm_bool, traj_scs); 

                }, {
                    arm_prefix+"x_sv", 
                    arm_prefix+"y_sv", 
                    arm_prefix+"dxdz_sv", 
                    arm_prefix+"dydz_sv", 
                    arm_prefix+"dpp_sv", 
                    arm_prefix+"position_vtx"
                }); 

                rna.Define(arm_prefix+"traj_hcs_smeared_MS", [&rand, &hole_and_RMS, &targets, arm, arm_bool, &find_bestfit_hole, V1_zcut, V3_zcut, this](
                    double x, 
                    double y, 
                    double dxdz, 
                    double dydz, 
                    double dpp
                ) {
                    Trajectory_t traj_scs{
                        x, 
                        y, 
                        dxdz + rand.Gaus()*kRMS_multiple_scattering, 
                        dydz + rand.Gaus()*kRMS_multiple_scattering,
                        dpp 
                    };

                    return ApexOptics::SCS_to_HCS(arm_bool, traj_scs); 

                }, {
                    arm_prefix+"x_sv", 
                    arm_prefix+"y_sv", 
                    arm_prefix+"dxdz_sv", 
                    arm_prefix+"dydz_sv", 
                    arm_prefix+"dpp_sv", 
                }); 
                
                rna.Define(arm_prefix+"traj_hcs_smeared_vdc", [&rand, arm_bool, vdc_smearing_dxdz_sv, vdc_smearing_dydz_sv, this](
                    double x, 
                    double y, 
                    double dxdz, 
                    double dydz, 
                    double dpp
                ) {
                    Trajectory_t traj_scs{
                        x, 
                        y, 
                        dxdz + rand.Gaus()*vdc_smearing_dxdz_sv, 
                        dydz + rand.Gaus()*vdc_smearing_dydz_sv,
                        dpp 
                    };

                    return ApexOptics::SCS_to_HCS(arm_bool, traj_scs); 

                }, {
                    arm_prefix+"x_sv", 
                    arm_prefix+"y_sv", 
                    arm_prefix+"dxdz_sv", 
                    arm_prefix+"dydz_sv", 
                    arm_prefix+"dpp_sv"
                }); 
            
            }

            /// @return invariant mass, given target trajectories
            auto compute_invariant_mass = [hrs_momentum](Trajectory_t R_traj, Trajectory_t L_traj){
                double R_p = hrs_momentum*(1. + R_traj.dpp); 
                double L_p = hrs_momentum*(1. + L_traj.dpp); 
                
                TVector3 R_momentum( R_traj.dxdz, R_traj.dydz, 1. ); R_momentum = R_momentum.Unit() * R_p; 
                TVector3 L_momentum( L_traj.dxdz, L_traj.dydz, 1. ); L_momentum = L_momentum.Unit() * L_p; 
                
                double E = R_p + L_p; 
                double P2 = (R_momentum + L_momentum).Mag2(); 
                //ignoring correction from electron/positron mass
                return E*E - P2; 
            };

            rna.Define("m2_reco_full",   compute_invariant_mass,    {"R_traj_hcs_smeared_full",     "L_traj_hcs_smeared_full"}); 
            rna.Define("m2_reco_optics", compute_invariant_mass,    {"R_traj_hcs_smeared_optics",   "L_traj_hcs_smeared_optics"}); 
            rna.Define("m2_reco_MS",     compute_invariant_mass,    {"R_traj_hcs_smeared_MS",       "L_traj_hcs_smeared_MS"}); 
            rna.Define("m2_reco_vdc",    compute_invariant_mass,    {"R_traj_hcs_smeared_vdc",      "L_traj_hcs_smeared_vdc"}); 

            /// @param m2_name name of invariant mass^2 column 
            /// @return stddev of this reconstruction about the mean 
            auto Get_inv_mass_stddev = [&rna, &mA](const char* m2_name) {
                const auto err_name = Form("err_%s",m2_name);
                auto my_df = rna.Get(); 
                return *my_df.Define(err_name, [&mA](double m2_reco, double m)
                    {
                        mA = m; 
                        return sqrt(m2_reco) - m;   
                    }, {m2_name, "invariant_mass"}).StdDev(err_name); 
            }; 

            /*rna.Define("invariant_mass_error", [&mA](double m_reco, double m)
            { 
                mA = m; 
                return m_reco - m; 
            }, {"m_reco", "invariant_mass"}); 

            hist_error = rna.Get().Histo1D<double>(
                {"h_err", "Invariant Mass error.;reco-m_{A} - m_{A};", 200, -8., +8.}, 
                "invariant_mass_error"
            );*/

            pts_error_total .push_back( Get_inv_mass_stddev("m2_reco_full") );
            pts_error_optics.push_back( Get_inv_mass_stddev("m2_reco_optics") );
            pts_error_MS    .push_back( Get_inv_mass_stddev("m2_reco_MS") );
            pts_error_vdc   .push_back( Get_inv_mass_stddev("m2_reco_vdc") );     

            printf("A' mass: %.0f MeV,  RMS %.3f MeV\n", mA, pts_error_total.back()); cout << flush; 
            
            pts_mass.push_back( mA );
        
        } catch (const std::exception& e) {
            Error(__func__, "Something went wrong trying to define the A'-prime mass error hist.\n what(): %s\n", e.what()); 
            return -1; 
        }

        ///double rms = hist_error->GetStdDev(); 
        
        
        
        /* 
        printf(" ~~ total rms, mrad (dxdz/dydz) %.3f / %.3f\n", 
            1.e3*sqrt(dxdz_wire_square_error/dxdz_wire_statweight), 
            1.e3*sqrt(dydz_wire_square_error/dydz_wire_statweight)
        );
 
        dxdz_total_square_error += dxdz_wire_square_error; 
        dydz_total_square_error += dydz_wire_square_error; 
    
        dxdz_total_statweight += dxdz_wire_statweight; 
        dydz_total_statweight += dydz_wire_statweight;*/  
    }

    //draw the graph of John Williamson's reconstructed A' resolution
    constexpr int n_pts_JW = 20; 

    //John Williamson's reported A' resolution 
    constexpr double RMS_JW[n_pts_JW] = {
        0.8456,
        0.8684,
        0.9027,
        0.9376,
        0.9773,
        1.0236,
        1.0653,
        1.1066,
        1.1383,
        1.1592,
        1.1800,
        1.1880,
        1.1946,
        1.1965,
        1.1877,
        1.1837,
        1.1546,
        1.129,
        1.1353,
        1.0862
    }; 
    double mass_JW[n_pts_JW];
    for (int i=0; i<n_pts_JW; i++) mass_JW[i] = 127.8 + 5.*((double)i);  

    //maximum extent graph = 
    double x_max = max<double>( pts_mass.back(),  mass_JW[n_pts_JW-1] );
    double x_min = min<double>( pts_mass.front(), mass_JW[0] ); 

    auto get_array_max = [](const double* _x, int _n)
    {
        double ret=_x[0]; for (int i=1; i<_n; i++) { ret = max( _x[i], ret ); } return ret; 
    }; 

    const double max_mass = get_array_max(pts_error_total.data(), pts_error_total.size());  

    //create a new TGraph which is meant to draw the 'frame' on which the other graphs will be drawn
    double pts_frame_M[]    = { x_min, x_max }; 
    double pts_frame_RMS[]  = { 0.,    0. }; 
    
    new TCanvas; 
    auto legend = new TLegend(); 
    auto g_frame = new TGraph(2, pts_frame_M, pts_frame_RMS ); 

    g_frame->SetMaximum( max_mass*1.1 );
    g_frame->SetMinimum( 0. );

    g_frame->SetTitle("m_{A} vs #sigma_{m_{A}} (MeV/c^{2})"); 
    g_frame->Draw(); 

    auto g_total    = new TGraph(pts_mass.size(), pts_mass.data(), pts_error_total.data()); 
    auto g_optics   = new TGraph(pts_mass.size(), pts_mass.data(), pts_error_optics.data()); 
    auto g_MS       = new TGraph(pts_mass.size(), pts_mass.data(), pts_error_MS.data()); 
    auto g_vdc      = new TGraph(pts_mass.size(), pts_mass.data(), pts_error_vdc.data()); 
    
    
    g_total->SetMarkerStyle(kFullSquare); 
    g_total->Draw("SAME PL"); 
    legend->AddEntry(g_total, "Total"); 

    g_optics->SetMarkerStyle(kOpenCircle); 
    g_optics->SetMarkerColor(kRed); 
    g_optics->SetLineColor(kRed); 
    g_optics->SetLineStyle(kDashed);
    g_optics->Draw("SAME PL");
    legend->AddEntry(g_optics, "Optics"); 
    
    g_MS->SetMarkerStyle(kPlus); 
    g_MS->SetMarkerColor(kBlue); 
    g_MS->SetLineColor(kBlue); 
    g_MS->SetLineStyle(kDashed);
    g_MS->Draw("SAME PL");
    legend->AddEntry(g_MS, "M.S."); 
    
    g_vdc->SetMarkerStyle(kMultiply); 
    g_vdc->SetMarkerColor(kBlack); 
    g_vdc->SetLineColor(kBlack); 
    g_vdc->SetLineStyle(kBlack);
    g_vdc->Draw("SAME PL");
    legend->AddEntry(g_vdc, "VDC Res."); 

    legend->Draw(); 

    return 0; 
}

#endif
