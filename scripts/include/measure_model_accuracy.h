#ifndef measure_model_accuracy_H
#define measure_model_accuracy_H

//APEX headers
#include <ApexOptics.h>
#include <AngleRecoTester.h>
#include "RDFNodeAccumulator.h"
#include <SieveHole.h> 
#include "ArmMode.h"
//ROOT headers
#include <ROOT/RDataFrame.hxx>  
#include <TRandom3.h> 
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h> 
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

//very holder for quantities that have one copy for each arm. 
template<typename T> struct BothArmValue_t {
    T R; 
    T L; 
};

/**
 * @brief measures & reports reconstruction accuracy for the left or right arm
 */
class MeasureModelAccuracy {
public: 
    
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

    std::string fPath_R_V1{"data/replay/real_R_V1.root"}; 
    std::string fPath_R_V2{"data/replay/real_R_V2.root"}; 
    std::string fPath_R_V3{"data/replay/real_R_V3-v01.root"}; 

    std::string fPath_L_V1{"data/replay/real_L_V1.root"}; 
    std::string fPath_L_V2{"data/replay/real_L_V2.root"}; 
    std::string fPath_L_V3{"data/replay/real_L_V3.root"}; 

    std::string fBranch_R_dxdz{"R_dxdz_sv"}; 
    std::string fBranch_R_dydz{"R_dydz_sv"}; 
    
    std::string fBranch_L_dxdz{"R_dxdz_sv"}; 
    std::string fBranch_L_dydz{"R_dydz_sv"}; 
    
    ArmMode::Bit fArm_mode = ArmMode::kBoth; 

    MeasureModelAccuracy() {}; 

    std::function<ROOT::RDF::RNode(bool, ROOT::RDF::RNode)> fModel; 

    int EvaluateModel ( 
        std::function<ROOT::RDF::RNode(bool,ROOT::RDF::RNode)> optics_model
    ); 

    int Evaluate_singleMass(const std::string& path_mass); 
}; 


int MeasureModelAccuracy::EvaluateModel(std::function<ROOT::RDF::RNode(bool, ROOT::RDF::RNode)> model)    
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

    BothArmQuantity_t<ROOT::RDataFrame*> 
        df_V1_in{0, 0}, 
        df_V2_in{0, 0}, 
        df_V2_in{0, 0}; 

    ROOT::RDataFrame 
        *df_R_V1_in{nullptr}, 
        *df_R_V2_in{nullptr}, 
        *df_R_V3_in{nullptr}, 
        *df_L_V1_in{nullptr}, 
        *df_L_V2_in{nullptr}, 
        *df_L_V3_in{nullptr};

    //try to construct the RDataFrame's 
    try { 
        df_R_V1_in = new ROOT::RDataFrame(tree_name, fPath_R_V1.c_str()); 
        df_R_V2_in = new ROOT::RDataFrame(tree_name, fPath_R_V2.c_str()); 
        df_R_V3_in = new ROOT::RDataFrame(tree_name, fPath_R_V3.c_str()); 
        df_L_V1_in = new ROOT::RDataFrame(tree_name, fPath_L_V1.c_str()); 
        df_L_V2_in = new ROOT::RDataFrame(tree_name, fPath_L_V2.c_str()); 
        df_L_V3_in = new ROOT::RDataFrame(tree_name, fPath_L_V3.c_str()); 
    }
    catch (const std::invalid_argument &e) { 
        Error(__func__, "Error trying to open one of the V-wire files.\n what(): %s", e.what()); 
        return -1; 
    }

    //get a list of sieve_holes

    constexpr double NaN_double = std::numeric_limits<double>::quiet_NaN(); 
    
    struct SieveHoleAndWire_t {
        SieveHole hole; 
        Trajectory_t traj; 
        double stats=0.; 
        double rms_dxdz=NaN_double;
        double rms_dydz=NaN_double; 
    }; 
    std::vector<SieveHole> R_sieve_holes = ApexOptics::ConstructSieveHoles(true); 
    std::vector<SieveHole> L_sieve_holes = ApexOptics::ConstructSieveHoles(false); 

    //const double average_y_hcs = *mc_Vwire_nodes[0].Define("y_hcs", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}).Mean("y_hcs");

    //a list of sieve-holes for each wire
    vector<SieveHoleAndWire_t> R_wire_holes[] = { {}, {}, {} }; 
    vector<SieveHoleAndWire_t> L_wire_holes[] = { {}, {}, {} }; 

    //now, we compute the get ready to actually perform the measurements 

    //these are flags for use with the angle tester
    const auto kDxdz   = AngleRecoTester::kDxdz; 
    const auto kDydz   = AngleRecoTester::kDydz; 
    const auto kSlopes = AngleRecoTester::kSlopes; 
    
    //get dataframes with the info we need
    ROOT::RDF::RNode df_R_Vwire[] = {
        model(true, *df_R_V1_in),
        model(true, *df_R_V2_in),
        model(true, *df_R_V3_in)
    };
    ROOT::RDF::RNode df_L_Vwire[] = {
        model(false, *df_L_V1_in),
        model(false, *df_L_V2_in),
        model(false, *df_L_V3_in)
    };
    
    constexpr double range_R_dxdz[] = {-0.05, +0.05}; 
    constexpr double range_R_dydz[] = {-0.04, +0.04}; 
    constexpr double range_L_dxdz[] = {-0.05, +0.05}; 
    constexpr double range_L_dydz[] = {-0.04, +0.04}; 
    
    AngleFitResult_t angles_dxdz[3]; 
    AngleFitResult_t angles_dydz[3]; 


    //each of these values are for wires 1, 2, 3.
    //                             V1      V2      V3 
    int row_minx[]           = {    4,      4,      5   };  //min row to fit
    int row_miny[]           = {    4,      4,      5   };
    int row_maxx[]           = {   12,     12,     11   };  //max row to fit
    int row_maxy[]           = {   12,     12,     11   };
    int col_minx[]           = {    3,      1,      0   };
    int col_miny[]           = {    3,      1,      0   };
    int col_maxx[]           = {    7,      6,      3   };
    int col_maxy[]           = {    7,      6,      3   };
    double row_cut_widthx[]  = { 1.35,   1.35,   1.35   };
    double row_cut_widthy[]  = { 1.35,   1.35,   1.35   };
    double col_cut_widthx[]  = { 0.50,   0.50,   0.65   };
    double col_cut_widthy[]  = { 0.50,   0.50,   0.65   };
    double bg_cut_widthx[]   = { 1.00,   1.00,   1.00   };
    double bg_cut_widthy[]   = { 1.00,   1.00,   0.70   };

    constexpr double vdc_smearing_dxdz_sv = 0.426e-3; //rad
    constexpr double vdc_smearing_dydz_sv = 0.138e-3; //rad
    
    
    //this counter will keep track of the overall stat-weighting used for each hole
    double dxdz_total_statweight = 0.; 
    double dydz_total_statweight = 0.; 

    double dxdz_total_square_error = 0.; 
    double dydz_total_square_error = 0.; 

    const string wire_names[] = {"V1", "V2", "V3"}; 

    AngleFitResult_t 
        R_angleFits_dxdz[3],
        R_angleFits_dydz[3],
        L_angleFits_dxdz[3],
        L_angleFits_dydz[3]; 

    //for each sieve-hole which has a dx/dz AND dy/dz RMS measured, keep the sieve hole, along with the RMS for both angles. 
    vector<HoleAndAngles_t> R_hole_and_RMS[3], L_hole_and_RMS[3]; 

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
        auto& df_R_wire = df_R_Vwire[iwire]; 
        auto& df_L_wire = df_L_Vwire[iwire]; 

        //now, go thru each sieve-hole and designate a wire
        AngleRecoTester R_angle_tester(true, target, df_R_wire); 

        R_angle_tester.SetDrawing(false); 

        R_angle_tester.SetRange_dxdz(range_R_dxdz[0], range_R_dxdz[1]); 
        R_angle_tester.SetRange_dydz(range_R_dydz[0], range_R_dydz[1]); 

        R_angle_tester.SetName_dxdz("reco_dxdz_sv");
        R_angle_tester.SetName_dydz("reco_dydz_sv"); 

        AngleRecoTester L_angle_tester(false, target, df_L_wire); 

        L_angle_tester.SetDrawing(false); 

        L_angle_tester.SetRange_dxdz(range_L_dxdz[0], range_L_dxdz[1]); 
        L_angle_tester.SetRange_dydz(range_L_dydz[0], range_L_dydz[1]); 

        L_angle_tester.SetName_dxdz("reco_dxdz_sv");
        L_angle_tester.SetName_dydz("reco_dydz_sv"); 
        
        //have verbose fitting of holes (for now)
        //angle_tester.SetQuiet(false); 
        //angle_tester.SetFlag(AngleRecoTester::kVerbose_fitting); 

        cout << "------------ RHRS --------------------------------------------------\n" << flush; 

        R_angleFits_dxdz[iwire] = R_angle_tester.Measure(kDxdz,     //measure dxdz
            row_minx[iwire],row_maxx[iwire],               //min row, max row 
            col_minx[iwire],col_maxx[iwire],               //min col, max col 
            row_cut_widthx[iwire], col_cut_widthx[iwire], bg_cut_widthx[iwire] ).value(); //row cut width, col cut width, background cut width 
        
        R_angleFits_dydz[iwire] = R_angle_tester.Measure(kDydz | kSlopes,     //measure dydz
            row_miny[iwire],row_maxy[iwire],                        //min row, max row 
            col_miny[iwire],col_maxy[iwire],                        //min col, max col 
            row_cut_widthy[iwire], col_cut_widthy[iwire], bg_cut_widthy[iwire] ).value(); //row cut width, col cut width, background cut width         
            
        
        cout << "------------ LHRS --------------------------------------------------\n" << flush; 
        
        L_angleFits_dxdz[iwire] = L_angle_tester.Measure(kDxdz,     //measure dxdz
            row_minx[iwire],row_maxx[iwire],               //min row, max row 
            col_minx[iwire],col_maxx[iwire],               //min col, max col 
            row_cut_widthx[iwire], col_cut_widthx[iwire], bg_cut_widthx[iwire] ).value(); //row cut width, col cut width, background cut width 
        
        L_angleFits_dydz[iwire] = L_angle_tester.Measure(kDydz | kSlopes,     //measure dydz
            row_miny[iwire],row_maxy[iwire],                        //min row, max row 
            col_miny[iwire],col_maxy[iwire],                        //min col, max col 
            row_cut_widthy[iwire], col_cut_widthy[iwire], bg_cut_widthy[iwire] ).value(); //row cut width, col cut width, background cut width         

        //if (iwire!=2) continue;
            
        /// @brief Matches angles which have both dx/dz and dy/dz measured
        /// @param hole_and_RMS the vector of 'HoleAndAngles_t' structs to save the results in
        /// @param fits_dxdz fits for dx/dz
        /// @param fits_dydz fits for dy/dz 
        auto Match_angles = [](vector<HoleAndAngles_t>& hole_and_RMS, vector<AngleFit_t> fits_dxdz, vector<AngleFit_t> fits_dydz)
        {   
            hole_and_RMS.clear(); 

            for (const auto& fit_dxdz : fits_dxdz) {
                for (const auto& fit_dydz : fits_dydz) {
                    if (fit_dxdz.hole == fit_dydz.hole) {
                        hole_and_RMS.push_back({
                            .hole = fit_dxdz.hole, 
                            .dxdz = fit_dxdz.overall_RMS, 
                            .dydz = fit_dydz.overall_RMS
                        }); 
                    } 
                }
            }            
            return; 
        }; 
        //_____________________________________________________________________________________
            
        Match_angles(R_hole_and_RMS[iwire], R_angleFits_dxdz[iwire].fits, R_angleFits_dydz[iwire].fits); 
        Match_angles(L_hole_and_RMS[iwire], L_angleFits_dxdz[iwire].fits, L_angleFits_dydz[iwire].fits); 

        //for each angle, add contributions of vdc smearing (approximately!)
        for (auto& hRMS : R_hole_and_RMS[iwire]) {
            hRMS.dxdz = sqrt( hRMS.dxdz*hRMS.dxdz + vdc_smearing_dxdz_sv*vdc_smearing_dxdz_sv ); 
            hRMS.dydz = sqrt( hRMS.dydz*hRMS.dydz + vdc_smearing_dydz_sv*vdc_smearing_dydz_sv ); 
        }

        printf(
            " ~~ Number of holes fit:\n"
            "       R - (dxdz/dydz/both): %zi / %zi / %zi\n"
            "       L - (dxdz/dydz/both): %zi / %zi / %zi\n", 
            R_angleFits_dxdz[iwire].fits.size(), 
            R_angleFits_dydz[iwire].fits.size(), 
            R_hole_and_RMS[iwire].size(), 

            L_angleFits_dxdz[iwire].fits.size(), 
            L_angleFits_dydz[iwire].fits.size(),
            L_hole_and_RMS[iwire].size() 
        );

    }

    vector<double> pts_mass, pts_error; 

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
            hist_error = df_mc 
                
                .Define("v_wire_id", [V1_zcut,V3_zcut](TVector3& vtx)
                {
                    if (vtx.z() > V3_zcut) return 2; //V3 
                    if (vtx.z() < V1_zcut) return 0; //V1
                    return 1;                        //V2 
                }, {"R_position_vtx"})

                //find which hole each wire is closest to each wire (R) 
                .Define("R_bestfit_hole_id", [&R_hole_and_RMS, &targets, &find_bestfit_hole](int iwire, double dxdz, double dydz, const TVector3& vtx_hcs)
                {   
                    return find_bestfit_hole(true, targets[iwire], R_hole_and_RMS[iwire], dxdz, dydz, vtx_hcs); 
                }, {"v_wire_id", "R_dxdz_sv", "R_dydz_sv", "R_position_vtx"})

                //find which hole each wire is closest to each wire (L) 
                .Define("L_bestfit_hole_id", [&L_hole_and_RMS, &targets, &find_bestfit_hole](int iwire, double dxdz, double dydz, const TVector3& vtx_hcs)
                {   
                    return find_bestfit_hole(false, targets[iwire], L_hole_and_RMS[iwire], dxdz, dydz, vtx_hcs); 
                }, {"v_wire_id", "L_dxdz_sv", "L_dydz_sv", "L_position_vtx"}) 

                .Define("R_traj_hcs_smeared", [&rand, &R_hole_and_RMS](double x, double y, double dxdz, double dydz, double dpp, int iwire, int ihole)
                {
                    const auto& hRMS = R_hole_and_RMS[iwire].at(ihole); 

                    Trajectory_t traj_scs{
                        x, 
                        y, 
                        dxdz + rand.Gaus()*hRMS.dxdz, 
                        dydz + rand.Gaus()*hRMS.dydz,
                        dpp 
                    };

                    return ApexOptics::SCS_to_HCS(true, traj_scs); 

                }, {"R_x_sv", "R_y_sv", "R_dxdz_sv", "R_dydz_sv", "R_dpp_sv", "v_wire_id", "R_bestfit_hole_id"})

                .Define("L_traj_hcs_smeared", [&rand, &L_hole_and_RMS](double x, double y, double dxdz, double dydz, double dpp, int iwire, int ihole)
                {
                    const auto& hRMS = L_hole_and_RMS[iwire].at(ihole); 

                    Trajectory_t traj_scs{
                        x, 
                        y, 
                        dxdz + rand.Gaus()*hRMS.dxdz, 
                        dydz + rand.Gaus()*hRMS.dydz,
                        dpp  + rand.Gaus()*1e-4
                    };

                    return ApexOptics::SCS_to_HCS(false, traj_scs); 
                    
                }, {"L_x_sv", "L_y_sv", "L_dxdz_sv", "L_dydz_sv", "L_dpp_sv", "v_wire_id", "L_bestfit_hole_id"})

                .Define("m2_reco", [hrs_momentum](Trajectory_t R_traj, Trajectory_t L_traj)
                {
                    double R_p = hrs_momentum*(1. + R_traj.dpp); 
                    double L_p = hrs_momentum*(1. + L_traj.dpp); 
                    
                    TVector3 R_momentum( R_traj.dxdz, R_traj.dydz, 1. ); R_momentum = R_momentum.Unit() * R_p; 
                    TVector3 L_momentum( L_traj.dxdz, L_traj.dydz, 1. ); L_momentum = L_momentum.Unit() * L_p; 
                    
                    double E = R_p + L_p; 
                    double P2 = (R_momentum + L_momentum).Mag2(); 
                    //ignoring correction from electron/positron mass
                    return E*E - P2; 
                }, {"R_traj_hcs_smeared", "L_traj_hcs_smeared"})
                
                .Define("m_reco", [hrs_momentum](double m2){ return sqrt(m2); }, {"m2_reco"})

                .Define("invariant_mass_error", [&mA](double m_reco, double m)
                { 
                    mA = m; 
                    return m_reco - m; 
                }, {"m_reco", "invariant_mass"})

                .Histo1D<double>({"h_err", "Invariant Mass error.;reco-m_{A} - m_{A};", 200, -8., +8.}, "invariant_mass_error"); 
        
        } catch (const std::exception& e) {
            Error(__func__, "Something went wrong trying to define the A'-prime mass error hist.\n what(): %s\n", e.what()); 
            return -1; 
        }

        double rms = hist_error->GetStdDev(); 
        printf("A' mass: %.0f MeV,  RMS %.3f MeV\n", mA, rms); cout << flush; 
        
        pts_error.push_back( rms );
        pts_mass.push_back( mA ); 
        
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

    new TCanvas; 
    auto g = new TGraph(pts_mass.size(), pts_mass.data(), pts_error.data()); 
    g->SetTitle("m_{A} vs #sigma_{m_{A}} (MeV/c^{2})"); 
    g->Draw(); 

    g->SetMarkerStyle(kOpenCircle); 
    g->Draw("SAME P"); 

    /*/overall weighted error 
    printf("total RMS, mrad (dxdz/dydz): %.3f / %.3f\n", 
        1.e3*sqrt(dxdz_total_square_error/dxdz_total_statweight), 
        1.e3*sqrt(dydz_total_square_error/dydz_total_statweight)
    );*/ 
    
    return 0; 
}

#endif
