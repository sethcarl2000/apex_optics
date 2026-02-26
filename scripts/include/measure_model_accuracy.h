#ifndef measure_model_accuracy_H
#define measure_model_accuracy_H

//APEX headers
#include <ApexOptics.h>
#include <AngleRecoTester.h>
#include "RDFNodeAccumulator.h"
#include <SieveHole.h> 
#include <ArmMode.h>
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

    T& operator() (ArmMode::Bit flag) { 
        if (flag & ArmMode::kRHRS) { return R; } else { return L; }
    }
};

//give this a general function, and it will perform some operation on both arms 


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

    std::function<ROOT::RDF::RNode(ArmMode::Bit, ROOT::RDF::RNode)> fModel; 

    int EvaluateModel ( 
        std::function<ROOT::RDF::RNode(ArmMode::Bit,ROOT::RDF::RNode)> optics_model
    ); 

    int Evaluate_singleMass(const std::string& path_mass); 
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

    BothArmValue_t<const double[2]> range_dxdz = { {-0.05, +0.05}, {-0.05, +0.05} }; 
    BothArmValue_t<const double[2]> range_dydz = { {-0.04, +0.04}, {-0.04, +0.04} }; 
    
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

    constexpr double vdc_smearing_dxdz_sv = 0.426e-3; //rad
    constexpr double vdc_smearing_dydz_sv = 0.138e-3; //rad
    
    
    //this counter will keep track of the overall stat-weighting used for each hole
    double dxdz_total_statweight = 0.; 
    double dydz_total_statweight = 0.; 

    double dxdz_total_square_error = 0.; 
    double dydz_total_square_error = 0.; 

    const string wire_names[] = {"V1", "V2", "V3"}; 

    BothArmValue_t<AngleFitResult_t>  
        angleFits_dxdz[3],
        angleFits_dydz[3]; 
    
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
            
            bool is_RHRS = (ArmMode::kRHRS & arm); 
            AngleRecoTester angle_tester(is_RHRS, target, df_wire(arm)); 

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
        
            //now, collect a list of holes for which we have measured **both** dxdz and dydz. 
            hole_and_RMS[iwire](arm).clear(); 

            for (const auto& fit_dxdz : angleFits_dxdz[iwire](arm).fits) {
                for (const auto& fit_dydz : angleFits_dydz[iwire](arm).fits) {
                    if (fit_dxdz.hole == fit_dydz.hole) {
                        hole_and_RMS[iwire](arm).push_back({
                            .hole = fit_dxdz.hole, 
                            .dxdz = fit_dxdz.overall_RMS, 
                            .dydz = fit_dydz.overall_RMS
                        }); 
                    } 
                }
            }    
                        
            //now, add an estimate for VDC smearing
            for (auto& hRMS : hole_and_RMS[iwire](arm)) {
                hRMS.dxdz = sqrt( hRMS.dxdz*hRMS.dxdz + vdc_smearing_dxdz_sv*vdc_smearing_dxdz_sv ); 
                hRMS.dydz = sqrt( hRMS.dydz*hRMS.dydz + vdc_smearing_dydz_sv*vdc_smearing_dydz_sv ); 
            }

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

    vector<double> pts_mass, pts_error; 

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

                rna.Define(arm_prefix+"traj_hcs_smeared", [&rand, &hole_and_RMS, &targets, arm, arm_bool, &find_bestfit_hole, V1_zcut, V3_zcut](
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
            }

            rna.Define("m2_reco", [hrs_momentum](Trajectory_t R_traj, Trajectory_t L_traj)
            {
                double R_p = hrs_momentum*(1. + R_traj.dpp); 
                double L_p = hrs_momentum*(1. + L_traj.dpp); 
                
                TVector3 R_momentum( R_traj.dxdz, R_traj.dydz, 1. ); R_momentum = R_momentum.Unit() * R_p; 
                TVector3 L_momentum( L_traj.dxdz, L_traj.dydz, 1. ); L_momentum = L_momentum.Unit() * L_p; 
                
                double E = R_p + L_p; 
                double P2 = (R_momentum + L_momentum).Mag2(); 
                //ignoring correction from electron/positron mass
                return E*E - P2; 
            }, {"R_traj_hcs_smeared", "L_traj_hcs_smeared"}); 
                
            rna.Define("m_reco", [hrs_momentum](double m2){ return sqrt(m2); }, {"m2_reco"}); 

            rna.Define("invariant_mass_error", [&mA](double m_reco, double m)
            { 
                mA = m; 
                return m_reco - m; 
            }, {"m_reco", "invariant_mass"}); 

            hist_error = rna.Get().Histo1D<double>(
                {"h_err", "Invariant Mass error.;reco-m_{A} - m_{A};", 200, -8., +8.}, 
                "invariant_mass_error"
            ); 
        
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
