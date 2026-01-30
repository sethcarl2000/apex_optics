#ifndef measure_model_accuracy_H
#define measure_model_accuracy_H

//APEX headers
#include <ApexOptics.h>
#include <AngleRecoTester.h>
#include "RDFNodeAccumulator.h"
//ROOT headers
#include <ROOT/RDataFrame.hxx>  
//stdlib headers
#include <vector>
#include <string> 
#include <limits> 
#include <functional> 
#include <stdexcept> 


/**
 * @brief measures & reports reconstruction accuracy for the left or right arm
 */
int measure_model_accuracy(
    const bool is_RHRS, 
    std::function<ROOT::RDF::RNode(bool,ROOT::RDF::RNode)> optics_model,
    const std::string path_V1="data/replay/real_L_V1.root",
    const std::string path_V2="data/replay/real_L_V2.root",
    const std::string path_V3="data/replay/real_L_V3.root",
    const std::string path_montecarlo="data/mc/mc_L_production.root",
    const std::string branch_dxdz="reco_dxdz_sv",
    const std::string branch_dydz="reco_dydz_sv"
) 
{        
    struct TargetAndPath_t { ApexOptics::OpticsTarget_t target; std::string path; };

    const char* tree_name="tracks_fp";

    //check to make sure we have the correct branches
    using namespace ApexOptics; 
    using namespace std; 

    auto is_nan = [](double x) { return x==x ? false : true; };

    //vertical wiress
    const OpticsTarget_t targets[] = {
        GetTarget("V1"),
        GetTarget("V2"),
        GetTarget("V3")
    };

    ROOT::RDataFrame *df_V1_in{nullptr}, *df_V2_in{nullptr}, *df_V3_in{nullptr}; 

    //try to construct the RDataFrame's 
    try { 
        df_V1_in = new ROOT::RDataFrame(tree_name, path_V1.c_str()); 
        df_V2_in = new ROOT::RDataFrame(tree_name, path_V2.c_str()); 
        df_V3_in = new ROOT::RDataFrame(tree_name, path_V3.c_str()); 
    }
    catch (const std::invalid_argument &e) { 
        Error(__func__, "Error trying to open one of the V-wire files.\n what(): %s", e.what()); 
        return -1; 
    }

    ROOT::RDataFrame *df_mc_in{nullptr}; 
    try {
        df_mc_in = new ROOT::RDataFrame(tree_name, path_montecarlo.c_str());
    } 
    catch (const std::invalid_argument &e) {
        Error(__func__, "Error trying to open the monte-carlo file.\n what(): %s", e.what()); 
        return -1; 
    }

    auto df_mc = *df_mc_in; delete df_mc_in; 

    //now, we compute the get ready to actually perform the measurements 

    //these are flags for use with the angle tester
    const auto dxdz   = AngleRecoTester::kDxdz; 
    const auto dydz   = AngleRecoTester::kDydz; 
    const auto slopes = AngleRecoTester::kSlopes; 
    
    //get dataframes with the info we need
    ROOT::RDF::RNode df_Vwire[] = {
        optics_model(is_RHRS, *df_V1_in),
        optics_model(is_RHRS, *df_V2_in),
        optics_model(is_RHRS, *df_V3_in)
    };
    
    constexpr double range_dxdz[] = {-0.05, +0.05}; 
    constexpr double range_dydz[] = {-0.04, +0.04}; 
    
    AngleFitResult_t angles_dxdz[3]; 
    AngleFitResult_t angles_dydz[3]; 

    //each of these values are for wires 1, 2, 3.
    //                            V1      V2      V3 
    int row_min[]           = {    4,      4,      4   };
    int row_max[]           = {   12,     12,     12   };
    int col_min[]           = {    1,      1,      1   };
    int col_max[]           = {    7,      8,      8   };
    double row_cut_width[]  = { 1.35,   1.35,   1.35   };
    double col_cut_width[]  = { 0.50,   0.50,   0.50   };
    double bg_cut_width[]   = { 1.00,   1.00,   1.00   };

    for (int iwire=0; iwire<3; iwire++) {

        AngleRecoTester angle_tester(is_RHRS, targets[iwire], df_Vwire[iwire]); 

        //angle_tester.SetDrawing(false); 

        angle_tester.SetRange_dxdz(range_dxdz[0], range_dxdz[1]); 
        angle_tester.SetRange_dydz(range_dydz[0], range_dydz[1]); 

        angle_tester.SetName_dxdz(branch_dxdz);
        angle_tester.SetName_dydz(branch_dydz); 

        if (iwire!=0) continue;

        angles_dxdz[iwire] = angle_tester.Measure(dxdz,              //measure dxdz
            row_min[iwire],row_max[iwire],               //min row, max row 
            col_min[iwire],col_max[iwire],                //min col, max col 
            row_cut_width[iwire], col_cut_width[iwire], bg_cut_width[iwire] ).value(); //row cut width, col cut width, background cut width 

        angles_dydz[iwire] = angle_tester.Measure(dydz | slopes,     //measure dydz
            row_min[iwire],row_max[iwire],               //min row, max row 
            col_min[iwire],col_max[iwire],                //min col, max col 
            row_cut_width[iwire], col_cut_width[iwire], bg_cut_width[iwire] ).value(); //row cut width, col cut width, background cut width 
    }
    
    return 0; 
}

#endif
