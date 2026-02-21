
//APEX headers
#include "include/RDFNodeAccumulator.h"
#include "include/TestAngleReco.h"
#include "include/ChainedOpticsModel.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include "include/measure_model_accuracy.h"
#include "include/Get_TParameter_from_TFile.h"
#include <ApexOptics.h> 
#include <ArmMode.h> 
//ROOT headers
#include <TVector3.h> 
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
//stdlib headers
#include <functional> 
#include <vector>
#include <string> 
#include <map> 

/// @brief Implement the 'chained' optics model, and then submit it to 'measure_model_accuracy.h' to be tested
/// @param arm valid options: 'RHRS', 'LHRS' or 'both'  
int opticsModel_chained(const std::string& arm) 
{
    using namespace std; 

    const map<string,ArmMode::Bit> arm_options{
        {"both", ArmMode::kBoth}, 
        {"RHRS", ArmMode::kRHRS}, 
        {"LHRS", ArmMode::kLHRS}
    }; 

    auto it_mode = arm_options.find(arm); 
    if (it_mode == arm_options.end()) {
        Error(__func__, "Invalid arm option given: '%s'. Valid options are: 'RHRS', 'LHRS', 'both'", arm.c_str());
        return -1;  
    }
    const auto mode = it_mode->second; 
    
    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

    const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
    const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

    const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
    const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};

    //create right-arm model_____________________________________________________________________________________
    ChainedOpticsModel* model_R = new ChainedOpticsModel(true); 

    model_R->CreateChainRev({

        // sv <= [Poly] <= fp
        {"data/poly/fits_30Dec/V123-v01_fp_sv_R_4ord.dat", branches_sv, 4} 
    }); 

    model_R->CreateChainFwd({
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_R_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_30Dec/V123_fp-fwd_fp_R_3ord.dat", branches_fp, 4} 
    }); 

    //create left-arm model______________________________________________________________________________________
    ChainedOpticsModel* model_L = new ChainedOpticsModel(false); 

    model_L->CreateChainRev({

        // sv <= [Poly] <= fp
        {"data/poly/fits_21Dec/V123_fp_sv_L_4ord.dat", branches_sv, 4} 
    }); 

    model_L->CreateChainFwd({
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_L_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_21Dec/V123_fp-fwd_fp_L_3ord.dat", branches_fp, 4} 
    }); 

    //const ptr to model, to make sure we don't modify it. 
    const ChainedOpticsModel* const_model_R = model_R; 
    const ChainedOpticsModel* const_model_L = model_L; 

    //___________________________________________________________________________________________________________
    auto optics_model = [const_model_L,const_model_R,&branches_fp,mode](ArmMode::Bit arm_mode, ROOT::RDF::RNode df) 
    {
        using namespace ApexOptics; 
        RDFNodeAccumulator rna(df); 
           
        rna.DefineIfMissing("position_vtx", [](TVector3 vtx){ return vtx; }, {"R_position_vtx"}); 

        if (arm_mode & ArmMode::kRHRS) { 

            rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
                {
                    return Trajectory_t{x,y,dxdz,dydz}; 
                }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}); 

            rna.Define("Xsv_reco", [const_model_R](const Trajectory_t& Xfp, const TVector3& vtx_hcs)
                {
                    return const_model_R->Compute_Xsv(Xfp, vtx_hcs); 
                }, {"Xfp", "position_vtx"}); 

            rna.Define("reco_dxdz_sv", [](const Trajectory_t& Xsv){ return Xsv.dxdz; }, {"Xsv_reco"}); 
            rna.Define("reco_dydz_sv", [](const Trajectory_t& Xsv){ return Xsv.dydz; }, {"Xsv_reco"}); 
        }
        
        if (arm_mode & ArmMode::kLHRS) {
            
            rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
                {
                    return Trajectory_t{x,y,dxdz,dydz}; 
                }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}); 

            rna.Define("Xsv_reco", [const_model_L](const Trajectory_t& Xfp, const TVector3& vtx_hcs)
                {
                    return const_model_L->Compute_Xsv(Xfp, vtx_hcs); 
                }, {"Xfp", "position_vtx"}); 

            rna.Define("reco_dxdz_sv", [](const Trajectory_t& Xsv){ return Xsv.dxdz; }, {"Xsv_reco"}); 
            rna.Define("reco_dydz_sv", [](const Trajectory_t& Xsv){ return Xsv.dydz; }, {"Xsv_reco"}); 
        }

        return rna.Get(); 
    };
    //___________________________________________________________________________________________________________

    MeasureModelAccuracy model_tester; 

    //now, measure the angles
    int ret_code = model_tester.EvaluateModel(optics_model); 
    
    delete model_R; delete model_L; 

    if (ret_code != 0) {
        Error(__func__, "Something went wrong trying to measure the accuracy. ret code: %i", ret_code); 
        return -1; 
    }
    return 0; 
}