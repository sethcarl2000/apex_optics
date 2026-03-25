
//APEX headers
#include "../include/RDFNodeAccumulator.h"
#include "../include/TestAngleReco.h"
#include "../include/ChainedOpticsModel.h"
#include "../include/Add_branch_from_Trajectory_t.h"
#include "ModularOpticsModel.h"
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
#include <stdexcept> 

/// @brief Implement the 'chained' optics model, and then submit it to 'measure_model_accuracy.h' to be tested
class ChainedModel : public ModularOpticsModel {
private: 

    /// The optics model to be used (R or L arm)
    ChainedOpticsModel *fModel_R, *fModel_L; 

    const std::vector<std::string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    const std::vector<std::string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const std::vector<std::string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

    const std::vector<std::string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
    const std::vector<std::string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

    const std::vector<std::string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
    const std::vector<std::string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};

    ArmMode::Bit fArmMode{ ArmMode::kNone }; 

public: 

    ChainedModel(ArmMode::Bit arm_mode=ArmMode::kBoth); 
    ~ChainedModel() {}; 

    ROOT::RDF::RNode DefineOutputs(ROOT::RDF::RNode node_in) const; 
};


ChainedModel::ChainedModel(ArmMode::Bit arm_mode)
    : fArmMode{arm_mode}
{   
    
    //create right-arm model_____________________________________________________________________________________
    fModel_R = new ChainedOpticsModel(true); 

    fModel_R->CreateChainRev({

        // sv <= [Poly] <= fp
        {"data/poly/fits_30Dec/V123-v01_fp_sv_R_4ord.dat", branches_sv, 4} 
    }); 

    fModel_R->CreateChainFwd({
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_R_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_30Dec/V123_fp-fwd_fp_R_3ord.dat", branches_fp, 4} 
    }); 

    //create left-arm model______________________________________________________________________________________
    fModel_L = new ChainedOpticsModel(false); 

    fModel_L->CreateChainRev({

        // sv <= [Poly] <= fp
        {"data/poly/fits_21Dec/V123_fp_sv_L_4ord.dat", branches_sv, 4} 
    }); 

    fModel_L->CreateChainFwd({
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_L_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_21Dec/V123_fp-fwd_fp_L_3ord.dat", branches_fp, 4} 
    }); 
} 

ROOT::RDF::RNode ChainedModel::DefineOutputs(ROOT::RDF::RNode node_in) const 
{
    using namespace std; 
    using namespace ApexOptics; 

    RDFNodeAccumulator rna(node_in); 
           
    rna.DefineIfMissing("position_vtx", [](TVector3 vtx){ return vtx; }, {"R_position_vtx"}); 

    switch (fArmMode) {
        case ArmMode::kBoth : {
            rna.Overwrite("R_Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, {"R_x_fp","R_y_fp","R_dxdz_fp","R_dydz_fp"}); 
        
            rna.Overwrite("L_Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, {"L_x_fp","L_y_fp","L_dxdz_fp","L_dydz_fp"}); 

            rna = Add_branches_from_Trajectory_t(rna.Get(), "R_Xsv_reco", {
                "R_x_sv", "R_y_sv", "R_dxdz_sv", "R_dydz_sv", "R_dpp_sv"
            }); 
            rna = Add_branches_from_Trajectory_t(rna.Get(), "L_Xsv_reco", {
                "L_x_sv", "L_y_sv", "L_dxdz_sv", "L_dydz_sv", "L_dpp_sv"
            }); 
            break; 
        }
        case ArmMode::kRHRS : {
            rna.Overwrite("Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}); 
            rna.Define("Xsv_reco", [this](const Trajectory_t& Xfp, const TVector3& vtx_hcs)
            {
                return fModel_R->Compute_Xsv(Xfp, vtx_hcs); 
            }, {"Xfp", "position_vtx"}); 
            rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_reco", {
                "reco_x_sv", "reco_y_sv", "reco_dxdz_sv", "reco_dydz_sv", "reco_dpp_sv"
            });
            break;  
        }
        case ArmMode::kLHRS : {
            rna.Overwrite("Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}); 
            rna.Define("Xsv_reco", [this](const Trajectory_t& Xfp, const TVector3& vtx_hcs)
            {
                return fModel_L->Compute_Xsv(Xfp, vtx_hcs); 
            }, {"Xfp", "position_vtx"}); 
            rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_reco", {
                "reco_x_sv", "reco_y_sv", "reco_dxdz_sv", "reco_dydz_sv", "reco_dpp_sv"
            });
            break;  
        }
        default : {
            throw invalid_argument("in <ChainedModel::DefineOutputs>: Model's arm mode specification is not valid; must be kLHRS, kRHRS, or kBoth.");
            break;  
        }
    }

    return rna.Get(); 
}