
//APEX headers
#include "../include/RDFNodeAccumulator.h"
#include "../include/TestAngleReco.h"
#include "../include/ChainedOpticsModel.h"
#include "../include/Add_branch_from_Trajectory_t.h"
#include "ModularOpticsModel.h"
#include <ApexOptics.h> 
#include <ArmMode.h> 
#include <RMatrix.h> 
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
class FwdProdModel : public ModularOpticsModel {
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

    int fPivotElement=0; 

public: 

    struct TrajAndZ_t {
        ApexOptics::Trajectory_t R_Xsv, L_Xsv; 
        double z_hcs; 
    };

    FwdProdModel(); 
    ~FwdProdModel() {}; 

    ROOT::RDF::RNode DefineOutputs(ROOT::RDF::RNode node_in) const; 
};


FwdProdModel::FwdProdModel()
{   
    //create right-arm model_____________________________________________________________________________________
    fModel_R = new ChainedOpticsModel(true); 

    fModel_R->CreateChainRev({

        // sv <= [Poly] <= fp
        {"data/poly/fits_29Mar/V123_fp_sv_R_4ord.dat", branches_sv, 4} 
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
        {"data/poly/fits_22Mar/V123_fp_sv_L_4ord.dat", branches_sv, 4} 
    }); 

    fModel_L->CreateChainFwd({
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_L_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_21Dec/V123_fp-fwd_fp_L_3ord.dat", branches_fp, 4} 
    }); 
} 

ROOT::RDF::RNode FwdProdModel::DefineOutputs(ROOT::RDF::RNode node_in) const 
{
    using namespace std; 
    using namespace ApexOptics; 
    using RVecD = ROOT::RVec<double>; 

    RDFNodeAccumulator rna(node_in); 
           
    rna.DefineIfMissing("position_vtx", [](TVector3 vtx){ return vtx; }, {"L_position_vtx"}); 

    rna.Overwrite("R_Xfp", [](const RVecD& x, const RVecD& y, const RVecD& dxdz, const RVecD& dydz)
    {
        return Trajectory_t{x[0],y[0],dxdz[0]-x[0]/6.,dydz[0]}; 
    }, {"R_x_fp","R_y_fp","R_dxdz_fp","R_dydz_fp"}); 

    rna.Overwrite("L_Xfp", [](const RVecD& x, const RVecD& y, const RVecD& dxdz, const RVecD& dydz)
    {
        return Trajectory_t{x[0],y[0],dxdz[0]-x[0]/6.,dydz[0]}; 
    }, {"L_x_fp","L_y_fp","L_dxdz_fp","L_dydz_fp"}); 

    rna.Define("R_Xsv_reco", [this](const Trajectory_t& Xfp)
    {
        return fModel_R->Compute_Xsv_first_guess(Xfp); 
    }, {"R_Xfp"});

    rna.Define("L_Xsv_reco", [this](const Trajectory_t& Xfp)
    {
        return fModel_L->Compute_Xsv_first_guess(Xfp); 
    }, {"L_Xfp"});

    rna.Define("reco_position_vtx", [](const Trajectory_t& R_Xsv_scs, const Trajectory_t& L_Xsv_scs)
    {   
        auto R_Xsv = SCS_to_HCS(true,  R_Xsv_scs); 
        auto L_Xsv = SCS_to_HCS(false, L_Xsv_scs); 
        
        TVector3 s1( R_Xsv.dxdz, R_Xsv.dydz, 1. );         
        TVector3 s2( L_Xsv.dxdz, L_Xsv.dydz, 1. ); 

        TVector3 r1( R_Xsv.x, R_Xsv.y, 0. ); 
        TVector3 r2( L_Xsv.x, L_Xsv.y, 0. ); 

        TVector3 R = r1 - r2;

        //dot products
        double s1s1 = s1.Mag2(); 
        double s1s2 = s1 * s2; 
        double s2s2 = s2.Mag2(); 

        RVecD B{ -R*s1, R*s2 }; 

        double det = s1s1*s2s2 - (s1s2*s1s2); 

        RMatrix A(2,2, {
            s2s2/det, s1s2/det,
            s1s2/det, s1s1/det  
        });

        RVecD t = A*B; 

        TVector3 closest_approach_R = s1*t[0] + r1; 
        TVector3 closest_approach_L = s2*t[1] + r2; 

        return 0.5*(closest_approach_L + closest_approach_R); 

    }, {"R_Xsv_reco", "L_Xsv_reco"}); 

    rna = Add_branches_from_Trajectory_t(rna.Get(), "R_Xsv_reco", {
        "R_x_sv", "R_y_sv", "R_dxdz_sv", "R_dydz_sv", "R_dpp_sv"
    }); 
    rna = Add_branches_from_Trajectory_t(rna.Get(), "L_Xsv_reco", {
        "L_x_sv", "L_y_sv", "L_dxdz_sv", "L_dydz_sv", "L_dpp_sv"
    }); 

    return rna.Get(); 
}