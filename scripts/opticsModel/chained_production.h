
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

    int fPivotElement=0; 

public: 

    struct TrajAndZ_t {
        ApexOptics::Trajectory_t R_Xsv, L_Xsv; 
        double z_hcs; 
    };

    ChainedModel(); 
    ~ChainedModel() {}; 

    ROOT::RDF::RNode DefineOutputs(ROOT::RDF::RNode node_in) const; 
};


ChainedModel::ChainedModel()
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

ROOT::RDF::RNode ChainedModel::DefineOutputs(ROOT::RDF::RNode node_in) const 
{
    using namespace std; 
    using namespace ApexOptics; 
    using RVecD = ROOT::RVec<double>; 

    RDFNodeAccumulator rna(node_in); 
           
    rna.DefineIfMissing("position_vtx", [](TVector3 vtx){ return vtx; }, {"L_position_vtx"}); 

    rna.Overwrite("R_Xfp", [](const RVecD& x, const RVecD& y, const RVecD& dxdz, const RVecD& dydz)
    {
        return Trajectory_t{x[0],y[0],dxdz[0],dydz[0]}; 
    }, {"R_x_fp","R_y_fp","R_dxdz_fp","R_dydz_fp"}); 

    rna.Overwrite("L_Xfp", [](const RVecD& x, const RVecD& y, const RVecD& dxdz, const RVecD& dydz)
    {
        return Trajectory_t{x[0],y[0],dxdz[0],dydz[0]}; 
    }, {"L_x_fp","L_y_fp","L_dxdz_fp","L_dydz_fp"}); 

    rna.Define("output", [this](const Trajectory_t& R_Xfp, const Trajectory_t& L_Xfp, const TVector3& vtx_hcs)
    {   
        const auto& R_rev = fModel_R->GetChain_rev(); 
        const auto& R_fwd = fModel_R->GetChain_fwd();
        
        const auto& L_rev = fModel_L->GetChain_rev(); 
        const auto& L_fwd = fModel_L->GetChain_fwd();
        
        auto get_dXsv = [this](const NPolyArrayChain* fwd, const ROOT::RVec<double>& Xsv_start, const ROOT::RVec<double>& Xfp)
        {                
            auto&& J = fwd->Jacobian(Xsv_start);

            RMatrix Ji(4,4); 
            ROOT::RVec<double> J_pivot; J_pivot.reserve(4);

            int j_mat=0; 
            for (int j=0; j<5; j++) {//loop over columns 

                if (j==fPivotElement) {
                    for (int i=0; i<4; i++) J_pivot.push_back( J.get(i,j) );
                }  else {
                    for (int i=0; i<4; i++) Ji.get(i,j_mat) = J.get(i,j); 
                    j_mat++; 
                }
            }

            Ji.Set_report_singular(false);

            auto dX = Ji.Solve( J_pivot*(-1.) ); 
            
            //check for NaN / invalid result 
            if (dX.size() != 4)                return Trajectory_t{std::numeric_limits<double>::quiet_NaN()}; 
            for (double& x : dX) { if (x != x) return Trajectory_t{std::numeric_limits<double>::quiet_NaN()}; } 

            dX.insert(dX.begin()+fPivotElement, 1. ); 

            return ApexOptics::RVec_to_Trajectory_t(dX); 
        }; 

        //create 'first guess' 
        auto R_Xfp_rvec = ApexOptics::Trajectory_t_to_RVec(R_Xfp);
        auto R_Xsv_fg = R_rev->Eval(R_Xfp_rvec);
        auto R_dXsv   = get_dXsv(R_fwd, R_Xsv_fg, R_Xfp_rvec);

        auto R_Xsv_1 = SCS_to_HCS( true, RVec_to_Trajectory_t(R_Xsv_fg) ); 
        auto R_Xsv_2 = R_Xsv_1 + SCS_to_HCS( true, R_dXsv * 1e-3 ); 

        //create 'first guess' 
        auto L_Xfp_rvec = ApexOptics::Trajectory_t_to_RVec(L_Xfp);
        auto L_Xsv_fg = L_rev->Eval(L_Xfp_rvec);
        auto L_dXsv   = get_dXsv(L_fwd, L_Xsv_fg, L_Xfp_rvec);

        auto L_Xsv_1 = SCS_to_HCS( false, RVec_to_Trajectory_t(L_Xsv_fg) );
        auto L_Xsv_2 = L_Xsv_1 + SCS_to_HCS( false, L_dXsv * 1e-3 );
        
        auto get_z_x = [&vtx_hcs](Trajectory_t& X, double& z, double& x) { 
            z = (vtx_hcs.y() - X.y)/X.dydz; 
            x = X.x + z*X.dxdz; 
        };
        
        double Rx1, Rz1, Rx2, Rz2, Rm, Rb; 
        
        get_z_x(R_Xsv_1, Rx1, Rz1); 
        get_z_x(R_Xsv_2, Rx2, Rz2); 
        
        Rm = (Rx2 - Rx1)/(Rz2 - Rz1); 
        Rb = Rx1 - Rm*Rz1;  

        double Lx1, Lz1, Lx2, Lz2, Lm, Lb; 

        get_z_x(L_Xsv_1, Lx1, Lz1); 
        get_z_x(L_Xsv_2, Lx2, Lz2); 
        
        Lm = (Lx2 - Lx1)/(Lz2 - Lz1); 
        Lb = Lx1 - Lm*Lz1; 

        //their z-intercept 
        double z = (Lb - Rb)/(Rm - Lm); 

        auto R_Xsv = R_Xsv_1 + (R_Xsv_2 - R_Xsv_1)*((z - Rz1)/(Rz2 - Rz1)); 
        auto L_Xsv = L_Xsv_1 + (L_Xsv_2 - L_Xsv_1)*((z - Lz1)/(Lz2 - Lz1)); 
        
        R_Xsv = HCS_to_SCS( true, R_Xsv );
        L_Xsv = HCS_to_SCS( false, L_Xsv );

        return std::pair<Trajectory_t, Trajectory_t>{ R_Xsv, L_Xsv };

    }, {"R_Xfp", "L_Xfp", "L_position_vtx"});

    rna.Define("R_Xsv_reco", [](const std::pair<Trajectory_t,Trajectory_t>& output)
    {
        return output.first; 
    }, {"output"});

    rna.Define("reco_position_vtx", [](const Trajectory_t& R_Xsv, const TVector3& vtx)
    {   
        auto X = SCS_to_HCS(true, R_Xsv); 
        double z_hcs = (vtx.y() - X.y)/X.dydz; 
        double x_hcs = vtx.x() + z_hcs*X.dxdz; 
        return TVector3( x_hcs, vtx.y(), z_hcs ); 
    }, {"R_Xsv_reco", "L_position_vtx"}); 
    
    rna.Define("L_Xsv_reco", [](const std::pair<Trajectory_t,Trajectory_t>& output)
    {
        return output.second; 
    }, {"output"});

    rna = Add_branches_from_Trajectory_t(rna.Get(), "R_Xsv_reco", {
        "R_x_sv", "R_y_sv", "R_dxdz_sv", "R_dydz_sv", "R_dpp_sv"
    }); 
    rna = Add_branches_from_Trajectory_t(rna.Get(), "L_Xsv_reco", {
        "L_x_sv", "L_y_sv", "L_dxdz_sv", "L_dydz_sv", "L_dpp_sv"
    }); 

    return rna.Get(); 
}