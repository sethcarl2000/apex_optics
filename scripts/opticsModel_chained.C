
//APEX headers
#include "include/RDFNodeAccumulator.h"
#include "include/TestAngleReco.h"
#include "include/ChainedOpticsModel.h"
#include "include/Add_branch_from_Trajectory_t.h"
#include "include/measure_model_accuracy.h"
#include "include/Get_TParameter_from_TFile.h"
#include <ApexOptics.h> 
//ROOT headers
#include <TVector3.h> 
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
//stdlib headers
#include <functional> 
#include <vector>
#include <string> 

int opticsModel_chained(bool _is_RHRS) 
{
    using namespace std; 

    const vector<string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};
    const vector<string> branches_q1{"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"};
    const vector<string> branches_fp{"x_fp","y_fp","dxdz_fp","dydz_fp"};

    const vector<string> branches_fwd_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};
    const vector<string> branches_fwd_fp{"fwd_x_fp","fwd_y_fp","fwd_dxdz_fp","fwd_dydz_fp"};

    const vector<string> branches_rev_sv{"fwd_x_sv","fwd_y_sv","fwd_dxdz_sv","fwd_dydz_sv","fwd_dpp_sv"};
    const vector<string> branches_rev_q1{"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"};

    ChainedOpticsModel* model = new ChainedOpticsModel(_is_RHRS); 

    model->CreateChainRev({

        // sv <= [Poly] <= fp
        {"data/poly/fits_21Dec/V123_fp_sv_L_4ord.dat", branches_sv, 4} 
    }); 

    model->CreateChainFwd({
    
        // sv => [Poly] => fp-fwd => _Poly_ => fp 
        {"data/poly/mc_sv_fp_L_4ord.dat",            branches_fp, 5},
        {"data/poly/fits_21Dec/V123_fp-fwd_fp_L_3ord.dat", branches_fp, 4} 
    }); 

    //const ptr to model, to make sure we don't modify it. 
    const ChainedOpticsModel* const_model = model; 

    //___________________________________________________________________________________________________________
    auto optics_model = [const_model,&branches_fp](bool is_RHRS, ROOT::RDF::RNode df) 
    {
        using namespace ApexOptics; 
        auto df_out = df
            
            .Define("Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, branches_fp) 

            .Define("Xsv_reco", [const_model](const Trajectory_t& Xfp, const TVector3& vtx_hcs)
            {
                return const_model->Compute_Xsv(Xfp, vtx_hcs); 
            }, {"Xfp", "position_vtx"})

            .Define("reco_dxdz_sv", [](const Trajectory_t& Xsv){ return Xsv.dxdz; }, {"Xsv_reco"})
            .Define("reco_dydz_sv", [](const Trajectory_t& Xsv){ return Xsv.dydz; }, {"Xsv_reco"}); 

        return df_out; 
    };
    //___________________________________________________________________________________________________________

    //now, measure the angles
    int ret_code = measure_model_accuracy(
        _is_RHRS, 
        optics_model
    ); 
    
    delete model; 

    if (ret_code != 0) {
        Error(__func__, "Something went wrong trying to measure the accuracy. ret code: %i", ret_code); 
        return -1; 
    }
    return 0; 
}