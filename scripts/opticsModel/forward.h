
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
class ForwardModel : public ModularOpticsModel {
private: 

    /// The optics model to be used (R or L arm)
    NPolyArray fParr;

    ArmMode::Bit fArmMode{ ArmMode::kNone }; 

    const std::vector<std::string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};


public: 

    ForwardModel(const char* path_polyarray); 
    ~ForwardModel() {}; 

    ROOT::RDF::RNode DefineOutputs(ROOT::RDF::RNode node_in) const; 
};


ForwardModel::ForwardModel(const char* path_polyarray)
{   
    
    //create right-arm model_____________________________________________________________________________________
    fParr = ApexOptics::Parse_NPolyArray_from_file(path_polyarray, {branches_sv}, 4); 
} 

ROOT::RDF::RNode ForwardModel::DefineOutputs(ROOT::RDF::RNode node_in) const 
{
    using namespace std; 
    using namespace ApexOptics; 

    RDFNodeAccumulator rna(node_in); 
           
    rna.Define("Xsv_reco", [this](double x, double y, double dxdz, double dydz)
    {
        return ApexOptics::RVec_to_Trajectory_t( fParr.Eval({x,y,dxdz,dydz}) ); 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 
    
    /*rna = Add_branches_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        "reco_x_sv", "reco_y_sv", "reco_dxdz_sv", "reco_dydz_sv", "reco_dpp_sv"
    });*/ 

    return Add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"reco_x_sv",    &Trajectory_t::x},
        {"reco_y_sv",    &Trajectory_t::y},
        {"reco_dxdz_sv", &Trajectory_t::dxdz},
        {"reco_dydz_sv", &Trajectory_t::dydz}, 
        {"reco_dpp_sv",  &Trajectory_t::dpp}
    }); 
}