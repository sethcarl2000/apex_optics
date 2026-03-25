
//APEX headers
#include "../include/RDFNodeAccumulator.h"
#include "../include/TestAngleReco.h"
#include "../include/ChainedOpticsModel.h"
#include "../include/Add_branch_from_Trajectory_t.h"
#include "ModularOpticsModel.h"
#include <ApexOptics.h>
#include <NPolyArray.h>
#include <NPoly.h> 
#include <ArmMode.h> 
//ROOT headers
#include <TVector3.h> 
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <RMatrix.h>
//stdlib headers
#include <functional> 
#include <vector>
#include <string> 
#include <map> 
#include <stdexcept> 
#include <cmath> 

/// @brief Implement the 'chained' optics model, and then submit it to 'measure_model_accuracy.h' to be tested
class PiecewiseModel : public ModularOpticsModel {
private: 

    /// The optics model to be used (R or L arm)
    NPoly fPoly_zhcs; 
    NPolyArray fParr_v1, fParr_v2, fParr_v3;

    ArmMode::Bit fArmMode{ ArmMode::kNone }; 

    const std::vector<std::string> branches_sv{"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"};

    /// positions of the 3 v-wires
    const double fZ_wire[3] = { -196.20e-3,   +3.80e-3, +203.80e-3 }; 

    /// matrix which helps us determine the coefficients of the polynomial 
    RMatrix fCoeff_matrix; 

public: 

    std::string fPath_poly_v1{"data/poly/fits_22Mar/V1_fp_sv_L_2ord.dat"}; 
    std::string fPath_poly_v2{"data/poly/fits_22Mar/V2_fp_sv_L_2ord.dat"}; 
    std::string fPath_poly_v3{"data/poly/fits_22Mar/V3_fp_sv_L_2ord.dat"}; 
    
    PiecewiseModel(const char* path_zhcs_poly="data/poly/fits_22Mar/V123_fp_z-hcs_L_5ord.dat"); 
    ~PiecewiseModel() {}; 

    ROOT::RDF::RNode DefineOutputs(ROOT::RDF::RNode node_in) const; 
};


PiecewiseModel::PiecewiseModel(const char* path_zhcs_poly)
{   
    
    //parse z-hcs polynomial 
    fPoly_zhcs = NPoly(4); 
    ApexOptics::Parse_NPoly_from_file(path_zhcs_poly, "z_hcs", &fPoly_zhcs); 

    
    //this is a very simple interpolation 
    double sum_z{0.}, sum_z2{0.}, sum_z3{0.}, sum_z4{0.}; 
    for (int i=0; i<3; i++) {
        sum_z  += fZ_wire[i]; 
        sum_z2 += std::pow( fZ_wire[i], 2 ); 
        sum_z3 += std::pow( fZ_wire[i], 3 ); 
        sum_z4 += std::pow( fZ_wire[i], 4 ); 
    }
    fCoeff_matrix = RMatrix(3,3, { 
        sum_z4, sum_z3, sum_z2, 
        sum_z3, sum_z2, sum_z, 
        sum_z2, sum_z,  3. 
    });

    //parse arrays for each wire
    fParr_v1 = ApexOptics::Parse_NPolyArray_from_file(fPath_poly_v1.data(), branches_sv, 4); 
    fParr_v2 = ApexOptics::Parse_NPolyArray_from_file(fPath_poly_v2.data(), branches_sv, 4); 
    fParr_v3 = ApexOptics::Parse_NPolyArray_from_file(fPath_poly_v3.data(), branches_sv, 4); 
} 

ROOT::RDF::RNode PiecewiseModel::DefineOutputs(ROOT::RDF::RNode node_in) const 
{
    using namespace std; 
    using namespace ApexOptics; 

    RDFNodeAccumulator rna(node_in); 
           
    rna.Define("z_hcs", [this](double x, double y, double dxdz, double dydz)
    {
        return fPoly_zhcs.Eval({x,y,dxdz,dydz}); 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 

    rna.Define("Xfp", [](double x, double y, double dxdz, double dydz)
    {
        return Trajectory_t{x, y, dxdz, dydz}; 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"});

    rna.Define("Xsv_reco", [this](const Trajectory_t& Xfp, double z_hcs)
    {
        //define a reconstruction for each arm 
        ROOT::RVec<double> X_wire[3]; 
        X_wire[0] = fParr_v1.Eval(ApexOptics::Trajectory_t_to_RVec(Xfp)); 
        X_wire[1] = fParr_v2.Eval(ApexOptics::Trajectory_t_to_RVec(Xfp)); 
        X_wire[2] = fParr_v3.Eval(ApexOptics::Trajectory_t_to_RVec(Xfp)); 

        ROOT::RVec<double> X_sv; X_sv.reserve(5); 

        for (int i=0; i<5; i++) {

            ROOT::RVec<double> B{0., 0., 0.}; 

            for (int iwire=0; iwire<3; iwire++) {

                double z_wire = fZ_wire[iwire]; 

                double y = X_wire[iwire][i];
        
                B[0] += y * z_wire * z_wire; 
                B[1] += y * z_wire; 
                B[2] += y; 
            }
            
            auto coeff = fCoeff_matrix.Solve( B ); 

            //now, our answer should be a 2nd-order polynomial in 
            X_sv.push_back(
                coeff[0] * (z_hcs*z_hcs) +
                coeff[1] * (z_hcs) +
                coeff[2]
            ); 
        }

        return ApexOptics::RVec_to_Trajectory_t(X_sv); 

    }, {"Xfp", "z_hcs"}); 

    return Add_branch_from_Trajectory_t(rna.Get(), "Xsv_reco", {
        {"reco_x_sv",    &Trajectory_t::x},
        {"reco_y_sv",    &Trajectory_t::y},
        {"reco_dxdz_sv", &Trajectory_t::dxdz},
        {"reco_dydz_sv", &Trajectory_t::dydz}, 
        {"reco_dpp_sv",  &Trajectory_t::dpp}
    }); 
}