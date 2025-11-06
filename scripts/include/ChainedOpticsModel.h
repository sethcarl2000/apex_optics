#ifndef ChainedOpticsModel_H
#define ChainedOpticsModel_H

//std headers
#include <vector>
#include <string> 
#include <limits> 

//apex-specific headers
#include <RMatrix.h>
#include <NPolyArrayChain.h>
#include <NPolyArray.h> 
#include <NPoly.h> 
#include <ApexOptics.h> 

//ROOT headers
#include <ROOT/RVec.hxx>

using ApexOptics::Trajectory_t; 

class ChainedOpticsModel {
private: 
    bool f_is_RHRS; 

    NPolyArrayChain fChainRev;
    NPolyArrayChain fChainFwd;  

    //can be from [0-4], representing which element of the vector we 'pivot' around. 
    int fPivotElement{0}; 

public:

    struct ChainConstructor_t {
        std::string path{}; 
        std::vector<std::string> coords{}; 
        int input_DoF{0}; 
    };

    ChainedOpticsModel(bool is_RHRS) : f_is_RHRS{is_RHRS} {/*noop*/}; 
    ~ChainedOpticsModel() {/*noop*/}; 

    //initialize the 'forward' chain
    void CreateChainFwd(const std::vector<ChainConstructor_t>& chain) {
        
        for (const auto& path_and_coord : chain) {

            const char* path  = path_and_coord.path.c_str(); 
            const auto  coord = path_and_coord.coords;
            const int   DoF   = path_and_coord.input_DoF;  

            fChainFwd.AppendArray( ApexOptics::Parse_NPolyArray_from_file(path, coord, DoF) ); 
        }
    }

    //initialize the 'reverse' chain 
    void CreateChainRev(const std::vector<ChainConstructor_t>& chain) {
        
        for (const auto& path_and_coord : chain) {

            const char* path  = path_and_coord.path.c_str(); 
            const auto  coord = path_and_coord.coords;
            const int   DoF   = path_and_coord.input_DoF;  

            fChainRev.AppendArray( ApexOptics::Parse_NPolyArray_from_file(path, coord, DoF) ); 
        }
    }

    Trajectory_t Compute_Xsv_first_guess(Trajectory_t Xfp) const {
        return ApexOptics::RVec_to_Trajectory_t(
            fChainRev.Eval(ApexOptics::Trajectory_t_to_RVec(Xfp)) 
        ); 
    }; 

    Trajectory_t Compute_Xsv(Trajectory_t Xfp, TVector3 vtx_hcs) const {

        using rvecd = ROOT::VecOps::RVec<double>; 

        //create 'first guess' 
        auto&& Xfp_rvec = ApexOptics::Trajectory_t_to_RVec(Xfp);

        auto&& Xsv_fg_rvec = fChainRev.Eval(Xfp_rvec); 
        
        auto&& Xsv_fg = ApexOptics::RVec_to_Trajectory_t(Xsv_fg_rvec);

        auto&& J = fChainFwd.Jacobian(Xsv_fg_rvec);

        RMatrix Ji(4,4); 
        rvecd J_pivot; J_pivot.reserve(4);

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

        auto&& dX = Ji.Solve( J_pivot*(-1.) ); 
        
        //check for NaN / invalid result 
        if (dX.size() != 4)                return Trajectory_t{std::numeric_limits<double>::quiet_NaN()}; 
        for (double& x : dX) { if (x != x) return Trajectory_t{std::numeric_limits<double>::quiet_NaN()}; } 

        dX.insert(dX.begin()+fPivotElement, 1. ); 

        Trajectory_t dXsv = ApexOptics::RVec_to_Trajectory_t(dX); 

        //now, project onto the target 
        Trajectory_t Xhcs_fg    = ApexOptics::SCS_to_HCS(f_is_RHRS, Xsv_fg); 
        Trajectory_t Xhcs_fg_dp = ApexOptics::SCS_to_HCS(f_is_RHRS, Xsv_fg + dXsv); 
        
        double y0 = Xhcs_fg.y       + (Xhcs_fg.dydz     * vtx_hcs.z()); 
        double y1 = Xhcs_fg_dp.y    + (Xhcs_fg_dp.dydz  * vtx_hcs.z()); 
            
        double d_pivot = ( vtx_hcs.y() - y0 )/( y1 - y0 ); 

        return Xsv_fg + (dXsv * d_pivot); 
    }   


};



#endif