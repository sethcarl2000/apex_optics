#ifndef ApexOptics_h_
#define ApexOptics_h_

//////////////////////////////////////////////////////////////////////////
//
// ApexOptics
//
// This is a namespace defintion to give us some helper functions which 
// can be used in conjuction with the ApexOptics library. 
//
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>  
#include "NPoly.h"
#include "NPolyArray.h"
#include "MultiLayerPerceptron.h"
#include <TVector3.h>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <ROOT/RVec.hxx>
#include <limits>

namespace ApexOptics {

    //given an array of input branches, and a 'target' output branch, will attempt to 
    // fit a polynomial (of order poly_order) from the inputs to the output. 
    std::map<std::string, NPoly*> Create_NPoly_fit( ROOT::RDF::RNode df, 
                                                    const int poly_order, 
                                                    const std::vector<std::string> &inputs, 
                                                    const std::vector<std::string> &outputs ); 

    //Given an std::map<string,NPoly*>, where the 'string' key is the name of the polynomial to be written to the 
    //file, this fcn will create (and truncate!) a '.dat' file which contains all information needed to 
    //reconstruct an NPoly's elements. 
    int Create_dbfile_from_polymap(bool is_RHRS, std::string path_outfile, std::map<std::string, NPoly*> polymap); 

    //Given a path to a db-file to use, the polynomial name, and a ptr to a polynomial to use, this will fill 
    // all NPoly elements found in the dbfile into the polynomial. returns the number of elements found. 
    // return code < 0 means that something has failed in the opening of the dbfile. 
    int Parse_NPoly_from_file(const char* path_dbfile, const char* poly_name, NPoly *poly);
    
    NPolyArray Parse_NPolyArray_from_file(const char* path_dbfile, const std::vector<std::string>& output_names, const int DoF); 

    //Given a MultiLayerPerceptron, create a dbifle output, so that it can be read later. 
    int Create_dbfile_from_mlp(const char* path_dbfile, const MultiLayerPerceptron* mlp); 

    //parse an NPoly from file
    MultiLayerPerceptron* Parse_mlp_from_file(const char* path_dbfile); 

    //these functions are meant to provide information about the APEX target geometry
    inline double Get_sieve_angle(bool _is_RHRS) {
        return ( _is_RHRS ? -5.372 : 5.366 ) * 0.0174532925199; //pi/180
    }
  
    inline double Get_HRS_angle(bool _is_RHRS) {
        return ( _is_RHRS ? -12.50 : 12.50 ) * 0.0174532925199; //pi/180
    }

    //units in mm. 
    inline TVector3 Get_APEX_Target_center() {
        return TVector3( 0., 0., -1053.7952e-3 );
    }
        
    //units in mm. These are in Target coordinates (TCS), obtained by rotating hall
    // coordinates first by '-sieve_angle' about the y-axis, then by pi/2 about the
    // z-axis. 
    inline TVector3 Get_sieve_pos(bool _is_RHRS) {
        return TVector3( _is_RHRS ?  -1.101e-3 :  -1.301e-3,
                         _is_RHRS ?  -3.885e-3 :   6.672e-3,
                         _is_RHRS ? 794.609e-3 : 795.766e-3 ); 
    }

    //a simple 'trajectory' object. the four 'geometric' coordinates are default-initialzied to 0., but the
    // 'dp/p' coordinate is default-initialized to 'nan', to indicate that this field is invalid if it is not defined. 
    // (for example, in VDC coordinate systems, in which the dp/p momentum coordinate cannot be defined). 
    struct Trajectory_t {
        double x{0.}, y{0.}, dxdz{0.}, dydz{0.}; 
        double dpp{std::numeric_limits<double>::quiet_NaN()}; 

        //addition of two 'Trajectory_t' objects
        Trajectory_t operator+(const Trajectory_t& rhs) const {
            return Trajectory_t{
                x    + rhs.x, 
                y    + rhs.y, 
                dxdz + rhs.dxdz, 
                dydz + rhs.dydz, 
                dpp  + rhs.dpp
            }; 
        }

        //multiplication of two 'Trajectory_t' objects
        Trajectory_t operator*(double a) const {
            return Trajectory_t{
                x       * a, 
                y       * a, 
                dxdz    * a, 
                dydz    * a, 
                dpp     * a
            }; 
        }

        //subtraction of two 'Trajectory_t' objects
        Trajectory_t operator-(const Trajectory_t& rhs) const {
            return Trajectory_t{*this} + ( rhs * -1. ); 
        }
    }; 

    //converts from the 'Hall coordinate system' to the 'Sieve coordinate system'. 
    Trajectory_t HCS_to_SCS(const bool is_RHRS, const Trajectory_t traj_hcs); 

    //converts from the 'Seive coordiante system' to the 'Hall coordinate system' 
    Trajectory_t SCS_to_HCS(const bool is_RHRS, const Trajectory_t traj_scs); 

    
    //this translates a DISPLACEMENT in HCS to SCS
    TVector3 HCS_to_SCS(const bool is_RHRS, TVector3 pos);
    
    //this translates a DISPLACEMENT in SCS to HCS
    TVector3 SCS_to_HCS(const bool is_RHRS, TVector3 dir); 


    //quick and dirty (slow) way to convert a Trajectory_t struct to an RVec
    ROOT::RVec<double> Trajectory_t_to_RVec(Trajectory_t X) noexcept;

    //quick and dirty (slow) way to convert an RVec<double> to a Trajectory_t struct
    Trajectory_t RVec_to_Trajectory_t(const ROOT::RVec<double>& V); 


    
    struct OpticsTarget_t {
        std::string name{"none"}; 

        //initialize all these with default value of 'Nan'. that way, we can check wich coordinates (if any) 
        // are precisely fixed by this particular target's geometry. 
        double 
            x_hcs{std::numeric_limits<double>::quiet_NaN()},
            y_hcs{std::numeric_limits<double>::quiet_NaN()}, 
            z_hcs{std::numeric_limits<double>::quiet_NaN()}; 

        bool operator==(const OpticsTarget_t& rhs) const { return rhs.name == name; } 
    };  

    //return a copy of a vector which is a list of all optics targets. 
    const std::vector<OpticsTarget_t> GetTargetList(); 

    //find a target in the list of target names. throws a std::inavlid_argument exception if name is invalid. 
    const OpticsTarget_t GetTarget(std::string target_name); 

};

#endif