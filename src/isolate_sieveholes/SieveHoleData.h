#ifndef SieveHoleData_h_
#define SieveHoleData_h_

///////////////////////////////////////////////
//
//  This defines a few structs: 
//
//  SieveHoleData
//  FPCoordPolynomial
//
///////////////////////////////////////////////

#include <ROOT/RVec.hxx> 
#include <TEllipse.h> 
#include <limits> 
#include <SieveHole.h> 
#include <ApexOptics.h>
#include <TVector3.h>  
#include <NPoly.h> 

//this just lists a few possible return values from the 'TCanvas::GetEvent()' method. 
enum ECanvasEventType { 
    kMouseButton1_down=1, 
    kMouseButton1_up=11, 
    kEnterObj=52,   //don't take my word on these last two; i'm not sure that's exactly waht they are. 
    kLeaveObj=53 
};  

const double xsv_draw_range[] = { -0.050, 0.050 }; 
const double ysv_draw_range[] = { -0.035, 0.035 }; 


struct FPcoordPolynomial {
    ROOT::RVec<double> poly{};
    
    //evaluate the polynomial 
    double Eval(double x_fp) const { 
        double val=0.; 
        for (int i=0; i<poly.size(); i++) val += pow( x_fp, i ) * poly[i]; 
        return val; 
    }; 
}; 

#define DOUBLE_NAN std::numeric_limits<double>::quiet_NaN()

struct EventData {
    //This will be coordinates in the sieve & focal plane 
    ApexOptics::Trajectory_t Xsv, Xfp; 

    //this is the 'raster index', which is defined as =0 for the minimum y-raster event, and =1 for the maximum y-raster event
    double raster_index;

    //this will be the react vertex, in the Sieve Coordinate System (SCS). 
    TVector3 vtx_scs;  
}; 

//stores information about 
struct SieveHoleData {
    
    SieveHoleData() {}; 
    SieveHoleData(const SieveHole& _hole) : hole{_hole} {}; 

    //copy constructor
    SieveHoleData(const SieveHoleData& shd) 
        : hole{shd.hole}, 
        cut_x{shd.cut_x}, cut_y{shd.cut_y},
        cut_width{shd.cut_width}, cut_height{shd.cut_height}, 
        y_fp{shd.y_fp}, dxdz_fp{shd.dxdz_fp}, dydz_fp{shd.dydz_fp},
        x_fp_min{shd.x_fp_min}, x_fp_max{shd.x_fp_max}, 
        is_evaluated{shd.is_evaluated}, 
        draw_circ{nullptr}, hole_cut{nullptr} 
        {}; 

    ~SieveHoleData() {}; 

    //clear all data for a clean slate
    void Clear() {
        cut_x=DOUBLE_NAN; 
        cut_y=DOUBLE_NAN; 
        cut_width=DOUBLE_NAN; 
        cut_height=DOUBLE_NAN; 
        y_fp=NPoly(2); dxdz_fp=NPoly(2); dydz_fp=NPoly(2); 
        x_fp_min=DOUBLE_NAN; 
        x_fp_max=DOUBLE_NAN; 
        is_evaluated=false;
    }; 

    SieveHole hole; 

    double cut_x     {DOUBLE_NAN};
    double cut_y     {DOUBLE_NAN};
    double cut_width {DOUBLE_NAN};
    double cut_height{DOUBLE_NAN}; 
    
    //These polynomials will be used to fit x_fp and rast_param to the fp-coordinates. 
    //the first coordinate is x_fp, the second coordinate is rast_param. 
    NPoly y_fp{2}, dxdz_fp{2}, dydz_fp{2}; 

    double x_fp_min{DOUBLE_NAN}; 
    double x_fp_max{DOUBLE_NAN};

    bool is_evaluated{false};

    TEllipse* draw_circ{nullptr}; 
    TEllipse* hole_cut{nullptr}; 

    //fraction of a sigma for which we will keeps events 
    double focalplane_cut_width = 1.8; 

    bool operator==(const SieveHoleData& rhs) const { return hole == rhs.hole; }
};

#endif 