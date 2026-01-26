
/***
 * @class AngleRecoTester
 * 
 * @brief Class which tests positions and slopes of reconstructed sieve-holes 
 * 
 ***/

#ifndef AngleRecoTester_H
#define AngleRecoTester_H

//APEX headers
#include <ApexOptics.h> 
#include <ROOT/RDataFrame.hxx> 
//ROOT headers
#include <TH1D.h> 
#include <TH2D.h> 
#include <TVirtualPad.h> 
//std-lib headers
#include <vector> 
#include <optional>
#include <string> 
#include <functional> 
#include <array> 

struct AngleFit_t {

    //the sieve hole associated with this fit
    SieveHole hole;
    //the position that this angle should be, based on target used and sieve-hole position 
    double angle_real;
    //the position of the peak, from the gaussian fit
    double angle_fit; 
    //the sigma of the gaussian fit
    double angle_sigma; 
    //the slope of the hole's tilt: d[phi]/d[theta], where phi=dy/dz theta=dx/dz
    double angle_slope; 
    //the amplitude of the gaussian fit (with background subtracted). 
    double amplitude; 
    //statistical error associated with hole-fit
    double angle_fit_staterr; 
}; 

struct AngleFitResult_t {
    
    double RMS_position;
    double RMS_smearing; 
    double RMS_overall; 

    std::vector<AngleFit_t> fits{}; 
}; 

class AngleRecoTester {
private: 

    const bool fIs_RHRS; 
    
    ApexOptics::OpticsTarget_t fTarget;

    ROOT::RDF::RNode fNode; 

    std::string fBranch_dxdz{"reco_dxdz_sv"}, fBranch_dydz{"reco_dydz_sv"}; 

    //do drawings along with fitting
    bool fDo_drawing=true; 

    //reduce the amount of printing to the stdout stream 
    bool fQuiet=false; 

    //drawing / fitting range for each coordinate
    double fDxdz_min{-0.04}, fDxdz_max{+0.04}, fDydz_min{-0.04}, fDydz_max{+0.04}; 

    //some miscellaneous state flags 
    int32_t fStateFlag;
 
public: 

    struct SlopeFit_t { double m, m_err, b; }; 

    enum EFlag : int32_t {
        kDxdz    = 1 << 0,  //measure dxdz
        kDydz    = 1 << 1,  //measure dydz 
        kSlopes  = 1 << 2,  //measure slopes
        kDraw_slopes = 1 << 3,
        kDraw_boxes  = 1 << 4, 
        kDraw_slope_points = 1 << 5 
    }; 

    AngleRecoTester(
        bool _is_RHRS, 
        ApexOptics::OpticsTarget_t _target, 
        ROOT::RDF::RNode _dataframe
    ); 

    ~AngleRecoTester() {};

    std::optional<AngleFitResult_t> Measure(
        int flag, 
        const int first_row,                                    
        const int last_row,      
        const int first_col,
        const int last_col,
        const double row_cut_width=1.2,
        const double col_cut_width=0.80,
        const double bg_cut_width=1.00
    ); 

    void SetName_dxdz(std::string dxdz) { fBranch_dxdz=dxdz; }
    void SetName_dydz(std::string dydz) { fBranch_dydz=dydz; }

    void SetDrawing(bool _do_drawing) { fDo_drawing=_do_drawing; }
    void SetQuiet(bool _quiet) { fQuiet=_quiet; }

    void SetRange_dxdz(double _min, double _max) { fDxdz_min=_min; fDxdz_max=_max; }; 
    void SetRange_dydz(double _min, double _max) { fDydz_min=_min; fDydz_max=_max; };  

    //turn on a state flag
    void SetFlag(int32_t flag)      { fStateFlag = fStateFlag | flag; } 
    //turn off a state flag
    void UnsetFlag(int32_t flag)    { fStateFlag = fStateFlag & (~flag); }  

private: 

    //fit the background (gaps between holes) with a polynomial
    std::function<double(double*,double*)> 
        FitBackground(TH1D* hist, const std::vector<std::array<double,2>>& exclusions, const int order);

    //measure the 'tilting' of a single hole 
    SlopeFit_t MeasureSlope(
        TH2D* h_data, 
        const double x0, 
        const double y0, 
        const std::array<double,2>& xlim, 
        const std::array<double,2>& ylim, 
        const double hole_sigma,
        const std::function<double(double*,double*)>& bg_fcn, 
        TVirtualPad* pad_2d=nullptr
    );

    ClassDef(AngleRecoTester,1); 
}; 

#endif 