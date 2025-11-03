#ifndef EvaluateCutFrame_h_
#define EvaluateCutFrame_h_

#include <TGFrame.h>
#include <TGWindow.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButton.h>
#include <TH2D.h>
#include <TLine.h>
#include "PickSieveHoleApp.h"
#include "SieveHoleData.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx> 
#include <vector> 

class EvaluateCutFrame : public TGMainFrame {
private: 
    TGHorizontalFrame *fFrame_canv; 
    TRootEmbeddedCanvas *fEcanvas; 

    TGHorizontalFrame *fFrame_buttons; 

    TGTextButton *fButton_Save; 
    TGTextButton *fButton_Reject; 

    struct HistAndLimit {
        
        TH2D *hist{nullptr}; 
        TLine *lim_low{nullptr}, *lim_high{nullptr}; 

        HistAndLimit(TH2D* _hist=nullptr) : hist{_hist} {}; 

        ~HistAndLimit() {
            if (lim_low)  delete lim_low; 
            if (lim_high) delete lim_high; 
        }

        void DrawLim_low (double val);
        void DrawLim_high(double val);  
    }; 
    TH2D *fHist_holes;
    //*fHist_yfp, *fHist_dxdzfp, *fHist_dydzfp; 
    HistAndLimit fY_fp, fDxdz_fp, fDydz_fp; 

    const double fFpcoord_cut_width; 

    PickSieveHoleApp *fParent; 

    ROOT::RDF::RNode *fRDF{nullptr}; 

    SieveHoleData *fSelectedSieveHole{nullptr}; 

    //type of canvas event passed to 'fEventType' 
    int fEventType{-1}; 

    double fX_min{DOUBLE_NAN}, fX_max{DOUBLE_NAN}; 

public: 
    EvaluateCutFrame(   const TGWindow *p, 
                        PickSieveHoleApp *_parent, 
                        const std::vector<EventData>& data, 
                        SieveHoleData *_hd,
                        const double fp_cut_width, 
                        const int n_raster_partitions, 
                        const char* branch_x="dxdz_sv", 
                        const char* branch_y="dydz_sv", 
                        const char* draw_option="col2",
                        const unsigned int palette=kBird    );
     
    ~EvaluateCutFrame(); 

    //slots for button signals
    void DoSave(); 
    void DoReject(); 
    
    void HandleCanvasClicked(); 
    
    //create fit points by fitting the profile of each x-bin of the histogram
    struct FitPoint_t { double x,y,sigma,N; };
    std::vector<FitPoint_t> CreatePointsFromHist(TH2D* hist); 

    //create a polynomial fit of the given degree from the points given
    ROOT::RVec<double> FitPolynomialToPoints(const std::vector<FitPoint_t>& points, const int poly_degree); 
      
    void Draw_Hist_Points_Poly(TH2D* hist, const std::vector<FitPoint_t>& points, const ROOT::RVec<double>& poly, const char* draw_option); 

    void UpdateButtons(); 
    void DrawLimits();     

    ClassDef(EvaluateCutFrame, 1); 
};

#endif 