#ifndef EvaluateCutFrame_h_
#define EvaluateCutFrame_h_

#include <TGFrame.h>
#include <TGWindow.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButton.h>
#include <TH2D.h>
#include <TLine.h>
#include <TObject.h> 
#include "PickSieveHoleApp.h"
#include "SieveHoleData.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx> 
#include <vector> 
#include <NPoly.h> 
#include <ApexOptics.h>
#include <functional> 

class EvaluateCutFrame : public TGMainFrame {
private: 

    PickSieveHoleApp* fParent; 

    //valid states of the app
    enum AppState {
        kNone = 0,      //startup of the app / intermediate state (disables buttons)
        kPickLimits,    //picking limits for fitting on the focal-plane plots
        kEvaluated,     //fits evaluated, ready to accept / reject 
    };
    AppState fAppState{kNone}; 

    TGHorizontalFrame *fFrame_canv; 
    TRootEmbeddedCanvas *fEcanvas; 

    TGHorizontalFrame *fFrame_buttons; 

    TGTextButton *fButton_Save; 
    TGTextButton *fButton_Reject; 
    TGTextButton *fButton_Fit; 
    TGTextButton *fButton_ResetFits; 

    struct HistAndLimit {
        
        TH2D *hist{nullptr}; 
        TLine *lim_low{nullptr}, *lim_high{nullptr}; 
        
        TF1 *fit_lo; 
        TF1 *fit_hi; 

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

    ROOT::RDF::RNode *fRDF{nullptr}; 

    SieveHoleData *fSelectedSieveHole{nullptr}; 

    //the function this points to should return 'true' if a given event is inside the cut, and 'false' otherwise
    std::function<bool(const EventData&)> fCutFcn; 

    //data for all events 
    const std::vector<EventData>* fData;  

    //type of canvas event passed to 'fEventType' 
    int fEventType{-1}; 

    //order of polynomial we're fitting FP-coordinates with (the higest exponent in the polynomial)
    int fPolynomialOrder{3}; 
    int fRasterPolyOrder{2}; 

    //minimum draw range of x_fp
    const double fXfp_draw_min{-0.65}; 
    //maximum draw range of x_fp
    const double fXfp_draw_max{+0.65}; 

    double fX_min{DOUBLE_NAN}, fX_max{DOUBLE_NAN}; 

    NPoly FitPolynomialToFP(
        const std::function<bool(const EventData&)>& is_inside_cut,     // the cut for events on this particular hole
        double ApexOptics::Trajectory_t::*coord                         // the coordinate we're fitting to (y, dxdz, dydz)
    ) const; 

    /// @brief  Attempts to delete 'obj' from all drawn sub-pads of this application (noop when passed a nullptr -- it's fine to do.)  
    void DeleteDrawnObject(TObject* obj); 

public: 
    EvaluateCutFrame(   const TGWindow *p, 
                        const std::vector<EventData>&data, 
                        SieveHoleData *_hd,
                        const double fp_cut_width, 
                        const char* branch_x="dxdz_sv", 
                        const char* branch_y="dydz_sv", 
                        const char* draw_option="col2",
                        const unsigned int palette=kBird    );
     
    ~EvaluateCutFrame(); 

    //slots for button signals
    void DoSave(); 
    void DoReject(); 
    void DoFit();
    void DoResetFits();  
    
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