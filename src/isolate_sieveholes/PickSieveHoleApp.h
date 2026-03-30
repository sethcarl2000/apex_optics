#ifndef PickSieveHoleApp_h_
#define PickSieveHoleApp_h_
////////////////////////////////////////////////////////////
//  
//  The interactive app which picks sieve-holes, and generates polynomials with each sieve-hole =
//  reconstruct sieve-hole coordinates. 
//
////////////////////////////////////////////////////////////

// APEX headers
#include "SieveHoleData.h"
#include <ApexOptics.h> 
#include <ModularOpticsModel.h>
// TGUI headers
#include <TGFrame.h>
#include <TGWindow.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGSlider.h>
#include <TGButton.h>
// ROOT headers
#include <TH2D.h>
#include <TVector3.h>
#include <ROOT/RDataFrame.hxx>
#include <TColor.h>
// stdlib headers
#include <vector>
#include <string> 

class EvaluateCutFrame; 
class SaveOutputFrame; 

class PickSieveHoleApp : public TGMainFrame {
private:

    //pointer to the (one) singleton instance of this app 
    static PickSieveHoleApp *fInstance; 

    enum AppState : int32_t {

        kNone = 0, 
        
        kAppDisabled  = 1 << 0, // if this bitflag is active, then user interactivity is disabled
        
        kHoleSelected   = 1 << 1, // if this bitflag is active, then a sieve-hole has been selected 

        kEvaluateCutFrame = (kAppDisabled | kHoleSelected | (1 << 2)), // The 'EvaluateCutFrame' app is active; this app is disabled 

        kSelectHole             = (1 << 3), // Selecting new sieve-hole (not yet selected)
        kAlreadyEvaluated       = (kHoleSelected | 1 << 4), // selected a hole which has already been evaluated 
        kPickCutOnPlot          = (kHoleSelected | 1 << 5), // selected a sieve hole which has not been evaluated yet
        kReadyToEvaluate        = (kHoleSelected | 1 << 6)
    }; 
    AppState fAppState{kNone}; 
    

    UInt_t fWindow_width, fWindow_height; 
    const TGWindow* fptr_TGWindow{nullptr}; 
    std::string fPath_infile; 

    TGHorizontalFrame*   fFrame_canvPick;
    
    TRootEmbeddedCanvas* fEcanvas_data;
    TRootEmbeddedCanvas* fEcanvas_drawing; 
    int fEventType; 

    TGHorizontalFrame *fFrame_PickHoleButtons; 
    TGHorizontalFrame *fFrame_numbers; 
    
    
    ApexOptics::OpticsTarget_t fTarget; //the optics target we're looking at 
    
    
    //buttons, sorted by the 'status' in which they appear
    TGTextButton*  fButton_Exit;     //Exit

    //additional buttons when a hole is picked on the plot
    TGTextButton*  fButton_Evaluate; //evaluate this hole

    //button which appears when a hole is picked which has already been evaluated
    TGTextButton*  fButton_Delete;   //delete this hole which has already been evaluated

    //button to reset drawn plots
    TGTextButton*  fButton_ResetPlots; 
    

    //label for the currently selected hole
    TGLabel* fLabel_selectedHole; 



    TGNumberEntry* fNumber_cutWidth; 
    TGLabel*       fLabel_cutWidth;  
    const double fCutWidth_default =  1.75; //units are in mm 
    const double fCutWidth_max     = 10.00; //units are in mm 

    TGNumberEntry* fNumber_cutHeight;
    TGLabel*       fLabel_cutHeight; 
    const double fCutHeight_default =  1.75; //units are in mm 
    const double fCutHeight_max     = 10.00; //units are in mm 

    TGSlider*      fSlider_cutAngle; 
    TGLabel*       fLabel_cutAngle; 
    const double fCutAngle_default =  0.; //units: degrees
    const double fCutAngle_max     = 30.; //units: degrees

    
    //entering a file path in a text window*/ 
    TGTextButton* fButton_Save_Output; 
    //Status of sieve-hole picking window
    enum EWindow { 
        kWindow_PickSieveHole=0, //the default window, evaluate the sieve-holes
        kWindow_EvalSieveHole=1  //the window which is used to evaluate the sieve-holes, once one is picked. 
    }; 
    EWindow fCurrentWindow{kWindow_PickSieveHole}; 

    //possible status in the sieve-hole picking window
    //this is a guide of which buttons should be visible in each of these states: 
    //
    //      EPickStatus     |   Buttons
    //      ________________|___________________________________________
    //      kNoneSelected   |   Exit
    //      kPickHoleOnPlot |   Exit
    //      kReadyToEval    |   Exit Evaluate
    //      kOldSelected    |   Exit Delete

    enum EPickStatus { 
        kNoneSelected=0,    //no sieve hole is selected. the default state
        kPickHoleOnPlot=1,  //a new hole is selected; we are ready to pick it on the plot
        kReadyToEval=2,     //a sieve-hole has been picked on the plot, we are ready to evaluate it  
        kOldSelected=3      //a hole has been selected which is already evaluated. we will ask the user if they want to delete it.
    }; 
    EPickStatus fCurrentPickStatus{kNoneSelected}; 

    TH2D *fSieveHist;   //histogram to draw actual data
    TH2D *fHoleDrawingHist; //histogram which has schematic drawing of all sieve-holes

    //size of sieve holes in drawing
    const double fDrawSize_big   = 3.0e-3; 
    const double fDrawSize_small = 1.5e-3;

    std::vector<SieveHoleData> fSieveHoleData; 
    SieveHoleData *fSelectedSieveHole{nullptr}; 

    TVector3 fReactVertex; 

    //the optics model for this particular wire/run 
    const ModularOpticsModel* fOpticsModel{nullptr}; 

    //max allowable dist between the fp-coord value and the actual value
    const double fFpcoord_cut_width{0.0055}; 

    //RDF node, with which we will do analysis 
    ROOT::RDF::RNode* fRDF{nullptr};  

    std::vector<EventData> fEventData{}; 

    const std::string fTreeName = "tracks_fp"; 
    const bool f_is_RHRS; 

    double fRastMin, fRastMax; 

    const char* fDrawingOption; 
    unsigned int fPalette;

    //drawing ranges
    double xsv_draw_range[2] = { -0.0675, +0.0675 }; 
    double ysv_draw_range[2] = { -0.0350, +0.0350 }; 
    double xfp_draw_range[2] = { -0.725, +0.725 };  

    /// @brief  Attempts to delete 'obj' from all drawn sub-pads of this application (noop when passed a nullptr -- it's fine to do.)  
    void DeleteDrawnObject(TObject* obj); 

public:
    
    PickSieveHoleApp(
        const TGWindow* p, 
        UInt_t w, 
        UInt_t h, 
        const bool is_RHRS, 
        const char* path_infile,
        const char* target_name="none",
        const char* branch_x="dxdz_sv",
        const char* branch_y="dydz_sv", 
        const char* drawing_option="col2", 
        unsigned int palette=kSunset
    );


    

    //launch the application 
    void LaunchApplication(); 

    //add an optics model to use with the reconstruction
    void SetOpticsModel(const ModularOpticsModel* model) { fOpticsModel=model; }; 

    /// @return ptr to singleton instance of PickSieveHoleApp
    static PickSieveHoleApp* Instance() { return fInstance; } 
    
    ~PickSieveHoleApp();
    
    void CloseWindow();     //close the window
    void DoExit();          //exit button (exits)
    void WriteOutput(); 
    void DoneWriteOutput(); 


    void DoEvaluate();        //execute the evaluation loop for the selected sievehole
    void DoneEvaluate();      //this is executed upon exiting the EvaluateCutApp 

    void DoDelete();          //delete the sieve-hole which is currently selected 

    void DoResetPlots();      //reset and redraw all plots
    
    void HandleCanvasClick_data();    //handle the canvas being clicked (data histogram)
    void HandleCanvasClick_drawing(); //handle the canvas beign clicked (drawing histogram)

    /// @brief Sets draw-range in x-sv 
    void SetDrawRange_x(double min, double max) { xsv_draw_range[0]=min; xsv_draw_range[1]=max; }
    /// @brief Sets draw-range in y-sv 
    void SetDrawRange_y(double min, double max) { ysv_draw_range[0]=min; ysv_draw_range[1]=max; }
    /// @brief Sets draw-range in x-fp
    void SetDrawRange_xfp(double min, double max) { xfp_draw_range[0]=min; xfp_draw_range[1]=max; }

    /// @return min draw range in x-sv coordinate
    double GetDrawRange_x_min() const { return xsv_draw_range[0]; }
    /// @return max draw range in x-sv coordinate
    double GetDrawRange_x_max() const { return xsv_draw_range[1]; }
    
    /// @return min draw range in x-sv coordinate
    double GetDrawRange_y_min() const { return ysv_draw_range[0]; }
    /// @return max draw range in x-sv coordinate
    double GetDrawRange_y_max() const { return ysv_draw_range[1]; }
    
    /// @return min draw range in x-sv coordinate
    double GetDrawRange_xfp_min() const { return xfp_draw_range[0]; }
    /// @return max draw range in x-sv coordinate
    double GetDrawRange_xfp_max() const { return xfp_draw_range[1]; }
        

    void UpdateButtons();     //update buttons to reflect current state

    //public methods needed by SaveOutputFrame
    TVector3 GetReactVertex() const { return fReactVertex; }

private: 
    std::string fPathInfile, fBranchX, fBranchY; 
public: 

    bool   Get_IsRHRS()     const { return f_is_RHRS; }
    std::string Get_PathInfile() const { return fPathInfile; }
    std::string Get_BranchX()    const { return fBranchX; }
    std::string Get_BranchY()    const { return fBranchY; }

    //this is called when the size / orientation of the hole cut is updated. 
    void SetCutSizeAngle(); 
    //this is called any time the drawing of all hole-cuts needs to be refreshed. 
    void DrawHoleCuts(); 

    void DrawSieveHoles(); 
    void DeselectSieveHole(); 

    ClassDef(PickSieveHoleApp, 1)
};

#endif 