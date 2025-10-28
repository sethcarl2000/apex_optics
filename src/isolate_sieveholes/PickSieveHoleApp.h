#ifndef PickSieveHoleApp_h_
#define PickSieveHoleApp_h_
////////////////////////////////////////////////////////////
//  
//  The interactive app which picks sieve-holes, and generates polynomials with each sieve-hole =
//  reconstruct sieve-hole coordinates. 
//
////////////////////////////////////////////////////////////

#include <TGFrame.h>
#include <TGWindow.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TH2D.h>
#include <TVector3.h>
#include <ApexOptics.h> 
#include "SieveHoleData.h"
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <string> 
#include <TColor.h>

class PickSieveHoleApp : public TGMainFrame {
private:

    TGHorizontalFrame*   fFrame_canvPick;
    
    TRootEmbeddedCanvas* fEcanvas_data;
    TRootEmbeddedCanvas* fEcanvas_drawing; 
    int fEventType; 

    TGHorizontalFrame *fFrame_PickHoleButtons; 
    TGHorizontalFrame *fFrame_numbers; 
    //buttons, sorted by the 'status' in which they appear
    TGTextButton*  fButton_Exit;     //Exit

    ApexOptics::OpticsTarget_t fTarget; //the optics target we're looking at 

    //total number of different 'fits' we're going to do 
    int fNRastPartitions{1}; 
    
    //additional buttons when a hole is picked on the plot
    TGTextButton*  fButton_Evaluate; //evaluate this hole

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

    //button which appears when a hole is picked which has already been evaluated
    TGTextButton* fButton_Delete;   //delete this hole which has already been evaluated
    
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

    //max allowable dist between the fp-coord value and the actual value
    const double fFpcoord_cut_width{0.0055}; 

    //RDF node, with which we will do analysis 
    ROOT::RDF::RNode* fRDF{nullptr};  

    const std::string fTreeName = "tracks_fp"; 
    const bool f_is_RHRS; 

    const char* fDrawingOption; 
    unsigned int fPalette;

public:
    PickSieveHoleApp(const TGWindow* p, 
                     UInt_t w, 
                     UInt_t h, 
                     const bool is_RHRS, 
                     const char* path_infile,
                     const char* target_name="none",
                     const char* branch_x="dxdz_sv",
                     const char* branch_y="dydz_sv", 
                     const int n_rast_partitions=3,
                     const char* drawing_option="col2", 
                     unsigned int palette=kSunset);
    
    ~PickSieveHoleApp();
    
    void CloseWindow();     //close the window
    void DoExit();          //exit button (exits)
    void WriteOutput(); 

    void DoEvaluate();        //execute the evaluation loop for the selected sievehole
    void DoneEvaluate();      //this is executed upon exiting the EvaluateCutApp 

    void DoDelete();          //delete the sieve-hole which is currently selected 
    
    void HandleCanvasClick_data();    //handle the canvas being clicked (data histogram)
    void HandleCanvasClick_drawing(); //handle the canvas beign clicked (drawing histogram)

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


    //this is called when the size of the hole cut is updated. 
    void SetCutSize(); 
    //this is called any time the drawing of all hole-cuts needs to be refreshed. 
    void DrawHoleCuts(); 

    void DrawSieveHoles(); 
    void DeselectSieveHole(); 

    ClassDef(PickSieveHoleApp, 1)
};

#endif 