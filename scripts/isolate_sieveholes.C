////////////////////////////////////////////////////////////////////////////////
//
//  This is a macro which creates an application which allows the creation of training data 
//  by selecting & isolating sieve-hole data 
//
//  This is the interactive app which creates a PolynomialCut, given a TH2.   
//
////////////////////////////////////////////////////////////////////////////////


#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h> 
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TApplication.h>
#include "TRandom3.h"
#include <iostream>
#include <vector>
#include <string> 
#include "TEllipse.h"
#include "TH2D.h"
#include "TStyle.h"


using namespace std; 

//this struct contains info about sieve holes
struct SieveHole {

    int row,col; 
    double x,y,radius_front,radius_back; 
    bool is_big; 

    //defining (overloading) this operator lets us use the std::find() function on a vector<SieveHole> 
    bool operator==(const SieveHole& rhs) const { return ((row==rhs.row) && (col==rhs.col)); }

    //defining this operator lets us use so-called 'ordered sets', like std::map, which can employ clever algorithms 
    // to search for a given element very qucikly
    bool operator<(const SieveHole& rhs) const { 
        if (row < rhs.row) return true;
        return col < rhs.col; 
    }; 

}; 

//this will contian all the data we need to make a sieve-hole output. 
struct FPcoordPolynomial {
    ROOT::RVec<double> poly{};
    
    //evaluate the polynomial 
    double Eval(double x) const { 
        double val=0.; 
        for (int i=0; i<poly.size(); i++) val += pow( x, i ) * poly[i]; 
        return val; 
    }; 
}; 

#define DOUBLE_NAN std::numeric_limits<double>::quiet_NaN()

//stores information about 
struct SieveHoleData {
    
    SieveHoleData(const SieveHole& _hole) : hole{_hole} {}; 

    ~SieveHoleData() { 
        if (draw_circ) delete draw_circ; 
        if (hole_cut)  delete hole_cut; 
    }; 

    SieveHole hole; 

    double cut_x     {DOUBLE_NAN};
    double cut_y     {DOUBLE_NAN};
    double cut_width {DOUBLE_NAN};
    double cut_height{DOUBLE_NAN}; 
    
    FPcoordPolynomial y_fp{}, dxdz_fp{}, dydz_fp{}; 

    double x_fp_min{DOUBLE_NAN}; 
    double x_fp_max{DOUBLE_NAN};

    bool is_evaluated{false};

    TEllipse* draw_circ{nullptr}; 
    TEllipse* hole_cut{nullptr}; 

    bool operator==(const SieveHoleData& rhs) const { return hole == rhs.hole; }
};

//this class is a mini-window which pops up when the user clicks the 'save' button. 
//_________________________________________________________________________________________________________________________________
class SaveOutputFrame : public TGMainFrame {
private: 
    TGTextButton *fButton_Save; 
    TGTextButton *fButton_Exit; 

    TGTextEntry *fTextEntry; 

    const std::vector<SieveHoleData> *fSieveHoleData; 

public: 
    SaveOutputFrame(const TGWindow *p, UInt_t w, UInt_t h, const std::vector<SieveHoleData>* _shd); 
    ~SaveOutputFrame() { Cleanup(); } 

    void DoSave(); 
    void DoExit(); 

    ClassDef(SaveOutputFrame, 1); 
}; 
//_________________________________________________________________________________________________________________________________
SaveOutputFrame::SaveOutputFrame(const TGWindow *p, UInt_t w, UInt_t h, const std::vector<SieveHoleData>* _shd)
    : TGMainFrame(p, w, h), 
      fSieveHoleData{_shd} 
{
    if (!fSieveHoleData) {
        throw logic_error("in <SaveOutputFrame::SaveOutputFrame>: ptr to SieveHoleData vector passed is null"); 
        return; 
    }

    // Set up the main frame
    SetCleanup(kDeepCleanup);

    auto bframe = new TGHorizontalFrame(this, 900, 500); 

    auto Add_button = [this](   TGHorizontalFrame *frame, 
                                TGTextButton* button, 
                                string button_label, 
                                string method, 
                                TGLayoutHints* hints  ) 
    {
        button_label = "&" + button_label; 
        button = new TGTextButton(frame, button_label.c_str(), 1); 
        button->Connect("Clicked()", "SaveOutputFrame", this, method.c_str()); 
        frame->AddFrame(button, hints); 
    };

    Add_button(bframe, fButton_Save, "Save", "DoSave()", new TGLayoutHints(kLHintsLeft  | kLHintsExpandX , 10, 10, 10, 5));
    Add_button(bframe, fButton_Exit, "Exit", "DoExit()", new TGLayoutHints(kLHintsRight | kLHintsExpandX , 10, 10, 10, 5));
    
    AddFrame(bframe, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 10, 10, 10, 10)); 

    fTextEntry = new TGTextEntry(this);
    AddFrame(fTextEntry, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 20, 20, 20, 20));  

    char buff[200]; sprintf(buff, "Save %zi sieve holes?", fSieveHoleData->size()); 

    SetWindowName("Save holes to output?");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}
//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoSave() 
{ 
    if (!fTextEntry) return; 
    cout << "saved file: " << fTextEntry->GetBuffer()->GetString() << endl; 

    DoExit(); 
    /*noop*/ 
}
//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoExit() { CloseWindow(); }
//_________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________

enum ECanvasEventType { 
    kMouseButton1_down=1, 
    kMouseButton1_up=11, 
    kEnterObj=52, 
    kLeaveObj=53 
};  
    

//_________________________________________________________________________________________________________________________________

//constructs a container of SieveHole structs with accurate positions, row/column indices, and sizes. 
// all units in mm for sieve hole positions and sizes. 
vector<SieveHole> Construct_sieve_holes(bool is_RHRS) 
{
    vector<SieveHole> sieve_holes; 

    const double dx = 5.842e-3; //all units in meters

    //the sign is different here, because the L & R-sieves are mirror-images of
    // of one another
    const double dy = is_RHRS ? 4.826e-3 : -4.826e-3;
    
    //the first row is 8 rows above (-x) the center hole
    const double x0 = -dx * 8.; 
    //the first column is 7-spaces +y from the center
    const double y0 =  dy * 7.; 

    //sieve thickness
    const double sieve_thickness = 12.7e-3; 

    //two possible hole radii
    const double holeR_small = 0.6985e-3; 
    const double holeR_big   = 1.3462e-3; 

    const double holeR_widened = 0.9525e-3; //radius of the holes whose exits are widened. 

    const int nRows=17; 
    for (int row=0; row<nRows; row++) { 
        //17 rows in total, numbered 0-16.
        // row 0 is the highest in HCS, so lowest-X in TCS.
        // Recall that in TCS, the central 'big hole' is the origin of x & y.
        
        //even rows are shifted in +y half a col-spacing
        bool evenRow = (row % 2==0); 
        int nCols;
        if (evenRow) { nCols = 15; }
        else  {
            if (row==1 || row==15) {nCols = 12;} //each odd row has a 'gap' missing, except these two
            else                   {nCols = 11;}
        } 
        
        for (int col=0; col<nCols; col++) { 
        
            //check wheter this hole is a 'big' hole
            bool bigHole ((row==8 && col==7) || (row==12 && col==3)); 

            //get this hole's x-position
            double x = x0 + ((double)row)*dx;
            double y = y0 - ((double)col)*dy;

            if (!evenRow) { 
                y += -dy/2.; //odd-holes are shifted in y a bit
                //skip over the gap which happens in some rows, but not these. 
                if ( (row!=1 && row!=15) && col>5 ) y += -dy; 
            }
            
            //we're ready to define our new hole. 
            SieveHole new_hole{
                .row    = row,
                .col    = col,
                .x      = x,
                .y      = y,
                .radius_front   = (bigHole ? holeR_big : holeR_small), 
                .radius_back    = (bigHole ? holeR_big : holeR_small),
                .is_big = bigHole
            }; 

            //check to see if this is one of the holes where the exit-hole is wider than the entrance-hole. 
            //this is true of the top-3 and bottom-3 rows, but on those rows, it is NOT so for the lateral-most 3 holes. 
            if (row <= 2 || row >= 14) {
                if (col <= 12) {
                    new_hole.radius_back = holeR_widened; 
                }
            }

            sieve_holes.push_back(new_hole); 
 
        }//for (int col=0; col<nCols; col++) 
    }//for (int row=0; row<nRows; row++) 
    
    return sieve_holes; 
}

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
    EWindow fCurrentWindow; 

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
    EPickStatus fCurrentPickStatus; 

    TH2D *fSieveHist;   //histogram to draw actual data
    TH2D *fHoleDrawingHist; //histogram which has schematic drawing of all sieve-holes

    //size of sieve holes in drawing
    const double fDrawSize_big   = 3.0e-3; 
    const double fDrawSize_small = 1.5e-3;

    std::vector<SieveHoleData> fSieveHoleData; 
    SieveHoleData *fSelectedSieveHole{nullptr}; 

    //RDF node, with which we will do analysis 
    ROOT::RDataFrame* fRDF{nullptr};  

    static constexpr std::string fTreeName{"tracks_fp"}; 
    const bool f_is_RHRS; 

    const char* fDrawingOption; 
    unsigned int fPalette;

public:
    PickSieveHoleApp(const TGWindow* p, 
                     UInt_t w, 
                     UInt_t h, 
                     const bool is_RHRS, 
                     const char* path_infile,
                     const char* branch_x="dxdz_sv",
                     const char* branch_y="dydz_sv", 
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

    void HandleCanvasClick_eval();    //handle the 'eval' canvas being clicked. 

    void UpdateButtons(); //update buttons to reflect current state

    void DoEvalSave()   {/*noop*/}; 
    void DoEvalReject() {/*noop*/}; 

    //this is called when the size of the hole cut is updated. 
    void SetCutSize(); 
    //this is called any time the drawing of all hole-cuts needs to be refreshed. 
    void DrawHoleCuts(); 

    void DrawSieveHoles(); 
    void DeselectSieveHole(); 

    //this is called to draw the 'DrawWindow_pickHole()' method
    void DrawWindow_pickHole(); 


    ClassDef(PickSieveHoleApp, 1)
};

//_________________________________________________________________________________________________________________________________
class EvaluateCutFrame : public TGMainFrame {
private: 
    TGHorizontalFrame *fFrame_canv; 
    TRootEmbeddedCanvas *fEcanvas; 

    TGHorizontalFrame *fFrame_buttons; 

    TGTextButton *fButton_Save; 
    TGTextButton *fButton_Reject; 

    TH2D *fHist_holes, *fHist_yfp, *fHist_dxdzfp, *fHist_dydzfp; 

    PickSieveHoleApp *fParent; 

    ROOT::RDataFrame *fRDF{nullptr}; 

    SieveHoleData *fSelectedSieveHole{nullptr}; 

    //type of canvas event passed to 'fEventType' 
    int fEventType{-1}; 

public: 
    EvaluateCutFrame(   const TGWindow *p, 
                        PickSieveHoleApp *_parent, 
                        ROOT::RDataFrame *_rdf, 
                        SieveHoleData *_hd,
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

    void DrawCuts() {/*noop*/}; 

    ClassDef(EvaluateCutFrame, 1); 
};
//_________________________________________________________________________________________________________________________________


using namespace std; 

const double xsv_draw_range[] = { -0.050, 0.050 }; 
const double ysv_draw_range[] = { -0.035, 0.035 }; 

//______________________________________________________________________________________________________________________
PickSieveHoleApp::PickSieveHoleApp( const TGWindow* p, 
                                    UInt_t w, 
                                    UInt_t h, 
                                    const bool is_RHRS, 
                                    const char* path_infile, 
                                    const char* coordname_x,
                                    const char* coordname_y,
                                    const char* drawing_option, 
                                    unsigned int palette) 
    : TGMainFrame(p, w, h),
    f_is_RHRS(is_RHRS), 
    fDrawingOption{drawing_option},
    fPalette{palette}
{
    fEventType=-1; 

    //to begin with, we are in the 'pick sieve hole' window
    fCurrentWindow     = kWindow_PickSieveHole;  
    fCurrentPickStatus = kNoneSelected; 

    // Set up the main frame
    SetCleanup(kDeepCleanup);
    
    //setup the canvas, and split it 
    // Create embedded canvas
    //fGFrame_pick    = new TGHorizontalFrame(this, 1400, 700); 

    fFrame_canvPick = new TGHorizontalFrame(this, 1400, 700); 

    fEcanvas_data   = new TRootEmbeddedCanvas("ECanvas_data", fFrame_canvPick, 700, 700);
    fFrame_canvPick->AddFrame(fEcanvas_data, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 5)); 
    
    TCanvas* canvas = fEcanvas_data->GetCanvas();
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(fPalette); 

    canvas->cd(); 
    //first things first, lets set up & draw the histogram: 
    ROOT::EnableImplicitMT(); 

    fRDF = new ROOT::RDataFrame(fTreeName.c_str(), path_infile); 

    fSieveHist = (TH2D*)fRDF->Histo2D<double>({
        "h_sieve_cpy", 
        "Sieve-coordinates;dx/dz_{sv};dy/dz_{sv}", 
        200, xsv_draw_range[0], xsv_draw_range[1], 
        200, ysv_draw_range[0], ysv_draw_range[1]
    }, 
        coordname_x, 
        coordname_y
    )->Clone("h_sieve_cpy"); 

    //make sure this histogram is not 'attached' to anything (won't be deleted by anything else). 
    fSieveHist->SetDirectory(0); 
    
    //draw the histogram and update the canvas. 
    fSieveHist->DrawCopy(fDrawingOption); 
    /*canvas->Modified(); 
    canvas->Update();*/  

    DrawHoleCuts(); 

    //connect the canvas to the 'HandleCanvasClick_data()' event
    canvas->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 
        "PickSieveHoleApp",
        this,
        "HandleCanvasClick_data()"
    ); 

    //construct the vector of sieve holes...
    for (const SieveHole& hole : Construct_sieve_holes(f_is_RHRS)) fSieveHoleData.emplace_back(hole); 

    fEcanvas_drawing = new TRootEmbeddedCanvas("ECanvas_drawing", fFrame_canvPick, 700, 700);
    fFrame_canvPick->AddFrame(fEcanvas_drawing, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0)); 

    canvas = fEcanvas_drawing->GetCanvas(); 
    canvas->cd(); 

    //now, create the 'sieve hole drawing' histogram, and draw it. 
    fHoleDrawingHist = new TH2D(
        "h_hole_drawing", 
        "Sieve-hole positions;x_{sv};y_{sv}", 
        200, -0.050, 0.050, 
        200, -0.045, 0.045
    ); 
    //draw thie histogram, draw the canvas
    fHoleDrawingHist->SetDirectory(0); 
    fHoleDrawingHist->DrawCopy(); 

    //draw the sieve-holes
    DrawSieveHoles(); 

    //connect the canvas to the 'HandleCanvasClick_drawing()' method 
    canvas->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 
        "PickSieveHoleApp",
        this,
        "HandleCanvasClick_drawing()"
    );

    AddFrame(fFrame_canvPick, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0)); 

    //now, we can start adding buttons
    fFrame_PickHoleButtons = new TGHorizontalFrame(this, 1400, 50);
    
    //Evaluate button
    fButton_Evaluate = new TGTextButton(fFrame_PickHoleButtons, "&Evaluate", 1); 
    fButton_Evaluate->Connect("Clicked()", "PickSieveHoleApp", this, "DoEvaluate()"); 
    fFrame_PickHoleButtons->AddFrame(fButton_Evaluate, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    
    //this frame contains all the hole size parameters
    fFrame_numbers = new TGHorizontalFrame(fFrame_PickHoleButtons, 500, 50); 
    
    //number entry label
    fLabel_cutWidth = new TGLabel(fFrame_numbers, "Cut width (mm): "); 
    fFrame_numbers->AddFrame(fLabel_cutWidth, new TGLayoutHints( kLHintsLeft | kLHintsCenterY, 20, 0, 5, 5)); 
    //Number entry - width of cut
    fNumber_cutWidth  = new TGNumberEntry(fFrame_numbers, fCutWidth_default );
    fNumber_cutWidth->SetLimits(TGNumberFormat::kNELLimitMinMax, 0., fCutWidth_max); 
    fNumber_cutWidth->SetNumAttr(TGNumberFormat::kNEAPositive); 
    fNumber_cutWidth->SetNumStyle(TGNumberFormat::kNESReal); 
    fNumber_cutWidth->Connect("ValueSet(Long_t)", "PickSieveHoleApp", this, "SetCutSize()"); 
    fNumber_cutWidth->GetNumberEntry()->Connect("ReturnPressed()", "PickSieveHoleApp", this, "SetCutSize()"); 
    fNumber_cutWidth->SetButtonToNum(false); //setting this to 'false' ensures that the 'ValueChanged' signal is sent to 'SetCutSize'
    fFrame_numbers->AddFrame(fNumber_cutWidth, new TGLayoutHints( kLHintsLeft | kLHintsCenterY, 20, 10, 5, 5)); 

    //number entry label
    fLabel_cutHeight = new TGLabel(fFrame_numbers, "Cut height (mm): "); 
    fFrame_numbers->AddFrame(fLabel_cutHeight, new TGLayoutHints( kLHintsLeft | kLHintsCenterY, 20, 0, 5, 5)); 
    //number entry - height of cut
    fNumber_cutHeight = new TGNumberEntry(fFrame_numbers, fCutHeight_default);  
    fNumber_cutHeight->SetLimits(TGNumberFormat::kNELLimitMinMax, 0., fCutHeight_max); 
    fNumber_cutHeight->SetNumAttr(TGNumberFormat::kNEAPositive); 
    fNumber_cutHeight->SetNumStyle(TGNumberFormat::kNESReal); 
    fNumber_cutHeight->Connect("ValueSet(Long_t)", "PickSieveHoleApp", this, "SetCutSize()"); 
    fNumber_cutHeight->GetNumberEntry()->Connect("ReturnPressed()", "PickSieveHoleApp", this, "SetCutSize()"); 
    fNumber_cutHeight->SetButtonToNum(false); //setting this to 'false' ensures that the 'ValueChanged' signal is sent to 'SetCutSize'
    fFrame_numbers->AddFrame(fNumber_cutHeight, new TGLayoutHints( kLHintsLeft | kLHintsCenterY, 20, 10, 5, 5)); 

    fLabel_selectedHole = new TGLabel(fFrame_numbers, "No hole selected"); 
    fFrame_numbers->AddFrame(fLabel_selectedHole, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 20,0,0,0)); 
    
    //add the number frame 
    fFrame_PickHoleButtons->AddFrame(fFrame_numbers, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 0,0,0,0)); 

    //Delete button
    fButton_Delete      = new TGTextButton(fFrame_PickHoleButtons, "&Delete", 1); 
    fButton_Delete->Connect("Clicked()", "PickSieveHoleApp", this, "DoDelete()"); 
    fFrame_PickHoleButtons->AddFrame(fButton_Delete, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //Save Output button
    fButton_Save_Output = new TGTextButton(fFrame_PickHoleButtons, "&Save Output", 1); 
    fButton_Save_Output->Connect("Clicked()", "PickSieveHoleApp", this, "WriteOutput()"); 
    fFrame_PickHoleButtons->AddFrame(fButton_Save_Output, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //Exit button
    fButton_Exit        = new TGTextButton(fFrame_PickHoleButtons, "&Exit", 1); 
    fButton_Exit->Connect("Clicked()", "PickSieveHoleApp", this, "DoExit()"); 
    fFrame_PickHoleButtons->AddFrame(fButton_Exit, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    AddFrame(fFrame_PickHoleButtons, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 5, 10));
    
    UpdateButtons(); 

    SetWindowName("Hole-isolation interactive tool");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

}
//_____________________________________________________________________________________________________________________________________
PickSieveHoleApp::~PickSieveHoleApp() {
    // Cleanup
    if (fSieveHist)       delete fSieveHist; 
    if (fHoleDrawingHist) delete fHoleDrawingHist; 
    Cleanup();
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DrawWindow_pickHole()
{

} 
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DrawSieveHoles()
{
    //basic struct to set style & color of circle
    struct CircleStyle { 
        unsigned int 
            line_color,
            line_style, 
            fill_color, 
            fill_style;  
    }; 

    auto set_circle_style = [](TEllipse *circ, const CircleStyle& style) {
        circ->SetLineColor(style.line_color);
        circ->SetLineStyle(style.line_style);
        circ->SetFillColor(style.fill_color);
        circ->SetFillStyle(style.fill_style); 
    }; 

    //hole is not selected and not evaluated
    const CircleStyle style_default {
        .line_color = 1, //black   
        .line_style = 1, //solid
        .fill_color = 0, //no fill 
        .fill_style = 0  //no fill 
    }; 

    //hole is currently selected
    const CircleStyle style_selected {
        .line_color = 1,   
        .line_style = 1, 
        .fill_color = 2,    //red 
        .fill_style = 1001, //dense dotted style
    }; 

    //hole hole has already been evaluated
    const CircleStyle style_evaluated {
        .line_color = 1,   
        .line_style = 1, 
        .fill_color = 1,   //grey 
        .fill_style = 3001, //dense dotted style
    }; 

    if (!fEcanvas_drawing) {
        Warning("PickSieveHolesApp::DrawSieveHoles", "Tried to execute when ptr for 'drawing'"
                " embedded canvas (fEcanvas_drawing) is null"); 
        return; 
    }

    TCanvas* canv = fEcanvas_drawing->GetCanvas(); 
    if (!canv) {
        Warning("PickSieveHolesApp::DrawSieveHoles", "Tried to execute when ptr for 'drawing'"
                " canvas (from: fEcanvas_drawing->GetCanvas()) is null"); 
        return; 
    }
    canv->cd(); 

    //draw all sieve-holes
    for (SieveHoleData& hole_data : fSieveHoleData) {
        auto& circ = hole_data.draw_circ;

        double draw_rad = hole_data.hole.is_big ? fDrawSize_big : fDrawSize_small; 

        //if circ already exists, delete it to prevent a memory leak
        if (circ) delete circ; 
        
        //draw the the new circle with the given hole-position and radius
        circ = new TEllipse(
            hole_data.hole.x, 
            hole_data.hole.y, 
            draw_rad,
            draw_rad
        ); 

        //set line/fill style based on the status of the hole (has it been selected / drawn / etc.)
        if ((fSelectedSieveHole != nullptr) && (hole_data == *fSelectedSieveHole)) {
            set_circle_style( circ, style_selected ); 
        } else {
            if (hole_data.is_evaluated) { set_circle_style( circ, style_evaluated ); }
            else                        { set_circle_style( circ, style_default ); } 
        }
        
        circ->Draw();  
    }
    canv->Modified(); 
    canv->Update(); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::SetCutSize()
{
    if (!fSelectedSieveHole) return; 
    if (!fNumber_cutHeight || !fNumber_cutHeight) {
        throw logic_error("in <PickSieveHoleApp::SetCutSize>: One or both of number-entry objects are null."); 
        return; 
    }

    if (!fSelectedSieveHole->is_evaluated) {
        fSelectedSieveHole->cut_width  = fNumber_cutWidth ->GetNumber() * 1e-3;
        fSelectedSieveHole->cut_height = fNumber_cutHeight->GetNumber() * 1e-3;
    }

    DrawHoleCuts(); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DrawHoleCuts() 
{
    const unsigned int line_color         = kBlue; 
    const unsigned int line_size          = 1; 
    const unsigned int line_size_selected = 3; 

    if (!fEcanvas_data) {
        Warning("PickSieveHolesApp::DrawHoleCuts", "Tried to execute when ptr for 'drawing'"
                " embedded canvas (fEcanvas_drawing) is null"); 
        return; 
    }

    TCanvas* canv = fEcanvas_data->GetCanvas(); 
    if (!canv) {
        Warning("PickSieveHolesApp::DrawHoleCuts", "Tried to execute when ptr for 'drawing'"
                " canvas (from: fEcanvas_drawing->GetCanvas()) is null"); 
        return; 
    }
    canv->cd(); 

    //draw all sieve-holes
    for (SieveHoleData& hole_data : fSieveHoleData) {

        auto& circ = hole_data.hole_cut;

        //if circ already exists, delete it to prevent a memory leak
        if (circ) delete circ; 

        //check if the cuts have been defiend for this hole
        if (hole_data.cut_x      != hole_data.cut_x ||
            hole_data.cut_y      != hole_data.cut_y ||
            hole_data.cut_width  != hole_data.cut_width ||
            hole_data.cut_height != hole_data.cut_height) continue; 

        circ = new TEllipse(
            hole_data.cut_x, 
            hole_data.cut_y, 
            hole_data.cut_width,
            hole_data.cut_height
        );  

        circ->SetLineColor(line_color); 

        //set line/fill style based on the status of the hole (has it been selected / eval'd / etc.)
        if ((fSelectedSieveHole != nullptr) && (hole_data == *fSelectedSieveHole)) { 
            circ->SetLineWidth(line_size_selected); 
        } else {
            circ->SetLineWidth(line_size); 
        }

        circ->SetFillStyle(0); //no fill (transparent)
        circ->Draw(); 
    }
    canv->Modified(); 
    canv->Update(); 
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoEvaluate() 
{
    //check if there is a hole picked. 
    if (!fSelectedSieveHole) {
        Error("DoEvaluate", "No sieve hole is selected. It shoul not have been possible to get here..."); 
        return; 
    }

    //check if our status is 'kReadyToEval' 
    /*if (fCurrentPickStatus != kReadyToEval) {
        Error("DoEvaluate", "Tried to eval when fCurrentPickStatus is not 'kRreadyToEval'. It should not have been possible to get here..."); 
        return; 
    }*/ 

    //call something like: 
    // DrawWindow_EvalSieveHole() 
    {
        fCurrentWindow = kWindow_EvalSieveHole; 
    }
    
    new EvaluateCutFrame(gClient->GetRoot(), this, fRDF, fSelectedSieveHole, "dxdz_sv", "dydz_sv", "col", kBird); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoneEvaluate()
{
    //once the EvaluateCutFrame object is done working, then we reutrn here. 
    //now, we've evaluated this sievehole

    if (!fSelectedSieveHole) {
        throw logic_error("in <PickSieveHoleApp::DoneEvaluate>: no sieve hole is selected."); 
        return; 
    }

    //this might be changed by the EvaulateCutFrame app
    gStyle->SetPalette(fPalette); 

    if (fSelectedSieveHole->is_evaluated) {

        cout << "Evaluated hole: row " << fSelectedSieveHole->hole.row << ", col " << fSelectedSieveHole->hole.col << endl; 
    
    } else { 

        cout << "Rejected hole: row " << fSelectedSieveHole->hole.row << ", col " << fSelectedSieveHole->hole.col << endl; 
    }

    //do re-drawing of buttons
    fCurrentWindow     = kWindow_PickSieveHole; 
    DeselectSieveHole(); 

    //redraw the sieve holes, and update the button layout
    DrawSieveHoles(); 
    DrawHoleCuts(); 
    UpdateButtons(); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoDelete() 
{
    //check if there is a hole picked. 
    if (!fSelectedSieveHole) {
        Error("DoDelete", "No sieve hole is selected. It shoul not have been possible to get here..."); 
        return; 
    }

    //check if our status is 'kReadyToEval' 
    if (fCurrentPickStatus != kOldSelected) {
        Error("DoDelete", "Tried to eval when fCurrentPickStatus is not 'kOldSelected'. It should not have been possible to get here..."); 
        return; 
    }

    //create a new, 'fresh' SieveHoleData, in which everything is default-initialized except the SieveHole struct. 
    SieveHoleData clean_slate(fSelectedSieveHole->hole); 
    fSelectedSieveHole->~SieveHoleData(); 
    *fSelectedSieveHole = clean_slate; 
    
    //and, finally, deselect the hole. 
    DeselectSieveHole(); 

    //redraw all the sieveholes, and update the button layout
    DrawSieveHoles(); 
    DrawHoleCuts(); 
    UpdateButtons(); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::UpdateButtons()
{
    //first, unmap all buttons
    fButton_Delete      ->UnmapWindow(); 
    fButton_Evaluate    ->UnmapWindow(); 
    
    if (!fFrame_PickHoleButtons) { Error("UpdateButtons()", "fFrame_PickHoleButtons is null"); return; }
    
    auto redraw_button = [this](TGTextButton* button) {
        this->fFrame_PickHoleButtons->AddFrame(button, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 
        button->MapWindow(); 
    };  
    
    //now, handle specific cases. we're only redrawing buttons which are specific to these cases
    switch (fCurrentPickStatus) {
        case kNoneSelected   : break;  
        case kPickHoleOnPlot : break;  
        case kOldSelected    : fButton_Delete->MapWindow();   break; 
        case kReadyToEval    : fButton_Evaluate->MapWindow(); break;  
    }
    
    fFrame_PickHoleButtons->MapWindow();
    MapWindow(); 
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::CloseWindow() {
    // Called when window is closed via window manager
    gApplication->Terminate(0);
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoExit() {
    // Called when exit button is pressed
    gApplication->Terminate(0);
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::WriteOutput() {
    new SaveOutputFrame(gClient->GetRoot(), 900, 500, &fSieveHoleData); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_data() {
    
    //get canvas event information
    TCanvas* canvas = fEcanvas_data->GetCanvas();
    if (!canvas) return;

    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType =  canvas->GetEvent(); 

        //check if the status is valid
        if (fCurrentPickStatus != kPickHoleOnPlot) return; 
        if (!fSelectedSieveHole) return;     

        //if the mouse-button was just released
        if (fEventType == kMouseButton1_up) {

            fSelectedSieveHole->cut_width  = fNumber_cutWidth ->GetNumber() * 1e-3; 
            fSelectedSieveHole->cut_height = fNumber_cutHeight->GetNumber() * 1e-3; 

             // Get pixel coordinates
            Int_t px = canvas->GetEventX();
            Int_t py = canvas->GetEventY();
            
            // Convert pixel coordinates to histogram coordinates
            Double_t x = canvas->AbsPixeltoX(px);
            Double_t y = canvas->AbsPixeltoY(py);
            
            fSelectedSieveHole->cut_x = x;
            fSelectedSieveHole->cut_y = y;              

            printf("Data histogram clicked: %+.6f, %+.6f \n", x, y); cout << endl; 

            fCurrentPickStatus = kReadyToEval; 

            DrawSieveHoles(); 
            DrawHoleCuts(); 
            UpdateButtons();
         }// if (fEventType == kMouseButton1_up) 
    }// if (fEventType != canvas->GetEvent()) 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_drawing() {
    // Get canvas and event information
    TCanvas* canvas = fEcanvas_drawing->GetCanvas();
    if (!canvas) return;
    
    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType = canvas->GetEvent(); 

        //if the mouse-button was just released
        if (fEventType == kMouseButton1_up) {

            //deselect any selected sievehole 
            DeselectSieveHole(); 

             // Get pixel coordinates
            Int_t px = canvas->GetEventX();
            Int_t py = canvas->GetEventY();
            
            // Convert pixel coordinates to histogram coordinates
            Double_t x = canvas->AbsPixeltoX(px);
            Double_t y = canvas->AbsPixeltoY(py);
            
            printf("Drawing histogram clicked: %+.6f, %+.6f ", x, y); 

            for (auto& hole_data : fSieveHoleData) {
                const auto& hole = hole_data.hole; 
                
                const double rad2 = pow( hole.is_big ? fDrawSize_big : fDrawSize_small, 2 );
                
                if ( pow(hole.x - x, 2) + pow(hole.y - y, 2) < rad2 ) { //this hole was clicked. 

                    fSelectedSieveHole = &hole_data; 
                    
                    //check to see if this hole has been evaluated or not 
                    if (fSelectedSieveHole->is_evaluated) {
                        fCurrentPickStatus = kOldSelected;      //the user has selected a hole which has already been evaluated. 
                                                                //we will ask them if they want to delete the data from this evaluation.
                    } else {
                        fCurrentPickStatus = kPickHoleOnPlot;   //this hole has not been selected. we are ready to let the user 
                                                                //select a position on the data plot to make a cicular cut on events. 
                    }

                    char buff[200]; 
                    sprintf(buff, "selected: row %i, col %i. %s", 
                        fSelectedSieveHole->hole.row, 
                        fSelectedSieveHole->hole.col,
                        fSelectedSieveHole->is_evaluated ? "Already evaluated. Delete?" : "Not yet evaluated." 
                    ); 
                    fLabel_selectedHole->SetText(buff);  

                    break; 
                }
            }//for (auto& hole_data : fSieveHoleData)

            //if no sieve hole was picked, then return. 
            if (fSelectedSieveHole == nullptr) { cout << endl; }
            else {
                printf("; hole selected, row %i, col %i", 
                fSelectedSieveHole->hole.row, 
                fSelectedSieveHole->hole.col); cout << endl; 
            } 

            DrawSieveHoles(); 
            DrawHoleCuts(); 
            UpdateButtons(); 
        }// if (fEventType == kMouseButton1_up) 
    }// if (fEventType != canvas->GetEvent()) 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DeselectSieveHole()
{
    //deselect any selected sieve holes
    //check if we've already loaded a bit of data into this hole
    if (fSelectedSieveHole && !fSelectedSieveHole->is_evaluated) {

        SieveHoleData clean_slate = SieveHoleData(fSelectedSieveHole->hole); 
        
        //clean up the sieve hole by deleting this one
        fSelectedSieveHole->~SieveHoleData(); 
        
        //this hole has not been evaluated. therefore, we need to make sure that we 'clean it up' when we deselect. 
        //but if it *has* been evaluated, we don't want to wipe the data that has already been recorded. 
        
        *fSelectedSieveHole = clean_slate; 
    }
    fLabel_selectedHole->SetText("No hole selected"); 

    fSelectedSieveHole = nullptr; 
    fCurrentPickStatus = kNoneSelected; 
}

EvaluateCutFrame::EvaluateCutFrame( const TGWindow *p, 
                                    PickSieveHoleApp *_parent, 
                                    ROOT::RDataFrame *_rdf, 
                                    SieveHoleData *_hd,
                                    const char* branch_x, 
                                    const char* branch_y, 
                                    const char* draw_option,
                                    const unsigned int palette )
    : TGMainFrame( p, 1400, 800 ), 
    fParent{_parent},
    fRDF{_rdf}, 
    fSelectedSieveHole{_hd}
{   
    const int polynomial_degree =3; 

    //setup the window
    SetCleanup(kDeepCleanup); 

    fFrame_canv = new TGHorizontalFrame(this, 1400, 800); 

    fEcanvas = new TRootEmbeddedCanvas("ECanvas_eval", fFrame_canv, 1400, 800); 

    const double cut_x      = fSelectedSieveHole->cut_x; 
    const double cut_y      = fSelectedSieveHole->cut_y; 
    const double cut_width  = fSelectedSieveHole->cut_width; 
    const double cut_height = fSelectedSieveHole->cut_height; 

    auto hist_xy = (*fRDF).Histo2D<double>({"h_xy_temp", "", 
            200, xsv_draw_range[0], xsv_draw_range[1],
            200, ysv_draw_range[0], ysv_draw_range[1] }, branch_x, branch_y);     

    fHist_holes = (TH2D*)hist_xy->Clone("h_xy"); 

    auto df_output = (*fRDF)

        .Filter([cut_x, cut_y, cut_width, cut_height](double x, double y)
        {
            return pow( (x - cut_x)/cut_width, 2 ) + pow( (y - cut_y)/cut_height, 2 ) < 1.; 
        }, {branch_x, branch_y}); 
    
    
    auto hist_y_fp    = df_output.Histo2D<double>({"h_y_temp",    "y_fp vs x_fp",     30, -0.65, 0.65, 75, -0.070, 0.055}, "x_fp", "y_fp"); 
    auto hist_dxdz_fp = df_output.Histo2D<double>({"h_dxdz_temp", "dx/dz_fp vs x_fp", 30, -0.65, 0.65, 75, -0.035, 0.025}, "x_fp", "y_fp"); 
    auto hist_dydz_fp = df_output.Histo2D<double>({"h_dydz_temp", "dy/dz_fp vs x_fp", 30, -0.65, 0.65, 75, -0.060, 0.040}, "x_fp", "y_fp"); 

    TCanvas *canv = fEcanvas->GetCanvas(); 
   
    gStyle->SetPalette(palette); 
    
    canv->cd(); 
    canv->Divide(2,2, 0.001, 0.001); 

    canv->cd(1); fHist_holes->Draw("col2"); 
    auto circ = new TEllipse(
        fSelectedSieveHole->cut_x, 
        fSelectedSieveHole->cut_y, 
        fSelectedSieveHole->cut_width,
        fSelectedSieveHole->cut_height
    ); 
    circ->SetFillStyle(0); 
    circ->SetLineColor(kRed); 
    circ->SetLineWidth(2); 
    circ->Draw(); 

    //handle y_fp drawing 
    fHist_yfp    = (TH2D*)hist_y_fp    ->Clone("h_yfp");
    fHist_dxdzfp = (TH2D*)hist_dxdz_fp ->Clone("h_dxdzfp"); 
    fHist_dydzfp = (TH2D*)hist_dydz_fp ->Clone("h_dydzfp"); 

    auto Draw_and_fit = [&](TH2D* hist, ROOT::RVec<double>& poly)
    {
        auto points = CreatePointsFromHist(hist); 
        poly = FitPolynomialToPoints(points, polynomial_degree); 
        Draw_Hist_Points_Poly(hist, points, poly, draw_option); 
    }; 

    canv->cd(2); Draw_and_fit(fHist_yfp,    fSelectedSieveHole->y_fp.poly); 
    
    canv->cd(3); Draw_and_fit(fHist_dxdzfp, fSelectedSieveHole->dxdz_fp.poly); 
    
    canv->cd(4); Draw_and_fit(fHist_dydzfp, fSelectedSieveHole->dydz_fp.poly); 
    
    canv->Modified(); 
    canv->Update(); 
    canv->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "EvaluateCutFrame", this, "HandleCanvasClicked()"); 

    fFrame_canv->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,10)); 

    AddFrame(fFrame_canv, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10)); 

    fFrame_buttons = new TGHorizontalFrame(this, 1600, 50); 

    //Reject button
    fButton_Save = new TGTextButton(fFrame_buttons, "&Save", 1); 
    fButton_Save->Connect("Clicked()", "EvaluateCutFrame", this, "DoSave()"); 
    fFrame_buttons->AddFrame(fButton_Save, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //Reject button
    fButton_Reject = new TGTextButton(fFrame_buttons, "&Reject", 1); 
    fButton_Reject->Connect("Clicked()", "EvaluateCutFrame", this, "DoReject()"); 
    fFrame_buttons->AddFrame(fButton_Reject, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //connect the exit of this app to the 'DoneEvaluate' 
    Connect("CloseWindow()", "PickSieveHoleApp", fParent, "DoneEvaluate()"); 

    AddFrame(fFrame_buttons, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 10,10,10,5 )); 

    SetWindowName("Evaluating sieve-hole cut"); 
    MapWindow();
    Resize(GetDefaultSize()); 
    MapSubwindows(); 
}   
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::HandleCanvasClicked()
{
    //get canvas event information
    TCanvas* canvas = fEcanvas->GetCanvas();
    if (!canvas) return;

    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType =  canvas->GetEvent(); 

        //mouse button 1 was just released
        if (fEventType == kMouseButton1_up) {

            //check if the selected object is a TH2D. if not, then quit. 
            auto selected_ptr = canvas->GetSelected(); 
            if (!selected_ptr || selected_ptr->IsA() != TClass::GetClass<TH2D>()) return; 

            TH2D* selected_hist = (TH2D*)selected_ptr; 

            //check if its one of the three histograms of ours.
            
            //this hist is the hole drawing; it can't be used to reset the x-range
            if (string(selected_hist->GetName()) == "h_xy") return; 

            //so, wev've selected one of our 3 fp-coord histograms. 
            //let's reset the x-drawing range
        }
    }
}   
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoSave() 
{
    if (fSelectedSieveHole) {
        fSelectedSieveHole->is_evaluated = true; 
        //save the data of this sieve-hole cut
    } else {
        throw logic_error("in <EvaluateCutFrame::DoSave>: no sieve hole is selected (ptr is null)"); 
        return; 
    }
    //exit this window 
    CloseWindow();
    fParent->DoneEvaluate(); 
}
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoReject()
{
    if (fSelectedSieveHole) {
        fSelectedSieveHole->is_evaluated = false; 
        //reject the data of this sieve-hole (reset it)
    } else {
        throw logic_error("in <EvaluateCutFrame::DoSave>: no sieve hole is selected (ptr is null)"); 
        return; 
    }
    //exit this window 
    CloseWindow(); 
    fParent->DoneEvaluate(); 
} 
//_____________________________________________________________________________________________________________________________________
vector<EvaluateCutFrame::FitPoint_t> EvaluateCutFrame::CreatePointsFromHist(TH2D* hist) {

    if (!hist) {
        throw logic_error("in <EvaluateCutFrame::CreatePointsFromHist>: hist passed is null"); 
        return {}; 
    }

    //fit a graph and polynomial to each point
    const int min_stats = 10; 
    const double min_frac  = 0.10; 
    const double fit_radius = 8e-3; 
    auto x_axis = hist->GetXaxis(); 

    vector<FitPoint_t> pts; 

    double max_proj_integral=0.; 

    for (int b=1; b<x_axis->GetNbins(); b++) { 
    
        auto proj = hist->ProjectionY("proj",b,b); 

        if (proj->Integral()<min_stats) continue; 

        double fit_center = proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() ); 
        
        max_proj_integral = max<double>( max_proj_integral, proj->Integral() ); 
    
        auto gaus_fit = new TF1("gausFit", "gaus(0)", fit_center - fit_radius, fit_center + fit_radius); 
    
        gaus_fit->SetParameter(0, proj->GetMaximum());
        gaus_fit->SetParameter(1, fit_center);
        gaus_fit->SetParameter(2, 2.5e-3); 
    
        auto fit_result = proj->Fit("gausFit", "N Q S L R"); 

        if (!fit_result->IsValid()) continue; //skip if the fit failed

        pts.push_back({
            .x      = x_axis->GetBinCenter(b),              //center of bin
            .y      = fit_result->Parameter(1),             //center of gaussian, according to fit_result
            .sigma  = pow( fit_result->ParError(1), 2 ),    //error of this value from fit_result
            .N      = proj->Integral()                      //total stats for this profile 
        }); 
    }
    
    //prune the vectors. delete points which are less than  maxHeight * minStat_frac        
    for(auto it = pts.begin(); it != pts.end();) {
        if ( it->N < min_frac * max_proj_integral ) { pts.erase(it); }
        else                                        { it++; }
    }
    return pts; 
}
//_____________________________________________________________________________________________________________________________________
ROOT::RVec<double> EvaluateCutFrame::FitPolynomialToPoints(const vector<EvaluateCutFrame::FitPoint_t>& points, const int polynomial_degree)
{
    //minimum number of points
    if (points.size() < polynomial_degree+1) return {}; 

    RMatrix A(polynomial_degree+1, polynomial_degree+1, 0.); 

    ROOT::RVec<double> B(polynomial_degree+1, 0.); 

    for (const FitPoint_t& pt : points) {
        
        for (int i=0; i<=polynomial_degree; i++) { 

            B[i] += pt.y * pow( pt.x, i ) / pt.sigma; 
            
            for (int j=0; j<=polynomial_degree; j++) A.get(i,j) += pow( pt.x, i ) * pow( pt.x, j ) / pt.sigma; 
        }
    }
    auto coeffs = A.Solve( B ); 
    
    //check coeffs for NaN 
    for (double x : coeffs) if (x != x) { return{}; }

    return coeffs; 
}; 
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::Draw_Hist_Points_Poly(TH2D* hist, const vector<FitPoint_t>& points, const ROOT::RVec<double>& poly, const char* draw_option)
{
    hist->DrawCopy(draw_option); 

    if (points.empty()) return; 

    //draw the points & poylnomials
    for (const auto& pt : points) {
        auto box = new TBox( 
            pt.x - 2e-3, pt.y - 10. * sqrt(pt.sigma), 
            pt.x + 2e-3, pt.y + 10. * sqrt(pt.sigma)
        ); 
        box->SetLineColor(1); 
        box->SetLineWidth(1); 
        box->SetFillStyle(0); 
        box->Draw("SAME"); 
    }

    if (poly.empty()) return; 

    char buff[50]; sprintf(buff, "pol%zi", poly.size()-1); 
    
    char name[200]; sprintf(name, "%s_poly", hist->GetName()); 

    if (!points.empty()) {
        auto f1_poly = new TF1(name, buff, -0.65, 0.65); 
        f1_poly->SetLineColor(kRed); 
        f1_poly->SetLineWidth(2); 
        int i=0; for (double coeff : poly) f1_poly->SetParameter(i++, coeff); 
        f1_poly->Draw("SAME");
    } 
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
EvaluateCutFrame::~EvaluateCutFrame() { 
    
    //delete histograms
    if (fHist_yfp)    delete fHist_yfp;    
    if (fHist_dxdzfp) delete fHist_dxdzfp;
    if (fHist_dydzfp) delete fHist_dydzfp; 
    
    Cleanup(); 
}
//_____________________________________________________________________________________________________________________________________


//_____________________________________________________________________________________________________________________________________
int isolate_sieveholes( const bool is_RHRS, 
                        const char* path_infile,
                        const char* coord_x ="dxdz_sv",
                        const char* coord_y ="dydz_sv",
                        const char* drawing_option ="col2", 
                        const unsigned int palette =kSunset)
{
    new PickSieveHoleApp(gClient->GetRoot(), 
                         1400, 
                         700, 
                         is_RHRS,
                         path_infile, 
                         coord_x, 
                         coord_y, 
                         drawing_option,
                         palette);
    return 0; 
}

ClassImp(EvaluateCutFrame); 
ClassImp(PickSieveHoleApp); 
ClassImp(SaveOutputFrame); 