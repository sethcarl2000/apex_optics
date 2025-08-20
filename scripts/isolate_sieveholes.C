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

constexpr double NaN{std::numeric_limits<double>::quiet_NaN()}; 

//this will contian all the data we need to make a sieve-hole output. 
struct FPcoordPolynomial {
    double 
        x_min, 
        x_max; 
    ROOT::RVec<double> poly{}; 
}; 

//stores information about 
struct SieveHoleData {
    
    SieveHole hole; 

    double cut_x    {std::numeric_limits<double>::quiet_NaN()};
    double cut_y    {std::numeric_limits<double>::quiet_NaN()};
    double cut_width {std::numeric_limits<double>::quiet_NaN()};
    double cut_height{std::numeric_limits<double>::quiet_NaN()}; 
    
    FPcoordPolynomial fpcoord_poly; 

    bool is_evaluated{false};

    TEllipse* draw_circ{nullptr}; 

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
    TGHorizontalFrame*   fFrame_canv; 
    
    TRootEmbeddedCanvas* fEcanvas_data;
    TRootEmbeddedCanvas* fEcanvas_drawing; 
    int fEventType; 

    TGHorizontalFrame *fFrame_PickHoleButtons; 
    TGHorizontalFrame *fFrame_numbers; 
    //buttons, sorted by the 'status' in which they appear
    TGTextButton* fButton_Exit;     //Exit
    
    //additional buttons when a hole is picked on the plot
    TGTextButton*  fButton_Evaluate; //evaluate this hole

    TGLabel* fLabel_holename; 

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
    //      kOldSelected    |   Exit Reset (Delete)  

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


    static constexpr std::string fTreeName{"tracks_fp"}; 
    const bool f_is_RHRS; 

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
    void DoDelete();          //delete the sieve-hole which is currently selected 
    
    void HandleCanvasClick_data();    //handle the canvas being clicked (data histogram)
    void HandleCanvasClick_drawing(); //handle the canvas beign clicked (drawing histogram)

    void UpdateButtons(); 

    //this is called when the size of the hole cut is updated. 
    void SetCutSize(); 
    //this is called any time the drawing of all hole-cuts needs to be refreshed. 
    void DrawHoleCuts(); 

    void DrawSieveHoles(); 
    void DeselectSieveHole(); 

    enum ECanvasEventType { kMouseButton1_down=1, kMouseButton1_up=11, kEnterObj=52, kLeaveObj=53 };  
    
    ClassDef(PickSieveHoleApp, 1)
};


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
    f_is_RHRS(is_RHRS)
{
    fEventType=-1; 

    //to begin with, we are in the 'pick sieve hole' window
    fCurrentWindow     = kWindow_PickSieveHole;  
    fCurrentPickStatus = kNoneSelected; 

    // Set up the main frame
    SetCleanup(kDeepCleanup);
    
    //setup the canvas, and split it 
    // Create embedded canvas

    fFrame_canv   = new TGHorizontalFrame(this, 1400, 700); 

    fEcanvas_data = new TRootEmbeddedCanvas("ECanvas_data", this, 700, 700);
    fFrame_canv->AddFrame(fEcanvas_data, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 5)); 
    
    TCanvas* canvas = fEcanvas_data->GetCanvas();
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(palette); 

    canvas->cd(); 
    //first things first, lets set up & draw the histogram: 
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df(fTreeName.c_str(), path_infile); 

    fSieveHist = (TH2D*)df.Histo2D<double>({
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
    fSieveHist->Draw(drawing_option); 
    canvas->Update(); 

    //connect the canvas to the 'HandleCanvasClick_data()' event
    canvas->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 
        "PickSieveHoleApp",
        this,
        "HandleCanvasClick_data()"
    ); 

    //construct the vector of sieve holes...
    for (const SieveHole& hole : Construct_sieve_holes(f_is_RHRS)) fSieveHoleData.push_back({ .hole = hole }); 

    fEcanvas_drawing = new TRootEmbeddedCanvas("ECanvas_drawing", this, 700, 700);
    fFrame_canv->AddFrame(fEcanvas_drawing, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 5)); 

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
    fHoleDrawingHist->Draw(); 

    //draw the sieve-holes
    DrawSieveHoles(); 

    //connect the canvas to the 'HandleCanvasClick_drawing()' method 
    canvas->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 
        "PickSieveHoleApp",
        this,
        "HandleCanvasClick_drawing()"
    );

    AddFrame(fFrame_canv, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 5)); 

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
    
    
    //add the number frame 
    fFrame_PickHoleButtons->AddFrame(fFrame_numbers, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0,0,0,0)); 

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

    fSelectedSieveHole->cut_width  = fNumber_cutWidth ->GetNumber() * 1e-3;
    fSelectedSieveHole->cut_height = fNumber_cutHeight->GetNumber() * 1e-3;


    //DrawSieveHoleCuts(); 
}
//_____________________________________________________________________________________________________________________________________
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
    
    //... 
    //... compute the sieve-hole data ...
    //... 

    //now, we've evaluated this sievehole
    fSelectedSieveHole->is_evaluated = true; 

    cout << "Evaluated hole: row " << fSelectedSieveHole->hole.row << ", col " << fSelectedSieveHole->hole.col << endl; 

    //return to the 'pick sieve-hole' window, and deselect the sieve-hole. 
    //call something like: 
    // DrawWindow_PickSieveHole(); 
    {
        //do re-drawing of buttons
        fCurrentWindow     = kWindow_PickSieveHole; 
        DeselectSieveHole(); 
    }

    //redraw the sieve holes, and update the button layout
    DrawSieveHoles(); 
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
    SieveHoleData new_hole{ .hole = fSelectedSieveHole->hole }; 
    *fSelectedSieveHole = new_hole; 
    
    //and, finally, deselect the hole. 
    DeselectSieveHole(); 

    //redraw all the sieveholes, and update the button layout
    DrawSieveHoles(); 
    UpdateButtons(); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::UpdateButtons()
{
    //first, unmap all buttons
    fButton_Delete      ->UnmapWindow(); 
    fButton_Evaluate    ->UnmapWindow(); 
    fFrame_numbers      ->UnmapWindow(); 

    if (!fFrame_PickHoleButtons) { Error("UpdateButtons()", "fFrame_PickHoleButtons is null"); return; }
    
    auto redraw_button = [this](TGTextButton* button) {
        this->fFrame_PickHoleButtons->AddFrame(button, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 
        button->MapWindow(); 
    };  
    
    //now, handle specific cases. we're only redrawing buttons which are specific to these cases
    switch (fCurrentPickStatus) {
        case kNoneSelected   : break;  
        case kPickHoleOnPlot : fFrame_numbers->MapWindow(); break;  
        case kOldSelected    : fButton_Delete->MapWindow(); break; 
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
    
    //check if the status is valid
    if (fCurrentPickStatus != kPickHoleOnPlot) return; 

    double cut_radius = 1.75e-3; 

    //get canvas event information
    TCanvas* canvas = fEcanvas_data->GetCanvas();
    if (!canvas) return;

    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType = canvas->GetEvent(); 

        //if the mouse-button was just released
        if (fEventType == kMouseButton1_up) {

             // Get pixel coordinates
            Int_t px = canvas->GetEventX();
            Int_t py = canvas->GetEventY();
            
            // Convert pixel coordinates to histogram coordinates
            Double_t x = canvas->AbsPixeltoX(px);
            Double_t y = canvas->AbsPixeltoY(py);
            


            printf("Data histogram clicked: %+.6f, %+.6f \n", x, y); cout << endl; 
        }
    }
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

        //delete the ellipse drawn for this sieve hole 
        if (fSelectedSieveHole->draw_circ) delete fSelectedSieveHole->draw_circ;  
        
        //this hole has not been evaluated. therefore, we need to make sure that we 'clean it up' when we deselect. 
        //but if it *has* been evaluated, we don't want to wipe the data that has already been recorded. 
        SieveHoleData clean_slate{ .hole = fSelectedSieveHole->hole }; 

        *fSelectedSieveHole = clean_slate; 
    }
    fSelectedSieveHole = nullptr; 
    fCurrentPickStatus = kNoneSelected; 
}
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




ClassImp(PickSieveHoleApp); 
ClassImp(SaveOutputFrame); 