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

    FPcoordPolynomial fpcoord_poly; 

    bool is_evaluated{false};
    bool selected{false};  

    TEllipse* draw_circ{nullptr}; 
};

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
    
    //buttons, sorted by the 'status' in which they appear
    TGTextButton* fButton_Exit;     //Exit

    //buttons which appear when a NEW sieve hole is selected (kPickHoleOnPlot)
    /*TGTextButton* fButton_Deselect; //deselect this hole
    
    //additional buttons when a hole is picked on the plot
    TGTextButton* fButton_Evaluate; //evaluate this hole

    //button which appears when a hole is picked which has already been evaluated
    TGTextButton* fButton_Delete;   //delete this hole which has already been evaluated
    
    
    
    //entering a file path in a text window
    TGTextEntry*  fFilePathEntry;  //window to enter file path
    TGTextButton* fButton_write_output; */ 
    //Status of sieve-hole picking window
    enum EWindow { 
        kPickSieveHole=0, //the default window, evaluate the sieve-holes
        kEvalSieveHole=1  //the window which is used to evaluate the sieve-holes, once one is picked. 
    }; 
    EWindow fState_window; 

    //possible status in the sieve-hole picking window
    //this is a guide of which buttons should be visible in each of these states: 
    //
    //      EPickStatus     |   Buttons
    //      ________________|___________________________________________
    //      kNoneSelected   |   Exit
    //      kPickHoleOnPlot |   Exit Deselect
    //      kReadyToEval    |   Exit Deselect Evaluate
    //      kOldSelected    |   Exit Deselect Delete   

    enum EPickStatus { 
        kNoneSelected=0,    //no sieve hole is selected. the default state
        kPickHoleOnPlot=1,  //a new hole is selected; we are ready to pick it on the plot
        kReadyToEval=2,     //a sieve-hole has been picked on the plot, we are ready to evaluate it  
        kOldSelected=3      //a hole has been selected which is already evaluated. we will ask the user if they want to delete it.
    }; 
    EPickStatus fState_pickStatus; 

    TH2D *fSieveHist;   //histogram to draw actual data
    TH2D *fHoleDrawingHist; //histogram which has schematic drawing of all sieve-holes

    std::vector<SieveHoleData> fSieveHoleData; 

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
/*    void Evaluate();        //execute the evaluation loop for the selected sievehole
    void Deselect();        //de-select the current sieve hole */ 
    void HandleCanvasClick_data();    //handle the canvas being clicked (data histogram)
    void HandleCanvasClick_drawing(); //handle the canvas beign clicked (drawing histogram)

    void DrawSieveHoles(); 

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

    // Set up the main frame
    SetCleanup(kDeepCleanup);
    
    //setup the canvas, and split it 
    // Create embedded canvas

    fFrame_canv = new TGHorizontalFrame(this, 1400, 700); 

    fEcanvas_data    = new TRootEmbeddedCanvas("ECanvas_data", this, 700, 700);
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
    TGHorizontalFrame* bframe = new TGHorizontalFrame(this, 1400, 50);
    
    auto Add_button = [this](   TGHorizontalFrame *frame, 
                                TGTextButton* button, 
                                string button_label, 
                                string method, 
                                TGLayoutHints* hints    ) 
    {
        button_label = "&" + button_label; 
        button = new TGTextButton(frame, button_label.c_str(), 1); 
        button->Connect("Clicked()", "PickSieveHoleApp", this, method.c_str()); 
        frame->AddFrame(button, hints); 
    };

    Add_button(bframe, fButton_Exit,  "Exit",  "DoExit()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    ); 

    AddFrame(bframe, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 5, 10));
    
    SetWindowName("Hole-isolation interactive tool");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

#if 0
    //create horizontal frame for writing the output file
    // Create horizontal frame for labels and button
    TGHorizontalFrame* bframe = new TGHorizontalFrame(this, 1200, 50);

    Add_button(wframe, fWriteButton, "Write Output",  "WriteOutput()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    ); 

    fFilePathEntry = new TGTextEntry(wframe);
    wframe->AddFrame(fFilePathEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 20, 20, 5, 5)); 

    AddFrame(wframe, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 5, 10)); 

    // Create horizontal frame for labels and button
    TGHorizontalFrame* bframe = new TGHorizontalFrame(this, 1200, 50);
    
    // Create coordinate labels
    fXLabel = new TGLabel(bframe, "X: ---");
    fXLabel->SetTextJustify(kTextLeft);
    bframe->AddFrame(fXLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 10, 20, 5, 5));
    
    fYLabel = new TGLabel(bframe, "Y: ---");
    fYLabel->SetTextJustify(kTextLeft);
    bframe->AddFrame(fYLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 20, 20, 5, 5));

    Add_button(bframe, fExitButton,  "Exit",          "DoExit()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    ); 
    Add_button(bframe, fUndoButton,  "Undo",          "DeleteLast()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    );
    Add_button(bframe, fClearButton, "Clear",         "DeleteAll()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    );

    AddFrame(bframe, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 5, 10));
    
    // Draw histogram on canvas
    if (gStyle) {
        gStyle->SetOptStat(0); 
        gStyle->SetPalette(palette); 
    }
    TCanvas* canvas = fEcanvas->GetCanvas();
    canvas->cd();
    fHistogram->Draw(drawing_option);
    canvas->Update();
    
    // Connect canvas click event
    canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "PickSieveHoleApp", this, 
                    "HandleCanvasClick()");
    
    // Set window properties
#endif 
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
        .fill_style = 3001, //dense dotted style
    }; 

    //hole hole has already been evaluated
    const CircleStyle style_evaluated {
        .line_color = 1,   
        .line_style = 1, 
        .fill_color = 17,   //grey 
        .fill_style = 3004, //diagonal hatching style
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

    //assumes that the status of all sieve holes is recorded correctly 
    const double draw_size_big   = 3.0e-3; 
    const double draw_size_small = 1.5e-3;

    //draw all sieve-holes
    for (SieveHoleData& hole_data : fSieveHoleData) {
        auto& circ = hole_data.draw_circ;

        double draw_rad = hole_data.hole.is_big ? draw_size_big : draw_size_small; 

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
        if (hole_data.selected) {
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
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_data() {
    // Get canvas and event information
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

             // Get pixel coordinates
            Int_t px = canvas->GetEventX();
            Int_t py = canvas->GetEventY();
            
            // Convert pixel coordinates to histogram coordinates
            Double_t x = canvas->AbsPixeltoX(px);
            Double_t y = canvas->AbsPixeltoY(py);
            
            printf("Drawing histogram clicked: %+.6f, %+.6f \n", x, y); cout << endl; 
        }
    }
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