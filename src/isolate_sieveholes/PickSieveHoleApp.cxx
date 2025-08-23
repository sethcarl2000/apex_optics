#include "PickSieveHoleApp.h"
#include <TGClient.h>
#include <TCanvas.h> 
#include <TStyle.h>
#include <cstdio> 
#include "EvaluateCutFrame.h"
#include <TApplication.h>
#include "SaveOutputFrame.h"

using namespace std; 

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
    fBranchX(coordname_x),
    fBranchY(coordname_y),
    fPathInfile(path_infile),
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
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 

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

    //compute the react-vertex
    fReactVertex = TVector3(
        *(fRDF->Define("x", [](const TVector3& v){ return v.x(); }, {"position_vtx_scs"}).Mean("x")),
        *(fRDF->Define("y", [](const TVector3& v){ return v.y(); }, {"position_vtx_scs"}).Mean("y")),
        *(fRDF->Define("z", [](const TVector3& v){ return v.z(); }, {"position_vtx_scs"}).Mean("z"))
    ); 

    printf("avg. react vertex: (% 4.1f, % 4.1f, % 4.1f) (mm)", fReactVertex.x(), fReactVertex.y(), fReactVertex.z()); cout << endl; 

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
    for (const SieveHole& hole : SieveHole::ConstructSieveHoles(f_is_RHRS)) fSieveHoleData.emplace_back(hole); 

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

    SetWindowName(Form("Sieve-hole selection tool. data: '%s'", path_infile));
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
    if (fCurrentWindow != kWindow_PickSieveHole) return; 
    
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
    if (fCurrentWindow != kWindow_PickSieveHole) return; 

    //check if there is a hole picked. 
    if (!fSelectedSieveHole) {
        Error("DoEvaluate", "No sieve hole is selected. It shoul not have been possible to get here..."); 
        return; 
    }

    //check if our status is 'kReadyToEval' 
    if (fCurrentPickStatus != kReadyToEval) {
        Error("DoEvaluate", "Tried to eval when fCurrentPickStatus is not 'kRreadyToEval'. It should not have been possible to get here..."); 
        return; 
    } 

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
    if (fCurrentWindow != kWindow_PickSieveHole) return;  

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
    if (fCurrentWindow != kWindow_PickSieveHole) return; 
    
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

    if (fCurrentWindow != kWindow_PickSieveHole) return; 

    new SaveOutputFrame(
            gClient->GetRoot(), 
            900, 500,  
            fSieveHoleData,
            Get_IsRHRS(),
            GetReactVertex()
        );
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_data() {
    
    if (fCurrentWindow != kWindow_PickSieveHole) return; 

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
            fCurrentPickStatus = kReadyToEval; 

            DrawSieveHoles(); 
            DrawHoleCuts(); 
            UpdateButtons();
         }// if (fEventType == kMouseButton1_up) 
    }// if (fEventType != canvas->GetEvent()) 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_drawing() {
    
    if (fCurrentWindow != kWindow_PickSieveHole) return; 
    
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
            DrawSieveHoles(); 
            DrawHoleCuts(); 
            UpdateButtons(); 
        }// if (fEventType == kMouseButton1_up) 
    }// if (fEventType != canvas->GetEvent()) 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DeselectSieveHole()
{
    if (fCurrentWindow != kWindow_PickSieveHole) return; 
    
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

ClassImp(PickSieveHoleApp); 