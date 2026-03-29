#include "PickSieveHoleApp.h"
#include <TGClient.h>
#include <TCanvas.h> 
#include <TStyle.h>
#include <cstdio> 
#include "EvaluateCutFrame.h"
#include <TApplication.h>
#include "SaveOutputFrame.h"
#include <SieveHole.h> 
#include <isolate_sieveholes/SieveHoleData.h> 
#include <ApexOptics.h> 

using namespace std; 
using ApexOptics::Trajectory_t; 

//ptr to singleton instance of PickSieveHoleApp
PickSieveHoleApp* PickSieveHoleApp::fInstance = nullptr; 

namespace {

    // base unit, from which we define all others
    constexpr double meter = 1.; 

    // value of 1 mm, expressed in m
    constexpr double mm = 1/1000.; 

    // the value of one degree, expressed in radians 
    constexpr double deg = 3.14159265359 / 180.; 
}

//______________________________________________________________________________________________________________________
PickSieveHoleApp::PickSieveHoleApp( const TGWindow* p, 
                                    UInt_t w, 
                                    UInt_t h, 
                                    const bool is_RHRS, 
                                    const char* path_infile,
                                    const char* target_name,  
                                    const char* coordname_x,
                                    const char* coordname_y,
                                    const char* drawing_option, 
                                    unsigned int palette) 
    : fptr_TGWindow{p},
    fWindow_width{w},
    fWindow_height{h}, 
    f_is_RHRS(is_RHRS), 
    fPath_infile{path_infile},
    fBranchX(coordname_x),
    fBranchY(coordname_y),
    fPathInfile(fPath_infile.c_str()),
    fDrawingOption(drawing_option),
    fPalette{palette}
{
    //label ourselves as the one singleton instance of this class 
    fInstance = this; 
}

void PickSieveHoleApp::LaunchApplication()
{
    TGMainFrame(fptr_TGWindow, fWindow_width, fWindow_height); 
    fEventType=-1; 

    //to begin with, we are in the 'pick sieve hole' window
    fCurrentWindow     = kWindow_PickSieveHole;  
    fCurrentPickStatus = kNoneSelected; 

    // Set up the main frame    
    SetCleanup(kDeepCleanup);
    
    //setup the canvas, and split it 
    // Create embedded canvas
    // fGFrame_pick    = new TGHorizontalFrame(this, 1400, 700); 

    fFrame_canvPick = new TGHorizontalFrame(this, 1400, 700); 

    fEcanvas_data   = new TRootEmbeddedCanvas("ECanvas_data", fFrame_canvPick, 700, 700);
    fFrame_canvPick->AddFrame(fEcanvas_data, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 5)); 
    
    TCanvas* canvas = fEcanvas_data->GetCanvas();
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(fPalette); 

    canvas->cd(); 
    //first things first, lets set up & draw the histogram: 
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
    //ROOT::EnableImplicitMT(); 
    fRDF = (ROOT::RDF::RNode*)(new ROOT::RDataFrame(fTreeName.c_str(), fPath_infile.c_str())); 
    ROOT::RDataFrame df(fTreeName.data(), fPath_infile.c_str());
    
    const double rast_min = *df.Define("y", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}).Min("y");
    const double rast_max = *df.Define("y", [](TVector3 vtx){ return vtx.y(); }, {"position_vtx"}).Max("y");
    
    fRastMin = rast_min;
    fRastMax = rast_max; 

    // if 'true', then we'll use the optics model to define sieve-coordinates anew. if not, then we wont. 
    const bool use_optics_model = (fOpticsModel != nullptr); 

    if (use_optics_model) {

        auto df_with_sv_data = fOpticsModel->DefineOutputs(df); 

        fEventData = *df_with_sv_data

            .Define("Xsv",  [](double x, double y, double dxdz, double dydz, double dpp)
            {
                return Trajectory_t{x,y,dxdz,dydz,dpp}; 
            }, {"reco_x_sv","reco_y_sv","reco_dxdz_sv","reco_dydz_sv","reco_dpp_sv"}) 

            .Define("Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}) 

            .Define("rast_index", [rast_min,rast_max](TVector3 vtx)
            {
                return  ( vtx.y() - 0.5*(rast_max+rast_min) )/( 0.5*(rast_max-rast_min) ); 
            }, {"position_vtx"})

            .Define("EventData", [](Trajectory_t Xsv, Trajectory_t Xfp, double rast_index, TVector3 vtx_scs)
            {
                return EventData{ .Xsv=Xsv, .Xfp=Xfp, .raster_index=rast_index, .vtx_scs=vtx_scs }; 
            }, {"Xsv", "Xfp", "rast_index", "position_vtx_scs"})

            .Take<EventData>("EventData"); 

    } else {
        
        fEventData = *df 

            .Define("Xsv", [this, use_optics_model](double x, double y, double dxdz, double dydz, double dpp)
            {   
                return Trajectory_t{x,y,dxdz,dydz,dpp}; 
            }, {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"})

            .Define("Xfp", [](double x, double y, double dxdz, double dydz)
            {
                return Trajectory_t{x,y,dxdz,dydz}; 
            }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}) 

            .Define("rast_index", [rast_min,rast_max](TVector3 vtx)
            {
                return  ( vtx.y() - 0.5*(rast_max+rast_min) )/( 0.5*(rast_max-rast_min) ); 
            }, {"position_vtx"})

            .Define("EventData", [](Trajectory_t Xsv, Trajectory_t Xfp, double rast_index, TVector3 vtx_scs)
            {
                return EventData{ .Xsv=Xsv, .Xfp=Xfp, .raster_index=rast_index, .vtx_scs=vtx_scs }; 
            }, {"Xsv", "Xfp", "rast_index", "position_vtx_scs"})

            .Take<EventData>("EventData"); 
    }
    

    fSieveHist = new TH2D(
        "h_sieve_cpy", 
        "Sieve-coordinates;dx/dz_{sv};dy/dz_{sv}", 
        200, GetDrawRange_x_min(), GetDrawRange_x_max(), 
        200, GetDrawRange_y_min(), GetDrawRange_y_max()
    ); 
    
    for (const auto& ev : fEventData) fSieveHist->Fill( ev.Xsv.dxdz, ev.Xsv.dydz ); 

    //compute the react-vertex
    fReactVertex = TVector3(
        *(fRDF->Define("x", [](const TVector3& v){ return v.x(); }, {"position_vtx"}).Mean("x")),
        *(fRDF->Define("y", [](const TVector3& v){ return v.y(); }, {"position_vtx"}).Mean("y")),
        *(fRDF->Define("z", [](const TVector3& v){ return v.z(); }, {"position_vtx"}).Mean("z"))
    ); 

    printf("avg. react vertex: (% 4.1f, % 4.1f, % 4.1f) (mm)", 
        fReactVertex.x() *1e3, 
        fReactVertex.y() *1e3, 
        fReactVertex.z() *1e3
    ); cout << endl; 

    //make sure this histogram is not 'attached' to anything (won't be deleted by anything else). 
    fSieveHist->SetDirectory(0); 
    
    //connect the canvas to the 'HandleCanvasClick_data()' event
    canvas->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 
        "PickSieveHoleApp",
        this,
        "HandleCanvasClick_data()"
    ); 

    //construct the vector of sieve holes...
    for (const SieveHole& hole : ApexOptics::ConstructSieveHoles(f_is_RHRS)) fSieveHoleData.emplace_back(hole);
    

    fEcanvas_drawing = new TRootEmbeddedCanvas("ECanvas_drawing", fFrame_canvPick, 700, 700);
    fFrame_canvPick->AddFrame(fEcanvas_drawing, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0)); 

    canvas = fEcanvas_drawing->GetCanvas(); 
    canvas->cd(); 

    //now, create the 'sieve hole drawing' histogram, and draw it. 
    fHoleDrawingHist = new TH2D(
        "h_hole_drawing", 
        "Sieve-hole positions;x_{sv};y_{sv}", 
        200, -0.050, +0.050, 
        200, -0.045, +0.045
    ); 
    //draw thie histogram, draw the canvas
    fHoleDrawingHist->SetDirectory(0); 

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
    
    //Evaluate button
    fButton_Evaluate    = new TGTextButton(fFrame_PickHoleButtons, "&Evaluate", 1); 
    fButton_Evaluate->Connect("Clicked()", "PickSieveHoleApp", this, "DoEvaluate()"); 
    fFrame_PickHoleButtons->AddFrame(fButton_Evaluate, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //Reset Plots button
    fButton_ResetPlots  = new TGTextButton(fFrame_PickHoleButtons, "&Reset Plots", 1); 
    fButton_ResetPlots->Connect("Clicked()", "PickSieveHoleApp", this, "DoResetPlots()"); 
    fFrame_PickHoleButtons->AddFrame(fButton_ResetPlots, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 


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
    fNumber_cutWidth->Connect("ValueSet(Long_t)", "PickSieveHoleApp", this, "SetCutSizeAngle()"); 
    fNumber_cutWidth->GetNumberEntry()->Connect("ReturnPressed()", "PickSieveHoleApp", this, "SetCutSizeAngle()"); 
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
    fNumber_cutHeight->Connect("ValueSet(Long_t)", "PickSieveHoleApp", this, "SetCutSizeAngle()"); 
    fNumber_cutHeight->GetNumberEntry()->Connect("ReturnPressed()", "PickSieveHoleApp", this, "SetCutSizeAngle()"); 
    fNumber_cutHeight->SetButtonToNum(false); //setting this to 'false' ensures that the 'ValueChanged' signal is sent to 'SetCutSize'
    fFrame_numbers->AddFrame(fNumber_cutHeight, new TGLayoutHints( kLHintsLeft | kLHintsCenterY, 20, 10, 5, 5)); 

    // angle slider ----------------------------------------------------------------------------------------------------
    fLabel_cutAngle = new TGLabel(fFrame_numbers, Form("Cut angle: +00",0.)); 
    fFrame_numbers->AddFrame(fLabel_cutAngle, new TGLayoutHints( kLHintsLeft | kLHintsCenterY, 20, 5, 5, 5)); 
    // slider
    fSlider_cutAngle = new TGHSlider(fFrame_numbers, 20); 
    fSlider_cutAngle->Connect("Released()", "PickSieveHoleApp", this, "SetCutSizeAngle()"); 
    fSlider_cutAngle->SetPosition( (fSlider_cutAngle->GetMinPosition() + fSlider_cutAngle->GetMaxPosition())/2. ); 
    fFrame_numbers->AddFrame(fSlider_cutAngle, new TGLayoutHints( kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 20, 10, 5, 5)); 
    
    


    fLabel_selectedHole = new TGLabel(fFrame_numbers, "No hole selected"); 
    fFrame_numbers->AddFrame(fLabel_selectedHole, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 20,0,0,0)); 
    
    //add the number frame 
    fFrame_PickHoleButtons->AddFrame(fFrame_numbers, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 0,0,0,0)); 


    AddFrame(fFrame_PickHoleButtons, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 0, 0, 5, 10));
    
    
    SetWindowName(Form("Sieve-hole selection tool. data: '%s'", fPath_infile.c_str()));
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    fAppState = kSelectHole; 
    DoResetPlots(); 
    UpdateButtons(); 

}
//_____________________________________________________________________________________________________________________________________
PickSieveHoleApp::~PickSieveHoleApp() 
{
    // Cleanup
    if (fSieveHist)       delete fSieveHist; 
    if (fHoleDrawingHist) delete fHoleDrawingHist; 
    Cleanup();
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoResetPlots() 
{
    if (fAppState & kAppDisabled) return; 

    //disable user actions first
    const auto app_state = fAppState; 
    fAppState = AppState{fAppState | kAppDisabled};   

    //re-draw the data histogram
    auto canv_data = fEcanvas_data->GetCanvas(); 
    
    if (!canv_data) {
        Error(__func__, "Canvas for data-histogram is null here... something is wrong.");
        return; 
    }
    canv_data->Clear(); 
    canv_data->cd(); 
    fSieveHist->DrawCopy(fDrawingOption); 
    DrawHoleCuts(); 
    canv_data->Modified(); canv_data->Update(); 

    //re-draw the drawing histogram 
    auto canv_drawing = fEcanvas_drawing->GetCanvas(); 

    if (!canv_drawing) {
        Error(__func__, "Canvas for drawing-histogram is null here... something is wrong.");
        return; 
    }
    canv_drawing->Clear(); 
    canv_drawing->cd(); 
    fHoleDrawingHist->DrawCopy(); 
    DrawSieveHoles(); 
    canv_drawing->Modified(); canv_drawing->Update(); 

    //return the app to its previous state
    fAppState = app_state; 
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DrawSieveHoles()
{
    const char* const here = "DrawSieveHoles"; 
#ifdef DEBUG
    Info(here, "In %s body...", here); 
#endif

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
#ifdef DEBUG
    cout << "Entering SieveHoleData Loop..." << endl;  
#endif


    for (SieveHoleData& hole_data : fSieveHoleData) {

        auto& circ = hole_data.draw_circ;

        double draw_rad = hole_data.hole.is_big ? fDrawSize_big : fDrawSize_small;

        //draw the the new circle with the given hole-position and radius
        if (circ==nullptr) {
            circ = new TEllipse(
                hole_data.hole.x, 
                hole_data.hole.y, 
                draw_rad,
                draw_rad
            ); 
        } 

       
        //set line/fill style based on the status of the hole (has it been selected / drawn / etc.)
        if ((fSelectedSieveHole != nullptr) && (hole_data == *fSelectedSieveHole)) {
            set_circle_style( circ, style_selected ); 
        } else {
            if (hole_data.is_evaluated) { set_circle_style( circ, style_evaluated ); }
            else                        { set_circle_style( circ, style_default ); } 
        }
        
        circ->Draw();  
        //we use this line to tell root that WE will manage this object, so it should not try to delete it itself.
        //circ->SetBit(kCanDelete, false); 
    }
#ifdef DEBUG
    cout << "Done with loop." << endl; 
#endif 

#ifdef DEBUG
    Info(here, "Calling canv->Modified()...");
#endif 
    canv->Modified(); 
#ifdef DEBUG
    Info(here, "Calling canv->Update()...");
#endif 
    canv->Update(); 
#ifdef DEBUG
    Info(here, "Exiting %s body.", here); 
#endif
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::SetCutSizeAngle()
{
    //igore button command if the app is disabled
    if (fAppState & kAppDisabled) return; 

    if (fCurrentWindow != kWindow_PickSieveHole) return; 
    
    if (!fNumber_cutHeight || !fNumber_cutHeight || !fSlider_cutAngle) {
        throw logic_error("in <PickSieveHoleApp::SetCutSize>: One or both of number-entry objects are null."); 
        return; 
    }
    
    // get the updated cut angle 
    double slider_amplitude = (double)(fSlider_cutAngle->GetMaxPosition() - fSlider_cutAngle->GetMinPosition());
    double slider_val       = (double)(fSlider_cutAngle->GetPosition()    - fSlider_cutAngle->GetMinPosition());
    double index = slider_val/slider_amplitude; 

    // angle can be anywhere on the interval -90 => +90
    double angle = (index - 0.5) * fCutAngle_max * 2; 

    fLabel_cutAngle->ChangeText(Form("Cut angle: %+.0f",angle)); 

    //if no sieve-hole is selected, we don't need to do anything right now. 
    if (!fSelectedSieveHole) return; 

    if (!fSelectedSieveHole->is_evaluated) {
        
        // get the updated cut width 
        fSelectedSieveHole->cut_width  = fNumber_cutWidth ->GetNumber()*mm;
        fSelectedSieveHole->cut_height = fNumber_cutHeight->GetNumber()*mm;

        // min angle is 
        fSelectedSieveHole->cut_angle = angle; 
        
    }

    DrawHoleCuts(); 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DrawHoleCuts() 
{
    const char* const here = "DrawHoleCuts"; 

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

    //a few helper functions we need 
    //returns true if the passed double is NaN, false otherwise 
    auto is_nan = [](double x){ return bool(x!=x); }; 

    //returns true if the passed TObject is in the canvas' list of primitives, false otherwise. 
    auto is_owned_by_TCanvas = [canv](TObject* obj){
        for (const auto prim : *canv->GetListOfPrimitives()) { if (prim == obj) return true; }
        return false; 
    }; 

    //delete all ellipses from the pad, so we can re-draw them 
    
#ifdef DEBUG
    Info(here, "About to loop through all sieve-holes"); 
#endif 

    //draw all sieve-holes
    for (SieveHoleData& hole_data : fSieveHoleData) {

        auto& circ = hole_data.hole_cut;
#ifdef DEBUG
        //Info(here, "Checking to see if cuts are null... circ ptr: %p", circ); 
#endif  
        //check if the cuts have been defiend for this hole
        if (is_nan(hole_data.cut_x) ||
            is_nan(hole_data.cut_y) ||
            is_nan(hole_data.cut_width) ||
            is_nan(hole_data.cut_height)) {

            //check to see if this circle's been drawn. if it has, we need to remove it. 
            if (circ) {

#ifdef DEBUG 
                Info(here, "Found null hole with extant circ (%p). Removing from canvas list of primitives..."); 
#endif  
                DeleteDrawnObject(circ); 
            }
            continue; //move on to the next hole
        }

#ifdef DEBUG
        Info(here, "Found non-null hole. Checking to see if 'circ' ptr exists..."); 
#endif 
        //if circ has not been drawn, then make a new one. 
        if (circ==nullptr) {
#ifdef DEBUG
            Info(here, "Creating new Circ object..."); 
#endif 
            circ = new TEllipse(
                hole_data.cut_x, 
                hole_data.cut_y, 
                hole_data.cut_width,
                hole_data.cut_height
            ); 
        } 
#ifdef DEBUG
        else {
            Info(here, "Circ already exists.");
        } 
#endif 
        //position of the cut
        circ->SetX1(hole_data.cut_x); 
        circ->SetY1(hole_data.cut_y); 

        //size of the cut's major axes
        circ->SetR1(hole_data.cut_width); 
        circ->SetR2(hole_data.cut_height); 

        //the hole's orientation from the horizontal (angle of 0 means the 'cut_width' axis will be horizontal)
        circ->SetTheta(hole_data.cut_angle); 

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
#ifdef DEBUG
    Info(here, "Modifying & updating canvas."); 
#endif 
    canv->Modified(); 
    canv->Update(); 
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoEvaluate() 
{
    //igore button command if the app is disabled
    if (fAppState & kAppDisabled) return; 

    //ignore this command if we aren't in the right state 
    if (!(fAppState & kReadyToEvaluate)) return; 

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
    
    new EvaluateCutFrame(
        gClient->GetRoot(), 
        fEventData, 
        fSelectedSieveHole, 
        fFpcoord_cut_width, 
        "dxdz_sv", "dydz_sv", "col", 
        kBird
    );
    
    fAppState = kEvaluateCutFrame; 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoneEvaluate()
{
    const char* const here = "PickSieveHoleApp::DoneEvaluate"; 
#ifdef DEBUG
    Info(here, "Enter DoneEvaluate() body"); 
#endif 
    //once the EvaluateCutFrame object is done working, then we reutrn here. 
    //now, we've evaluated this sievehole

    if (!fSelectedSieveHole) {
        Warning(__func__, "no sieve hole is selected; something must have gone wrong."); 
        DrawSieveHoles(); 
#ifdef DEBUG
        Info(here, "Calling DrawHoleCuts...");
#endif
        DrawHoleCuts(); 
#ifdef DEBUG
        Info(here, "Calling UpdateButtons...");
#endif
        UpdateButtons(); 
        return; 
    } 
#ifdef DEBUG
    else { cout << "fSelectedSieveHole : " << fSelectedSieveHole << endl; }
#endif 

    //this might be changed by the EvaulateCutFrame app
    gStyle->SetPalette(fPalette); 

    if (fSelectedSieveHole->is_evaluated) {

        cout << "Evaluated hole: row " << fSelectedSieveHole->hole.row << ", col " << fSelectedSieveHole->hole.col << endl; 
        
    } else { 

        cout << "Rejected hole: row " << fSelectedSieveHole->hole.row << ", col " << fSelectedSieveHole->hole.col << endl; 
    }

    //do re-drawing of buttons
    fCurrentWindow     = kWindow_PickSieveHole; 
    
#ifdef DEBUG 
    Info(here, "Calling DeselectSieveHole..."); 
#endif 
    DeselectSieveHole(); 
#ifdef DEBUG
    Info(here, "Returned."); 
#endif

    //redraw the sieve holes, and update the button layout
#ifdef DEBUG
    Info(here, "Calling DrawSieveHoles...");
#endif
    DrawSieveHoles(); 
#ifdef DEBUG
    Info(here, "Calling DrawHoleCuts...");
#endif
    DrawHoleCuts(); 
#ifdef DEBUG
    Info(here, "Calling UpdateButtons...");
#endif
    UpdateButtons(); 
#ifdef DEBUG
    Info("DoneEvaluate", "Leaving DoneEvaluate() body."); 
#endif 
    fAppState = kSelectHole; 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoDelete() 
{
    //igore button command if the app is disabled
    if (fAppState & kAppDisabled) return; 

    if (fCurrentWindow != kWindow_PickSieveHole) return;  

    //check if there is a hole picked. 
    if (!fSelectedSieveHole) {
        Error("DoDelete", "No sieve hole is selected. It shoul not have been possible to get here..."); 
        return; 
    }

    //check if our status is 'kReadyToEval' 
    if (!(fAppState & kAlreadyEvaluated)) {
        Error(__func__, "Tried to eval when fCurrentPickStatus is not 'kOldSelected'. It should not have been possible to get here..."); 
        return; 
    }

    //create a new, 'fresh' SieveHoleData, in which everything is default-initialized except the SieveHole struct. 
    fSelectedSieveHole->Clear(); 

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
    fButton_ResetPlots  ->UnmapWindow(); 
    fButton_Save_Output ->UnmapWindow(); 
    fButton_Exit        ->UnmapWindow(); 

    //this button will be 'valid' in any state: 
    fButton_Exit->MapWindow(); 
    fButton_ResetPlots->MapWindow(); 
    fButton_Save_Output->MapWindow(); 
    
    //now, figure out which buttons should be drawn
    switch (fAppState) {
        case kSelectHole : break; 
        case kAlreadyEvaluated : { fButton_Delete->MapWindow(); break; }
        case kPickCutOnPlot : break; 
        case kReadyToEvaluate : { fButton_Evaluate->MapWindow(); break; }
        default : break; 
    }
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DeleteDrawnObject(TObject* obj) {       

#ifdef DEBUG
    Info(__func__, "Enter fcn body"); 
#endif 
    auto pad_data    = (TVirtualPad*)fEcanvas_data->GetCanvas();
    auto pad_drawing = (TVirtualPad*)fEcanvas_drawing->GetCanvas();  
    if (!pad_data || !pad_drawing) {  
        Warning(__func__, "ptrs to Canvases are null... something isn't right here."); 
        return;  
    }  

    if (!obj) {
#ifdef DEBUG
        printf("Info in <%s>: Ptr to TObject passed is null.\n", __func__);          
#endif
        return; 
    }

#ifdef DEBUG
    printf("Info in <%s>: Attepmting to delete TObject with name '%s'\n", __func__, obj->GetName());          
#endif 
    
    //remove the object from the TVirtualPad (and any sub-pads recursively)
    //we don't need to delete it, because (according to claude) ROOT's TPad system owns the TObjects drawn on the pads 
    pad_data   ->RecursiveRemove(obj); 
    pad_drawing->RecursiveRemove(obj);

#ifdef DEBUG
    Info(__func__, "Exit fcn body"); 
#endif
}
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
void PickSieveHoleApp::WriteOutput() 
{
    //igore button command if the app is disabled
    if (fAppState & kAppDisabled) return; 

    if (fCurrentWindow != kWindow_PickSieveHole) return; 

    //disable the app
    fAppState = AppState{fAppState | kAppDisabled}; 

    new SaveOutputFrame(
            gClient->GetRoot(), 
            900, 500,  
            fSieveHoleData,
            fRDF,
            Get_IsRHRS(),
            GetReactVertex(),
            fRastMin, fRastMax, 
            fFpcoord_cut_width
        );
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DoneWriteOutput()
{
    //re-enable interactive elements of the app 
    fAppState = AppState{fAppState & (~kAppDisabled)}; 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_data() 
{
    //igore button command if the app is disabled
    if (fAppState & kAppDisabled) return; 

    //get canvas event information
    TCanvas* canvas = fEcanvas_data->GetCanvas();
    if (!canvas) return;

    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType =  canvas->GetEvent(); 

        //check if the status is valid
        if (!(fAppState & (kPickCutOnPlot | kReadyToEvaluate))) return; 
        
        if (!fSelectedSieveHole) return;     

        //if the mouse-button was just released
        if (fEventType == kMouseButton1_up) {

            //fSelectedSieveHole->cut_width  = fNumber_cutWidth ->GetNumber() * 1e-3; 
            //fSelectedSieveHole->cut_height = fNumber_cutHeight->GetNumber() * 1e-3; 
            //update the selected sieve-hole cut size & cut orientation 
            SetCutSizeAngle(); 

            // Get pixel coordinates
            Int_t px = canvas->GetEventX();
            Int_t py = canvas->GetEventY();
            
            // Convert pixel coordinates to histogram coordinates
            Double_t x = canvas->AbsPixeltoX(px);
            Double_t y = canvas->AbsPixeltoY(py);
            
            fSelectedSieveHole->cut_x = x;
            fSelectedSieveHole->cut_y = y;              
            fCurrentPickStatus = kReadyToEval; 

            //mark the state as ready to evaluate 
            fAppState = kReadyToEvaluate; 

            DrawSieveHoles(); 
            DrawHoleCuts(); 
            UpdateButtons();
            return; 
         }// if (fEventType == kMouseButton1_up) 
    }// if (fEventType != canvas->GetEvent()) 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::HandleCanvasClick_drawing() 
{
    //igore button command if the app is disabled
    if (fAppState & kAppDisabled) return; 

    const char* const here = "HandleCanvasClick_drawing"; 

    // Get canvas and event information
    TCanvas* canvas = fEcanvas_drawing->GetCanvas();
    if (!canvas) return;

    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType = canvas->GetEvent(); 

        //if the mouse-button was just released
        if (fEventType == kMouseButton1_up) {

#ifdef DEBUG
            Info(here, "Handling kMouseButton1_up event."); 
#endif 

            //deselect any selected sievehole 
#ifdef DEBUG
            cout << "Calling DeselectSieveHole..." << endl; 
#endif 
            DeselectSieveHole(); 
            fAppState = kSelectHole; 
#ifdef DEBUG
            cout << "Done." << endl; 
#endif 

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
                        fAppState = kAlreadyEvaluated; 
                    } else {
                        fAppState = kPickCutOnPlot;     
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
#ifdef DEBUG
            Info(here, "Calling 'DrawSieveHoles'"); 
#endif
            DrawSieveHoles(); 
#ifdef DEBUG
            Info(here, "Calling 'DrawHoleCuts'"); 
#endif
            DrawHoleCuts(); 
#ifdef DEBUG
            Info(here, "Calling 'UpdateButtons'"); 
#endif
            UpdateButtons(); 
#ifdef DEBUG
            Info(here, "Done, exiting.'"); 
#endif
        }// if (fEventType == kMouseButton1_up) 
    }// if (fEventType != canvas->GetEvent()) 
}
//_____________________________________________________________________________________________________________________________________
void PickSieveHoleApp::DeselectSieveHole()
{
    //deselect any selected sieve holes
    //check if we've already loaded a bit of data into this hole
    if (fSelectedSieveHole && !fSelectedSieveHole->is_evaluated) {

        fSelectedSieveHole->Clear(); 

        //clean up the sieve hole by deleting this one
        //fSelectedSieveHole->~SieveHoleData(); 
        
        //this hole has not been evaluated. therefore, we need to make sure that we 'clean it up' when we deselect. 
        //but if it *has* been evaluated, we don't want to wipe the data that has already been recorded. 
        
        //*fSelectedSieveHole = clean_slate; 
    }
    fLabel_selectedHole->SetText("No hole selected"); 

    fSelectedSieveHole = nullptr; 
    fCurrentPickStatus = kNoneSelected; 
}

ClassImp(PickSieveHoleApp); 