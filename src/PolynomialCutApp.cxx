#include "PolynomialCutApp.h"
#include "TStyle.h"

using namespace std; 

//______________________________________________________________________________________________________________________
PolynomialCutApp::PolynomialCutApp( const TGWindow* p, 
                                    UInt_t w, 
                                    UInt_t h, 
                                    TH2 *hist, 
                                    const char* drawing_option, 
                                    unsigned int palette) 
    : TGMainFrame(p, w, h) 
{
    
    fEventType = -1; 
    fLineColor = kRed; 

    // Set up the main frame
    SetCleanup(kDeepCleanup);
    
    fHistogram = (TH2D*)hist->Clone("histogram_clone");
    fHistogram->SetDirectory(0); // Detach from file directory
    
    // Create embedded canvas
    fEcanvas = new TRootEmbeddedCanvas("ECanvas", this, 1200, 800);
    AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 5));
    
    auto Add_button = [this](TGHorizontalFrame *frame, 
                                     TGTextButton* button, 
                                     string button_label, 
                                     string method, 
                                     TGLayoutHints* hints) 
    {
        button_label = "&" + button_label; 
        button = new TGTextButton(frame, button_label.c_str(), 1); 
        button->Connect("Clicked()", "PolynomialCutApp", this, method.c_str()); 
        frame->AddFrame(button, hints); 
    };

    //create horizontal frame for writing the output file
    // Create horizontal frame for labels and button
    TGHorizontalFrame* wframe = new TGHorizontalFrame(this, 1200, 50);

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
    canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "PolynomialCutApp", this, 
                    "HandleCanvasClick()");
    
    // Set window properties
    SetWindowName("PolynomialCut creation tool");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}
//_____________________________________________________________________________________________________________________________________
PolynomialCutApp::~PolynomialCutApp() {
    // Cleanup
    if (fHistogram) delete fHistogram;
    Cleanup();
}
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::CloseWindow() {
    // Called when window is closed via window manager
    gApplication->Terminate(0);
}
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::DoExit() {
    // Called when exit button is pressed
    gApplication->Terminate(0);
}
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::WriteOutput() 
{
    const char* const here = "PolynomialCutApp::WriteOutput"; 
    //must have at least 3 vertices to write the output
    if (fVertices.size()<3) return; 
    
    TGTextBuffer *buffer = fFilePathEntry->GetBuffer(); 
    if (!buffer) return; 

    string path_outfile(buffer->GetString());
    
    //let's try to create a new PolynomialCut object from the given output-file name
    PolynomialCut polycut; 

    //we'll do a try-catch block here, as an exception will be thrown in this step if any of the segments overlap 
    try {   
        
        polycut.AddVertices(fVertices); 

    } catch (const PolynomialCut::InvalidVertexException& e) {

        cerr << "\nError in <" << here << ">:\n"
                " -- One or more of the segments in the PolynomialCut intersect, which is not allowed.\n"
                " -- Try removing bad vertices, and click 'Write Output' again.\n"
                " -- Exception thrown: '" << e.what() << "'" << endl; 
        return; 
    }

    //if we got here, it means that all vertices were added without exception (PolynomialCut checks each new vertex to make sure
    // that the segments generated don't intersect with any others). 

    //now, try and create the output file. 
    try {

        polycut.Create_dbfile(path_outfile.c_str()); 
    
    } catch (const std::exception& e) {

        cerr << "\nError in <" << here << ">:\n"
                " -- An issue was encountered when attempting to write the output db-file.\n"
                " -- Exception thrown: '" << e.what() << "'" << endl; 
        return;
    } 

    cout << "Info in <" << here << ">:\n" 
            " -- wrote PolynomialCut db file to '" << path_outfile << "'" << endl; 
    return; 
};
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::DeleteLast() 
{
    //delete the last point
    if (!fVertices.empty()) fVertices.erase( fVertices.end()-1 ); 
    RedrawLines(); 
}; 
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::DeleteAll() 
{
    //delete all points
    fVertices.clear(); 
    RedrawLines(); 
}; 
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::RedrawLines() 
{
    TCanvas* canvas = fEcanvas->GetCanvas();
    if (!canvas) return;
    canvas->cd(); 

    //erase all lines which previously existed
    for (auto line : fLines) delete line; 
    fLines.clear();  
    fLines.reserve(fVertices.size()); 

    auto Draw_line = [this](const PolynomialCut::Vertex_t& vtx_now, const PolynomialCut::Vertex_t& vtx_next)
    {
        auto line = new TLine(
            vtx_now.x,  vtx_now.y,
            vtx_next.x, vtx_next.y
        ); 
        line->SetLineColor(fLineColor); 
        line->Draw(); 
        
        this->fLines.push_back(line); 
    };

    if (fVertices.size() >= 2) {
        for (size_t i=0; i<fVertices.size()-1; i++) Draw_line( fVertices[i], fVertices[i+1] ); 

        //now, connect the last -> first vertices
        Draw_line( fVertices.back(), fVertices.front() ); 
    }

    canvas->Update(); 
    canvas->Modified(); 
}; 
//_____________________________________________________________________________________________________________________________________
void PolynomialCutApp::HandleCanvasClick() {
    // Get canvas and event information
    TCanvas* canvas = fEcanvas->GetCanvas();
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
            
            // Update labels with coordinates
            fXLabel->SetText(Form("X: %.3f", x));
            fYLabel->SetText(Form("Y: %.3f", y));

            fVertices.push_back({ .x=x, .y=y }); 

            RedrawLines(); 
        }
    }
}

ClassImp(PolynomialCutApp); 