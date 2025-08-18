#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TApplication.h>
#include <ROOT/RDataFrame.hxx>
#include "TRandom3.h"
#include <iostream>
#include <vector>
#include <string> 
#include "TLine.h"
#include "TH2.h"

using namespace std; 

class PolynomialCutApp : public TGMainFrame {
private:
    TRootEmbeddedCanvas* fEcanvas;
    
    //buttons
    TGTextButton* fExitButton;  //Exit
    TGTextButton* fWriteButton; //Write output
    TGTextButton* fClearButton; //Clear
    TGTextButton* fUndoButton;  //Undo 
    
    TGLabel* fXLabel;
    TGLabel* fYLabel;
    TH2D* fHistogram;

    vector<TLine*> fLines;
    vector<PolynomialCut::Vertex_t> fVertices;  
    
public:
    PolynomialCutApp(const TGWindow* p, UInt_t w, UInt_t h, TH2* hist);
    virtual ~PolynomialCutApp();
    
    void CloseWindow();     //close the window
    void DoExit();          //exit button (exits)
    void WriteOutput();     //
    void DeleteLast(); 
    void DeleteAll(); 
    void HandleCanvasClick();
    
    ClassDef(PolynomialCutApp, 0)
};

PolynomialCutApp::PolynomialCutApp(const TGWindow* p, UInt_t w, UInt_t h, TH2 *hist) 
    : TGMainFrame(p, w, h) {
    
    // Set up the main frame
    SetCleanup(kDeepCleanup);
    
    fHistogram = (TH2D*)hist->Clone("histogram_clone");
    fHistogram->SetDirectory(0); // Detach from file directory
    
    // Create embedded canvas
    fEcanvas = new TRootEmbeddedCanvas("ECanvas", this, 600, 400);
    AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 5));
    
    // Create horizontal frame for labels and button
    TGHorizontalFrame* hframe = new TGHorizontalFrame(this, 600, 50);
    
    // Create coordinate labels
    fXLabel = new TGLabel(hframe, "X: ---");
    fXLabel->SetTextJustify(kTextLeft);
    hframe->AddFrame(fXLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 10, 20, 5, 5));
    
    fYLabel = new TGLabel(hframe, "Y: ---");
    fYLabel->SetTextJustify(kTextLeft);
    hframe->AddFrame(fYLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX, 20, 20, 5, 5));
    

    auto Add_button = [hframe, this](TGTextButton* button, string button_label, string method, TGLayoutHints* hints) 
    {
        button_label = "&" + button_label; 
        button = new TGTextButton(hframe, button_label.c_str(), 1); 
        button->Connect("Clicked()", "PolynomialCutApp", this, method.c_str()); 
        hframe->AddFrame(button, hints); 
    };

    Add_button(fExitButton,  "Exit",            "DoExit()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    ); 
    Add_button(fWriteButton, "Write Output",    "WriteOutput()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    ); 
    Add_button(fUndoButton,  "Undo",            "DeleteLast()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    );
    Add_button(fClearButton, "Clear",           "DeleteAll()", 
        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)
    );

    
    AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 10));
    
    // Draw histogram on canvas
    TCanvas* canvas = fEcanvas->GetCanvas();
    canvas->cd();
    fHistogram->Draw("COLZ");
    canvas->Update();
    
    // Connect canvas click event
    canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "PolynomialCutApp", this, 
                   "HandleCanvasClick()");
    
    // Set window properties
    SetWindowName("RDataFrame Analysis");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}

PolynomialCutApp::~PolynomialCutApp() {
    // Cleanup
    if (fHistogram) delete fHistogram;
    Cleanup();
}

void PolynomialCutApp::CloseWindow() {
    // Called when window is closed via window manager
    gApplication->Terminate(0);
}

void PolynomialCutApp::DoExit() {
    // Called when exit button is pressed
    gApplication->Terminate(0);
}

void PolynomialCutApp::WriteOutput() {};
void PolynomialCutApp::DeleteLast() {}; 
void PolynomialCutApp::DeleteAll() {}; 


void PolynomialCutApp::HandleCanvasClick() {
    // Get canvas and event information
    TCanvas* canvas = fEcanvas->GetCanvas();
    if (!canvas) return;
    
    // Get pixel coordinates
    Int_t px = canvas->GetEventX();
    Int_t py = canvas->GetEventY();
    
    // Convert pixel coordinates to histogram coordinates
    Double_t x = canvas->AbsPixeltoX(px);
    Double_t y = canvas->AbsPixeltoY(py);
    
    // Update labels with coordinates
    fXLabel->SetText(Form("X: %.3f", x));
    fYLabel->SetText(Form("Y: %.3f", y));
}

ClassImp(PolynomialCutApp); 

// Example usage function
void create_polynomial_cut() {
    // Create a sample RDataFrame for demonstration
    // In practice, you would pass your actual RDataFrame
    ROOT::RDataFrame df(1e5);

    TRandom3 rand; 

    auto hist = (TH2*)df
        .Define("x", [&rand](){ return rand.Gaus(); }, {})
        .Define("y", [&rand](){ return rand.Gaus(); }, {})
        .Define("z", [](double x, double y){ return x*x - 0.5 * y*y; }, {"x", "y"})
        .Histo2D<double>({"h_test_cpy", "Test histogram", 200, -5., 5., 200, -5., 5.}, "x", "y")->Clone("h_test"); 

    // Create and run the GUI
    new PolynomialCutApp(gClient->GetRoot(), 800, 600, hist);
}
