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

class AnalysisFrame : public TGMainFrame {
private:
    TRootEmbeddedCanvas* fEcanvas;
    TGTextButton* fExitButton;
    TGLabel* fXLabel;
    TGLabel* fYLabel;
    TH2D* fHistogram;
    
public:
    AnalysisFrame(const TGWindow* p, UInt_t w, UInt_t h, ROOT::RDF::RNode df);
    virtual ~AnalysisFrame();
    
    void CloseWindow();
    void DoExit();
    void HandleCanvasClick();
    
    ClassDef(AnalysisFrame, 0)
};

AnalysisFrame::AnalysisFrame(const TGWindow* p, UInt_t w, UInt_t h, ROOT::RDF::RNode df) 
    : TGMainFrame(p, w, h) {
    
    // Set up the main frame
    SetCleanup(kDeepCleanup);
    
    // Perform analysis on the RDataFrame
    // Example analysis: create a 2D histogram from two columns
    // Assuming the DataFrame has columns "x" and "y" - adjust as needed
    auto hist = df.Histo2D({"analysis_hist", "Analysis Result;X;Y", 100, -5, 5, 100, -5, 5}, "x", "y");
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
    
    // Create exit button
    fExitButton = new TGTextButton(hframe, "&Exit", 1);
    fExitButton->Connect("Clicked()", "AnalysisFrame", this, "DoExit()");
    hframe->AddFrame(fExitButton, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5));
    
    AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 10));
    
    // Draw histogram on canvas
    TCanvas* canvas = fEcanvas->GetCanvas();
    canvas->cd();
    fHistogram->Draw("COLZ");
    canvas->Update();
    
    // Connect canvas click event
    canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "AnalysisFrame", this, 
                   "HandleCanvasClick()");
    
    // Set window properties
    SetWindowName("RDataFrame Analysis");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}

AnalysisFrame::~AnalysisFrame() {
    // Cleanup
    if (fHistogram) delete fHistogram;
    Cleanup();
}

void AnalysisFrame::CloseWindow() {
    // Called when window is closed via window manager
    gApplication->Terminate(0);
}

void AnalysisFrame::DoExit() {
    // Called when exit button is pressed
    gApplication->Terminate(0);
}

void AnalysisFrame::HandleCanvasClick() {
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

ClassImp(AnalysisFrame); 

// Example usage function
void test_application() {
    // Create a sample RDataFrame for demonstration
    // In practice, you would pass your actual RDataFrame
    ROOT::RDataFrame df(1e5);

    TRandom3 rand; 

    auto df_with_data = df
        .Define("x", [&rand](){ return rand.Gaus(); }, {})
        .Define("x", [&rand](){ return rand.Gaus(); }, {})
        .Define("z", [](double x, double y){ return x*x - 0.5 * y*y; }, {"x", "y"}); 

    // Create and run the GUI
    new AnalysisFrame(gClient->GetRoot(), 800, 600, df_with_data);
}

// Main function for standalone compilation
#ifdef STANDALONE
int main(int argc, char** argv) {
    TApplication app("AnalysisApp", &argc, argv);
    
    RunAnalysisGUI();
    
    app.Run();
    return 0;
}
#endif