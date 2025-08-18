////////////////////////////////////////////////////////////////////////////////
//
//  PolynomialCut
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
#include "TLine.h"
#include "TH2.h"
#include "PolynomialCut.h"


class PolynomialCutApp : public TGMainFrame {
private:
    TRootEmbeddedCanvas* fEcanvas;
    int fEventType; 
    
    //buttons
    TGTextButton* fExitButton;  //Exit
    TGTextButton* fWriteButton; //Write output
    TGTextButton* fClearButton; //Clear
    TGTextButton* fUndoButton;  //Undo 

    //text entry
    TGTextEntry*  fFilePathEntry;   //window to enter file path
    
    TGLabel* fXLabel;
    TGLabel* fYLabel;
    TH2D* fHistogram;

    std::vector<TLine*> fLines;
    unsigned int fLineColor; 
    std::vector<PolynomialCut::Vertex_t> fVertices;  
    
public:
    PolynomialCutApp(const TGWindow* p, UInt_t w, UInt_t h, TH2* hist, const char* drawing_option="col2", unsigned int palette=102);
    virtual ~PolynomialCutApp();
    
    void CloseWindow();     //close the window
    void DoExit();          //exit button (exits)
    void WriteOutput();     //Write cut to output file
    void DeleteLast();      //delete last vertex
    void DeleteAll();       //delete all vertices
    void HandleCanvasClick();   //handle the canvas being clicked

    void RedrawLines();     //force redraw of all lines

    enum ECanvasEventType { kMouseButton1_down=1, kMouseButton1_up=11, kEnterObj=52, kLeaveObj=53 };  
    
    ClassDef(PolynomialCutApp, 0)
};
