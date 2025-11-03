#ifndef SaveOutputFrame_h_
#define SaveOutputFrame_h_

//////////////////////////////////////////////////////////////
//
//  A small window which appears when the output is saved. 
//
//////////////////////////////////////////////////////////////

//ROOT GUI headers
#include <TGFrame.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h> 
#include <TGWindow.h>
#include <TGLabel.h> 
//misc. ROOT headers
#include <TVector3.h>
#include <ROOT/RDataFrame.hxx>
//std headers
#include <vector>
#include <string> 
//APEX headers
#include "SieveHoleData.h"


class SaveOutputFrame : public TGMainFrame {
private: 

    //vertical frame to contain the path-entry frame, event-entry frame, and button frames 
    //TGVerticalFrame     *fVFrame_all; 

    //Horizontal frame for the text-entry field to put the output file path in 
    TGHorizontalFrame   *fHFrame_outPath;  
    TGTextEntry         *fTextEntry_outPath; 
    TGLabel             *fLabel_outPath; 

    //Horizontal frame to enter the number of events in the output file. 
    TGHorizontalFrame   *fHFrame_nEvents; 
    TGNumberEntry       *fNumEntry_nEvents; 
    TGLabel             *fLabel_nEvents; 

    const unsigned long int fnEvents_default = 1e5; 
    const unsigned long int fnEvents_max     = 1e8; 
    //if the user tries to make an output file with less than this many events per hole, 
    // they will be warned that the results may not be reasonable. 
    const double fnEvents_min_events_per_hole_warning = 50.; 

    //'Save' and 'Exit' buttons to save the output
    TGHorizontalFrame   *fHFrame_buttons; 
    TGTextButton        *fButton_Save; 
    TGTextButton        *fButton_Exit; 

    std::vector<SieveHoleData> fSavedHoles; 

    ROOT::RDF::RNode *fNode{nullptr}; 

    std::string fBranch_vert, fBranch_horiz; 

    const double fFpcoord_cut_width; 

    TVector3 fReactVertex; 
    bool fIsRHRS; 

public: 
    SaveOutputFrame(const TGWindow *p, 
                    UInt_t w, 
                    UInt_t h, 
                    const std::vector<SieveHoleData>& shd,
                    ROOT::RDF::RNode *df, 
                    bool is_RHRS,
                    TVector3 vtx, 
                    const double fpcoord_cut_width = 0.0055,
                    const char* branch_horiztonal="dxdz_sv",
                    const char* branch_vertical="dydz_sv"     
                    ); 

    ~SaveOutputFrame() { Cleanup(); } 

    void DoSave(); 
    void DoExit(); 

    // a new number of events is specified. 
    void SetNEvents() {}; 

    ClassDef(SaveOutputFrame, 1); 
}; 

#endif