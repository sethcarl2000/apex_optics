#ifndef SaveOutputFrame_h_
#define SaveOutputFrame_h_

//////////////////////////////////////////////////////////////
//
//  A small window which appears when the output is saved. 
//
//////////////////////////////////////////////////////////////

#include <TGFrame.h>
#include <TGButton.h>
#include <vector>
#include <TGTextEntry.h>
#include <TGWindow.h>
#include "SieveHoleData.h"
#include <TVector3.h>
#include <ROOT/RDataFrame.hxx>
#include <string> 

class SaveOutputFrame : public TGMainFrame {
private: 
    TGTextButton *fButton_Save; 
    TGTextButton *fButton_Exit; 

    TGTextEntry *fTextEntry; 

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

    ClassDef(SaveOutputFrame, 1); 
}; 

#endif