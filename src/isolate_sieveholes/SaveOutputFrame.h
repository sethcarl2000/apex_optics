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

class SaveOutputFrame : public TGMainFrame {
private: 
    TGTextButton *fButton_Save; 
    TGTextButton *fButton_Exit; 

    TGTextEntry *fTextEntry; 

    std::vector<SieveHoleData> fSavedHoles; 

    TVector3 fReactVertex; 
    bool fIsRHRS; 

public: 
    SaveOutputFrame(const TGWindow *p, 
                    UInt_t w, 
                    UInt_t h, 
                    const std::vector<SieveHoleData>& shd,
                    bool is_RHRS,
                    TVector3 vtx ); 

    ~SaveOutputFrame() { Cleanup(); } 

    void DoSave(); 
    void DoExit(); 

    ClassDef(SaveOutputFrame, 1); 
}; 

#endif