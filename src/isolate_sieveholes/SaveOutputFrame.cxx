#include "SaveOutputFrame.h"
#include <ROOT/RDataFrame.hxx>
#include <TParameter.h>

using namespace std; 

//_______________________________________________________________________________________________________________
SaveOutputFrame::SaveOutputFrame(   const TGWindow *p, 
                                    UInt_t w, 
                                    UInt_t h, 
                                    const std::vector<SieveHoleData>& shd,
                                    bool is_RHRS,
                                    TVector3 vtx )
    : TGMainFrame(p, w, h)
{
    if (shd.empty()) CloseWindow(); 

    

    fIsRHRS = is_RHRS; 
    fReactVertex = vtx; 

    // Set up the main frame
    SetCleanup(kDeepCleanup);

    auto bframe = new TGHorizontalFrame(this, 900, 500); 
    
    fButton_Save = new TGTextButton(bframe, "&Save", 1); 
    fButton_Save->Connect("Clicked()", "SaveOutputFrame", this, "DoSave()"); 
    bframe->AddFrame(fButton_Save, new TGLayoutHints(kLHintsLeft  | kLHintsExpandX , 10, 10, 10, 5)); 

    fButton_Exit = new TGTextButton(bframe, "&Exit", 1); 
    fButton_Exit->Connect("Clicked()", "SaveOutputFrame", this, "DoExit()"); 
    bframe->AddFrame(fButton_Exit, new TGLayoutHints(kLHintsLeft  | kLHintsExpandX , 10, 10, 10, 5)); 
    
    AddFrame(bframe, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 10, 10, 10, 10)); 

    fTextEntry = new TGTextEntry(this);
    AddFrame(fTextEntry, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 20, 20, 20, 20));  

    //make copies of all the sieve holes which were sucessfully evaluated
    for (const SieveHoleData& dat : shd) {
        if (dat.is_evaluated) fSavedHoles.push_back(dat); 
    }

    char buff[200]; sprintf(buff, "Save %zi sieve holes?", fSavedHoles.size()); 

    SetWindowName(buff);
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}
//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoSave() 
{ 
    if (!fTextEntry) return; 
    
    const char* path_outfile = fTextEntry->GetBuffer()->GetString(); 

    //get the react vertex from the parent
    const TVector3 rvtx = fReactVertex; 

    //now, we save the output file 
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 

    ROOT::RDataFrame df(fSavedHoles.size()); 

    int i_elem =0; 

    auto snapshot = df 

        .Define("hole_data", [this, &i_elem]()
        { 
            printf("element %4i/%zi", i_elem, this->fSavedHoles.size()); cout << endl; 
            return this->fSavedHoles.at(i_elem++); 
        }, {})

        .Define("x_sv",     [](const SieveHoleData& hd){ return hd.hole.x; }, {"hole_data"})
        .Define("y_sv",     [](const SieveHoleData& hd){ return hd.hole.y; }, {"hole_data"}) 
        .Define("dxdz_sv",  [rvtx](const SieveHoleData& hd){ return ( hd.hole.x - rvtx.x() ) / ( 0. - rvtx.z() ); }, {"hole_data"})
        .Define("dydz_sv",  [rvtx](const SieveHoleData& hd){ return ( hd.hole.y - rvtx.y() ) / ( 0. - rvtx.z() ); }, {"hole_data"})

        .Define("a_y_fp",      [](const SieveHoleData& hd){ return hd.y_fp.poly; },    {"hole_data"})        
        .Define("a_dxdz_fp",   [](const SieveHoleData& hd){ return hd.dxdz_fp.poly; }, {"hole_data"})
        .Define("a_dydz_fp",   [](const SieveHoleData& hd){ return hd.dydz_fp.poly; }, {"hole_data"})

        .Define("hole_row",   [](const SieveHoleData& hd){ return hd.hole.row; },  {"hole_data"})
        .Define("hole_col",   [](const SieveHoleData& hd){ return hd.hole.col; },  {"hole_data"})

        .Define("x_fp_min",   [](const SieveHoleData& hd){ return hd.x_fp_min; },  {"hole_data"})
        .Define("x_fp_max",   [](const SieveHoleData& hd){ return hd.x_fp_max; },  {"hole_data"})

        .Define("hole_cut_x",   [](const SieveHoleData& hd){ return hd.cut_x; },  {"hole_data"})
        .Define("hole_cut_y",   [](const SieveHoleData& hd){ return hd.cut_y; },  {"hole_data"})

        .Define("hole_cut_width",  [](const SieveHoleData& hd){ return hd.cut_width;  },  {"hole_data"})
        .Define("hole_cut_height", [](const SieveHoleData& hd){ return hd.cut_height; },  {"hole_data"})

        .Define("position_vtx_scs", [rvtx](){ return rvtx; }, {})

        .Snapshot("hole_data", path_outfile, {
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",

            "a_y_fp",
            "a_dxdz_fp",
            "a_dydz_fp",

            "x_fp_min",
            "x_fp_max",

            "hole_cut_x",
            "hole_cut_y",

            "hole_cut_width",
            "hole_cut_height",

            "hole_row",
            "hole_col",

            "position_vtx_scs"
        }); 

    //now, create the output parameters we want to make
    auto file = new TFile(path_outfile, "UPDATE"); 

    if (!file || file->IsZombie()) {
        throw logic_error("in <SaveOutputFrame::DoSave>: putput file file is null/zombie. check path."); 
        return; 
    }

    auto param_is_RHRS = new TParameter<bool>("is_RHRS", fIsRHRS);  
    param_is_RHRS->Write();

    file->Close();
    delete file; 

    cout << "saved file: " << path_outfile << endl; 

    DoExit(); 
}
//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoExit() { CloseWindow(); }
//_________________________________________________________________________________________________________________________________

ClassImp(SaveOutputFrame); 