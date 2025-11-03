#include "SaveOutputFrame.h"
#include <ROOT/RDataFrame.hxx>
#include <TParameter.h>
#include <string> 
#include <cmath> 
#include <optional>
#include <ApexOptics.h> 
#include <TRandom3.h> 
#include <TFile.h> 
#include <TParameter.h> 

using namespace std; 
using ApexOptics::Trajectory_t; 

//_______________________________________________________________________________________________________________
SaveOutputFrame::SaveOutputFrame(   const TGWindow *p, 
                                    UInt_t w, 
                                    UInt_t h, 
                                    const std::vector<SieveHoleData>& shd,
                                    ROOT::RDF::RNode *df, 
                                    bool is_RHRS,
                                    TVector3 vtx, 
                                    const double fpcoord_cut_width,
                                    const char* branch_horizontal, 
                                    const char* branch_vertical
                                )
    : TGMainFrame(p, w, h), 
    fFpcoord_cut_width{fpcoord_cut_width},
    fNode{df},
    fBranch_horiz{branch_horizontal}, 
    fBranch_vert{branch_vertical}
{
    if (shd.empty()) CloseWindow(); 

    fIsRHRS = is_RHRS; 
    fReactVertex = vtx; 

    //check if node is ok
    if (fNode == nullptr) {
        throw logic_error("in <SaveOutputFrame::DoSave>: passed a nullptr for the dataframe!"); 
        return; 
    }

    // Set up the main frame
    SetCleanup(kDeepCleanup);

    // ~~~~ Text entry for the user to input the output file path: 
    fHFrame_outPath = new TGHorizontalFrame(this, w, 60);
    //
    // label
    fLabel_outPath = new TGLabel(fHFrame_outPath, "Output file rel. path: "); 
    fHFrame_outPath->AddFrame(fLabel_outPath, new TGLayoutHints(kLHintsLeft)); 
    //
    // text entry 
    fTextEntry_outPath = new TGTextEntry(fHFrame_outPath); 
    fHFrame_outPath->AddFrame(fTextEntry_outPath, new TGLayoutHints(kLHintsRight | kLHintsExpandX)); 
    //
    // Add this frame to the parent frame
    AddFrame(fHFrame_outPath, new TGLayoutHints(kLHintsTop | kLHintsExpandX)); 


    // ~~~ Number entry for the user to input the desired number of input events
    fHFrame_nEvents = new TGHorizontalFrame(this, w, 60); 
    // 
    // label
    fLabel_nEvents = new TGLabel(fHFrame_nEvents, "Number of output events: "); 
    fHFrame_nEvents->AddFrame(fLabel_nEvents, new TGLayoutHints(kLHintsLeft)); 
    //
    // number entry
    fNumEntry_nEvents = new TGNumberEntry(fHFrame_nEvents, 1e5); // 1e5 is the default number of events generated. 
    fNumEntry_nEvents->SetNumStyle(TGNumberEntry::kNESInteger); //number entered must be an integer
    fNumEntry_nEvents->SetNumAttr(TGNumberEntry::kNEAPositive); //number entered must be positive  
    fNumEntry_nEvents->SetLimits(TGNumberEntry::kNELLimitMinMax, 1., (double)fnEvents_max); //set the min and max reasonable limits
    fNumEntry_nEvents->GetNumberEntry()->Connect("ReturnPressed()", "SaveOutputFrame", this, "SetNEvents()"); 
    fHFrame_nEvents->AddFrame(fNumEntry_nEvents, new TGLayoutHints(kLHintsRight | kLHintsExpandX)); 
    //
    // Add this frame to the parent frame
    AddFrame(fHFrame_nEvents, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX)); 


    // ~~~ Add frame for the 'save' and 'Exit' buttons
    fHFrame_buttons = new TGHorizontalFrame(this, w, 60); 
    //
    // 'save' button
    fButton_Save = new TGTextButton(fHFrame_buttons, "&Save", 1); 
    fButton_Save->Connect("Clicked()", "SaveOutputFrame", this, "DoSave()"); 
    fHFrame_buttons->AddFrame(fButton_Save, new TGLayoutHints(kLHintsLeft  | kLHintsExpandX)); 
    //
    // 'exit' button
    fButton_Exit = new TGTextButton(fHFrame_buttons, "&Exit", 1); 
    fButton_Exit->Connect("Clicked()", "SaveOutputFrame", this, "DoExit()"); 
    fHFrame_buttons->AddFrame(fButton_Exit, new TGLayoutHints(kLHintsLeft  | kLHintsExpandX)); 
    //
    // Add this frame to the parent frame
    AddFrame(fHFrame_buttons, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX)); 

    //Text entry for the user to input the desired number of simulated events:     

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
//simple helper function to define tparameters to TFiles 
template<typename T> void Add_TParameter_to_TFile(const char* name, T val)
{
    auto param = new TParameter<T>(name, val); 
    param->Write(); 
}

//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoSave() 
{ 
    const char* const here = "SaveOutputFrame::DoSave"; 

    if (!fTextEntry_outPath) {
        throw logic_error("in <SaveOutputFrame::DoSave>: fTextEntry_outPath ptr is null!"); 
        return; 
    }
    if (!fNumEntry_nEvents)  {
        throw logic_error("in <SaveOutputFrame::DoSave>: fNumberEntry_nEvents ptr is null!"); 
        return; 
    }
    
    const char* path_outfile = fTextEntry_outPath->GetBuffer()->GetString(); 

    //check to see if the '.root' file extension is present 
    if (string(path_outfile).find(".root") == string::npos) {
        fprintf(stderr, "\nWarning - the desired file path '%s' does NOT have a '.root' extension... Save anyway? [y/N]: ", path_outfile);  
        string ans; cin >> ans; 
        //if they answer anything but 'yes', then abord this save. 
        if (ans.find("y")==string::npos && ans.find("Y")==string::npos) {
            cout << ans << " - Ok, we won't save this time." << endl; 
            return; 
        } else {
            cout << ans << " - Ok, saving anyway..." << endl; 
        }
    }

    //get the react vertex from the parent
    const TVector3 rvtx = fReactVertex; 

    //now, we save the output file 
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 

    //ROOT::RDataFrame df(fSavedHoles.size()); 

    int i_elem =0; 

    const unsigned long int n_events = (unsigned long)fNumEntry_nEvents->GetNumber(); 

    cout << "\nCreating outfile of hole-fits: '" << path_outfile << "' with " << n_events <<  " events..." << flush; 
    
    ROOT::RDataFrame df(n_events); 

    const size_t n_holes = fSavedHoles.size(); 
    TRandom3 rand(0); 

    df 
        .Define("hole_data", [n_holes, this, &rand]()
        {   
            //pick a random sieve hole
            return this->fSavedHoles[ rand.Integer(n_holes) ];
        }, {})

        .Define("hole_save_data", [&rand](const SieveHoleData& hd)
        {   
            //pick a random raster partition
            return hd.hole_save_data[ rand.Integer(hd.hole_save_data.size()) ]; 
        }, {"hole_data"})

        .Define("Xfp", [&rand](const SieveHoleData& hd, const HoleSaveData& hsd)
        {
            const double x_fp = rand.Uniform( hd.x_fp_min, hd.x_fp_max ); 

            return Trajectory_t{
                x_fp, 
                hsd.y_fp.Eval(x_fp),
                hsd.dxdz_fp.Eval(x_fp),
                hsd.dydz_fp.Eval(x_fp)
            }; 
        }, {"hole_data", "hole_save_data"})

        .Define("Xsv", [](const HoleSaveData& hsd)
        {
            return hsd.Xsv; 
        }, {"hole_save_data"})

        .Define("position_vtx_scs", [](const HoleSaveData& hsd)
        {
            return hsd.position_vtx_scs; 
        }, {"hole_save_data"})

        .Define("x_fp",     [](Trajectory_t X){ return X.x; },      {"Xfp"})
        .Define("y_fp",     [](Trajectory_t X){ return X.y; },      {"Xfp"})
        .Define("dxdz_fp",  [](Trajectory_t X){ return X.dxdz; },   {"Xfp"})
        .Define("dydz_fp",  [](Trajectory_t X){ return X.dydz; },   {"Xfp"})
        
        .Define("x_sv",     [](Trajectory_t X){ return X.x; },      {"Xsv"})
        .Define("y_sv",     [](Trajectory_t X){ return X.y; },      {"Xsv"})
        .Define("dxdz_sv",  [](Trajectory_t X){ return X.dxdz; },   {"Xsv"})
        .Define("dydz_sv",  [](Trajectory_t X){ return X.dydz; },   {"Xsv"})
        .Define("dpp_sv",   [](Trajectory_t X){ return X.dpp; },    {"Xsv"})

        .Define("position_vtx", [this](TVector3 vtx_scs)
        {
            return ApexOptics::SCS_to_HCS(this->fIsRHRS, vtx_scs);
        }, {"position_vtx_scs"})

        .Snapshot("tracks_fp", path_outfile, {
            "x_fp",
            "y_fp",
            "dxdz_fp",
            "dydz_fp",

            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",

            "position_vtx",
            "position_vtx_scs"
        }); 

    auto file = new TFile(path_outfile, "UPDATE"); 

    Add_TParameter_to_TFile<bool>("is_RHRS", fIsRHRS); 

    file->Close(); 
    delete file; 

    cout << "done." << endl; 

#if 0 
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

    cout << "done" << endl; 

    cout << "saved (hole) file: " << path_outfile << endl; 

    const double focalplane_cut_width = 0.0055; 

    const char* path_outfile_cuts = Form("%s-cutdata.root", prefix_outfile); 

    cout << "Creating outfile of cuts: '" << path_outfile_cuts << "..." << flush; 

    auto df_cuts = (*fNode)

        //now save all the 'cut' data (actual events)
        .Define("holedata_selected_opt", [this](  double dxdz_sv, //sieve coords
                                                                        double dydz_sv, 
                                                                        double x_fp, //target coords
                                                                        double y_fp, 
                                                                        double dxdz_fp, 
                                                                        double dydz_fp, 
                                                                        TVector3 vtx_scs    
                                                                    )
        {   
            std::optional<SieveHoleData> holedata_selected(std::nullopt); 

            for (const auto & holedata : this->fSavedHoles) {

                double dxdz_hole = (holedata.hole.x - vtx_scs.x())/(0. - vtx_scs.z()); 
                double dydz_hole = (holedata.hole.y - vtx_scs.y())/(0. - vtx_scs.z()); 

                //make cut on hole coords in fp 
                double hole_rad2 = 
                    pow( (dxdz_hole - dxdz_sv)/holedata.cut_width,  2) + 
                    pow( (dydz_hole - dydz_sv)/holedata.cut_height, 2);
                    
                if (hole_rad2 > 1.) continue; 
                
                //now, make a cut on event in the focal plane 
                //printf(" y_fp (polynomial vs real): %+.4f  %+.4f\n", holedata.y_fp.Eval(x_fp), y_fp); 
                const double& fp_cut = this->fFpcoord_cut_width; 

                if (fabs( holedata.y_fp   .Eval(x_fp) - y_fp )    > fp_cut) continue; 
                if (fabs( holedata.dxdz_fp.Eval(x_fp) - dxdz_fp ) > fp_cut) continue; 
                if (fabs( holedata.dydz_fp.Eval(x_fp) - dydz_fp ) > fp_cut) continue; 
                
                //make a copy for us to keep
                holedata_selected = optional<SieveHoleData>(SieveHoleData{holedata});
            }

            return holedata_selected; 

        }, {fBranch_horiz.c_str(), 
            fBranch_vert.c_str(), 
            "x_fp", "y_fp", "dxdz_fp", "dydz_fp", 
            "position_vtx_scs"})
        
        //filter out events which were not selected for a particular hole
        .Filter([](std::optional<SieveHoleData> holedata_opt){ return holedata_opt.has_value(); }, {"holedata_selected_opt"})

        .Define("Xsv", [](std::optional<SieveHoleData> holedata_opt, TVector3 vtx_scs)
        {   
            SieveHole hole = holedata_opt.value().hole; 

            //take the 'average' bewtween the front and back centers of the sieve hole 
            double dxdz_front = (hole.x - vtx_scs.x())/(0. - vtx_scs.z()); 
            double dydz_front = (hole.y - vtx_scs.y())/(0. - vtx_scs.z());

            double dxdz_back  = (hole.x - vtx_scs.x())/(0.0125 - vtx_scs.z()); 
            double dydz_back  = (hole.y - vtx_scs.y())/(0.0125 - vtx_scs.z());

            return Trajectory_t{ 
                .x      = hole.x, 
                .y      = hole.y, 
                .dxdz   = 0.5*(dxdz_front + dxdz_back), 
                .dydz   = 0.5*(dydz_front + dydz_back)
            }; 
        }, {"holedata_selected_opt", "position_vtx_scs"})

        .Redefine("x_sv",     [](Trajectory_t Xsv){ return Xsv.x;     }, {"Xsv"})
        .Redefine("y_sv",     [](Trajectory_t Xsv){ return Xsv.y;     }, {"Xsv"})
        .Redefine("dxdz_sv",  [](Trajectory_t Xsv){ return Xsv.dxdz;  }, {"Xsv"})
        .Redefine("dydz_sv",  [](Trajectory_t Xsv){ return Xsv.dydz;  }, {"Xsv"})
        
        .Snapshot("tracks_fp_cut", path_outfile_cuts, {
            "x_sv", 
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",

            "x_fp",
            "y_fp",
            "dxdz_fp",
            "dydz_fp",

            "position_vtx_scs"
        }); 

    cout << "done." << endl; 

    //now, create the output parameters we want to make
    file = new TFile(path_outfile_cuts, "UPDATE"); 
    if (!file || file->IsZombie()) {
        throw logic_error("in <SaveOutputFrame::DoSave>: putput file file is null/zombie. check path."); 
        return; 
    }
    param_is_RHRS = new TParameter<bool>("is_RHRS", fIsRHRS);  
    param_is_RHRS->Write();
    file->Close();
    delete file; 
#endif 

    DoExit(); 
}
//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoExit() { CloseWindow(); }
//_________________________________________________________________________________________________________________________________

ClassImp(SaveOutputFrame); 