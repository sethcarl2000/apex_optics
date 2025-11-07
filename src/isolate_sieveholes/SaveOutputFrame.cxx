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
#include <stdexcept> 

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
                                    double rast_min, double rast_max, 
                                    const double fpcoord_cut_width,
                                    const char* branch_horizontal, 
                                    const char* branch_vertical
                                )
    : TGMainFrame(p, w, h), 
    fFpcoord_cut_width{fpcoord_cut_width},
    fNode{df},
    fBranch_horiz{branch_horizontal}, 
    fBranch_vert{branch_vertical},
    fRastMin{rast_min}, fRastMax{rast_max},
    fIsRHRS{is_RHRS},
    fReactVertex{vtx}
{
    if (shd.empty()) CloseWindow(); 

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
    const char* const here = "DoSave"; 

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

    //parse a polynomial to estimate dpp_sv 
    NPoly poly_dpp(4); 
    try {
        ApexOptics::Parse_NPoly_from_file(fPath_to_dppPoly.c_str(), "dpp_sv", &poly_dpp); 
    } catch (const std::exception& e) {
        throw logic_error(Form(
            "in <SaveOutputFrame::DoSave>: something went wrong trying to parse polynomial from file "
            "'%s'.\n"
            " what(): %s"
            ,
            fPath_to_dppPoly.c_str(),
            e.what()
        )); 
        return; 
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

    try {
        df 
            .Define("hole_data", [n_holes, this, &rand]()
            {   
                //pick a random sieve hole
                return this->fSavedHoles[ rand.Integer(n_holes) ];
            }, {})

            .Define("rast_param", [&rand]()
            {   
                //pick a random raster position
                return rand.Uniform(-1., +1.);  
            }, {})

            .Define("Xfp", [&rand](const SieveHoleData& hd, double rast_param)
            {
                const double x_fp = rand.Uniform( hd.x_fp_min, hd.x_fp_max ); 

                return Trajectory_t{
                    x_fp, 
                    hd.y_fp   .Eval({x_fp, rast_param}),
                    hd.dxdz_fp.Eval({x_fp, rast_param}),
                    hd.dydz_fp.Eval({x_fp, rast_param})
                }; 
            }, {"hole_data", "rast_param"})

            .Define("position_vtx", [this](double rast_param)
            {
                //so we used the average react-vertex, but let's put accurate y-raster information back into it. 
                TVector3 position_vtx = this->fReactVertex; 

                double y_hcs = 0.5*(this->fRastMax + this->fRastMin); //middle of raster span
                
                y_hcs += 0.5*(this->fRastMax - this->fRastMin) * rast_param;  

                return TVector3(
                    position_vtx.x(),
                    y_hcs,
                    position_vtx.z()
                ); 

            }, {"rast_param"})


            .Define("position_vtx_scs", [this](const TVector3& vtx_hcs)
            {
                return ApexOptics::HCS_to_SCS(this->fIsRHRS, vtx_hcs); 
            }, {"position_vtx"})

            .Define("Xsv", [this, &poly_dpp](const TVector3& vtx_scs, Trajectory_t Xfp, SieveHoleData& hd)
            {   
                //placeholder estimate of dp/p_sv
                double dpp_est = poly_dpp.Eval({Xfp.x, Xfp.y, Xfp.dxdz, Xfp.dydz}); 

                const SieveHole& hole = hd.hole; 

                double dxdz = (hole.x - vtx_scs.x()) / (0. - vtx_scs.z()); 
                double dydz = (hole.y - vtx_scs.y()) / (0. - vtx_scs.z()); 

                return Trajectory_t{
                    hd.hole.x, 
                    hd.hole.y, 
                    dxdz, 
                    dydz, 
                    dpp_est
                }; 
            }, {"position_vtx_scs", "Xfp", "hole_data"})

            .Define("x_fp",     [](Trajectory_t X){ return X.x; },      {"Xfp"})
            .Define("y_fp",     [](Trajectory_t X){ return X.y; },      {"Xfp"})
            .Define("dxdz_fp",  [](Trajectory_t X){ return X.dxdz; },   {"Xfp"})
            .Define("dydz_fp",  [](Trajectory_t X){ return X.dydz; },   {"Xfp"})
            
            .Define("x_sv",     [](Trajectory_t X){ return X.x; },      {"Xsv"})
            .Define("y_sv",     [](Trajectory_t X){ return X.y; },      {"Xsv"})
            .Define("dxdz_sv",  [](Trajectory_t X){ return X.dxdz; },   {"Xsv"})
            .Define("dydz_sv",  [](Trajectory_t X){ return X.dydz; },   {"Xsv"})
            .Define("dpp_sv",   [](Trajectory_t X){ return X.dpp; },    {"Xsv"})

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
    
    } catch (const std::exception& e) {
        cerr << endl; 
        Error(here, "Something went wrong trying to create the output file.\n what(): %s", e.what()); 
        DoExit(); 
        return;
    } 

    cout << "done." << endl; 
    DoExit(); 
}
//_________________________________________________________________________________________________________________________________
void SaveOutputFrame::DoExit() { CloseWindow(); }
//_________________________________________________________________________________________________________________________________

ClassImp(SaveOutputFrame); 