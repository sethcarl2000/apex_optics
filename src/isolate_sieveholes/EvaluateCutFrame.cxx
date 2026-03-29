// APEX headers
#include <isolate_sieveholes/EvaluateCutFrame.h>
#include <RMatrix.h>
#include <NPoly.h> 
#include <ApexOptics.h> 
// ROOT headers
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TPad.h> 
#include <TVirtualPad.h> 
#include <TCanvas.h>
#include <TObject.h> 
#include <TBox.h>
#include <TGraphErrors.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h> 
#include <Math/Factory.h> 
#include <ROOT/RVec.hxx>
// stdlib headers
#include <cmath> 
#include <thread> 

using namespace std; 
using ApexOptics::Trajectory_t; 

namespace {
    constexpr double pi = 3.14159265359; 

}

EvaluateCutFrame::EvaluateCutFrame( const TGWindow *p, 
                                    const vector<EventData>& data, 
                                    SieveHoleData *_hd,
                                    const double fp_cut_width, 
                                    const char* branch_x, 
                                    const char* branch_y, 
                                    const char* draw_option,
                                    const unsigned int palette )
    : TGMainFrame( p, 1400, 800 ), 
    fRDF{nullptr}, 
    fData{&data},
    fSelectedSieveHole{_hd}, 
    fFpcoord_cut_width{fp_cut_width}
{   
    const char* const here = "EvaluateCutFrame::constructor"; 

    fParent = PickSieveHoleApp::Instance(); 
    if (!fParent) {
        Error(__func__, "static ptr to singleton instance of 'PickSieveHoleApp' is null... something has gone terribly wrong here."); 
        return; 
    }

#ifdef DEBUG
    Info(here, "Entering constructor body"); 
#endif

    // in this 'app state', we are picking limits to evaluate the fits within  
    fAppState = kNone; 

    const int polynomial_degree =3; 

    //setup the window
    SetCleanup(kDeepCleanup); 

    fFrame_canv = new TGHorizontalFrame(this, 1400, 800); 

    fEcanvas = new TRootEmbeddedCanvas("ECanvas_eval", fFrame_canv, 1400, 800); 

    const double cut_x      = fSelectedSieveHole->cut_x; 
    const double cut_y      = fSelectedSieveHole->cut_y; 
    const double cut_width  = fSelectedSieveHole->cut_width; 
    const double cut_height = fSelectedSieveHole->cut_height; 
    const double cut_angle  = fSelectedSieveHole->cut_angle; 

    const double cos_theta  = cos( fSelectedSieveHole->cut_angle * pi / 180. );
    const double sin_theta  = sin( fSelectedSieveHole->cut_angle * pi / 180. );  

    fHist_holes = new TH2D("h_xy_temp", "", 
        75, fParent->GetDrawRange_x_min(), fParent->GetDrawRange_x_max(),
        75, fParent->GetDrawRange_y_min(), fParent->GetDrawRange_y_max() 
    ); 
    fY_fp    = HistAndLimit(new TH2D("h_y_fp",    "y_fp vs x_fp",     75, fParent->GetDrawRange_xfp_min(), fParent->GetDrawRange_xfp_max(), 75, -0.070, +0.055)); 
    fDxdz_fp = HistAndLimit(new TH2D("h_dxdz_fp", "dx/dz_fp vs x_fp", 75, fParent->GetDrawRange_xfp_min(), fParent->GetDrawRange_xfp_max(), 75, -0.035, +0.025)); 
    fDydz_fp = HistAndLimit(new TH2D("h_dydz_fp", "dy/dz_fp vs x_fp", 75, fParent->GetDrawRange_xfp_min(), fParent->GetDrawRange_xfp_max(), 75, -0.060, +0.040)); 
    

    /// @return 'true' if a given event is inside the sieve-coordinate cut, and 'false' if not. 
    fCutFcn = [cos_theta, sin_theta, cut_width, cut_height, cut_x, cut_y](const EventData& ev)
    {
        //now, only include events which are inside our 'cut ellipse' 
        double x = ev.Xsv.dxdz - cut_x; 
        double y = ev.Xsv.dydz - cut_y; 
        
        double u = ( cos_theta*x + sin_theta*y )/cut_width; 
        double v = ( cos_theta*y - sin_theta*x )/cut_height; 

        if ( u*u + v*v < 1. ) { return true; } else { return false; }
    };
    //_________________________________________________________________________________________

    for (const auto& ev : *fData) {

        fHist_holes->Fill( ev.Xsv.dxdz, ev.Xsv.dydz );

        if (fCutFcn(ev)) {
            fY_fp.hist   ->Fill( ev.Xfp.x, ev.Xfp.y );
            fDxdz_fp.hist->Fill( ev.Xfp.x, ev.Xfp.dxdz ); 
            fDydz_fp.hist->Fill( ev.Xfp.x, ev.Xfp.dydz ); 
        }
    }

#ifdef DEBUG
    printf("Number of events in cut: %zi\n", fData->size()); 
#endif 

    //create the histograms
    TCanvas *canv = fEcanvas->GetCanvas(); 

    gStyle->SetPalette(palette); 
    
    canv->cd(); 
    canv->Divide(2,2); 

    canv->cd(1); fHist_holes->Draw("col2"); 
    auto circ = new TEllipse(
        fSelectedSieveHole->cut_x, 
        fSelectedSieveHole->cut_y, 
        fSelectedSieveHole->cut_width,
        fSelectedSieveHole->cut_height, 
        0, 360,
        cut_angle
    ); 

    circ->SetFillStyle(0); 
    circ->SetLineColor(kRed); 
    circ->SetLineWidth(2); 
    circ->Draw(); 
    
    
#ifdef DEBUG
    cout << "Cut histogram drawn" << endl;  
#endif
    //gStyle->SetOptStat(0); 
    canv->cd(2); fY_fp.hist     ->Draw(draw_option); 
    canv->cd(3); fDxdz_fp.hist  ->Draw(draw_option);
    canv->cd(4); fDydz_fp.hist  ->Draw(draw_option);

#ifdef DEBUG
    cout << "FP-histograms drawn" << endl; 
#endif

    canv->Modified(); 
    canv->Update(); 
    canv->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "EvaluateCutFrame", this, "HandleCanvasClicked()"); 

    fFrame_canv->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,10)); 

    AddFrame(fFrame_canv, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10)); 

    fFrame_buttons = new TGHorizontalFrame(this, 1600, 50); 

    //Reject button
    // the Reject button is active in every state. 
    fButton_Reject = new TGTextButton(fFrame_buttons, "&Reject", 1); 
    fButton_Reject->Connect("Clicked()", "EvaluateCutFrame", this, "DoReject()"); 
    fFrame_buttons->AddFrame(fButton_Reject, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5));

    //Save button
    fButton_Save = new TGTextButton(fFrame_buttons, "&Save", 1); 
    fButton_Save->Connect("Clicked()", "EvaluateCutFrame", this, "DoSave()"); 
    fFrame_buttons->AddFrame(fButton_Save, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 
    
    //Fit button
    fButton_Fit = new TGTextButton(fFrame_buttons, "&Fit", 1); 
    fButton_Fit->Connect("Clicked()", "EvaluateCutFrame", this, "DoFit()"); 
    fFrame_buttons->AddFrame(fButton_Fit, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5));

    //Reset Fits button
    fButton_ResetFits = new TGTextButton(fFrame_buttons, "&Reset Fits", 1); 
    fButton_ResetFits->Connect("Clicked()", "EvaluateCutFrame", this, "DoResetFits()"); 
    fFrame_buttons->AddFrame(fButton_ResetFits, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5));

    //connect the exit of this app to the 'DoneEvaluate' 
    Connect("CloseWindow()", "PickSieveHoleApp", fParent, "DoneEvaluate()"); 

    AddFrame(fFrame_buttons, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 10,10,10,5 )); 


    SetWindowName("Evaluate hole cut"); 
    MapWindow();
    Resize(GetDefaultSize()); 
    MapSubwindows(); 
        
    // in this 'app state', we are picking limits to evaluate the fits within  
    fAppState = kPickLimits; 
    UpdateButtons(); 



#ifdef DEBUG 
    Info(here, "Exiting constructor"); 
#endif 
}   
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoResetFits()
{
#ifdef DEBUG
    Info(__func__, "Enter Fcn body"); 
#endif
    
    //Get a ptr to the canvas
    auto canv = fEcanvas->GetCanvas(); 
    if (!canv) {
        throw logic_error("in <EvaluateCutFrame::DoResetFits>: "
            "ptr to canvas is null"
        ); 
        return; 
    }

#ifdef DEBUG
    Info(__func__, "Wiping polynomials"); 
#endif
    //wipe the polynomials
    fSelectedSieveHole->y_fp    = NPoly(2); 
    fSelectedSieveHole->dxdz_fp = NPoly(2); 
    fSelectedSieveHole->dydz_fp = NPoly(2); 

#ifdef DEBUG
    Info(__func__, "Deleting drawn fits (if they exist)"); 
#endif
    //delete drawn fits
    DeleteDrawnObject(fY_fp.fit_hi); 
    DeleteDrawnObject(fY_fp.fit_lo); 
    DeleteDrawnObject(fDxdz_fp.fit_hi); 
    DeleteDrawnObject(fDxdz_fp.fit_lo); 
    DeleteDrawnObject(fDydz_fp.fit_hi); 
    DeleteDrawnObject(fDydz_fp.fit_lo); 

#ifdef DEBUG
    Info(__func__, "Flushing the canvas (and all sub-pads)"); 
#endif
    canv->Modified(); 
    canv->Update(); 

#ifdef DEBUG
    Info(__func__, "Setting fAppState = kPickLimits, and updating buttons"); 
#endif

    //update the buttons
    fAppState = kPickLimits; 
    UpdateButtons(); 
    DrawLimits(); 

#ifdef DEBUG
    Info(__func__, "End Fcn body"); 
#endif
}
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoFit()
{
#ifdef DEBUG
    Info(__func__, "Enter fcn body"); 
#endif  

    if (fAppState != kPickLimits) {
        throw logic_error("in <EvaluateCutFrame::DoFit>: "
            "Call to 'DoFit' made when application is not in state 'kPickLimits'; this should not be possible."
        ); 
        return; 
    }

    if (is_nan(fX_min) || is_nan(fX_max)) {
        throw logic_error("in <EvaluateCutFrame::DoFit>: "
            "One or both x_fp limits are NaN; it should not be possible to execute fit in this program state"
        );
        return;
    } 

    //now, we can do the fits.
    TCanvas *canv = fEcanvas->GetCanvas(); 

    if (!canv) return; 

    //disable other buttons while we're fitting
#ifdef DEBUG
    Info(__func__, "Setting fAppState to 'kNone"); 
#endif
    fAppState = kNone; 
    UpdateButtons(); 


    //draws the results of the fits on the tcanavs
    //_______________________________________________________________________________________________
    auto DrawFitsOnCanvas = [this](NPoly poly, HistAndLimit* hl, TVirtualPad* vpad)
    {   
        vpad->cd(); 

        auto fcn_poly = [poly](double *x, double *par){
            return poly.Eval({x[0], par[0]}); 
        };
        auto tf1_poly = new TF1("poly", fcn_poly, fX_min, fX_max, 1);

        tf1_poly->SetLineWidth(1); 

        auto& fit_lo = hl->fit_lo; 
        auto& fit_hi = hl->fit_hi;

        //draw the minimum raster value of our fit
        tf1_poly->SetParameter(0, -1.); 
        tf1_poly->SetLineColor(kRed); 

        DeleteDrawnObject(fit_lo); 
        fit_lo = (TF1*)tf1_poly->Clone(Form("%s_fit_lo",hl->hist->GetName()));
        fit_lo->Draw("SAME");

        //draw the maximum raster value of our fit
        tf1_poly->SetParameter(0, +1.);
        tf1_poly->SetLineColor(kBlack); 

        DeleteDrawnObject(fit_hi); 
        fit_hi = (TF1*)tf1_poly->Clone(Form("%s_fit_hi",hl->hist->GetName()));
        fit_hi->Draw("SAME");

        delete tf1_poly; 
    }; 
    //______________________________________________________________________________________________

    //these polynomials map from {x_fp, y_hcs} => {focal plane variables}
    NPoly poly_y_fp(2), poly_dxdz_fp(2), poly_dydz_fp(2); 

    //get the min/max of the raster, and split it up into 'n_raster_partitions'  
    cout << "Fitting polynomials___{ y... " << flush; 

    poly_y_fp    = FitPolynomialToFP(fCutFcn, &Trajectory_t::y); 
    if (poly_y_fp.Get_nDoF() != 2) {
        //something failed with this fit. tell the user to try again. 
        cout << "\n Fit for fp-coordinate 'y-fp' failed. please try again with different limits." << endl; 
        DoResetFits(); 
        return; 
    } 
    fSelectedSieveHole->y_fp    = poly_y_fp; 
    DrawFitsOnCanvas(poly_y_fp,    &fY_fp,      canv->cd(2));       
    canv->Modified(); canv->Update();  
    
    cout << "dx/dz... " << flush; 
    
    poly_dxdz_fp = FitPolynomialToFP(fCutFcn, &Trajectory_t::dxdz); 
    if (poly_dxdz_fp.Get_nDoF() != 2) {
        //something failed with this fit. tell the user to try again. 
        cout << "\n Fit for fp-coordinate 'dx/dz-fp' failed. please try again with different limits." << endl; 
        DoResetFits();
        return; 
    } 
    fSelectedSieveHole->dxdz_fp = poly_dxdz_fp; 
    DrawFitsOnCanvas(poly_dxdz_fp, &fDxdz_fp,   canv->cd(3)); 
    canv->Modified(); canv->Update();  
    
    cout << "dy/dz... " << flush; 
    
    poly_dydz_fp = FitPolynomialToFP(fCutFcn, &Trajectory_t::dydz); 
    if (poly_dydz_fp.Get_nDoF() != 2) {
        //something failed with this fit. tell the user to try again. 
        cout << "\n Fit for fp-coordinate 'dy/dz-fp' failed. please try again with different limits." << endl; 
        DoResetFits(); 
        return; 
    } 
    fSelectedSieveHole->dydz_fp = poly_dydz_fp; 
    DrawFitsOnCanvas(poly_dydz_fp, &fDydz_fp,   canv->cd(4)); 
    canv->Modified(); canv->Update();  
    
    cout << "}. done"<< endl; 

#ifdef DEBUG
    Info(__func__, "Setting fAppState to 'kEvaluated'"); 
#endif

    fAppState = kEvaluated; 
    UpdateButtons(); 
    DrawLimits(); 

#ifdef DEBUG
    Info(__func__, "Exit fcn body"); 
#endif
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::HandleCanvasClicked()
{
    if (fAppState != kPickLimits) return; 

    //get canvas event information
    TCanvas* canvas = fEcanvas->GetCanvas();
    if (!canvas) return;

    const map<string,int> hist_names_canv_map {
        {"h_y_fp",    2},
        {"h_dxdz_fp", 3},
        {"h_dydz_fp", 4}
    }; 
    
    //new event registered 
    if (fEventType != canvas->GetEvent()) {
        fEventType =  canvas->GetEvent(); 

        //mouse button 1 was just released
        if (fEventType == kMouseButton1_up) {

            //check if the selected object is a TH2D. if not, then quit. 
            auto selected_ptr = canvas->GetSelected(); 
            if (!selected_ptr || selected_ptr->IsA() != TClass::GetClass<TH2D>()) return; 

            TH2D* selected_hist = (TH2D*)selected_ptr; 

            //check if its one of the three histograms of ours.
            string selected_hist_name = string(selected_hist->GetName()); 

            if (selected_hist_name == "h_xy") return; //this is the 'hole drawing' hist; do nothing

            auto it = hist_names_canv_map.find(selected_hist_name); 
            if (it == hist_names_canv_map.end()) {
                ostringstream oss; 
                oss << "in <EvaluateCutFrame::HandleCanvasClicked>: invalid histogram name encountered. name: '" << selected_hist_name << "'"; 
                throw logic_error(oss.str()); 
                return; 
            }

            canvas->cd(0); 
            double x = canvas->cd(it->second)->AbsPixeltoX( canvas->GetEventX() );
            canvas->cd(0); 

            if (is_nan(fX_min)) { 
                fX_min = x; 
            } else {
                //in this case, only the first (min) has been selected. 
                if (is_nan(fX_max)) {

                    //the new limit picked is larger than the minimum the user already picked; fX_max = x; 
                    if (x > fX_min) { fX_max = x; }
                    //the new limit picked is smaller than the minimum the user already picked; fX_min = x; 
                    else            { fX_max = fX_min; fX_min = x; } 

                }
                //in this case, then both have been selected. wipe out both, and redraw them.     
                else                  { fX_min = x; fX_max = DOUBLE_NAN; }
            }
            
            DrawLimits(); 
            UpdateButtons(); 
        }
    }
}   
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DrawLimits()
{
    auto canv = fEcanvas->GetCanvas(); 
    if (!canv) {
        throw logic_error("in <EvaluateCutFrame::DrawLimits>: canvas from 'fEcanvas' is null.");
        return; 
    }
    
    if (is_nan(fX_min) || is_nan(fX_max)) return; 

    //determine whether or not to draw the limits
    //draw box limits on each histogram
    const float        fill_alpha = 0.30; 
    const unsigned int fill_color = kRed; 
    const unsigned int line_width = (fAppState == kEvaluated) ? 2 : 1; //draw thicker lines if the fits are evaluated
    const unsigned int fill_style = 3004; //diagonal hatching
    
    const unsigned int line_color = kRed;

    int i_canv =1;  

    for (HistAndLimit* hl : {&fY_fp, &fDxdz_fp, &fDydz_fp}) {

        canv->cd(++i_canv); 

        auto y_ax = hl->hist->GetYaxis(); 

        DeleteDrawnObject(hl->lim_low); 

        hl->lim_low = new TLine(
            fX_min, y_ax->GetXmin(), 
            fX_min, y_ax->GetXmax() 
        );  
         
        hl->lim_low->SetLineWidth(line_width); 
        hl->lim_low->SetLineColor(line_color); 
        hl->lim_low->Draw(); 
        
        DeleteDrawnObject(hl->lim_high); 

        hl->lim_high = new TLine(
            fX_max, y_ax->GetXmin(), 
            fX_max, y_ax->GetXmax() 
        );   
        
        hl->lim_high->SetLineWidth(line_width); 
        hl->lim_high->SetLineColor(line_color);
        hl->lim_high->Draw(); 
    }

    canv->cd(0); 
    canv->Modified();
    canv->Update(); 
}
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::UpdateButtons()
{
    //first -- unmap all buttons, *except* the reject button (the only button which ought to be usable in any state)
    fButton_Save     ->UnmapWindow(); 
    fButton_Fit      ->UnmapWindow(); 
    fButton_ResetFits->UnmapWindow(); 
    fButton_Reject   ->UnmapWindow();  
    
    fButton_Reject->MapWindow(); 

    switch (fAppState) {

        case kNone : break; 

        case kPickLimits : {
            
            //check if valid limits have been picked. if not, only the 'reject button should be available 
            if (is_nan(fX_min) || is_nan(fX_max)) break;
            
            //otherwise, we're ready to fit
            fButton_Fit->MapWindow(); 
            break; 
        }

        case kEvaluated : {

            //if valid limits have been picked. 
            fButton_ResetFits->MapWindow();
            fButton_Save->MapWindow(); 
            break; 
        }
    }
}
//__________________________________________________e___________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoSave() 
{
    if (fSelectedSieveHole) {
        fSelectedSieveHole->is_evaluated = true; 

        if (is_nan(fX_min) || is_nan(fX_max)) {
            throw logic_error("in <EvaluateCutFrame::DoSave>: fX_min and/or fX_max is NaN, which should not be possible.");
            return; 
        }
        fSelectedSieveHole->x_fp_min = fX_min; 
        fSelectedSieveHole->x_fp_max = fX_max; 

        //save the data of this sieve-hole cut
    } else {
        throw logic_error("in <EvaluateCutFrame::DoSave>: no sieve hole is selected (ptr is null)"); 
        return; 
    }
    //exit this window 
    CloseWindow();
    fParent->DoneEvaluate(); 
}
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoReject()
{
#ifdef DEBUG 
    cout << 
        "in <EvaluateCutFrame::DoReject>..."
        "\nptrs: "
        "\n   fSelectedSieveHole  : " << fSelectedSieveHole << 
        "\n   fParent             : " << fParent << endl;
#endif 


    if (fSelectedSieveHole) {
        fSelectedSieveHole->is_evaluated = false; 
        //reject the data of this sieve-hole (reset it)
    } else {
        throw logic_error("in <EvaluateCutFrame::DoSave>: no sieve hole is selected (ptr is null)"); 
        return; 
    }
    //exit this window 
    CloseWindow(); 
    fParent->DoneEvaluate(); 

#ifdef DEBUG 
    cout << "Exiting <EvaluateCutFrame::DoReject>" << endl; 
#endif 
} 
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DeleteDrawnObject(TObject* obj) {       

#ifdef DEBUG
    Info(__func__, "Enter fcn body"); 
#endif 
    auto pad = (TVirtualPad*)fEcanvas->GetCanvas(); 
    if (!pad) {  
        throw logic_error("in <EvaluateCutFrame::DeleteDrawnObject>: Ptr to canvas from fECanvas (2nd arg) is null");
        return;  
    }  

    if (!obj) {
#ifdef DEBUG
        printf("Info in <%s>: Ptr to TObject passed is null.\n", __func__);          
#endif
        return; 
    }

#ifdef DEBUG
    printf("Info in <%s>: Attepmting to delete TObject with name '%s'\n", __func__, obj->GetName());          
#endif 
    
    //remove the object from the TVirtualPad (and any sub-pads recursively)
    pad->RecursiveRemove(obj); 

#ifdef DEBUG
    Info(__func__, "Exit fcn body"); 
#endif
}
//_____________________________________________________________________________________________________________________________________
vector<EvaluateCutFrame::FitPoint_t> EvaluateCutFrame::CreatePointsFromHist(TH2D* hist) {

    if (!hist) {
        throw logic_error("in <EvaluateCutFrame::CreatePointsFromHist>: hist passed is null"); 
        return {}; 
    }

    //fit a graph and polynomial to each point
    const int min_stats = 10; 
    const double min_frac  = 0.10; 
    const double fit_radius = 8e-3; 
    auto x_axis = hist->GetXaxis(); 

    vector<FitPoint_t> pts; 

    double max_proj_integral=0.; 

    for (int b=1; b<x_axis->GetNbins(); b++) { 
    
        auto proj = hist->ProjectionY("proj",b,b); 

        if (proj->Integral()<min_stats) continue; 

        double fit_center = proj->GetXaxis()->GetBinCenter( proj->GetMaximumBin() ); 
        
        max_proj_integral = max<double>( max_proj_integral, proj->Integral() ); 
    
        auto gaus_fit = new TF1("gausFit", "gaus(0)", fit_center - fit_radius, fit_center + fit_radius); 
    
        gaus_fit->SetParameter(0, proj->GetMaximum());
        gaus_fit->SetParameter(1, fit_center);
        gaus_fit->SetParameter(2, 2.5e-3); 
    
        auto fit_result = proj->Fit("gausFit", "N Q S L R"); 

        if (!fit_result->IsValid()) continue; //skip if the fit failed

        pts.push_back({
            .x      = x_axis->GetBinCenter(b),              //center of bin
            .y      = fit_result->Parameter(1),             //center of gaussian, according to fit_result
            .sigma  = pow( fit_result->ParError(1), 2 ),    //error of this value from fit_result
            .N      = proj->Integral()                      //total stats for this profile 
        }); 
    }
    
    //prune the vectors. delete points which are less than  maxHeight * minStat_frac        
    for(auto it = pts.begin(); it != pts.end();) {
        if ( it->N < min_frac * max_proj_integral ) { pts.erase(it); }
        else                                        { it++; }
    }
    return pts; 
}
//_____________________________________________________________________________________________________________________________________
ROOT::RVec<double> EvaluateCutFrame::FitPolynomialToPoints(const vector<EvaluateCutFrame::FitPoint_t>& points, const int polynomial_degree)
{
    //minimum number of points
    if (points.size() < polynomial_degree+1) return {}; 

    RMatrix A(polynomial_degree+1, polynomial_degree+1, 0.); 

    ROOT::RVec<double> B(polynomial_degree+1, 0.); 

    for (const FitPoint_t& pt : points) {
        
        for (int i=0; i<=polynomial_degree; i++) { 

            B[i] += pt.y * pow( pt.x, i ) / pt.sigma; 
            
            for (int j=0; j<=polynomial_degree; j++) A.get(i,j) += pow( pt.x, i ) * pow( pt.x, j ) / pt.sigma; 
        }
    }
    auto coeffs = A.Solve( B ); 
    
    //check coeffs for NaN 
    for (double x : coeffs) if (x != x) { return{}; }

    return coeffs; 
}; 
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::Draw_Hist_Points_Poly(TH2D* hist, const vector<FitPoint_t>& points, const ROOT::RVec<double>& poly, const char* draw_option)
{
    hist->Draw(draw_option); 

    if (points.empty()) return; 

    //draw the points & poylnomials
    for (const auto& pt : points) {
        auto box = new TBox( 
            pt.x - 2e-3, pt.y - 10. * sqrt(pt.sigma), 
            pt.x + 2e-3, pt.y + 10. * sqrt(pt.sigma)
        ); 
        box->SetLineColor(1); 
        box->SetLineWidth(1); 
        box->SetFillStyle(0); 
        box->Draw("SAME"); 
    }

    if (poly.empty()) return; 

    char buff[50]; sprintf(buff, "pol%zi", poly.size()-1); 
    
    char name[200]; sprintf(name, "%s_poly", hist->GetName()); 

    if (!points.empty()) {

        auto f1_poly = new TF1(name, buff, fParent->GetDrawRange_xfp_min(), fParent->GetDrawRange_xfp_max()); 
        f1_poly->SetLineColor(kRed); 
        f1_poly->SetLineWidth(2); 
        int i=0; for (double coeff : poly) f1_poly->SetParameter(i++, coeff); 
        
        //draw the up and down limits of the cut
        f1_poly->SetParameter(0, poly[0] + fFpcoord_cut_width); 
        f1_poly->DrawCopy("SAME"); 

        f1_poly->SetParameter(0, poly[0] - fFpcoord_cut_width); 
        f1_poly->DrawCopy("SAME"); 
    } 
}
//_____________________________________________________________________________________________________________________________________
NPoly EvaluateCutFrame::FitPolynomialToFP(
    const function<bool(const EventData&)>& is_inside_cut, 
    double Trajectory_t::*coord
) const 
{
    using namespace ROOT::VecOps; 
#ifdef DEBUG
    Info(__func__, "In body..."); 
#endif

    //first, let's assemble our polynomial. 
    //create an empty polynomial (no elements) with 2 input DoF. 
    NPoly poly(2); 

    //add the elements which are functions of x-fp only.
    // the first element of 'powers' of each element is the exponent to which x_fp is raised. 
    // the second element is the power to which 'rast_param' is raised. 
    for (int i=0; i<=fPolynomialOrder; i++) poly.Add_element({.powers={i,0}, .coeff=1.}); 
    for (int i=0; i<=fRasterPolyOrder; i++) poly.Add_element({.powers={i,1}, .coeff=1.}); 
    
    //one quadratic element for 'rast_param' (i'm guessing here...). 
    poly.Add_element({.powers={0,2}, .coeff=1.}); 

    const unsigned int n_elems = poly.Get_nElems(); 

    //now, let's do a chi-square fit 
    RVec<double> B(n_elems, 0.);
    RMatrix A(n_elems,n_elems, 0.); 
    
    for (const auto& event : *fData) {

        //see if this event is inside the cut
        if (is_inside_cut(event)==false) continue; 

        //check to see if this even is inside the x_fp cut
        if (event.Xfp.x < fX_min || event.Xfp.x > fX_max ) continue; 

        auto&& Xmu = poly.Eval_noCoeff({event.Xfp.x, event.raster_index}); 
        
        double z = event.Xfp.*coord; 

        for (int i=0; i<n_elems; i++) {

            B[i] += z * Xmu[i]; 

            for (int j=0; j<n_elems; j++) A.get(i,j) += Xmu[i] * Xmu[j]; 
        }
    }   

    A.Set_report_singular(false); 
    auto&& coeffs = A.Solve(B); 

    //something went wrong with the chi-square 
    if (coeffs.size() != n_elems) { 
        
#ifdef DEBUG
        Info(__func__, "Chi-square it failed. returning null NPoly"); 
#endif
        return NPoly(0); 
    }
    //set the coeffs in our polynomial: 
    for (int i=0; i<n_elems; i++) { poly.Get_elem(i)->coeff = coeffs[i]; }

#ifdef DEBUG
    Info(__func__, "Chi-square Fit succeeded."); 
    Info("FitPolynomialToFP", "Printing polynomial..."); 
    poly.Print(); 
#endif

    //this chi2 fit is nice... but we are anticipating a lot of outlier data. 
    //thus, we need to find a way to cancel out this noise. 

    NPoly poly_objective(poly);     

    const double objective_sigma = 0.0025; 

    //____________________________________________________________________________________________
    auto objective_fcn = [n_elems,objective_sigma,&poly_objective,&is_inside_cut,coord,this](const double *par)
    {
        //set the parameters of the polynomal
        for (int i=0; i<n_elems; i++) poly_objective.Get_elem(i)->coeff = par[i]; 

        double val=0.; 

        for (const auto& event : *fData) {

            //check if this event is inside the cut
            if (is_inside_cut(event)==false) continue;  

            //if we just added this to the objective function sum, it would be a regular chi2-fit... 
            double arg = ( event.Xfp.*coord - poly_objective.Eval({event.Xfp.x, event.raster_index}) ) / objective_sigma;

            //but we stick it in the arg of a gaussian, so that we can mitigate the influence of outlier events. 
            val += -exp( -0.5*arg*arg ); 
        }

        return val; 
    };
    //____________________________________________________________________________________________
    ROOT::EnableImplicitMT(); 

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"); 
    
    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
#ifdef DEBUG
    minimizer->SetPrintLevel(2); 
#else 
    minimizer->SetPrintLevel(0); 
#endif 

    auto objective_functor = ROOT::Math::Functor(objective_fcn, n_elems);
    
    minimizer->SetFunction(objective_functor);

    //set the list of variables
    int i_var=0; 
    for (int i=0; i<n_elems; i++) {
        minimizer->SetVariable(i, Form("elem_%i",i), poly.Get_elem(i)->coeff, 1e-5);
    }
    
    bool fit_status = minimizer->Minimize();

    //if the fancy minuit fit fails, then just return the chi2-fit one. 
    if (!fit_status) return poly; 

    const double *coeffs_minuit = minimizer->X(); 

    for (int i=0; i<n_elems; i++) poly.Get_elem(i)->coeff = coeffs_minuit[i];

    return poly; 
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
EvaluateCutFrame::~EvaluateCutFrame() { 
    
    //delete histograms
    if (fHist_holes)   delete fHist_holes;
    if (fY_fp.hist)    delete fY_fp.hist; 
    if (fDxdz_fp.hist) delete fDxdz_fp.hist; 
    if (fDydz_fp.hist) delete fDydz_fp.hist; 

    Cleanup(); 
}
//_____________________________________________________________________________________________________________________________________
