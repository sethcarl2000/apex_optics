#include "EvaluateCutFrame.h"
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <RMatrix.h>
#include <TBox.h>

using namespace std; 


EvaluateCutFrame::EvaluateCutFrame( const TGWindow *p, 
                                    PickSieveHoleApp *_parent, 
                                    ROOT::RDF::RNode *_rdf, 
                                    SieveHoleData *_hd,
                                    const double fp_cut_width, 
                                    const char* branch_x, 
                                    const char* branch_y, 
                                    const char* draw_option,
                                    const unsigned int palette )
    : TGMainFrame( p, 1400, 800 ), 
    fParent{_parent},
    fRDF{_rdf}, 
    fSelectedSieveHole{_hd}, 
    fFpcoord_cut_width{fp_cut_width}
{   
    const int polynomial_degree =3; 

    //setup the window
    SetCleanup(kDeepCleanup); 

    fFrame_canv = new TGHorizontalFrame(this, 1400, 800); 

    fEcanvas = new TRootEmbeddedCanvas("ECanvas_eval", fFrame_canv, 1400, 800); 

    const double cut_x      = fSelectedSieveHole->cut_x; 
    const double cut_y      = fSelectedSieveHole->cut_y; 
    const double cut_width  = fSelectedSieveHole->cut_width; 
    const double cut_height = fSelectedSieveHole->cut_height; 

    auto hist_xy = (*fRDF).Histo2D<double>({"h_xy_temp", "", 
            200, xsv_draw_range[0], xsv_draw_range[1],
            200, ysv_draw_range[0], ysv_draw_range[1] }, branch_x, branch_y);     

    fHist_holes = (TH2D*)hist_xy->Clone("h_xy"); 

    auto df_output = (*fRDF)

        .Filter([cut_x, cut_y, cut_width, cut_height](double x, double y)
        {
            return pow( (x - cut_x)/cut_width, 2 ) + pow( (y - cut_y)/cut_height, 2 ) < 1.; 
        }, {branch_x, branch_y}); 
    
    
    auto hist_y_fp    = df_output.Histo2D<double>({"h_y_temp",    "y_fp vs x_fp",     30, -0.65, 0.65, 75, -0.070, 0.055}, "x_fp", "y_fp"); 
    auto hist_dxdz_fp = df_output.Histo2D<double>({"h_dxdz_temp", "dx/dz_fp vs x_fp", 30, -0.65, 0.65, 75, -0.035, 0.025}, "x_fp", "dxdz_fp"); 
    auto hist_dydz_fp = df_output.Histo2D<double>({"h_dydz_temp", "dy/dz_fp vs x_fp", 30, -0.65, 0.65, 75, -0.060, 0.040}, "x_fp", "dydz_fp"); 

    TCanvas *canv = fEcanvas->GetCanvas(); 
   
    gStyle->SetPalette(palette); 
    
    canv->cd(); 
    canv->Divide(2,2, 0.001, 0.001); 

    canv->cd(1); fHist_holes->Draw("col2"); 
    auto circ = new TEllipse(
        fSelectedSieveHole->cut_x, 
        fSelectedSieveHole->cut_y, 
        fSelectedSieveHole->cut_width,
        fSelectedSieveHole->cut_height
    ); 
    circ->SetFillStyle(0); 
    circ->SetLineColor(kRed); 
    circ->SetLineWidth(2); 
    circ->Draw(); 

    //create the histograms
    fY_fp    = HistAndLimit((TH2D*)hist_y_fp    ->Clone("h_yfp"));
    fDxdz_fp = HistAndLimit((TH2D*)hist_dxdz_fp ->Clone("h_dxdzfp")); 
    fDydz_fp = HistAndLimit((TH2D*)hist_dydz_fp ->Clone("h_dydzfp")); 

    auto Draw_and_fit = [&](TH2D* hist, ROOT::RVec<double>& poly)
    {
        auto points = CreatePointsFromHist(hist); 
        poly = FitPolynomialToPoints(points, polynomial_degree); 
        Draw_Hist_Points_Poly(hist, points, poly, draw_option); 
    }; 

    canv->cd(2); Draw_and_fit(fY_fp.hist,    fSelectedSieveHole->y_fp.poly); 
    
    canv->cd(3); Draw_and_fit(fDxdz_fp.hist, fSelectedSieveHole->dxdz_fp.poly); 
    
    canv->cd(4); Draw_and_fit(fDydz_fp.hist, fSelectedSieveHole->dydz_fp.poly); 
    
    canv->Modified(); 
    canv->Update(); 
    canv->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "EvaluateCutFrame", this, "HandleCanvasClicked()"); 

    fFrame_canv->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,10)); 

    AddFrame(fFrame_canv, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10)); 

    fFrame_buttons = new TGHorizontalFrame(this, 1600, 50); 

    //Reject button
    fButton_Save = new TGTextButton(fFrame_buttons, "&Save", 1); 
    fButton_Save->Connect("Clicked()", "EvaluateCutFrame", this, "DoSave()"); 
    fFrame_buttons->AddFrame(fButton_Save, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //Reject button
    fButton_Reject = new TGTextButton(fFrame_buttons, "&Reject", 1); 
    fButton_Reject->Connect("Clicked()", "EvaluateCutFrame", this, "DoReject()"); 
    fFrame_buttons->AddFrame(fButton_Reject, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 20, 10, 5, 5)); 

    //connect the exit of this app to the 'DoneEvaluate' 
    Connect("CloseWindow()", "PickSieveHoleApp", fParent, "DoneEvaluate()"); 

    AddFrame(fFrame_buttons, new TGLayoutHints(kLHintsBottom | kLHintsExpandX, 10,10,10,5 )); 


    SetWindowName("Evaluate hole cut"); 
    MapWindow();
    Resize(GetDefaultSize()); 
    MapSubwindows(); 

    UpdateButtons(); 
}   
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::HandleCanvasClicked()
{
    //get canvas event information
    TCanvas* canvas = fEcanvas->GetCanvas();
    if (!canvas) return;

    const map<string,int> hist_names_canv_map {
        {"h_yfp",    2},
        {"h_dxdzfp", 3},
        {"h_dydzfp", 4}
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
                throw logic_error("in <EvaluateCutFrame::HandleCanvasClicked>: invalid histogram name encountered."); 
                return; 
            }

            canvas->cd(0); 
            double x = canvas->cd(it->second)->AbsPixeltoX( canvas->GetEventX() );
            canvas->cd(0); 

            if (fX_min != fX_min) { 
                fX_min = x; 
            } else {
                 //in this case, only the first (min) has been selected. 
                if (fX_max != fX_max) { fX_max = x; }
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

    //determine whether or not to draw the limits
    //draw box limits on each histogram
    const float        fill_alpha = 0.30; 
    const unsigned int fill_color = kRed; 
    const unsigned int fill_style = 3004; //diagonal hatching
    
    const unsigned int line_color = kRed;

    int i_canv =1;  

    for (HistAndLimit* hl : {&fY_fp, &fDxdz_fp, &fDydz_fp}) {

        canv->cd(++i_canv); 

        if (fX_min != fX_min || fX_max != fX_max) return; 
            
        auto y_ax = hl->hist->GetYaxis(); 

        if (hl->lim_low) delete hl->lim_low; 
        hl->lim_low = new TLine(
            fX_min, 
            y_ax->GetXmin(), 
            fX_min,
            y_ax->GetXmax() 
        );  
         
        hl->lim_low->SetLineColor(line_color); 
        hl->lim_low->Draw(); 
            

        if (hl->lim_high) delete hl->lim_high; 
        hl->lim_high = new TLine(
            fX_max, 
            y_ax->GetXmin(), 
            fX_max,
            y_ax->GetXmax() 
        );   

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
    if (fX_min == fX_min && fX_max == fX_max) { 
        fButton_Save->MapWindow(); 
    } else { 
        fButton_Save->UnmapWindow(); 
    }
}
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::DoSave() 
{
    if (fSelectedSieveHole) {
        fSelectedSieveHole->is_evaluated = true; 

        if (fX_max != fX_max || fX_min != fX_min) {
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
        auto f1_poly = new TF1(name, buff, -0.65, 0.65); 
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
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________________________
EvaluateCutFrame::~EvaluateCutFrame() { 
    
    //delete histograms
    Cleanup(); 
}
//_____________________________________________________________________________________________________________________________________
