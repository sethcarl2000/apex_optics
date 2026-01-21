#include "EvaluateCutFrame.h"
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <RMatrix.h>
#include <TBox.h>
#include <ApexOptics.h> 
#include <TGraphErrors.h> 
#include <RMatrix.h> 
#include <Math/Minimizer.h>
#include <Math/Functor.h> 
#include <Math/Factory.h> 
#include <NPoly.h> 
#include <RMatrix.h> 
#include <ROOT/RVec.hxx>

using namespace std; 
using ApexOptics::Trajectory_t; 

EvaluateCutFrame::EvaluateCutFrame( const TGWindow *p, 
                                    PickSieveHoleApp *_parent, 
                                    const vector<EventData>& data, 
                                    SieveHoleData *_hd,
                                    const double fp_cut_width, 
                                    const int n_raster_partitions, 
                                    const char* branch_x, 
                                    const char* branch_y, 
                                    const char* draw_option,
                                    const unsigned int palette )
    : TGMainFrame( p, 1400, 800 ), 
    fParent{_parent},
    fRDF{nullptr}, 
    fSelectedSieveHole{_hd}, 
    fFpcoord_cut_width{fp_cut_width}
{   
    const char* const here ="EvaluateCutFrame(constructor)"; 

#ifdef DEBUG
    Info(here, "In constructor body"); 
#endif

    const int polynomial_degree =3; 

    //setup the window
    SetCleanup(kDeepCleanup); 

    fFrame_canv = new TGHorizontalFrame(this, 1400, 800); 

    fEcanvas = new TRootEmbeddedCanvas("ECanvas_eval", fFrame_canv, 1400, 800); 

    const double cut_x      = fSelectedSieveHole->cut_x; 
    const double cut_y      = fSelectedSieveHole->cut_y; 
    const double cut_width  = fSelectedSieveHole->cut_width; 
    const double cut_height = fSelectedSieveHole->cut_height; 

    fHist_holes = new TH2D("h_xy_temp", "", 
        75, xsv_draw_range[0], xsv_draw_range[1],
        75, ysv_draw_range[0], ysv_draw_range[1] 
    ); 
    fY_fp    = HistAndLimit(new TH2D("h_y_fp",    "y_fp vs x_fp",     75, xfp_draw_range[0], xfp_draw_range[1], 75, -0.070, +0.055)); 
    fDxdz_fp = HistAndLimit(new TH2D("h_dxdz_fp", "dx/dz_fp vs x_fp", 75, xfp_draw_range[0], xfp_draw_range[1], 75, -0.035, +0.025)); 
    fDydz_fp = HistAndLimit(new TH2D("h_dydz_fp", "dy/dz_fp vs x_fp", 75, xfp_draw_range[0], xfp_draw_range[1], 75, -0.060, +0.040)); 
    
    for (const auto& ev : data) {

        if (pow((ev.Xsv.dxdz-cut_x)/cut_width, 2) + pow((ev.Xsv.dydz-cut_y)/cut_height, 2) < 1.) { 
            fY_fp.hist   ->Fill( ev.Xfp.x, ev.Xfp.y );
            fDxdz_fp.hist->Fill( ev.Xfp.x, ev.Xfp.dxdz ); 
            fDydz_fp.hist->Fill( ev.Xfp.x, ev.Xfp.dydz ); 
        }
        
        fHist_holes->Fill( ev.Xsv.dxdz, ev.Xsv.dydz );
    }

#ifdef DEBUG
    printf("Number of events in cut: %zi\n", data.size()); 
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
        fSelectedSieveHole->cut_height
    ); 
    circ->SetFillStyle(0); 
    circ->SetLineColor(kRed); 
    circ->SetLineWidth(2); 
    circ->Draw(); 

    auto is_inside_cut = [=](const EventData& event)
    {
        return pow((event.Xsv.dxdz-cut_x)/cut_width, 2) + pow((event.Xsv.dydz-cut_y)/cut_height, 2) < 1.; 
    };
    
    //draws the results of the fits on the tcanavs
    auto DrawFitsOnCanvas = [this](NPoly poly)
    {   
        auto fcn_poly = [poly](double *x, double *par){
            return poly.Eval({x[0], par[0]}); 
        };
        auto tf1_poly = new TF1("poly", fcn_poly, this->fXfp_draw_min, this->fXfp_draw_max, 1);

        tf1_poly->SetLineWidth(1); 

        //draw the minimum raster value of our fit
        tf1_poly->SetParameter(0, -1.); 
        tf1_poly->SetLineColor(kRed); 
        tf1_poly->DrawCopy("SAME");

        //draw the maximum raster value of our fit
        tf1_poly->SetParameter(0, +1.);
        tf1_poly->SetLineColor(kBlack); 
        tf1_poly->DrawCopy("SAME");

        delete tf1_poly; 
    }; 

    //get the min/max of the raster, and split it up into 'n_raster_partitions'  
    cout << "Fitting polynomials___{ y... " << flush; 
    auto poly_y_fp    = FitPolynomialToFP(data, is_inside_cut, &Trajectory_t::y); 
    cout << "dx/dz... " << flush; 
    auto poly_dxdz_fp = FitPolynomialToFP(data, is_inside_cut, &Trajectory_t::dxdz); 
    cout << "dy/dz... " << flush; 
    auto poly_dydz_fp = FitPolynomialToFP(data, is_inside_cut, &Trajectory_t::dydz); 
    cout << "}. done"<< endl; 
    
    
#ifdef DEBUG
    cout << "Cut histogram drawn" << endl;  
#endif
    //gStyle->SetOptStat(0); 
    canv->cd(2); fY_fp.hist     ->Draw(draw_option); DrawFitsOnCanvas(poly_y_fp); 
    canv->cd(3); fDxdz_fp.hist  ->Draw(draw_option); DrawFitsOnCanvas(poly_dxdz_fp); 
    canv->cd(4); fDydz_fp.hist  ->Draw(draw_option); DrawFitsOnCanvas(poly_dydz_fp); 

#ifdef DEBUG
    cout << "FP-histograms drawn" << endl; 
#endif

    canv->Modified(); 
    canv->Update(); 
    canv->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "EvaluateCutFrame", this, "HandleCanvasClicked()"); 

    fSelectedSieveHole->y_fp    = poly_y_fp; 
    fSelectedSieveHole->dxdz_fp = poly_dxdz_fp; 
    fSelectedSieveHole->dydz_fp = poly_dydz_fp; 

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

#ifdef DEBUG 
    Info(here, "Exiting constructor"); 
#endif 
}   
//_____________________________________________________________________________________________________________________________________
void EvaluateCutFrame::HandleCanvasClicked()
{
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
        auto f1_poly = new TF1(name, buff, xfp_draw_range[0], xfp_draw_range[1]); 
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
NPoly EvaluateCutFrame::FitPolynomialToFP(const vector<EventData>& data, function<bool(const EventData&)> is_inside_cut, double Trajectory_t::*coord)
{
    using namespace ROOT::VecOps; 
#ifdef DEBUG
    Info("FitPolynomialToFP", "In body..."); 
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
    
    for (const auto& event : data) {

        //see if this event is inside the cut
        if (is_inside_cut(event)==false) continue; 

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
        return NPoly(0); 
#ifdef DEBUG
        cout << "Fit failed." << endl; 
#endif
    }
    //set the coeffs in our polynomial: 
    for (int i=0; i<n_elems; i++) { poly.Get_elem(i)->coeff = coeffs[i]; }

#ifdef DEBUG
    cout << "Fit succeeded." << endl; 
    Info("FitPolynomialToFP", "Printing polynomial..."); 
    poly.Print(); 
#endif

    //this chi2 fit is nice... but we are anticipating a lot of outlier data. 
    //thus, we need to find a way to cancel out this noise. 

    NPoly poly_objective(poly);     

    const double objective_sigma = 0.0025; 

    //____________________________________________________________________________________________
    auto objective_fcn = [n_elems,objective_sigma,&poly_objective,&data,&is_inside_cut,coord](const double *par)
    {
        //set the parameters of the polynomal
        for (int i=0; i<n_elems; i++) poly_objective.Get_elem(i)->coeff = par[i]; 

        double val=0.; 

        for (const auto& event : data) {

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
