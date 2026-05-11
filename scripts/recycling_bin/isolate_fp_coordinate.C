#include <TCanvas.h> 
#include <ROOT/RDataFrame.hxx>
#include <TStyle.h> 
#include <TColor.h> 
#include <TH2D.h> 
#include <TEllipse.h> 


int isolate_fp_coordinate() {

    //create the dataframe
    ROOT::RDataFrame df("tracks_fp", "data/mc/mc_L_production.root"); 

    const double cent_x_fp = 0.; 
    const double cent_y_fp = 0.02; 
    const double cent_dxdz_fp = 0.01;
    const double cent_dydz_fp = 0.; 
    
    double scale_x_fp     = 0.6;
    double scale_dxdz_fp  = 0.03; 

    double scale_y_fp     = 0.06; 
    double scale_dydz_fp  = 0.05; 

    const double decay = 0.95; 

    gStyle->SetPalette(kSunset);
    gStyle->SetOptStat(0);  

    auto c = new TCanvas("c", "cut", 1600, 800); 
    c->Divide(2,2, 0.01,0.01); 

    for (size_t i=0; i<50; i++) {
        
        auto df_out = df
            .Filter([=](double x, double y, double dxdz, double dydz)
            { 
                return 
                    pow((x - cent_x_fp)/scale_x_fp, 2) + 
                    pow((dxdz - cent_dxdz_fp)/scale_dxdz_fp, 2) < 1. 
                    &&
                    pow((y - cent_y_fp)/scale_y_fp, 2) +
                    pow((dydz - cent_dydz_fp)/scale_dydz_fp, 2) < 1.;   

            }, {"x_fp","y_fp","dxdz_fp","dydz_fp"}); 


        auto h_x_y_sieve = df_out.Histo2D<double>({"h_xy_sv", "y_{sv} vs x_{sv}", 200, -0.04, 0.04, 200, -0.05, 0.03}, "x_sv", "y_sv"); 

        auto h_dx_dy_sieve = df_out.Histo2D<double>({"h_xy_sv", "dy/dz_{sv} vs dx/dz_{sv}", 200, -0.04, 0.04, 200, -0.05, 0.03}, "dxdz_sv", "dydz_sv"); 
            
        auto h_x_dx = df_out.Histo2D<double>({"h_x_dx", "dx/dz_{fp} vs x_{fp}", 200, -0.6, 0.6, 200, -0.03, 0.03}, "x_fp", "dxdz_fp"); 

        auto h_y_dy = df_out.Histo2D<double>({"h_y_dy", "dy/dz_{fp} vs y_{fp}", 200, -0.06, 0.06, 200, -0.05, 0.05}, "y_fp", "dydz_fp"); 

        c->cd(1); h_x_y_sieve->DrawCopy("col2"); 
        c->cd(3); h_dx_dy_sieve->DrawCopy("col2"); 

        c->cd(2); h_x_dx->DrawCopy("col2"); 
        auto circ1 = new TEllipse(cent_x_fp, cent_dxdz_fp, scale_x_fp, scale_dxdz_fp);
        circ1->SetFillStyle(0); 
        circ1->SetLineColor(kRed); 
        circ1->Draw();

        c->cd(4); h_y_dy->DrawCopy("col2"); 
        auto circ2 = new TEllipse(cent_y_fp, cent_dydz_fp, scale_y_fp, scale_dydz_fp);
        circ2->SetFillStyle(0); 
        circ2->SetLineColor(kRed); 
        circ2->Draw();

        c->SaveAs("histos/constriction.gif+20"); 

        c->Modified(); 
        c->Update(); 

        scale_x_fp    *= decay; 
        scale_y_fp    *= decay; 
        scale_dxdz_fp *= decay; 
        scale_dydz_fp *= decay; 
    }

    return 0; 
}