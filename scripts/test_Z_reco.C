#include <ApexOptics.h>
#include <NPoly.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <ROOT/RDataFrame.hxx>
#include <TLegend.h>

using namespace std; 

int test_Z_reco(const char* path_poly)
{
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df_v1("tracks_fp", "data/sieve_holes/fits_29Mar/V1.root"); //"data/replay/real_R_V1-v01.root"); 
    ROOT::RDataFrame df_v2("tracks_fp", "data/sieve_holes/fits_29Mar/V2.root"); //"data/replay/real_R_V2-v01.root"); 
    ROOT::RDataFrame df_v3("tracks_fp", "data/sieve_holes/fits_29Mar/V3.root"); //"data/replay/real_R_V3-v01.root"); 

    const char* poly_name = "z_hcs_reco"; 
    //parse the polynomial
    NPoly pol_z(4); 
    ApexOptics::Parse_NPoly_from_file(path_poly, poly_name, &pol_z); 

    if (pol_z.Get_nElems() < 1) {
        Error(__func__, "Zero elements of polynomial '%s' parsed from file '%s'.", poly_name, path_poly); 
        return -1; 
    }

    auto hist_v1 = (TH1D*)df_v1
        
        .Define("z_hcs", [&pol_z](double x, double y, double dxdz, double dydz, double dpp)
        {
            return pol_z.Eval({x,y,dxdz,dydz}); 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp", "dpp_sv"})
        
        .Histo1D<double>({"h_z_v1", Form("z_{hcs} reco: %s;z_{hcs} (m);",path_poly), 200, -0.5,+0.5}, "z_hcs")->Clone(); 

    
    auto hist_v2 = (TH1D*)df_v2
        
        .Define("z_hcs", [&pol_z](double x, double y, double dxdz, double dydz, double dpp)
        {
            return pol_z.Eval({x,y,dxdz,dydz}); 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp", "dpp_sv"})
        
        .Histo1D<double>({"h_z_v2", Form("z_{hcs} reco: %s;z_{hcs} (m);",path_poly), 200, -0.5,+0.5}, "z_hcs")->Clone(); 

    auto hist_v3 = (TH1D*)df_v3
        
        .Define("z_hcs", [&pol_z](double x, double y, double dxdz, double dydz, double dpp)
        {
            return pol_z.Eval({x,y,dxdz,dydz}); 
        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp", "dpp_sv"})
        
        .Histo1D<double>({"h_z_v3", Form("z_{hcs} reco: %s;z_{hcs} (m);",path_poly), 200, -0.5,+0.5}, "z_hcs")->Clone(); 


    double histo_max = 0.; 
    histo_max = max<double>( histo_max, hist_v1->GetMaximum() ); 
    histo_max = max<double>( histo_max, hist_v2->GetMaximum() ); 
    histo_max = max<double>( histo_max, hist_v3->GetMaximum() ); 
        
    new TCanvas; 
    gStyle->SetOptStat(0); 

    hist_v2->SetFillStyle(0); 
    hist_v2->SetFillColor(kBlack); 
    hist_v2->SetLineColor(kBlack);
    hist_v2->GetYaxis()->SetRangeUser( 0., histo_max * 1.05 );  
    hist_v2->Draw(); 

    hist_v1->SetFillStyle(3003); 
    hist_v1->SetFillColor(kBlue); 
    hist_v1->SetLineColor(kBlue); 
    hist_v1->Draw("SAME"); 

    hist_v3->SetFillStyle(3004); 
    hist_v3->SetFillColor(kRed); 
    hist_v3->SetLineColor(kRed); 
    hist_v3->Draw("SAME"); 

    auto legend = new TLegend; 

    legend->AddEntry(hist_v1, "V1"); 
    legend->AddEntry(hist_v2, "V2"); 
    legend->AddEntry(hist_v3, "V3"); 
    legend->Draw(); 

    return 0; 
}