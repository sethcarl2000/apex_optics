#include <TCanvas.h> 
#include <TGraph.h> 
#include <TPad.h>
#include <TFile.h> 
#include <TTree.h> 
#include <TStyle.h> 
#include <vector>
#include <string> 

using namespace std; 

int plot_all_coordinates(const char* path_infile)
{
    auto file = new TFile(path_infile, "READ"); 

    auto tree = (TTree*)file->Get("tracks_fp"); 

    auto canv = new TCanvas("c", "All coordinate correlations", 1600, 600); 

    int i_canv=1; 

    gStyle->SetOptStat(0); 
    gStyle->SetPalette(kSunset); 

    const vector<string> coords_horiz{"x_q1", "dxdz_q1",           "y_q1", "dydz_q1", "dpp_q1"};  
    const vector<string> coords_vert {"x_fp", "dxdz_fp + x_fp/6.", "y_fp", "dydz_fp"};  
    
    canv->Divide(coords_horiz.size(),coords_vert.size(), 0,0); 

    for (const auto& str_v : coords_vert) {
        for (const auto& str_h : coords_horiz) {

            char* buff = Form("%s:%s>>h_%i(200,-1,-1,200,-1,-1)", str_v.c_str(), str_h.c_str(), i_canv); 
            
            canv->cd(i_canv); 
            
            tree->Draw(buff, "", "col2");
            
            i_canv++; 
        }
    }

    return 0; 
}