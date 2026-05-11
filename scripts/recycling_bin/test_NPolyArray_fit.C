#include <NPolyArray.h> 
#include <stdexcept> 
#include <ROOT/RDataFrame.hxx> 
#include <ApexOptics.h> 
#include <ROOT/RVec.hxx>
#include <ROOT/RResultPtr.hxx> 
#include "include/RDFNodeAccumulator.h"
#include <iostream> 
#include <string> 
#include <vector> 
#include <TCanvas.h> 
#include <TStyle.h> 

using namespace std; 
using namespace ROOT::VecOps; 

using ApexOptics::Trajectory_t; 


//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Trajectory_t struct. 
ROOT::RDF::RNode add_branch_from_Trajectory_t(ROOT::RDF::RNode df, const char* branch_in, map<string, double Trajectory_t::*> branches)
{
    const int n_nodes = branches.size() + 1; 
    RVec<ROOT::RDF::RNode> nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Trajectory_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Trajectory_t::coord' 
        nodes.push_back( nodes.back()
            
            .Define(branch_name, [coord](const Trajectory_t& track) { return track.*coord; }, {branch_in})
        ); 
    }
    
    return nodes.back(); 
}


int test_NPolyArray_fit(const char* path_infile, const char* path_NPolyArray, const vector<string>& inputs, const vector<string> outputs) 
{
    const char* const here = "test_NPolyArray_fit"; 

    //try to parse the polynomials from files. catch any exceptions, if they occur. 
    NPolyArray parr; 
    
    try {
        parr = ApexOptics::Parse_NPolyArray_from_file(path_NPolyArray, outputs, inputs.size()); 
    } catch (const std::exception& e) {
        
        Error(here, "Something went wrong trying to parse NPolyArray.\n what(): %s", e.what()); 
        return -1; 
    }

    
    //now, try to construct the RDataFrame. 

    ROOT::RDataFrame *df = nullptr; 

    try {        
        ROOT::EnableImplicitMT();         
        df = new ROOT::RDataFrame("tracks_fp", path_infile); 
    
    } catch (const std::exception& e) {

        Error(here, "Someting went wrong trying to construct the RDataFrame.\n what(): %s", e.what()); 
        return -1; 
    } 

    RDFNodeAccumulator rna(*df); 

    // 'X' is the input vector, 'Z' is the output vector
    rna.Define("X", [](){ return RVec<double>{}; }, {}); 
    rna.Define("Z", [](){ return RVec<double>{}; }, {}); 
    
    for (const auto& input : inputs) {
        rna.Overwrite("X", [](RVec<double>& V, double x){ V.push_back(x); return V; }, {"X", input.c_str()}); 
    }

    for (const auto& output : outputs) {
        rna.Overwrite("Z", [](RVec<double>& V, double x){ V.push_back(x); return V; }, {"Z", output.c_str()}); 
    }

    //now, perform the reconstruction. 
    rna.Define("reco_Z", [&parr](const RVec<double>& X){ return parr.Eval(X); }, {"X"}); 

    vector<ROOT::RDF::RResultPtr<TH2D>> histos{}; 

    int i_out=0; 
    for (const auto& output : outputs) {

        const char* cstr_output = output.c_str(); 
        const char* reco_output = Form("reco_%s", cstr_output); 

        rna.Overwrite(reco_output, [i_out](const RVec<double>& V){ return V.at(i_out); }, {"Z"}); 
        i_out++; 

        const char* h_name  = Form("h_%s", cstr_output); 
        const char* h_title = Form("Reco - %s;%s;%s-reco", cstr_output,cstr_output,cstr_output); 

        histos.push_back( 
            rna.Get().Histo2D<double>({h_name, h_title, 200, -1, -1, 200, -1, -1}, cstr_output, reco_output)
        ); 
    }

    auto canv = new TCanvas("c", "canv", 1600, 650);

    canv->Divide(2,2); 
    
    gStyle->SetOptStat(0); 
    gStyle->SetPalette(kSunset); 
    
    int i_canv=1; 
    for (auto& hist : histos) {

        canv->cd(i_canv); 
        hist->DrawCopy("col2"); 

        i_canv++; 
        //print only the first 4 histos
        if (i_canv>4) break; 
    }
    
    return 0; 
}