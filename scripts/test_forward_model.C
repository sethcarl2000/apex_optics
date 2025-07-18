
#include <fstream>
#include <iostream>
#include <sstream>
#include <string> 
#include <map>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

using namespace std; 
using namespace ROOT::VecOps; 

//this will be a temporary funct., that I may absorb into NPoly.cxx. 
int parse_poly_from_file(const char* path_dbfile, const char* poly_name, NPoly *poly) 
{
    //uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuer89999999999999999999999999999999999uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
    // -bear 

    //parse each line separately. in each line, each 'token' is separated by whitespace. 
    //each line must have the following format: 
    //
    //      dxdz_sv   0   0   2   0   -3.044778720e-01
    //
    //The first token is the name of the polynomial to which this element belongs.
    // the next 4 tokens are the power to which this polynomial must be raised. 
    // the last token is the coefficient of this element.
    
    //The file must be headed by the line: 
    //
    //      poly-DoF 4
    
    const char* const here = "parse_poly_from_file"; 

    //open the db file
    ifstream dbfile(path_dbfile); 
    
    if (!dbfile.is_open()) {
        Error(here, "unable to open db file '%s'", path_dbfile); 
        return -1; 
    }

    //now, read the file
    string line; 
    
    //get the DoF of this poly
    getline(dbfile, line); 
    istringstream iss_init(line);
    
    string token; 

    iss_init >> token; 
    if (token != "poly-DoF") {
        Error(here, "Missing 'poly-DoF [n]' header at top of dbfile '%s'", path_dbfile); 
        return -1; 
    }

    int poly_DoF; 
    iss_init >> poly_DoF; 

    if (poly_DoF != poly->Get_nDoF()) {
        Error(here, "Poly DoF in db file (%i) does not match DoF of passed Poly (%i)", 
            poly_DoF, poly->Get_nDoF()); 
        return -1; 
    }

    int start_nElems = poly->Get_nElems(); 

    //now, we can ready the rest of the file. 
    while (getline(dbfile, line)) {

        //parse the string into token (delimited by whitespace!)
        istringstream iss(line); 

        string elem_name; iss >> elem_name; 

        //this line is not the poly you're looking for
        if (elem_name != poly_name) continue; 
        
        //now, we can read the powers / coefficient
        RVec<int> powers(poly_DoF, 0); 
        double coeff; 

        //read the powers
        for (int &pow : powers) iss >> pow;

        //read the coefficient
        iss >> coeff; 

        poly->Add_element(powers, coeff); 
    }
    
    dbfile.close(); 

    return poly->Get_nElems() - start_nElems; //noop
}

struct Track_t { 
    double x,y,dxdz,dydz; 
}; 

//_______________________________________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Track_t struct. 
ROOT::RDF::RNode add_branch_from_Track_t(ROOT::RDF::RNode df, const char* branch_in, map<string, double Track_t::*> branches)
{
    const int n_nodes = branches.size() + 1; 
    RVec<ROOT::RDF::RNode> df_nodes; 

    df_nodes.push_back(df); 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Track_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Track_t::coord' 
        auto new_node = df_nodes.at(df_nodes.size()-1)

            .Define(branch_name, [coord](const RVec<Track_t>& tracks)
            {
                RVec<double> ret; ret.reserve(tracks.size()); 

                for (const Track_t& track : tracks) ret.push_back( track.*coord ); 

                return ret; 

            }, {branch_in}); 

        df_nodes.push_back(new_node); 
    }
    
    return df_nodes.at(df_nodes.size()-1); 
}

//_______________________________________________________________________________________________________________________________________________
int test_forward_model( const char* path_infile="data/replay/replay.4768.root",
                        const char* path_dbfile="data/csv/db_mc_fp_sv_V1_L_2ord.dat",  
                        const char* tree_name="track_data" ) 
{
    const char* const here = "test_forward_model"; 

    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }
    
    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(here, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the proper TTree
    TTree* tree = (TTree*)infile->Get(tree_name); 
    if (!tree) {
        Error(here, "could not find TTree '%s'", tree_name); 
        return 1; 
    }
       
    delete tree; 
    infile->Close(); 
    delete infile; 

    //try to get some meta data from the db file: 
    ifstream dbfile(path_dbfile); 
    
    //check if the file can be opened
    if (!dbfile.is_open()) {
        Error(here, "unable to open db file '%s'", path_dbfile); 
        return -1; 
    }

    //now, read the file
    string line;
    istringstream iss_init;
    string token; 
    int int_buffer; 
    
    //read the DoF of the poly
    getline(dbfile, line); iss_init = istringstream(line);
    
    iss_init >> token; 
    if (token != "poly-DoF") {
        Error(here, "Missing 'poly-DoF [n]' header at top of dbfile '%s'", path_dbfile); 
        return -1; 
    }
    iss_init >> int_buffer; 
    const int poly_DoF = int_buffer; 


    //read which arm to use
    getline(dbfile, line); iss_init = istringstream(line);
    
    iss_init >> token; 
    if (token != "is-RHRS") {
        Error(here, "Missing 'is-RHRS [1/0]' header at top of dbfile '%s'", path_dbfile); 
        return -1; 
    }
    iss_init >> int_buffer; 
    const bool is_RHRS = (int_buffer==1); 

    dbfile.close(); 



    //now, read the polynomials from the data file
    map<string, unique_ptr<NPoly>> pols; 

    //add the polynomials we want to parse
    pols["x_sv"]    = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["y_sv"]    = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dxdz_sv"] = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dydz_sv"] = unique_ptr<NPoly>(new NPoly(poly_DoF)); 

    //parse all elements for each of these
    cout << "parsing elements for all polynomials...\n" << flush; 
    
    for (auto it = pols.begin(); it != pols.end(); it++) { 
        
        string poly_name = it->first; 
        
        //try to parse all elements for this polynomial
        int elems_found = parse_poly_from_file(path_dbfile, poly_name.data(), it->second.get()); 

        if (elems_found < 0) {
            Error(here, "fatal error trying to parse polynomial '%s'. db file: '%s'", poly_name.data(), path_dbfile); 
            return 1; 
        }

        if (elems_found < 1) {
            Warning(here, "did not parse any elements for polynomial '%s', check name.", poly_name.data()); 
        }

        printf("-- %3i elems found for poly '%s'\n", it->second->Get_nElems(), poly_name.data()); 
    }
    cout << "parsing done." << endl; 

    //now, we're ready to deal with the data. 


    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 

    //probably not the most elegant way to do this, but here we are. 
    const NPoly *pol_x      = pols["x_sv"].get(); 
    const NPoly *pol_y      = pols["y_sv"].get(); 
    const NPoly *pol_dxdz   = pols["dxdz_sv"].get(); 
    const NPoly *pol_dydz   = pols["dydz_sv"].get(); 
    
    const char* bn_x_tra    = is_RHRS ? "R_tr_tra_x"  : "L_tr_tra_x"; 
    const char* bn_y_tra    = is_RHRS ? "R_tr_tra_y"  : "L_tr_tra_y"; 
    const char* bn_dxdz_tra = is_RHRS ? "R_tr_tra_th" : "L_tr_tra_th"; 
    const char* bn_dydz_tra = is_RHRS ? "R_tr_tra_ph" : "L_tr_tra_ph"; 
    
    auto df_reco = df 
        
        .Define("tracks_fp", 
            []
            (RVec<double> &v_x,
             RVec<double> &v_y,
             RVec<double> &v_dxdz,
             RVec<double> &v_dydz)
            {   
                //Here, we will collect all VDC trajectory information into a single struct, and convert from TRANSPORT Coordinates (tra)
                // to FOCAL PLANE Coordinates (fp). 

                //create the vector of our 'Track_t' struct. '.reserve()' tells the vector what size it should expect to eventually be. 
                RVec<Track_t> tracks; 
                tracks.reserve(v_x.size());

                //fill our new vector
                for (int i=0; i<v_x.size(); i++) {
                    tracks.push_back({
                        .x = v_x.at(i),
                        .y = v_y.at(i),
                        .dxdz = v_dxdz.at(i) - v_x.at(i)/6.,  //the only difference between TRANSPORT (tra) and FOCAL-PLANE (fp) coordinates.
                        .dydz = v_dydz.at(i) 
                    });    
                }
                
                return tracks; 

            }, {bn_x_tra, bn_y_tra, bn_dxdz_tra, bn_dydz_tra}) 
        
        .Define("tracks_sv", 
            [pol_x,pol_y,pol_dxdz,pol_dydz]
            (const RVec<Track_t> &tracks_fp)
            {
                //as before, create a vector of Track_t structs to store the sieve tracks. 
                RVec<Track_t> tracks_sv; 
                tracks_sv.reserve(tracks_fp.size());

                //loop through all tracks; use the 'forward' polynomials to reconstruct 'sieve' data. 
                for (const Track_t& trk_fp : tracks_fp) {

                    //convert the Track_t struct to an RVec<double>
                    RVec<double> X_fp = {
                        trk_fp.x,
                        trk_fp.y,
                        trk_fp.dxdz,
                        trk_fp.dydz
                    }; 
                    
                    //use ouqqr polynomials to compute the sieve coordinates
                    tracks_sv.push_back({
                        .x    = pol_x->Eval(X_fp),
                        .y    = pol_y->Eval(X_fp),
                        .dxdz = pol_dxdz->Eval(X_fp),
                        .dydz = pol_dydz->Eval(X_fp)
                    }); 
                }

                return tracks_sv;

            }, {"tracks_fp"});

    auto df_sv = add_branch_from_Track_t(df_reco, "tracks_sv", {
        {"x_sv",    &Track_t::x},
        {"y_sv",    &Track_t::y},
        {"dxdz_sv", &Track_t::dxdz},
        {"dydz_sv", &Track_t::dydz}
    }); 

    auto df_fp = add_branch_from_Track_t(df_sv,   "tracks_fp", {
        {"x_fp",    &Track_t::x},
        {"y_fp",    &Track_t::y},
        {"dxdz_fp", &Track_t::dxdz},
        {"dydz_fp", &Track_t::dydz}
    }); 

    //create both histograms
    auto hist_xy        
        = df_fp.Histo2D({"h_xy", "Sieve-plane projection;x_sv;y_sv", 200, -0.050, 0.055, 200, -0.035, 0.030}, "x_sv", "y_sv"); 
        
    auto hist_angles    
        = df_fp.Histo2D({"h_angles", "Sieve-plane projection;dx/dx_sv;dy/dz_sv", 200, -0.05, 0.06, 200, -0.04, 0.03}, "dxdz_sv", "dydz_sv"); 
    
    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile); 

    new TCanvas("c1", c_title); 
    hist_xy->DrawCopy("col2");

    new TCanvas("c2", c_title); 
    hist_angles->DrawCopy("col2"); 

    return 0; 
}