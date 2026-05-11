
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
        //Error(here, "Poly DoF in db file (%i) does not match DoF of passed Poly (%i)", 
        //    poly_DoF, poly->Get_nDoF()); 
        //return -1;
        //this polynomial cannot belong to this file, it has the wrong DoF! 
        return 0;  
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
    double x,y,dxdz,dydz,dpp; 
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
//if you want to use the 'fp-sv' polynomial models, then have path_dbfile_2="". otherwise, the program will assume that the *first* dbfile
// provided (path_dbfile_1) is the q1=>sv polynomials, and the *second* dbfile provided (path_dbfile_2) are the fp=>sv polynomials. 
int test_forward_model_tweak( const char* path_infile="data/replay/replay.4768.root",
                              const char* path_dbfile_1="data/csv/db_mc_V2_q1_sv_L_2ord.dat",  
                              const char* path_dbfile_2="data/csv/db_mc_V2_fp_q1_L_2ord.dat",
                              const char* path_outgif="histos/tweak_fwd_model.gif",
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

    bool use_fp_q1_sv_mode(false);
    if (string(path_dbfile_2)=="") {
        
        use_fp_q1_sv_mode=false; 
        Info(here, "Using the fp=>sv polynomial model."); 
    } else {

        use_fp_q1_sv_mode=true; 
        Info(here, "Using the fp=>q1=>sv polynomial model.");     
    }


    //try to get some meta data from the db file: 
    ifstream dbfile1(path_dbfile_1); 
    
    //check if the file can be opened
    if (!dbfile1.is_open()) {
        Error(here, "unable to open db file '%s'", path_dbfile_1); 
        return -1; 
    }

    //if we're going to use the second db file, check that one too. 
    ifstream dbfile2; 
    if (use_fp_q1_sv_mode) {

        dbfile2 = ifstream(path_dbfile_2);
        //check if the file can be opened
        if (!dbfile2.is_open()) {
            Error(here, "unable to open db file '%s'", path_dbfile_2); 
            return -1; 
        }  
        dbfile2.close(); 
    }
 
    //now, read the file
    string line;
    istringstream iss_init;
    string token; 
    int int_buffer; 
    
    //read the DoF of the poly
    getline(dbfile1, line); iss_init = istringstream(line);
    
    iss_init >> token; 
    if (token != "poly-DoF") {
        Error(here, "Missing 'poly-DoF [n]' header at top of dbfile '%s'", path_dbfile_1); 
        return -1; 
    }
    iss_init >> int_buffer; 
    int poly_DoF = int_buffer; 


    //read which arm to use
    getline(dbfile1, line); iss_init = istringstream(line);
    
    iss_init >> token; 
    if (token != "is-RHRS") {
        Error(here, "Missing 'is-RHRS [1/0]' header at top of dbfile '%s'", path_dbfile_1); 
        return -1; 
    }
    iss_init >> int_buffer; 
    const bool is_RHRS = (int_buffer==1); 

    dbfile1.close(); 

    //now, read the polynomials from the data file
    map<string, unique_ptr<NPoly>> pols; 

    //add the polynomials we want to parse
    // -- q1=>sv -or- fp=>sv
    poly_DoF = (use_fp_q1_sv_mode ? 5 : 4); 
    pols["x_sv"]    = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["y_sv"]    = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dxdz_sv"] = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dydz_sv"] = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    
    // -- fp=>q1
    poly_DoF = 4; 
    pols["x_q1"]    = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["y_q1"]    = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dxdz_q1"] = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dydz_q1"] = unique_ptr<NPoly>(new NPoly(poly_DoF)); 
    pols["dpp_q1"]  = unique_ptr<NPoly>(new NPoly(poly_DoF));
    
    //parse all elements for each of these
    cout << "parsing elements for all polynomials...\n" << flush; 
    
    for (auto it = pols.begin(); it != pols.end(); it++) { 
        
        string poly_name = it->first; 

        //try to parse all elements for this polynomial. search both db_files. 
        int elems_found = 0; 
        elems_found += parse_poly_from_file(path_dbfile_1, poly_name.data(), it->second.get()); 
        
        if (use_fp_q1_sv_mode) {
            elems_found += parse_poly_from_file(path_dbfile_2, poly_name.data(), it->second.get()); 
        } 

        if (elems_found < 0) {
            Error(here, "fatal error trying to parse polynomial '%s'.", poly_name.data()); 
            return 1; 
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

    const NPoly* pol_x_q1    = pols["x_q1"].get(); 
    const NPoly* pol_y_q1    = pols["y_q1"].get(); 
    const NPoly* pol_dxdz_q1 = pols["dxdz_q1"].get(); 
    const NPoly* pol_dydz_q1 = pols["dydz_q1"].get(); 
    const NPoly* pol_dpp_q1  = pols["dpp_q1"].get(); 

    
    const char* bn_x_tra    = is_RHRS ? "R_tr_tra_x"  : "L_tr_tra_x"; 
    const char* bn_y_tra    = is_RHRS ? "R_tr_tra_y"  : "L_tr_tra_y"; 
    const char* bn_dxdz_tra = is_RHRS ? "R_tr_tra_th" : "L_tr_tra_th"; 
    const char* bn_dydz_tra = is_RHRS ? "R_tr_tra_ph" : "L_tr_tra_ph"; 
    
    auto df_tracks_fp = df 
        
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
                        .dxdz = v_dxdz.at(i), 
                        .dydz = v_dydz.at(i), 
                        .dpp  = 0. 
                    });    
                }
                
                return tracks; 

            }, {bn_x_tra, bn_y_tra, bn_dxdz_tra, bn_dydz_tra}); 

    //now, create a 'frame' for each different tweak of the parameter we've chosen.
    const int n_frames = 20; 
    const char* tweak_coord = "dydz_fp"; 
    const double tweak_range = 15e-3; //+/- magnitude of perturbation range

#define MAKE_OUTGIF true
    
    //compute change in 'tweak' parameter
    const double dx = n_frames>0 ? 2.*tweak_range/((double)n_frames - 1) : 0.; 
    
    double perturb = -tweak_range; 
    
    //what fps do you want the gif to run at? 
    double gif_fps = 10.0; 

    //for some reason, the 'TCanvas::SaveAs' method takes centi-seconds (ms/10) as the argument for how long each frame is. 
    // for example, if you wanted your gif to me 10 fps, each frame would last 100 ms, so you would need to say: 
    // myCanvas->SaveAs("output.gif+10"). I don't know why they would do this. 
    //
    char saveas_outgif[300]; 
    sprintf(saveas_outgif, "%s+%i", path_outgif, (int)(100.*(1./gif_fps)) ); 
    
    char c_title[255]; 
    sprintf(c_title, "data:'%s', db:'%s'", path_infile, path_dbfile_1); 

    //set the color pallete
    gStyle->SetPalette(kSunset); 

    auto canv = new TCanvas("c1", c_title); 


    //now lets create all the frames
    printf("Printing %i frames to '%s' ...", n_frames, path_outgif); cout << endl; 
    for (int i=0; i<n_frames; i++) { 
        
        auto df_frame = df_tracks_fp 

            .Define("tracks_sv", 
                [perturb,
                use_fp_q1_sv_mode,
                pol_x,pol_y,pol_dxdz,pol_dydz, 
                pol_x_q1, pol_y_q1, pol_dxdz_q1, pol_dydz_q1, pol_dpp_q1]
                (const RVec<Track_t> &tracks_fp)
                {
                    //as before, create a vector of Track_t structs to store the sieve tracks. 
                    RVec<Track_t> tracks_sv; 
                    tracks_sv.reserve(tracks_fp.size());

                    //loop through all tracks; use the 'forward' polynomials to reconstruct 'sieve' data. 
                    for (const Track_t& trk_fp : tracks_fp) {

                        //add the 'tweak' to whichever coordinate you're trying to tweak

                        //convert the Track_t struct to an RVec<double>
                        RVec<double> X_fp = {
                            trk_fp.x,           //x_fp
                            trk_fp.y,           //y_fp
                            trk_fp.dxdz,        //dxdz_fp
                            trk_fp.dydz + perturb         //dydz_fp
                        }; 
                        
                        if (use_fp_q1_sv_mode) {
                            //frist evaluate fp=>q1, then q1=>sv
                            RVec<double> X_q1 = {
                                pol_x_q1->Eval(X_fp),     //x_q1
                                pol_y_q1->Eval(X_fp),     //y_q1
                                pol_dxdz_q1->Eval(X_fp),  //dxdz_q1
                                pol_dydz_q1->Eval(X_fp),  //dydz_q1
                                pol_dpp_q1->Eval(X_fp)    //dpp_q1
                            };

                            //use ouqqr polynomials to compute the sieve coordinates
                            tracks_sv.push_back({
                                .x    = pol_x->Eval(X_q1),      //x_sv
                                .y    = pol_y->Eval(X_q1),      //y_sv
                                .dxdz = pol_dxdz->Eval(X_q1),   //dxdz_sv
                                .dydz = pol_dydz->Eval(X_q1),   //dydz_sv
                                .dpp  = 0.                      //dpp_sv
                            });

                        } else {
                            //use ouqqr polynomials to compute the sieve coordinates
                            tracks_sv.push_back({
                                .x    = pol_x->Eval(X_fp),      //x_sv
                                .y    = pol_y->Eval(X_fp),      //y_sv
                                .dxdz = pol_dxdz->Eval(X_fp),   //dxdz_sv
                                .dydz = pol_dydz->Eval(X_fp),   //dydz_sv
                                .dpp  = 0.                      //dpp_sv
                            });

                        }//if (use_fp_q1_sv_mode)

                    }

                    return tracks_sv;

                }, {"tracks_fp"});


        auto df_sv = add_branch_from_Track_t(df_frame, "tracks_sv", {
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

        char histo_name[200]; 
        sprintf(histo_name, "Sieve-plane proj: %s %+.8f;x_sv;y_sv", tweak_coord, perturb); 

        //create both histograms
        auto hist_xy        
            = df_fp
            
            //correct for react-vertex position
            .Define("x_react_vtx_fix", [](RVec<double> x, TVector3 r){return x + r.x();}, {"x_sv", "react_vertex_TCS"})
            .Define("y_react_vtx_fix", [](RVec<double> y, TVector3 r){return y + r.y();}, {"y_sv", "react_vertex_TCS"})

            .Histo2D({"h_xy", histo_name, 200, -0.040, 0.045, 200, -0.045, 0.010}, 
                    "x_react_vtx_fix", "y_react_vtx_fix");

        hist_xy->DrawCopy("col2");
 
#if MAKE_OUTGIF
        canv->SaveAs(saveas_outgif); 
#endif 

        printf(" --- completed frame %-4i: %s %+.4e", i, tweak_coord, perturb); cout << endl;
        perturb += dx;  
    }

    
    return 0; 
}