#include "TROOT.h"

using namespace std; 

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


//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int fitpoints_mc_sv_fp(  bool is_RHRS=false,
                         const int poly_order=2,
                         const char* path_infile="",
                         const char* path_dbfile="data/csv/db_prod_mc_sv_fp_L_2ord.dat",  
                         const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fitpoints_mc_sv_fp"; 

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

    ROOT::EnableImplicitMT(); 
    Info(here, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame
    ROOT::RDataFrame df(tree_name, path_infile); 


    const double hrs_momentum = 1104.0; 

    vector<string> branches_fp = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    const int DoF_sv = 5; 

    vector<NPoly> poly_vec; 

    for (const string& str : branches_fp) {

        NPoly poly(DoF_sv); 

        parse_poly_from_file(path_dbfile, str.data(), &poly); 

        poly_vec.push_back(poly); 
    }

    //now that we have created the polys, we can create the NPolyModel object
    NPolyModel pmod(poly_vec); 
    
}