#include <memory>
#include <string> 
#include <map>
#include <fstream>
#include <iomanip> 
#include "TFile.h"
#include "TVector3.h"
#include "TParameter.h"


using namespace std; 

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int fitpoints_mc_fp_sv( bool is_RHRS=false,
                        const int poly_order=2,
                        const char* path_infile="",
                        const char* stem_outfile="data/csv/db_fwd",  
                        const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fit_points_mc_forward"; 

    auto infile = new TFile(path_infile, "READ");

    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(here, "libApexOptics could not be loaded."); 
        return 1; 
    }

    //check if we can load the apex optics lib
    if (gSystem->Load("libRMatrix") < 0) {
        Error(here, "libRMatrix could not be loaded."); 
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

    
    const int nDoF = 4; //the 4 DoF are: x_fp, y_fp, dxdz_fp, dydz_fp

    auto poly = new NPoly(nDoF,poly_order); 

    //get number of coefficients
    const int n_elems = poly->Get_nElems(); 


    //this will be used for the least-squares calculation to find the best coefficients
    RMatrix A_init(n_elems, n_elems, 0.); 


    auto df_output = df

        .Define("X_elems", [poly](double x, double y, double dxdz, double dydz)
    {
        return poly->Eval_noCoeff({x, y, dxdz - x/6., dydz}); 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})

        .Define("x_sv",      [](TVector3 v){ return v.x(); },        {"position_sieve"})
        .Define("y_sv",      [](TVector3 v){ return v.y(); },        {"position_sieve"})
        .Define("dxdz_sv",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_sieve"})
        .Define("dydz_sv",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_sieve"}); 
         
    //now, put these outputs into a vector (so we know to make a seperate polynomial for each of them). 
    vector<string> output_branches = {
        "x_sv", 
        "y_sv",
        "dxdz_sv",
        "dydz_sv"
    };  


    //these 'result pointers' will let us see the result for each element of the least-squares fit matrix
    ROOT::RDF::RResultPtr<double> A_elems[n_elems][n_elems]; 

    for (int i=0; i<n_elems; i++) {
        for (int j=0; j<n_elems; j++) {
            A_elems[i][j] = df_output.Define("a_val", [i,j](const ROOT::RVec<double> &X){ return X[i]*X[j]; }, {"X_elems"}).Sum("a_val"); 
        }
    }

    //these are the 'outputs'; what the polynomial maps onto. 
    //the 'outside' vector will have .size()=num_of_output_branches (output_branches.size()), and the 'inside' vector  
    // will have .size()=num_of_input_branches
    const int n_out_branches = output_branches.size(); 
    ROOT::RDF::RResultPtr<double> B_ptr[n_elems][n_out_branches];  
    
    int i_out=0; 
    for (const string& str : output_branches) {

        //'book' the calculations for each element
        for (int i=0; i<n_elems; i++) {
            //this RResultPtr<double> will be the sum of the branch 'str' (see the construction of the output_branches vector above)
            //we call '.data()' method to get a 'char*' that the string is wrapping
            B_ptr[i][i_out] = df_output
                .Define("val", [i](const ROOT::RVec<double> &X, double y){ return X[i] * y; }, {"X_elems", str.data()})
                .Sum("val"); 
        }
        i_out++; 
    }

    //this matrix will be used to find the best coefficients for each element 
    RMatrix A(n_elems, n_elems); 

    vector<double> b_vec[n_out_branches]; 

    cout << "--filling matrix..." << flush; 
    //now, fill the matrix, and the 'b' values
    for (int i=0; i<n_elems; i++) {

        for (int obr=0; obr<output_branches.size(); obr++) {
            b_vec[obr].push_back( *(B_ptr[i][obr]) ); 
        }

        for (int j=0; j<n_elems; j++) {
             A.get(i,j) = *(A_elems[i][j]);
        }  
    }
    
    cout << "done." << endl; 

    //now, we can actually solve the linear equation for the tensor coefficients. 
    //we put this in 'std::map' form, so that we can access any of the 'coefficient vectors' indexed by the name of the input branch. 
    map<string, vector<double>> poly_coeffs; 

    cout << "--Solving linear system(s)..." << flush; 
    int obr=0; 
    for (const string& out_branch : output_branches) {
        poly_coeffs[out_branch] = A.Solve( b_vec[obr++] ); 
    }
    cout << "done." << endl; 

    //construct the full path to the outfile from the stem provided
    string path_outfile(stem_outfile); 

    char buffer[25];
    //specify the arm to use 
    path_outfile += (is_RHRS ? "_R" : "_L"); 

    //specify the order of the polynomial
    sprintf(buffer, "_%iord", poly_order); 
    path_outfile += string(buffer); 

    path_outfile += ".dat";



    fstream outfile(path_outfile, ios::out | ios::trunc); 

    if (!outfile.is_open()) {
        Error(here, "unable to open output file: '%s'", path_outfile.data()); 
        return 1; 
    }

    printf("--writing output file '%s'...", path_outfile.data()); cout << flush; 
    //now, we make the output file. 
    string output_format = "%s "; for (int i=0; i<nDoF; i++) output_format += "%i "; 
    output_format += "%+.8e\n"; 

    //write the poly DoF
    outfile << "poly-DoF " << nDoF << endl; 
    //write the arm
    outfile << "is-RHRS " << (is_RHRS ? "1" : "0") << endl; 

    map<string, NPoly*> poly_models; 

    for (auto it = poly_coeffs.begin(); it != poly_coeffs.end(); it++) {
        
        //get the name of the polynomial
        const char* poly_name = it->first.data(); 
        
        //add this to our list of models
        NPoly *poly_model = new NPoly(nDoF); 
        poly_models[it->first] = poly_model; 

        const vector<double>& coeffs = it->second; 

        char buffer[25]; 
        //now, we will print all the elements.
        for (int i=0; i<n_elems; i++) {
            outfile << poly_name; 
            
            const NPoly::NPolyElem* elem = poly->Get_elem(i);            
            
            for (int pow : elem->powers) {
                sprintf(buffer, " %3i", pow);
                outfile << buffer;  
            }

            //the '%+.9e' format produces scientific-notation floating-point output with 10 sig figures. 
            sprintf(buffer, "   %+.9e", coeffs.at(i)); 
            outfile << buffer << endl; 

            //add this element to our model polynomial
            poly_model->Add_element(elem->powers, coeffs.at(i));
        }
    }   
    outfile.close(); 
    cout << "done." << endl; 

    
    //draw the reults of all models
    vector<ROOT::RDF::RNode> error_nodes{ df_output }; 
    for (auto it = poly_models.begin(); it != poly_models.end(); it++) {

        //get the polynomial and its name
        const char* poly_name = it->first.data(); 
        NPoly *poly           = it->second; 
        
        //name the new output branch 
        char br_output_name[50];
        sprintf(br_output_name, "error_%s", poly_name);

        auto new_node = error_nodes.at(error_nodes.size()-1) 

            .Define(br_output_name, [poly](double target, double x, double y, double dxdz, double dydz)
            {
                return (target - poly->Eval({x, y, dxdz - x/6., dydz}))*1e3; 

            }, {poly_name, "x_fp", "y_fp", "dxdz_fp", "dydz_fp"}); 

        error_nodes.push_back(new_node); 
    }

    auto df_error = error_nodes.at(error_nodes.size()-1); 

    char b_c_title[120]; 
    sprintf(b_c_title, "Errors of different coords: %s", path_outfile.data()); 
    auto c = new TCanvas("c", b_c_title, 1200, 800); 

    c->Divide(2,2, 0.005,0.005); 
    
    c->cd(1); 
    auto h_x = df_error.Histo1D({"h_x", "Error of x_sv;mm", 200, -10, 10}, "error_x_sv"); 
    h_x->DrawCopy(); 
    c->cd(2); 
    auto h_y = df_error.Histo1D({"h_y", "Error of y_sv;mm", 200, -10, 10}, "error_y_sv"); 
    h_y->DrawCopy(); 
    c->cd(3); 
    auto h_dxdz = df_error.Histo1D({"h_dxdz", "Error of dxdz_sv;mrad", 200, -2, 2}, "error_dxdz_sv"); 
    h_dxdz->DrawCopy(); 
    c->cd(4); 
    auto h_dydz = df_error.Histo1D({"h_dydz", "Error of dydz_sv;mrad", 200, -2, 2}, "error_dydz_sv"); 
    h_dydz->DrawCopy(); 


    //delete our template polynomial
    delete poly; 

    //delete our poly models
    for (auto it = poly_models.begin(); it != poly_models.end(); it++ ) delete it->second;

    return 0;
}