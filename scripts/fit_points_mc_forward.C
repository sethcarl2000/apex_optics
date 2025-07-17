#include <memory>
#include <string> 
#include "TFile.h"
#include "TVector3.h"


using namespace std; 

int fit_points_mc_forward(const char* path_infile="", const char* tree_name="tracks_fp") 
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

    //order of fitting polynomial to use
    const int poly_order = 2;     

    const int nDoF = 4; //the 4 DoF are: x_fp, y_fp, dxdz_fp, dydz_fp

    auto poly = new NPoly(nDoF,poly_order); 

    //get number of coefficients
    const int n_elems = poly->Get_nElems(); 


    //this will be used for the least-squares calculation to find the best coefficients
    RMatrix A_init(n_elems, n_elems, 0.); 

    //this function adds the data from each event to the least-square matrix
    auto aggregate_matrix = [n_elems](RMatrix A, ROOT::RVec<double> vec)
    {   
        for (int i=0; i<n_elems; i++)
            for (int j=0; j<n_elems; j++) A.get(i,j) += vec.at(i) * vec.at(j); 
        
        return A; 
    }; 

    //this function aggregates the matrices produces from each thread
    auto merge_matrix = [n_elems](vector<RMatrix> &A_vec) 
    {
        RMatrix A(n_elems, n_elems, 0.); 
        for (const auto& A_elem : A_vec) A = A + A_elem; 
        return A;  
    }; 

    double initval=0.;

    auto df_output = df

        .Define("X_elems", [poly](double x, double y, double dxdz, double dydz)
    {
        return poly->Eval_noCoeff({x,y,dxdz,dydz}); 
    }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})
        
        .Define("x_sieve",      [](TVector3 v){ return v.x(); },        {"position_sieve"})
        .Define("y_sieve",      [](TVector3 v){ return v.y(); },        {"position_sieve"})
        .Define("dxdz_sieve",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_sieve"})
        .Define("dydz_sieve",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_sieve"}); 
         
    //now, put these outputs into a vector (so we know to make a seperate polynomial for each of them). 
    vector<string> output_branches = {
        "x_sieve", 
        "y_sieve",
        "dxdz_sieve",
        "dydz_sieve"
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

    //now, we can actually solve the linear equation for the tensor coefficients 
    vector<double> poly_coeffs[n_out_branches]; 

    cout << "--Solving linear equation..." << flush; 
    for (int obr=0; obr<n_out_branches; obr++) {
        poly_coeffs[obr] = A.Solve( b_vec[obr] ); 
    }
    cout << "done." << endl; 
    
    delete poly; 
    return 0;
}