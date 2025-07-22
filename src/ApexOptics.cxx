
#include "ApexOptics.h"
#include <iostream>
#include <stdio.h>
#include <cstdio>
#include <cmath>
#include <ROOT/RVec.hxx> 
#include <ROOT/RResultPtr.hxx>
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std; 
using namespace ROOT::VecOps; 

//__________________________________________________________________________________________________________________
unique_ptr<NPoly> ApexOptics::Create_NPoly_fit( ROOT::RDF::RNode df, 
                                                const int poly_order, 
                                                const vector<string> &inputs, 
                                                const char* output )
{
    const char* const here = "ApexOptics::Create_NPoly_fit"; 
    
    const int nDoF = inputs.size();  

    if (nDoF < 1 ) {
        fprintf(stderr, "Error in <%s>: No inputs given.\n", here); 
        return unique_ptr<NPoly>(nullptr); 
    }

    //check for the existence of all required columns    
    vector<string> column_names = df.GetColumnNames(); 
    vector<string> missing_columns{}; 
    
    vector<string> all_columns_needed = inputs; inputs.push_back(string(output)); 

    //loop through all branches needed, check if any are missing.
    for ( const string& column_needed : all_columns_needed ) { 
        
        bool is_found(false);    
        
        for ( const string& column_found : column_names ) if (column_found == column_needed) is_found = true; 
        
        if (!is_found) missing_columns.push_back(column_needed); 
    }

    //some missing columns were found. 
    if (missing_columns.size()>0) {
        fprintf(stderr, "Error in <%s>: Missing columns: ", here); 
        for (const string& col : missing_columns) {
            fprintf(stderr, "'%s' ", col.data()); 
        }
        cerr << "\n"; 
        return unique_ptr<NPoly>(nullptr); 
    }

    return unique_ptr<NPoly>(nullptr); 
    

    //now that we know all the necessary columns exist, we can proceed: 
    //first, define the vector of inputs. 
    // since this function can take a variable number of inputs, we need to 'tack' each input onto a single vector. 
    vector<ROOT::RDF::RNode> df_vec{df.Define("inputs", [](double _inp){return RVec<double>{};}, {inputs[0].data()})}; 

    int i_branch=0; 
    for (const string& input_name : inputs) {

        //define the next node. run thru this loop to add all inputs into a single RVec<double>. 
        auto new_df = df_vec.at(df_vec.size()-1) 
        
            .Redefine("inputs",  [](RVec<double> &input_vec, double _inp) 
            {
                input_vec.push_back(_inp);
                return input_vec; 
            }, {"inputs", input_name.data()}); 

        df_vec.push_back(new_df); 
    }

    //The polynomial we output will follow the same template as this one:
    auto poly_template = unique_ptr<NPoly>(new NPoly(nDoF, poly_order)); 

    const int n_elems = poly_template->Get_nElems(); 

   
    //now that we have put all the inputs into a vector, we can continue. 
    auto df_output = df_vec.at(df_vec.size()-1)
    
        .Define("X_elems", [poly_template](const RVec<double> &inputs) 
        {
            return poly_template->Eval_noCoeff(inputs); 
        } {"inputs"}); 


    
    ROOT::RDF::RResultPtr<double> A_ptr[n_elems][n_elems]; 

    for (int i=0; i<n_elems; i++) {
        for (int j=0; j<n_elems; j++) {
            A_ptr[i][j] = df_output
                .Define("val", [i,j](const ROOT::RVec<double> &X){ return X[i]*X[j]; }, {"X_elems"}).Sum("val"); 
        }
    }

    //these are the 'outputs'; what the polynomial maps onto. 
    //the 'outside' vector will have .size()=num_of_output_branches (output_branches.size()), and the 'inside' vector  
    // will have .size()=num_of_input_branches
    ROOT::RDF::RResultPtr<double> B_ptr[n_elems];  
    
    //'book' the calculations for each element
    for (int i=0; i<n_elems; i++) {

        //this RResultPtr<double> will be the sum of the branch 'str' (see the construction of the output_branches vector above)
        B_ptr[i] = df_output
            .Define("val", [i](const ROOT::RVec<double> &X, double y){ return X[i] * y; }, {"X_elems", output}).Sum("val"); 
    }
    
    //this matrix will be used to find the best coefficients for each element 
    TMatrixD A(n_elems, n_elems); 
    TVectorD B(n_elems); 

    //now, fill the matrix, and the 'b' values. 
    // When we request access to the RResultPtr objects created above (which is what the B_ptr and A_ptr arrays are), 
    // this is the moment that the RDataFrame actually loops through all events in the dataframe.  
    cout << "--filling matrix..." << flush; 
    for (int i=0; i<n_elems; i++) {
        
        B(i) =  *(B_ptr[i]); 
        
        for (int j=0; j<n_elems; j++) {
            A(i,j) = *(A_ptr[i][j]);
        }  
    }
    cout << "done." << endl; 

    //check to see if TMatrixD considers this matrix to be singular 
    if (fabs(A.Determinant()) < 1e-10) {
        fprintf(stderr, "Error in <%s>: A-matrix is singular (|det| < 1e-16)", here); 
        return unique_ptr<NPoly>(nullptr); 
    }

    //now, we can actually solve the linear equation for the tensor coefficients. 
    cout << "--Solving linear system(s)..." << flush; 
    
    //Solve the linear system of equations to get our best-fit coefficients
    auto coeffs = A.Invert() * B; 

    auto poly = unique_ptr<NPoly>(new NPoly(nDoF)); 

    for (int i=0; i<n_elems; i++) { 
        
        //we're gonna make a copy of each element in the 'template' polynomial
        NPoly::NPolyElem *elem = poly_template->Get_elem(i); 

        //now, we just add our coefficient we just computed to it
        poly->Add_element( elem->powers, coeffs(i) ); 
    }

    cout << "done." << endl; 

    return poly; 
}
//__________________________________________________________________________________________________________________


//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
