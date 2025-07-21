
#include "ApexOptics.h"
#include <iostream>
#include <stdio.h>
#include <cstdio>
#include <ROOT/RVec.hxx> 
#include <ROOT/RResultPtr.hxx>
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std; 
using namespace ROOT::VecOps; 

//__________________________________________________________________________________________________________________
NPoly ApexOptics::Create_NPoly_fit(ROOT::RDF::RNode df, const int poly_order, const vector<string> &inputs, const char* output)
{
    const int nDoF = inputs.size();  

    if (nDoF < 1 ) return NPoly(0); 

    //first, define the vector of inputs. 
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

   

    //now that we have put all the inputs into a vector, we can continue. 
    auto df_output = df_vec.at(df_vec.size()-1)
    
        .Define("X_elems", [poly](const RVec<double> &inputs) 
        {
            return poly->Eval_noCoeff(inputs); 
        } {"inputs"}); 


    NPoly *poly = new NPoly(nDoF, poly_order); 

    const int n_elems = poly->Get_nElems(); 

    ROOT::RDF::RResultPtr<double> A_ptr[n_elems][n_elems]; 

    for (int i=0; i<n_elems; i++) {
        for (int j=0; j<n_elems; j++) {
            A_ptr[i][j] = df
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
        B_ptr[i] = df
            .Define("val", [i](const ROOT::RVec<double> &X, double y){ return X[i] * y; }, {"X_elems", output}).Sum("val"); 
    }
    
    //this matrix will be used to find the best coefficients for each element 
    RMatrix A(n_elems, n_elems); 

    vector<double> b_vec[n_outputs]; 

    cout << "--filling matrix..." << flush; 
    //now, fill the matrix, and the 'b' values
    for (int i=0; i<n_elems; i++) {

        for (int obr=0; obr<n_outputs; obr++) {
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
    //use the RMatrix::Solve() method to solve the system of lin. equations corresponding to each 'output'. Store the answer in our map. 
    int obr=0; 
    for (const string& out_branch : outputs) {
        poly_coeffs[out_branch] = A.Solve( b_vec[obr++] ); 
    }
    cout << "done." << endl; 





    return NPoly(nDoF); 
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
