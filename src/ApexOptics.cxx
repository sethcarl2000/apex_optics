
#include "ApexOptics.h"
#include <iostream>
#include <stdio.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <ROOT/RVec.hxx> 
#include <ROOT/RResultPtr.hxx>
#include "RMatrix.h"

using namespace std; 
using namespace ROOT::VecOps; 

//__________________________________________________________________________________________________________________
NPoly* ApexOptics::Create_NPoly_fit( ROOT::RDF::RNode df, 
                                     const int poly_order, 
                                     const vector<string> &inputs, 
                                     const char* output )
{
    const char* const here = "ApexOptics::Create_NPoly_fit"; 
    
    const int nDoF = inputs.size();  

    if (nDoF < 1 ) {
        fprintf(stderr, "Error in <%s>: No inputs given.\n", here); 
        return nullptr; 
    }

    //check for the existence of all required columns    
    vector<string> column_names = df.GetColumnNames(); 
    vector<string> missing_columns{}; 
    
    //make a vector of all the columns we need to find to proceed
    vector<string> all_columns_needed(inputs); 
    all_columns_needed.push_back(string(output)); 

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
        return nullptr; 
    }


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
    NPoly* poly_template = new NPoly(nDoF, poly_order); 

    const int n_elems = poly_template->Get_nElems(); 

   
    //now that we have put all the inputs into a vector, we can continue. 
    auto df_output = df_vec.at(df_vec.size()-1)
    
        .Define("X_elems", [poly_template](const RVec<double> &inputs) 
        {
            return poly_template->Eval_noCoeff(inputs); 
        }, {"inputs"}); 

    
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
    RMatrix A(n_elems, n_elems); 
    vector<double> B; 

    //now, fill the matrix, and the 'b' values. 
    // When we request access to the RResultPtr objects created above (which is what the B_ptr and A_ptr arrays are), 
    // this is the moment that the RDataFrame actually loops through all events in the dataframe.  
    //cout << "--filling matrix..." << flush; 
    for (int i=0; i<n_elems; i++) {
        
        B.push_back( *(B_ptr[i]) ); 
        
        for (int j=0; j<n_elems; j++) {
            A.get(i,j) = *(A_ptr[i][j]);
        }  
    }
    //cout << "done." << endl; 

    //now, we can actually solve the linear equation for the tensor coefficients. 
    //cout << "--Solving linear system(s)..." << flush; 
    
    //Solve the linear system of equations to get our best-fit coefficients
    auto coeffs = A.Solve( B ); 

    auto poly = new NPoly(nDoF); 

    for (int i=0; i<n_elems; i++) { 
        
        //we're gonna make a copy of each element in the 'template' polynomial
        NPoly::NPolyElem *elem = poly_template->Get_elem(i); 

        //now, we just add our coefficient we just computed to it
        poly->Add_element( elem->powers, coeffs.at(i) ); 
    }

    //cout << "done." << endl; 

    delete poly_template; 
    return poly; 
}
//__________________________________________________________________________________________________________________
int ApexOptics::Create_dbfile_from_polymap(bool is_RHRS, string path_outfile, map<string, NPoly*> polymap) 
{            
    const char* const here = "ApexOptics::Create_dbfile_from_polymap"; 

    fstream outfile(path_outfile, ios::out | ios::trunc); 

    if (!outfile.is_open()) {
        fprintf(stderr, "Error in <%s>: Unable to open output file: '%s'", here, path_outfile.data()); 
        return 1; 
    }

    printf("--writing output file '%s'...", path_outfile.data()); cout << flush; 
    //now, we make the output file. 

    //write the poly DoF
    const int nDoF = polymap.begin()->second->Get_nDoF(); 
    outfile << "poly-DoF " << nDoF << endl; 
    //write the arm
    outfile << "is-RHRS " << (is_RHRS ? "1" : "0") << endl; 
    
    //assumes dbfile was opened successfully! 
    for (auto it = polymap.begin(); it != polymap.end(); it++) {
    
        //get the name of the polynomial
        const char* poly_name   = it->first.data(); 
        NPoly      *poly        = it->second; 

        char buffer[25]; 
        //now, we will print all the elements.
        for (int i=0; i<poly->Get_nElems(); i++) {
            outfile << poly_name; 
            
            const NPoly::NPolyElem* elem = poly->Get_elem(i);            
            
            for (int pow : elem->powers) {
                sprintf(buffer, " %3i", pow);
                outfile << buffer;  
            }

            //the '%+.9e' format produces scientific-notation floating-point output with 10 sig figures. 
            sprintf(buffer, "   %+.9e", elem->coeff); 
            outfile << buffer << endl; 
        }
    }

    outfile.close();
    cout << "done." << endl; 
    return 0; 
}
//__________________________________________________________________________________________________________________
//this will be a temporary funct., that I may absorb into NPoly.cxx. 
int ApexOptics::Parse_NPoly_from_file(const char* path_dbfile, const char* poly_name, NPoly *poly) 
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
    
    const char* const here = "ApexOptics::Parse_NPoly_from_file"; 

    //open the db file
    ifstream dbfile(path_dbfile); 
    
    if (!dbfile.is_open()) {
        fprintf(stderr, "Error in <%s>: Unable to open db file '%s'\n", here, path_dbfile); 
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
        fprintf(stderr, "Error in <%s>: Missing 'poly-DoF [n]' header at top of dbfile '%s'\n", here, path_dbfile); 
        return -1; 
    }

    int poly_DoF; 
    iss_init >> poly_DoF; 


    //this file can't contain any elements of this polynomial; it has the wrong DoF. 
    if (poly_DoF != poly->Get_nDoF()) return 0;  
    

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

    return poly->Get_nElems() - start_nElems;
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
