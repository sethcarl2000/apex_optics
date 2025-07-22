#include <memory>
#include <string> 
#include <map>
#include <fstream>
#include <iomanip> 
#include "TFile.h"
#include "TVector3.h"
#include "TParameter.h"


using namespace std; 

//this is a helper function which will take the names of of the 'X_elems_*' branches as input, along with the target 
//output branches, and return a std::map<string, NPoly*>, where the key of each element is the name of the polynomial. 
//______________________________________________________________________________________________________________________________
map<string,NPoly*> find_bestfit_poly_coeffs( ROOT::RDF::RNode df, 
                                             NPoly *poly_template, 
                                             const char* X_elems_name, 
                                             vector<string> outputs )
{
    const char* const here = "find_bestfit_poly_coeffs";

    const int nDoF    = poly_template->Get_nDoF(); 
    const int n_elems = poly_template->Get_nElems();

    const int n_outputs = outputs.size(); 

    ROOT::RDF::RResultPtr<double> A_elems[n_elems][n_elems]; 

    for (int i=0; i<n_elems; i++) {
        for (int j=0; j<n_elems; j++) {
            A_elems[i][j] = df
                .Define("val", [i,j](const ROOT::RVec<double> &X){ return X[i]*X[j]; }, {X_elems_name}).Sum("val"); 
        }
    }

    //these are the 'outputs'; what the polynomial maps onto. 
    //the 'outside' vector will have .size()=num_of_output_branches (output_branches.size()), and the 'inside' vector  
    // will have .size()=num_of_input_branches
    const int n_out_branches = outputs.size(); 
    ROOT::RDF::RResultPtr<double> B_ptr[n_elems][n_outputs];  
    
    int i_out=0; 
    for (const string& str : outputs) {

        //'book' the calculations for each element
        for (int i=0; i<n_elems; i++) {
            //this RResultPtr<double> will be the sum of the branch 'str' (see the construction of the output_branches vector above)
            //we call '.data()' method to get a 'char*' that the string is wrapping
            B_ptr[i][i_out] = df
                .Define("val", [i](const ROOT::RVec<double> &X, double y){ return X[i] * y; }, {X_elems_name, str.data()}).Sum("val"); 
        }
        i_out++; 
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

    //now, create the std::map which will store our output polynomials, indexed by their names. 
    map<string, NPoly*> poly_map;
    
    for (auto it = poly_coeffs.begin(); it != poly_coeffs.end(); it++) {

        //get polynomial name and vector of coefficients
        string         poly_name = it->first; 
        vector<double> coeff_vec = it->second; 

        //create a new NPoly, and add it to our output map  
        NPoly *poly = new NPoly(nDoF); 
        poly_map[poly_name] = poly; 
            
        for (int i=0; i<poly_template->Get_nElems(); i++) { 
            
            NPoly::NPolyElem *elem = poly_template->Get_elem(i); 
            poly->Add_element( elem->powers, coeff_vec.at(i) );
        }
    }

    return poly_map; 
}
//______________________________________________________________________________________________________________________________

//if you hand this helper function an 'fstream' object, and a sdt::map<string, NPoly*>, it will output all the elements of each
// polynomial into the dbfile. 
//______________________________________________________________________________________________________________________________
int create_dbfile_from_polymap(bool is_RHRS, string path_outfile, map<string, NPoly*> polymap) 
{            
    const char* const here = "create_dbfile_from_polymap"; 

    fstream outfile(path_outfile, ios::out | ios::trunc); 

    if (!outfile.is_open()) {
        Error(here, "unable to open output file: '%s'", path_outfile.data()); 
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
//______________________________________________________________________________________________________________________________
    

//creates db '.dat' files for two sepearate polynomials:
// the fp_q1 polynomials map from FOCAL PLANE coordinates to Q1 FRONT coordinates.
// the q1_sv polynomials map from Q1 FRONT coordinates to SIEVE coordinates. 
int fitpoints_mc_fp_q1_sv(  bool is_RHRS=false,
                            const int poly_fpq1_order=2,
                            const int poly_q1sv_order=2,
                            const char* path_infile="",
                            const char* stem_outfile="data/csv/db_mc",  
                            const char* tree_name="tracks_fp" ) 
{
    const char* const here = "fitpoints_mc_fp_q1_sv"; 

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

    
    const int nDoF_fpq1 = 4; //the 4 DoF are: x_fp, y_fp, dxdz_fp, dydz_fp

    const int nDoF_q1sv = 5; //the 5 DoF are: x_q1, y_q1, dxdz_q1, dydz_q1, dpp_q1

    const int hrs_momentum = 1104.0; 

    auto poly_fpq1 = new NPoly(nDoF_fpq1, poly_fpq1_order); 
    auto poly_q1sv = new NPoly(nDoF_q1sv, poly_q1sv_order); 

    auto df_output = df

        .Define("X_elems_fpq1", [poly_fpq1](double x, double y, double dxdz, double dydz)
        {
            return poly_fpq1->Eval_noCoeff({
                x,              //x_fp 
                y,              //y_fp
                dxdz - x/6.,    //dxdz_fp
                dydz            //dydz_fp
            }); 

        }, {"x_fp", "y_fp", "dxdz_fp", "dydz_fp"})

        .Define("x_q1",      [](TVector3 v){ return v.x(); },        {"position_Q1"})
        .Define("y_q1",      [](TVector3 v){ return v.y(); },        {"position_Q1"})
        .Define("dxdz_q1",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_Q1"})
        .Define("dydz_q1",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_Q1"})
        .Define("dpp_q1",    [hrs_momentum](TVector3 v){ return (v.Mag()-hrs_momentum)/hrs_momentum; }, {"momentum_Q1"})

        .Define("X_elems_q1sv", [poly_q1sv, hrs_momentum](TVector3 pos, TVector3 mom)
        {
            return poly_q1sv->Eval_noCoeff({
                pos.x(),                                //x_q1 
                pos.y(),                                //y_q1
                mom.x()/mom.z(),                        //dxdz_q1
                mom.y()/mom.z(),                        //dydz_q1
                (mom.Mag() - hrs_momentum)/hrs_momentum //dpp_q1
            });

        }, {"position_Q1", "momentum_Q1"})

        .Define("x_sv",      [](TVector3 v){ return v.x(); },        {"position_sieve"})
        .Define("y_sv",      [](TVector3 v){ return v.y(); },        {"position_sieve"})
        .Define("dxdz_sv",   [](TVector3 v){ return v.x()/v.z(); },  {"momentum_sieve"})
        .Define("dydz_sv",   [](TVector3 v){ return v.y()/v.z(); },  {"momentum_sieve"}); 
         
    //now, put these outputs into a vector (so we know to make a seperate polynomial for each of them). 
    vector<string> branches_fp = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 
    
    vector<string> branches_q1 = {
        "x_q1", 
        "y_q1",
        "dxdz_q1",
        "dydz_q1",
        "dpp_q1"
    };  

    vector<string> branches_sv = {
        "x_sv", 
        "y_sv",
        "dxdz_sv",
        "dydz_sv" 
    }; 


    //test the ApexOptics::Create_NPoly_fit() fcn
    map<string, unique_ptr<NPoly>> polymap_test; 

    auto poly_x_q1_test = ApexOptics::Create_NPoly_fit(df_output, poly_fpq1_order, branches_fp, "x_q1"); 


    //these 'result pointers' will let us see the result for each element of the least-squares fit matrix
    //first, find the best coffeficients for the polynomials.
    cout << "Creating polynomials for fp => q1..." << endl; 
    map<string, NPoly*> polymap_fpq1 = find_bestfit_poly_coeffs(df_output, poly_fpq1, "X_elems_fpq1", branches_q1); 
    cout << "done." << endl;    

    //test the differences between the new and old methods
    for (int i=0; i<poly_x_q1_test->Get_nElems(); i++) {
        auto elem_new = poly_x_q1_test->Get_elem(i); 
        auto elem_old = polymap_fpq1["x_q1"]->Get_elem(i); 
        printf("elem %2i : pows(", i); 
        for (int j=0; j<poly_x_q1_test->Get_nDoF(); j++) {
            printf("%i ", elem_new->powers.at(j)); 
        }
        printf(") coeff: new/old/diff: %+.8e/%+.8e/%+.8e\n", elem_new->coeff,elem_old->coeff,elem_new->coeff - elem_old->coeff); 
    }
    cout << endl; 
    return 0; 

    cout << "Creating polynomials for q1 => sv..." << endl; 
    map<string, NPoly*> polymap_q1sv = find_bestfit_poly_coeffs(df_output, poly_q1sv, "X_elems_q1sv", branches_sv); 
    cout << "done." << endl;  


    //construct the full path to the outfile from the stem provided
    string path_outfile; 
    char buffer[25];
    
    //create the fp => q1 output file ___________________________________________
    path_outfile = string(stem_outfile); 

    //specify the arm to use 
    path_outfile += "_fp_q1"; 

    path_outfile += (is_RHRS ? "_R" : "_L"); 

    //specify the order of the polynomial
    sprintf(buffer, "_%iord", poly_fpq1_order); 
    path_outfile += string(buffer); 

    path_outfile += ".dat";

    //write all elements of the output file
    create_dbfile_from_polymap(is_RHRS, path_outfile, polymap_fpq1); 


    //create the q1 => sv output file ___________________________________________
    path_outfile = string(stem_outfile); 

    //specify the arm to use 
    path_outfile += "_q1_sv"; 

    path_outfile += (is_RHRS ? "_R" : "_L"); 

    //specify the order of the polynomial
    sprintf(buffer, "_%iord", poly_q1sv_order); 
    path_outfile += string(buffer); 

    path_outfile += ".dat";

    //write all elements of the output file
    create_dbfile_from_polymap(is_RHRS, path_outfile, polymap_q1sv); 

    //delete our template polynomial, and all model polynomials
    delete poly_fpq1; 
    for (auto it = polymap_fpq1.begin(); it != polymap_fpq1.end(); it++ ) delete it->second;

    delete poly_q1sv; 
    for (auto it = polymap_q1sv.begin(); it != polymap_q1sv.end(); it++ ) delete it->second;


    return 0;
}