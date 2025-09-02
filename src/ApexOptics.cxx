
#include "ApexOptics.h"
#include <iostream>
#include <stdio.h>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <algorithm> 
#include <thread> 
#include <stdexcept> 
#include <ROOT/RVec.hxx> 
#include <ROOT/RResultPtr.hxx>
#include "RMatrix.h"
#include <sstream> 
#include <stdexcept> 

using namespace std; 
using namespace ROOT::VecOps; 

//__________________________________________________________________________________________________________________
map<string, NPoly*> ApexOptics::Create_NPoly_fit( ROOT::RDF::RNode df, 
                                                  const int poly_order, 
                                                  const vector<string> &inputs, 
                                                  const vector<string> &outputs )
{
    const char* const here = "ApexOptics::Create_NPoly_fit"; 
    
    const int DoF_in  = inputs.size();  
    const int DoF_out = outputs.size(); 

    if (DoF_in < 1 ) {
        fprintf(stderr, "Error in <%s>: No inputs given.\n", here); 
        return {}; 
    }

    if (DoF_out < 1) {
        fprintf(stderr, "Error in <%s>: No outputs given.\n", here); 
        return {};
    }

    //check for the existence of all required columns    
    vector<string> column_names = df.GetColumnNames(); 
    vector<string> missing_columns{};  
    
    //make a vector of all the columns we need to find to proceed
    vector<string> all_columns_needed; 
    copy( inputs.begin(),  inputs.end(),  back_inserter(all_columns_needed) ); 
    copy( outputs.begin(), outputs.end(), back_inserter(all_columns_needed) ); 

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
        return {}; 
    }

    struct TrainingData_t {
        RVec<double> inputs, outputs; 
    }; 

    //now that we know all the necessary columns exist, we can proceed: 
    //first, define the vector of inputs. 
    // since this function can take a variable number of inputs, we need to 'tack' each input onto a single vector. 
    vector<ROOT::RDF::RNode> df_vec{ df
        .Define("inputs",  [DoF_in ](){ RVec<double> v; v.reserve(DoF_in);  return v; })
        .Define("outputs", [DoF_out](){ RVec<double> v; v.reserve(DoF_out); return v; })
    }; 

    for (const string& name : inputs) {   //add all the inputs
        //define the next node. run thru this loop to add all inputs into a single RVec<double>. 
        df_vec.push_back( df_vec.back() 
        
            .Redefine("inputs",  [](RVec<double> &vec, double _val) 
            {
                vec.push_back(_val);
                return vec; 
            }, {"inputs", name.c_str()})
        ); 
    }
    for (const string& name : outputs) {   //add all the outputs
        //define the next node. run thru this loop to add all inputs into a single RVec<double>. 
        df_vec.push_back( df_vec.back() 
        
            .Redefine("outputs",  [](RVec<double> &vec, double _val) 
            {
                vec.push_back(_val);
                return vec; 
            }, {"outputs", name.c_str()})
        ); 
    }
    
    //now, we get a vector of all the inputs/outputs
    vector<TrainingData_t> training_data = *df_vec.back()

        .Define("training_data", [](RVec<double>& inputs, RVec<double>& outputs)
        {
            return TrainingData_t{ 
                .inputs  = inputs, 
                .outputs = outputs 
            }; 
        }, {"inputs", "outputs"})

        .Take<TrainingData_t>("training_data"); 
    
    //now that we have our vector, we can go ahead and compute the matrix and 'B' vectors 
    //The polynomial we output will follow the same template as this one:
    NPoly* poly_template = new NPoly(DoF_in, poly_order); 

    //number of threads 
    const size_t n_events_train        = training_data.size(); 
    const size_t n_threads             = thread::hardware_concurrency(); 
    const size_t n_events_per_thread   = n_events_train / n_threads; 
    const size_t remainder             = n_events_train % n_threads; 


    const int n_elems = poly_template->Get_nElems(); 
    
    //initialize the 'partial' matrix, one for each thread
    RMatrix A_partial[n_threads]; 
    RVec<RVec<double>> B_partial[n_threads]; 
    
    for (int i=0; i<n_threads; i++) {
        A_partial[i] = RMatrix(n_elems, n_elems, 0.); 
        for (int j=0; j<DoF_out; j++) B_partial[i].emplace_back(n_elems, 0.);  
    }

    auto Process_event_range = [&training_data, poly_template, DoF_out, n_elems]
        (RMatrix& A, RVec<RVec<double>>& B, size_t start, size_t end) 
    {   
        //check to make sure we were not passed an out-of-range vector
        if (end > training_data.size()) {
            ostringstream oss; 
            throw logic_error("in <ApexOptics::Create_NPoly_fit::Process_event_range>: attempted to access element"
                   " in 'training-data' vector which is out-of-range."); 
            return; 
        }

        //construct the A-matrix, and B-vector
        for (size_t i=start; i<end; i++) { 
            const auto& data = training_data[i]; 

            auto X_elems = poly_template->Eval_noCoeff(data.inputs); 

            for (int j=0; j<n_elems; j++) {
                
                for (int i_out=0; i_out<DoF_out; i_out++) B.at(i_out).at(j) += data.outputs.at(i_out) * X_elems[j]; 

                for (int k=0; k<n_elems; k++) A.at(j,k) += X_elems[j] * X_elems[k]; 
            }
        }

        return; 
    }; 

    
    //create & launch each thread (thank you to claude for teaching me how to use std::thread!)
    vector<thread> threads; threads.reserve(n_threads); 

    size_t start = 0; 
    for (size_t t=0; t<n_threads; t++) {

        size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0); 
        
        threads.emplace_back([&Process_event_range, &A_partial, &B_partial, start, end, t]
        {
            Process_event_range(A_partial[t], B_partial[t], start, end); 
        });

        start = end; 
    }
    
    //synchronize all the threads once they're done
    for (auto& thread : threads) thread.join();

    //combine all partial results in one sum
    RMatrix& A = A_partial[0]; 
    RVec<RVec<double>>& B = B_partial[0]; 

    for (int i=1; i<n_threads; i++) {
        A += A_partial[i]; 
        B += B_partial[i]; 
    }

    //Solve the linear system of equations to get our best-fit coefficients
    map<string, NPoly*> polymap; 

    int i_br=0; 
    for (const string& out_branch : outputs) {
     
        auto coeffs = A.Solve( B[i_br++] ); 

        auto poly = new NPoly(DoF_in); 

        for (int i=0; i<n_elems; i++) { 
            
            //we're gonna make a copy of each element in the 'template' polynomial
            const auto *elem = poly_template->Get_elem(i); 

            //now, we just add our coefficient we just computed to it
            poly->Add_element( elem->powers, coeffs.at(i) ); 
        }

        polymap[out_branch] = poly; 
    }
    //cout << "done." << endl; 

    delete poly_template; 
    return polymap; 
}
//__________________________________________________________________________________________________________________
int ApexOptics::Create_dbfile_from_mlp(const char* path_dbfile, const MultiLayerPerceptron* mlp)
{
    fstream dbfile(path_dbfile, ios::out | ios::trunc);

    const char* here = "ApexOptics::Create_dbfile_from_mlp"; 

    if (!dbfile.is_open()) {
        ostringstream oss; 
        oss << "in <" << here << ">: Unable to open file" << path_dbfile; 
        throw invalid_argument(oss.str());  
        return 1; 
    }

    const int n_layers = mlp->Get_n_layers(); 

    dbfile << "n-layers " << mlp->Get_n_layers(); 
    dbfile << "\nstructure "; for (int l=0; l<n_layers; l++) dbfile << mlp->Get_layer_size(l) << " "; 
    
    for (int l=0; l<n_layers-1; l++) {
        
        dbfile << "\nlayer-weights " << l; 

        char buff[50]; 
        
        for (int j=0; j<mlp->Get_layer_size(l+1); j++) {

            dbfile << "\n";  
            
            for (int k=0; k<mlp->Get_layer_size(l)+1; k++) {    
                sprintf(buff, "%+.9e ", mlp->Get_weight(l, j, k));
                dbfile << buff; 
            }        
        }
        
    }
    dbfile.close();
    
    return 0; 
}
//__________________________________________________________________________________________________________________
int ApexOptics::Create_dbfile_from_polymap(bool is_RHRS, string path_outfile, map<string, NPoly*> polymap) 
{            
    const char* const here = "ApexOptics::Create_dbfile_from_polymap"; 

    fstream outfile(path_outfile, ios::out | ios::trunc); 

    if (!outfile.is_open()) {
        fprintf(stderr, "Error in <%s>: Unable to open output file: '%s'", here, path_outfile.data()); 
        return -1; 
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
MultiLayerPerceptron* ApexOptics::Parse_mlp_from_file(const char* path_dbfile)
{
    const char* const here = "ApexOptics::Parse_mlp_from_file"; 
    
    //Parse a MLP block-by-block. 
    //the first line which must start "n-layers", is the number of layers 
    ifstream dbfile(path_dbfile); 

    if (!dbfile.is_open()) {
        fprintf(stderr, "Error in <%s>: unable to open dbfile '%s'\n", here, path_dbfile); 
        return nullptr; 
    }
    
    istringstream iss; 
    string line, token; 

    getline(dbfile, line); iss = istringstream(line); 

    size_t n_layers;
    iss >> token >> n_layers; 
    if (token != "n-layers") {
        fprintf(stderr, "Error in <%s>: Missing 'n-layers [n]' header at top of dbfile '%s'\n", here, path_dbfile); 
        return nullptr; 
    }

    getline(dbfile, line); iss = istringstream(line); 

    iss >> token; 
    if (token != "structure") {
        fprintf(stderr, "Error in <%s>: Missing 'structure [n] [n] ...' header at top of dbfile '%s'\n", here, path_dbfile); 
        return nullptr; 
    }
    RVec<int> structure; 
    for (int i=0; i<n_layers; i++) { int layer_size; iss >> layer_size; structure.push_back(layer_size); }

    if (structure.size() != n_layers) {
        fprintf(stderr, "Error in <%s>: n-layers listed does not match structure in '%s'\n", here, path_dbfile); 
        return nullptr;
    }

    MultiLayerPerceptron* mlp = new MultiLayerPerceptron(structure); 

    //now, we're ready to parse the layers.
    for (size_t l=0; l<n_layers-1; l++) {

        getline(dbfile, line); //the header for this layer
            
        for (int j=0; j<structure[l+1]; j++) {

            getline(dbfile, line); iss = istringstream(line); 
            
            for (int k=0; k<structure[l]+1; k++) iss >> mlp->Weight(l, j, k); 
        }
    }

    return mlp; 
}

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
        ostringstream oss; 
        oss << "in <" << here << ">: Unable to open file '" << path_dbfile << "'"; 
        throw invalid_argument(oss.str()); 
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
NPolyArray ApexOptics::Parse_NPolyArray_from_file(const char* path_dbfile, const vector<string>& output_names, const int DoF) 
{
    const char* const here = "ApexOptics::Parse_NPolyArray_from_file";
 
    //parse all relevant polys from file
    vector<NPoly> poly_vec; 

    try {
        for (const string& output_branch : output_names) {

            NPoly poly(DoF); 
    
            ApexOptics::Parse_NPoly_from_file(path_dbfile, output_branch.data(), &poly);        

            poly_vec.push_back(poly); 
        }
    
    } catch (const std::exception& e) { //check if an exception has been caught when trying to parse the polynomial from file. 
        
        ostringstream oss;  
        oss << "in <" << here << ">: Exception caught trying to parse NPolyArray from file.\n what(): " << e.what(); 
        throw invalid_argument(oss.str());

        NPolyArray parr_error;
        parr_error.Set_status(NPolyArray::kError); 
        return parr_error;  
    }

    return NPolyArray(poly_vec); 
}  
//__________________________________________________________________________________________________________________
ApexOptics::Trajectory_t ApexOptics::HCS_to_SCS(const bool is_RHRS, const ApexOptics::Trajectory_t traj_hcs) 
{
    //converts from the 'Hall coordinate system' (HCS) to the 'Sieve coordinate system' (SCS).

    TVector3 dir( traj_hcs.dxdz, traj_hcs.dydz, 1. );

    TVector3 pos( traj_hcs.x, traj_hcs.y, 0. ); 

    pos.RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
    pos.RotateZ( TMath::Pi()/2. ); 
    
    dir.RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
    dir.RotateZ( TMath::Pi()/2. ); 
    
    //have the position be relative to the central sieve hole. 
    pos += -ApexOptics::Get_sieve_pos(is_RHRS);

    //this will be our new trajectory. 
    Trajectory_t traj_scs{}; 
    traj_scs.dpp = traj_hcs.dpp; 

    //compute the new slopes
    traj_scs.dxdz = dir.x() / dir.z(); 
    traj_scs.dydz = dir.y() / dir.z(); 

    //use these new slopes to project the track onto the z=0 plane in HCS 
    traj_scs.x    = pos.x()  -  traj_scs.dxdz * pos.z(); 
    traj_scs.y    = pos.y()  -  traj_scs.dydz * pos.z(); 

    return traj_scs; 
}
//__________________________________________________________________________________________________________________
ApexOptics::Trajectory_t ApexOptics::SCS_to_HCS(const bool is_RHRS, const ApexOptics::Trajectory_t traj_scs) 
{
    //converts from the 'Hall coordinate system' (HCS) to the 'Sieve coordinate system' (SCS).
    //direction (SCS)

    auto dir = TVector3( traj_scs.dxdz, traj_scs.dydz, 1. );

    auto pos = TVector3( traj_scs.x, traj_scs.y, 0. ) + ApexOptics::Get_sieve_pos(is_RHRS); 

    //rotate both the position and the direction
    dir.RotateZ( -TMath::Pi()/2. ); 
    dir.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

    pos.RotateZ( -TMath::Pi()/2. ); 
    pos.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 

    //this will be our new trajectory. 
    Trajectory_t traj_hcs{}; 

    traj_hcs.dpp = traj_scs.dpp; 
    
    //compute the new slopes
    traj_hcs.dxdz = dir.x() / dir.z(); 
    traj_hcs.dydz = dir.y() / dir.z(); 

    //use these new slopes to project the track onto the z=0 plane in HCS 
    traj_hcs.x    = pos.x()  -  traj_hcs.dxdz * pos.z(); 
    traj_hcs.y    = pos.y()  -  traj_hcs.dydz * pos.z(); 

    return traj_hcs; 
}
//__________________________________________________________________________________________________________________
//if X.dpp is NOT NaN, then it is defined; so add it to the vector. otherwise, don't add it. 
ROOT::RVec<double> ApexOptics::Trajectory_t_to_RVec(ApexOptics::Trajectory_t X) noexcept {
    if ( X.dpp == X.dpp ) {
        return ROOT::RVec<double>{ X.x, X.y, X.dxdz, X.dydz, X.dpp };
    } else {
        return ROOT::RVec<double>{ X.x, X.y, X.dxdz, X.dydz };
    }
}
//__________________________________________________________________________________________________________________
//quick and dirty (slow) way to convert an RVec<double> to a Trajectory_t struct
ApexOptics::Trajectory_t ApexOptics::RVec_to_Trajectory_t(const ROOT::RVec<double>& V) {
    switch (V.size()) {
        case 4 : return Trajectory_t{V[0], V[1], V[2], V[3]};           //dpp is not defined here (==NaN)
        case 5 : return Trajectory_t{V[0], V[1], V[2], V[3], V[4]};     //dpp IS defined here. 
        default : {
            ostringstream oss; 
            oss << "in <ApexOptics::RVec_to_Trajectory_t>: ROOT::RVec<double> container passed is invalid size (" 
                << V.size() << "), must be either 4 or 5."; 
            throw invalid_argument(oss.str());
            return Trajectory_t{}; 
        }
    }
}
//__________________________________________________________________________________________________________________
TVector3 ApexOptics::HCS_to_SCS(const bool is_RHRS, TVector3 pos) 
{
    //converts from the 'Hall coordinate system' (HCS) to the 'Sieve coordinate system' (SCS).

    pos.RotateY( -ApexOptics::Get_sieve_angle(is_RHRS) ); 
    pos.RotateZ( TMath::Pi()/2. ); 
    
    //have the position be relative to the central sieve hole. 
    pos += -ApexOptics::Get_sieve_pos(is_RHRS);

    return pos; 
}
//__________________________________________________________________________________________________________________
TVector3 ApexOptics::SCS_to_HCS(const bool is_RHRS, TVector3 pos) 
{
    //converts from the 'Hall coordinate system' (HCS) to the 'Sieve coordinate system' (SCS).
    pos += ApexOptics::Get_sieve_pos(is_RHRS);

    pos.RotateZ( -TMath::Pi()/2. ); 
    pos.RotateY( ApexOptics::Get_sieve_angle(is_RHRS) ); 
    
    return pos; 
}
//__________________________________________________________________________________________________________________
const std::vector<ApexOptics::OpticsTarget_t> ApexOptics::GetTargetList() 
{
    return std::vector<OpticsTarget_t>{
        { .name="none" }, 

        { .name="O1",  .z_hcs=-296.30e-3 },  //carbon foil targets (Z is known precisely)
        { .name="O2",  .z_hcs=-215.30e-3 },
        { .name="O3",  .z_hcs=-146.10e-3 },
        { .name="O4",  .z_hcs= -71.30e-3 },
        { .name="O5",  .z_hcs=  78.90e-3 },
        { .name="O6",  .z_hcs= 153.70e-3 },
        { .name="O7",  .z_hcs= 222.90e-3 },
        { .name="O8",  .z_hcs= 303.70e-3 },

        { .name="V1",  .x_hcs=-3.23e-3, .z_hcs=-296.30e-3 },  //vertical wire targets (X & Z are known precisely)
        { .name="V2",  .x_hcs=-0.72e-3, .z_hcs=-215.30e-3 },
        { .name="V3",  .x_hcs=+1.73e-3, .z_hcs=-146.10e-3 },

        { .name="H1",  .y_hcs=-5.60e-3, .z_hcs=-246.25e-3 },  //horiztonal wire targets (Y & Z are known precisely)
        { .name="H2",  .y_hcs=-0.48e-3, .z_hcs= -96.20e-3 },
        { .name="H3",  .y_hcs=+4.52e-3, .z_hcs=+103.75e-3 },
        { .name="H4",  .y_hcs=+9.53e-3, .z_hcs=+253.75e-3 }
    };
}
//__________________________________________________________________________________________________________________
const ApexOptics::OpticsTarget_t ApexOptics::GetTarget(string target_name) 
{
    const auto target_list = GetTargetList(); 

    auto it = std::find_if( target_list.begin(), target_list.end(), 
        [&target_name](const OpticsTarget_t& elem) {
            return elem.name == target_name; 
        }
    ); 

    //if target was found
    if (it != target_list.end()) return *it; 

    //if target was NOT found
    ostringstream oss; 
    oss << "in <ApexOptics::GetTarget(string)>: input arguement '"<< target_name <<"' is not a valid target name. See ApexOptics::GetTargetList() in ApexOptics.cxx"
            " for list of valid target names."; 

    throw invalid_argument(oss.str()); 
    
    return OpticsTarget_t{ .name="null" }; 
}
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________
