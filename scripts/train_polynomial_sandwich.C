#include "TROOT.h"
#include <iostream>
#include <chrono> 
#include <vector> 
#include <thread>
#include <fstream>
#include <sstream>
#include <limits> 
#include <ROOT/RResultPtr.hxx> 

using namespace std; 
using namespace ROOT::VecOps; 

struct TrainingData_t {
    ROOT::RVec<double> inputs, outputs; 
}; 

double rv_mag2(const ROOT::RVec<double>& v) {
    double ret=0.; 
    for (const double& x : v) ret += x * x; 
    return ret; 
}

struct NPolyArrayChain {

    vector<NPolyArray> arrays; 

    //evaluate all arrays, starting with arrays[0], feeding the input of each one into the next. 
    RVec<double> Eval(const RVec<double>& X) const { 
        RVec<double> out{X}; 
        for (const NPolyArray& array : arrays) out = array.Eval(out);
        return out;  
    }

    RMatrix Jacobian(const RVec<double>& X) const { 
        RMatrix J = RMatrix::Square_identity(arrays[0].Get_DoF_in()); 
        for (const NPolyArray& array : arrays) {
            auto Ji = std::move(array.Jacobian(X)); 
            J = Ji * J; 
        };
        return J; 
    }
}; 

//__________________________________________________________________________________________________________________
//Compute the gradient w/r/t each weight for the mlp (with respect to the given range of elements). 
void Process_event_range(   ROOT::RVec<double>& dW, 
                            const NPolyArrayChain* chain, 
                            const std::vector<TrainingData_t>& data_vec,
                            size_t start, 
                            size_t end  )
{
    const char* const here = "Process_event_range"; 

    //Check to make sure that the 'end' parameter given is not out-of-range. 
    if (end > data_vec.size()) {
        ostringstream oss; 
        oss << "Error in <" << here << ">: Argument 'size_t end' is invalid (" << end << "), input 'data_vec' size is:" << data_vec.size(); 
        throw std::logic_error(oss.str()); 
        return; 
    }

    const auto& output_array = chain->arrays[2]; 
    const auto& input_array  = chain->arrays[0];     

    const size_t DoF_out = output_array.Get_DoF_out(); 
    const size_t DoF_in  = input_array.Get_DoF_out(); 

    for (size_t i=start; i<end; i++) { 
        
        const auto& data = data_vec[i]; 
        
        //get the input to the last layer
        auto X_mid  = chain->arrays[0].Eval( data.inputs ); 
        auto X_last = chain->arrays[1].Eval( X_mid ); 
        auto Z_err  = chain->arrays[2].Eval( X_last ) - data.outputs; 

        
        //get the product of the jacobian of the last two layers
        RMatrix J2 = std::move(chain->arrays[2].Jacobian( X_last )); 
        RMatrix J1 = std::move(chain->arrays[1].Jacobian( X_mid )); 

        RMatrix J_21 = std::move( J2 * J1 );

        //check matrix for NaN. this can happen because of the way the jacobian is computed. 
        // if there is a nan-element, then skip this point.
        bool has_nan=false;  
        for (double &x : J_21.Data()) if (x != x) { has_nan=true; break; }
        if (has_nan) continue;  

        int i_elem =0;  

        //this is the gradient of the function w/r/t each of the input weights
        for (int j=0; j<DoF_in; j++) {

            for(double weight : input_array.Get_poly(j)->Eval_noCoeff(data.inputs)) {

                for (int i_out=0; i_out<DoF_out; i_out++) dW[i_elem] += Z_err[i_out] * J_21.get(i_out, j) * weight; 
                i_elem++; 
            }
        }

        //this is the gradient of the function w/r/t each of the output weights
        for (int i_out=0; i_out<DoF_out; i_out++) 
            for (double weight : output_array.Get_poly(i_out)->Eval_noCoeff(X_last)) 
                dW[i_elem++] += weight * Z_err[i_out]; 
        
    }//for (size_t i=start; i<end; i++)
}
//__________________________________________________________________________________________________________________

//__________________________________________________________________________________________________________________
// Compute the error w/r/t each element in the given range of the 'data_vec' for the current mlp.
double Evaluate_error_range(const std::vector<TrainingData_t>& data_vec, 
                            const NPolyArrayChain *chain,
                            const size_t start, 
                            const size_t end )
{
    const char* const here = "Evaluate_error_range"; 

    //Check to make sure that the 'end' parameter given is not out-of-range. 
    if (end > data_vec.size()) {
        ostringstream oss; 
        oss << "Error in <" << here << ">: Argument 'size_t end' is invalid (" << end << "), input 'data_vec' size is:" << data_vec.size(); 
        throw std::logic_error(oss.str()); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }
    double error=0.; 

    for (size_t i=start; i<end; i++) {
        auto& data = data_vec[i]; 
        error += rv_mag2( data.outputs - chain->Eval(data.inputs) ); 
    }

    return error; 
}


#define RUN_WITHOUT_MULTITHREADDING false

//____________________________________________________________________________________________________________________________________
int train_polynomial_sandwich(  const int n_grad_iterations = 10,
                                const char* path_infile = "",
                                const char* path_dbfile = "",
                                const char* path_outfile = "",
                                const int input_layer_order = 1, 
                                const int output_layer_order = 1, 
                                const char* tree_name="tracks_fp" )
{
    const char* const here = "train_polynomial_sandwich"; 
    
    //if there are more training events than this in the 'path_infile' root file, then only use this many. 
    const int max_events_train = 1e7; 

    //put in a number here which is orders of magnitude larger than the maximum reasonable error you expect. 
    //if the error of any particular iteration exceeds this, then quit. 
    const double max_error = 1e3; 


    vector<string> branches_output = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    };

    vector<string> branches_input  = {
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"    
    }; 

    //this is the structure of the hidden layers. the eventual network will append an input layer and output layer on either side
    // of this set of hidden layers. 
    //
    // For example: for a network with the structure (4 inputs) => 6 => 6 => (3 outputs), this vector should be: 
    //  RVec<int> mlp_structure{6,6}; 
    //

    //create the dataframe
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
    
    ROOT::RDataFrame df(tree_name, path_infile); 
    
    //check for existence of branches: 
    vector<string> column_names = df.GetColumnNames(); 
    vector<string> missing_columns{};  
    
    //make a vector of all the columns we need to find to proceed
    vector<string> all_columns_needed; 
    copy( branches_input.begin(),  branches_input.end(),  back_inserter(all_columns_needed) ); 
    copy( branches_output.begin(), branches_output.end(), back_inserter(all_columns_needed) ); 

    //loop through all branches needed, check if any are missing.
    for ( const string& column_needed : all_columns_needed ) { 
        
        bool is_found(false);    
        
        for ( const string& column_found : column_names ) if (column_found == column_needed) is_found = true; 
        
        if (!is_found) missing_columns.push_back(column_needed); 
    }

    //some missing columns were found. 
    if (missing_columns.size()>0) {
        fprintf(stderr, "Error in <%s>: From file '%s', missing columns: ", here, path_infile); 
        for (const string& col : missing_columns) {
            fprintf(stderr, "\n  '%s'", col.data()); 
        }
        cerr << "\n"; 
        return 1; 
    }

    
    //if more training events than this exist, cap them at 'max_events_train' 
    const int n_events_train = min<int>( *df.Count(), max_events_train ); 

    printf("Using %i training events\n", n_events_train);

    //now that we know all the columns exist, we can continue; 
    
    vector<ROOT::RDF::RNode> input_nodes{ df 
        .Range(0, n_events_train)
        .Define("inputs",  [](){ return RVec<double>{}; }, {})
        .Define("outputs", [](){ return RVec<double>{}; }, {})
    }; 

    //add the input branches
    for (const auto& str : branches_input ) {

        auto new_node = input_nodes.back() 

            .Redefine("inputs",  [](RVec<double>& V, double x) { V.push_back(x); return V; }, {"inputs", str.data()}); 

        input_nodes.push_back(new_node); 
    }
    //add the output branches
    for (const auto& str : branches_output ) {
        
        auto new_node = input_nodes.back() 

            .Redefine("outputs", [](RVec<double>& V, double x) { V.push_back(x); return V; }, {"outputs", str.data()}); 

        input_nodes.push_back(new_node); 
    }
    
    cout << "Compiling vector of " << n_events_train << " training events..." << flush; 

    vector<TrainingData_t> training_data = *input_nodes.back() 

        .Define("data", [](RVec<double>& inputs, RVec<double>& outputs)
        {   
            return TrainingData_t{
                .inputs  = inputs,
                .outputs = outputs
            }; 
        }, {"inputs", "outputs"})
        
        .Take<TrainingData_t>("data"); //the 'Take' command makes RDataFrame compile a vector for us, of the given type. 
    
    cout << "done." << endl; 
    //___________________________________________________________________________________________________________________________
    NPolyArray poly_array = ApexOptics::Parse_NPolyArray_from_file(path_dbfile, branches_output, (int)branches_input.size()); 

    if (poly_array.Get_status() != NPolyArray::kGood) {
        Error(here, "NPolyArray from file '%s' did not parse successfully", path_dbfile); 
        return -1; 
    }


    //Returns an NPolyArray of the given DoF and oder, with only the 'linear' elements nonzero (in other words, this poly array is the identity matrix). 
    auto Create_Identity_NPolyArray = [](const int DoF, const int order)
    {   
        vector<NPoly> polys; 
        for (int i=0; i<DoF; i++) {

            NPoly poly(DoF, order); 
            //set all elements to zero 
            for (int j=0; j<poly.Get_nElems(); j++) poly.Get_elem(j)->coeff = 0.; 

            //Now, set all 'linear' elments to 1.
            RVec<int> powers(DoF, 0); 
            powers[i] = 1; 

            poly.Find_element(powers)->coeff = 1.; 

            polys.push_back(poly); 
        }

        return NPolyArray(polys); 
    };

    //create the 'NPolyArrayChain' which will be the ~actual~ polynomial we want, as well as the 'bread' polynomials on either side
    NPolyArrayChain *sandwich = new NPolyArrayChain({
        Create_Identity_NPolyArray(poly_array.Get_DoF_in(),  input_layer_order), //the input layer 'bread' 
        poly_array, 
        Create_Identity_NPolyArray(poly_array.Get_DoF_out(), output_layer_order) //the output layer 'bread' 
    }); 

    auto& array_input  = sandwich->arrays[0]; 
    auto& array_output = sandwich->arrays[2]; 

    const size_t n_elems_input  = array_input .Get_poly(0)->Get_nElems() * array_input .Get_DoF_out(); 
    const size_t n_elems_output = array_output.Get_poly(0)->Get_nElems() * array_output.Get_DoF_out(); 
    
    const int DoF_in  = poly_array.Get_DoF_in(); 
    const int DoF_out = poly_array.Get_DoF_out(); 
            
    ///////////////////////////////////////////////////////////////////////////////////////
    //  
    //  HYPERPARAMETERS - 
    //  
    //the extent to which 
    double eta          = 1e-4;

    //the fraction of the 'existing' gradient which stays behind at the last step
    double momentum     = 1.0000 - 4e-3;
    
    //the number of epochs between graph updates: 
    const int update_period = 50; 
    //  
    //  
    //  
    ////////////////////////////////////////////////////////////////////////////////////////
    printf("Hyperparameters: ~~~\n"); 
    printf(" - Eta:               %.9f\n", eta);
    printf(" - Momentum:          %.9f\n", momentum);
    printf(" - (Drag-parameter):  %.3e\n\n", (1. - momentum)/eta); 
    //mlp->Print(); 

    double x_epoch[n_grad_iterations];
    double y_error[n_grad_iterations];

    for (int i=0; i<n_grad_iterations; i++) { x_epoch[i]= i + 1; y_error[i]=0.; }

    TGraph* graph = nullptr;  
    
    char canv_title[200]; 
    sprintf(canv_title, "Momentum: %.6f, Eta: %.6f", momentum, eta); 
    auto canvas = new TCanvas("c", canv_title); 

    

    //number of threads 
    const size_t n_threads             = ( RUN_WITHOUT_MULTITHREADDING ? 1 : std::thread::hardware_concurrency() ); 
    const size_t n_events_per_thread   = n_events_train / n_threads; 
    const size_t remainder             = n_events_train % n_threads; 
    
    //use all threads in parallel to evaluate the error of the passed mlp
    //______________________________________________________________________________________________________________
    auto Evaluate_error = [ &training_data, 
                            n_threads,
                            n_events_train,
                            n_events_per_thread,
                            remainder, 
                            DoF_out](const NPolyArrayChain* chain)
    {   
        RVec<double> error_threads(n_threads, 0.); 
        
        size_t start = 0; 
        std::vector<thread> threads; threads.reserve(n_threads); 

        //launch all threads to evaluate the error
        for (size_t t=0; t<n_threads; t++) {
            //establish event range to run over
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0);
            
            threads.emplace_back([&training_data, &error_threads, start, end, chain, t]
            {
                error_threads[t] = Evaluate_error_range(training_data, chain, start, end); 
            }); 
            start = end; 
        }

        //wait for the threads to finish
        for (auto &thread : threads) thread.join(); 

        //add the sum from all threads
        double error=0.; 
        for (auto &err : error_threads) error += err;

        //return the final error
        return sqrt( error / ((double)n_events_train * DoF_out)); 
    };
    //______________________________________________________________________________________________________________

    //the error for each generation. 
    double error = Evaluate_error(sandwich); 
    double last_error = error; 

    //initialize the update-vector
    RVec<double> dW(n_elems_input + n_elems_output, 0.); 
    RVec<double> best_weights;
    
    for (int i=0; i<DoF_in; i++) 
        for (int j=0; j<array_input.Get_poly(i)->Get_nElems(); j++) 
            best_weights.push_back( array_input.Get_poly(i)->Get_elemCoeff(j) ); 

    for (int i=0; i<DoF_out; i++) 
        for (int j=0; j<array_output.Get_poly(i)->Get_nElems(); j++) 
            best_weights.push_back( array_output.Get_poly(i)->Get_elemCoeff(j) ); 

    
    printf("Running in parallel with %i threads...\n\n", (int)n_threads); 
    
    double best_error  = 1e30; 
    double start_error = error; 

    auto Save_weights = [&best_weights, &best_error, DoF_in, DoF_out, n_elems_input, n_elems_output]
        (const NPolyArrayChain* chain, double error) 
    {   
        const auto& array_input  = chain->arrays[0]; 
        const auto& array_output = chain->arrays[2]; 

        if (error < best_error) { 

            int i_elem=0; 

            for (int i=0; i<DoF_in; i++) 
                for (int j=0; j<array_input.Get_poly(i)->Get_nElems(); j++) 
                    best_weights[i_elem++] = array_input.Get_poly(i)->Get_elemCoeff(j); 

            for (int i=0; i<DoF_out; i++) 
                for (int j=0; j<array_output.Get_poly(i)->Get_nElems(); j++) 
                    best_weights[i_elem++] = array_output.Get_poly(i)->Get_elemCoeff(j); 

            best_error = error; 
            return true; 
        }
        return false; 
    };

    //gradient descent trials. We will try to 'match' the inputs of the  
    for (int i=0; i<n_grad_iterations; i++) {

        //create a dW_partial vector for each thread
        RVec<double> dW_partial[n_threads]; 
        for (size_t t=0; t<n_threads; t++) dW_partial[t] = RVec<double>(n_elems_input + n_elems_output, 0.); 

        //create & launch each thread (thank you to claude for teaching me how to use std::thread!)
        vector<thread> threads; threads.reserve(n_threads); 

        size_t start = 0; 
        for (size_t t=0; t<n_threads; t++) {
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0); 
            
            threads.emplace_back([&training_data, start, end, sandwich, &dW_partial, t]{
               
                Process_event_range(dW_partial[t], sandwich, training_data, start, end);  
            });
            start = end; 
        }

        //synchronize all the threads once they're done
        for (auto& thread : threads) thread.join();
        
        //apply momentum
        dW *= momentum; 

        //combine all partial results from each thread 
        for (size_t t=0; t<n_threads; t++) dW += dW_partial[t] * eta; 

        int i_elem=0; 
        for (int j=0; j<DoF_in; j++)                    //input layer
            for (int k=0; k<array_input.Get_poly(j)->Get_nElems(); k++) 
                array_input .Get_poly(j)->Get_elem(k)->coeff += -dW[i_elem++] / ((double)n_events_train); 

        for (int j=0; j<DoF_out; j++)                   //output layer
            for (int k=0; k<array_output.Get_poly(j)->Get_nElems(); k++) 
                array_output.Get_poly(j)->Get_elem(k)->coeff += -dW[i_elem++] / ((double)n_events_train); 
        
        //compute error with new weights
        error = Evaluate_error(sandwich);

        //print information about this epoch to stdout
        printf("\r -- epoch %i, error: % .4e (%+.4e); current 'best' error saved: % .4e  -- (progress %3.1f)", 
            i, 
            error, 
            error - last_error, 
            best_error,
            100.*((double)i+1)/((double)n_grad_iterations) ); cout << flush; 
        
        if (error > max_error) {
            printf("\nMaximum error exceeded (%.4e > %.4e). iteration loop terminated.\n", error, max_error); 
            break; 
        }

        if (error != error) {
            printf("\nNan Error reached. iteration loop terminated.\n"); 
            break; 
        }

        //update the graph & redraw
        if (i % update_period == 0) {

            Save_weights(sandwich, error); 
            
            if (i==0) {
                for (int i=0; i<n_grad_iterations; i++) y_error[i] = log(error)/log(10.);
                graph = new TGraph(n_grad_iterations, x_epoch, y_error);
            }    
        
            char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.4e;Epoch;log_{10}(RMS)", i+1, n_grad_iterations, error ); 
            graph->SetTitle(g_title);
            
            graph ->Draw(); 
            canvas->Modified(); 
            canvas->Update(); 
        }
        graph->SetPointY(i, log(error)/log(10.));
        last_error = error; 
    
    }// for (int i=0; i<n_grad_iterations; i++) {
    
    graph->SetPointY(n_grad_iterations-1, log(error)/log(10.)); 
        
    char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.4e;Epoch;log_{10}(RMS)", n_grad_iterations, n_grad_iterations, error ); 
    graph->SetTitle(g_title);
    
    graph ->Draw(); 
    canvas->Modified(); 
    canvas->Update(); 

    if (best_error > start_error) {
        Error(here, "Best error found in interation loop is worse than starting error. no new MLP saved.");
        return 1; 
    }

    Save_weights(sandwich, Evaluate_error(sandwich)); 

    cout << endl; 
    printf("Best error recorded: %.5e\n", best_error); 

    int i_elem=0; 
    for (int i=0; i<DoF_in; i++)                    //input layer
        for (int j=0; j<array_input.Get_poly(i)->Get_nElems(); j++) 
            array_input.Get_poly(i)->Get_elem(j)->coeff = best_weights.at(i_elem++); 

    for (int i=0; i<DoF_out; i++)                   //output layer
        for (int j=0; j<array_output.Get_poly(i)->Get_nElems(); j++) 
            array_output.Get_poly(i)->Get_elem(j)->coeff = best_weights.at(i_elem++); 


    cout << "output wedge: " << endl; 
    for (int i=0; i<DoF_out; i++) {
        for (int j=0; j<array_output.Get_poly(i)->Get_nElems()-1; j++) {
            printf(" % .4f", array_output.Get_poly(i)->Get_elemCoeff(j));
        }
        printf(" | % .4f\n", array_output.Get_poly(i)->Get_elemCoeff( array_output.Get_poly(i)->Get_nElems()-1) ); 
    } 

    cout << "input wedge: " << endl; 
    for (int i=0; i<DoF_in; i++) {
        for (int j=0; j<array_input.Get_poly(i)->Get_nElems()-1; j++) {
            printf(" % .4f", array_input.Get_poly(i)->Get_elemCoeff(j));
        }
        printf(" | % .4f\n", array_input.Get_poly(i)->Get_elemCoeff( array_input.Get_poly(i)->Get_nElems()-1) ); 
    } 

    cout << "nesting arrays..." << flush; 

    //now, actually 'bake' this polynomial, by integrating the padding polynomials we've put on either side. 
    auto nested_array = NPolyArray::Nest( array_output, NPolyArray::Nest( poly_array, array_input ) ); 
    
    cout << "done." << endl; 

    //now, save it in a file 
    map<string, NPoly*> polymap; 

    int i_pol=0; 
    for (const auto& str : branches_output) {

        polymap[str] = nested_array.Get_poly(i_pol++); 
    }

    printf("\nCreating output file '%s'\n", path_outfile);

    ApexOptics::Create_dbfile_from_polymap(false, path_outfile, polymap);

    return 0; 
}