#include "TROOT.h"
#include <iostream>
#include <chrono> 
#include <vector> 
#include <thread>
#include <ROOT/RResultPtr.hxx> 

using namespace std; 
using namespace ROOT::VecOps; 

struct TrainingData_t {
    ROOT::RVec<double> inputs, outputs; 
}; 

double rv_mag2(const ROOT::RVec<double>& v) {
    double ret=0.; 
    for (const double& xx : v * v) ret += xx; 
    return ret; 
}

//__________________________________________________________________________________________________________________
//Process events from index 'start' to index 'end' in the 'data_vec' vector. 
void Process_event_range(   ROOT::RVec<ROOT::RVec<double>>& dW, 
                            const MultiLayerPerceptron* mlp, 
                            const std::vector<TrainingData_t>& data_vec, 
                            size_t start, 
                            size_t end  )
{
    for (size_t i=start; i<end; i++) { 
        
        const auto& data = data_vec.at(i); 
        
        RVec<double> Z_err = mlp->Eval(data.inputs) - data.outputs; 

        //this is the gradient of each output coordinate (i), w/r/t each weight in the network. 
        auto weight_gradient = mlp->Weight_gradient(data.inputs); 

        //now, compute the gradient of the **loss function** w/r/t each weight: 
        for (int l=0; l<mlp->Get_n_layers()-1; l++) {                   // (l) - index of network layer
            
            for (int j=0; j<mlp->Get_layer_size(l+1); j++) {        // (j) - index of current layer output
                for (int k=0; k<mlp->Get_layer_size(l)+1; k++) {    // (k) - index of current layer input

                    for (int i=0; i<mlp->Get_DoF_out(); i++) {                  // (i) - index of output layer
                
                        dW.at(l).at( j*(mlp->Get_layer_size(l)+1) + k ) += Z_err[i] * weight_gradient.get(i,l,j,k);   
                    }
                }
            }
        }

    }//for (size_t i=start; i<end; i++)
}
//__________________________________________________________________________________________________________________

//will store the minimum, maximum, and mean of each branch
struct BranchLimits_t {
    BranchLimits_t(string _name, 
                    ROOT::RDF::RResultPtr<double> _min_p, 
                    ROOT::RDF::RResultPtr<double> _max_p, 
                    ROOT::RDF::RResultPtr<double> _mean_p) 
        : name(_name), min_ptr(_min_p), max_ptr(_max_p), mean_ptr(_mean_p) {}; 
    

    string name; 
    ROOT::RDF::RResultPtr<double> min_ptr, max_ptr, mean_ptr; 
    double min, max, mean; 
    void Compute_limits() { min=*min_ptr; max=*max_ptr; mean=*mean_ptr; }
}; 

#define TOY_MLP_DATA true

int train_new_mlp(  const int n_grad_iterations = 10,
                    const char* path_infile = "",
                    const char* path_outfile = "",
                    const char* tree_name="tracks_fp" )
{
    const char* const here = "train_new_mlp"; 

    vector<string> branches_input   = {
        "x_inp", 
        "y_inp", 
        "z_inp"
    };

    vector<string> branches_output  = {
        "x_out",
        "y_out",
        "z_out"
    }; 

    //this is the structure of the hidden layers. the eventual network will append an input layer and output layer on either side
    // of this set of hidden layers. 
    //
    // For example: for a network with the structure (4 inputs) => 6 => 6 => (3 outputs), this vector should be: 
    //  RVec<int> mlp_structure{6,6}; 
    //
    RVec<int> mlp_structure{3}; 

    //add the input / output layers to this network. 
    mlp_structure.insert( mlp_structure.begin(), branches_input.size() ); 
    mlp_structure.push_back( branches_output.size() ); 

#if TOY_MLP_DATA
    const int n_events_train = 1.5e5;                 
    //___________________________________________________________________________________________________________________________
    //first, we create a 'dummy' mlp, and use it to create data: 
    auto mlp_target = new MultiLayerPerceptron(mlp_structure); 

    const int DoF_in  = mlp_target->Get_DoF_in(); 
    const int DoF_out = mlp_target->Get_DoF_out(); 

    //initialize with random, gaussian weights
    for (int l=0; l<mlp_target->Get_n_layers()-1; l++) {
        for (double& weight : mlp_target->Get_layer(l)) weight = gRandom->Gaus(); 
    }

    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df_create(n_events_train);

    auto df_output = df_create

        .Define("X", [DoF_in]()
        {   
            RVec<double> X; 
            X.reserve(DoF_in); 
            for (int i=0; i<DoF_in; i++) X.push_back( -1. + 2.*gRandom->Rndm() ); 
            return X; 
        }, {})

        .Define("x_inp", [](const RVec<double>& X){ return X[0]; }, {"X"})
        .Define("y_inp", [](const RVec<double>& X){ return X[1]; }, {"X"})
        .Define("z_inp", [](const RVec<double>& X){ return X[2]; }, {"X"})

        .Define("Z", [mlp_target](const RVec<double>& X)
        {
            return mlp_target->Eval(X); 
        }, {"X"})

        .Define("x_out", [](const RVec<double>& Z){ return Z[0]; }, {"Z"})
        .Define("y_out", [](const RVec<double>& Z){ return Z[1]; }, {"Z"})
        .Define("z_out", [](const RVec<double>& Z){ return Z[2]; }, {"Z"});
    
    printf("making snapshot of %i 'fake' training events...", n_events_train); cout << flush; 
    df_output.Snapshot("training_data", "data/misc/mlp_train_data.root"); 
    cout << "done." << endl;
    //___________________________________________________________________________________________________________________________
#endif 

    //to create this data, we want to run in sequential mode. 
    //create the dataframe
    if (ROOT::IsImplicitMTEnabled()) {
        ROOT::DisableImplicitMT(); 
    }
    
#if TOY_MLP_DATA 
    ROOT::RDataFrame df("training_data", "data/misc/mlp_train_data.root"); 
#else 
    ROOT::RDataFrame df(tree_name, path_infile); 
#endif 
    
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
        fprintf(stderr, "Error in <%s>: Missing columns: ", here); 
        for (const string& col : missing_columns) {
            fprintf(stderr, "'%s' ", col.data()); 
        }
        cerr << "\n"; 
        return 1; 
    }

    
    vector<BranchLimits_t> lim_inputs, lim_outputs; 
#if !TOY_MLP_DATA
    const int n_events_train = *df.Count(); 
#endif
    //now that we know all the columns exist, we can continue; 
    vector<TrainingData_t> training_data; training_data.reserve(n_events_train); 

    vector<ROOT::RDF::RNode> input_nodes{ df 
        .Define("inputs",  [](){ return RVec<double>{}; }, {})
        .Define("outputs", [](){ return RVec<double>{}; }, {})
    }; 

    //add the input branches
    for (const auto& str : branches_input ) {

        auto new_node = input_nodes.back() 

            .Redefine("inputs", [](RVec<double>& V, double x) { V.push_back(x); return V; }, {"inputs", str.data()}); 

        input_nodes.push_back(new_node); 
    }
    //add the output branches
    for (const auto& str : branches_output ) {
        
        auto new_node = input_nodes.back() 

            .Redefine("outputs", [](RVec<double>& V, double x) { V.push_back(x); return V; }, {"outputs", str.data()}); 

        input_nodes.push_back(new_node); 
    }

    //get limits of input branches
    for (auto& str : branches_input) {
        auto node = input_nodes.back();  
        lim_inputs.push_back(BranchLimits_t( str, node.Min(str.data()), node.Max(str.data()), node.Mean(str.data()) )); 
    }
    //get limits of output branches
    for (auto& str : branches_output) {
        auto node = input_nodes.back();  
        lim_outputs.push_back(BranchLimits_t( str, node.Min(str.data()), node.Max(str.data()), node.Mean(str.data()) )); 
    }

    //now, compute the limits. 
    cout << "Computing min/max/mean of all input & output branches..." << flush; 
    for (auto& limit : lim_inputs ) limit.Compute_limits(); 
    for (auto& limit : lim_outputs) limit.Compute_limits(); 
    cout << "done." << endl; 


    auto n_events_run = input_nodes.back() 

        .Define("data", [&training_data, &lim_inputs, &lim_outputs](RVec<double>& inputs, RVec<double>& outputs)
        {   
            //normalize the inputs/outputs so that they all lie on the range x = [-1, +1]
            for (int i=0; i<inputs.size(); i++) {
                double& val = inputs[i]; 
                auto& limit = lim_inputs[i]; 

                //val = 2.*( val - limit.mean )/( limit.max - limit.min ); 
            }
            for (int i=0; i<outputs.size(); i++) {
                double& val = outputs[i]; 
                auto& limit = lim_outputs[i]; 

                //val = 2.*( val - limit.mean )/( limit.max - limit.min ); 
            }

            training_data.push_back({
                .inputs  = inputs,
                .outputs = outputs
            }); 
            return 1.; 
        }, {"inputs", "outputs"}).Sum("data"); //we Force RDataFrame to compute a sum, so that this last part executes
    

    auto start  = chrono::high_resolution_clock::now(); 
    
    cout << "number of events run: " << (int)*n_events_run << endl;  
        
    auto end    = chrono::high_resolution_clock::now(); 

    double time_elapsed = chrono::duration_cast<chrono::microseconds>( end - start ).count(); 
    
    printf("done.\ntime elapsed: %f seconds (%f us/event)\n", time_elapsed * 1e-6, time_elapsed/((double)n_events_train) ); 
    //___________________________________________________________________________________________________________________________
    

    //we want this mlp to eventually match our target mlp 
    auto mlp = new MultiLayerPerceptron(mlp_structure); 

    //initialize with random, gaussian weights
    for (int l=0; l<mlp->Get_n_layers()-1; l++) {
        for (double& weight : mlp->Get_layer(l)) weight = gRandom->Gaus(); 
    }


    printf("\n\n~~~~~~~~~~~~~~~~~ Differences: (after)");
    for (int l=0; l<mlp->Get_n_layers()-1; l++) {

        RVec<double> diff = mlp->Get_layer(l) - mlp_target->Get_layer(l); 
        printf("\nLayer %i => %i", l, l+1); 
        int i_elem=0;
        double error=0;  
        for (int j=0; j<mlp->Get_layer_size(l+1); j++) { 

            error += pow( diff.at(i_elem), 2 ); 
            printf("\n  -  % .4e --- ", diff.at(i_elem++) );

            for (int k=0; k<mlp->Get_layer_size(l); k++) {
                error += pow( diff.at(i_elem), 2 ); 
                printf("% .4e ", diff.at(i_elem++));
            } 
        }
        printf( "   ---   rms of layer = %f\n ", sqrt(error/((double)diff.size())) ); 
        
    }
    cout << endl; 
            
    //the extent to which 
    double eta          = 0.200; 

    //the fraction of the 'existing' gradient which stays behind at the last step
    double momentum     = 0.850; 
 


    printf("~~~~~~~~~~~~~~~~~ New mlp:"); 
    mlp->Print(); 

    double x_epoch[n_grad_iterations];
    double y_error[n_grad_iterations];
    for (int i=0; i<n_grad_iterations; i++) { x_epoch[i]=i; y_error[i]=0.; }

    TGraph* graph = nullptr;  
    
    char canv_title[200]; 
    sprintf(canv_title, "Momentum: %.3f, Eta: %.3f", momentum, eta); 
    auto canvas = new TCanvas("c", canv_title); 

        
    //number of threads 
    const size_t n_threads             = thread::hardware_concurrency(); 
    const size_t n_events_per_thread   = n_events_train / n_threads; 
    const size_t remainder             = n_events_train % n_threads; 
    
    RVec<RVec<double>> dW(mlp->Get_n_layers(), {}); 
    
    //zero out the weight-update vector
    for (int l=0; l<mlp->Get_n_layers()-1; l++)   
        dW[l] = RVec<double>( (mlp->Get_layer_size(l)+1) * mlp->Get_layer_size(l+1), 0. ); 
    

    printf("Running in parallel with %i threads...\n\n", (int)n_threads); 

    //gradient descent trials. We will try to 'match' the inputs of the  
    for (int i=0; i<n_grad_iterations; i++) {

        //loop over all training data 

        //create a dW_partial vector for each thread
        RVec<RVec<double>> dW_partial[n_threads]; for (int t=0; t<n_threads; t++) dW_partial[t] = RVec<RVec<double>>(mlp->Get_n_layers(), {}); 
        for (size_t t=0; t<n_threads; t++) {
            for (int l=0; l<mlp->Get_n_layers()-1; l++) { 
                dW_partial[t][l] = RVec<double>( (mlp->Get_layer_size(l)+1) * mlp->Get_layer_size(l+1), 0. ); 
            }
        }

        //create & launch each thread (thank you to claude for teaching me how to use std::thread!)
        vector<thread> threads; threads.reserve(n_threads); 

        size_t start = 0; 
        for (size_t t=0; t<n_threads; t++) {
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0); 
            
            threads.emplace_back([&training_data, start, end, mlp, &dW_partial, t]{
                Process_event_range(dW_partial[t], mlp, training_data, start, end);  
            });
            start = end; 
        }

        //synchronize all the threads once they're done
        for (auto& thread : threads) thread.join();

        for (size_t t=0; t<n_threads; t++) dW += dW_partial[t] * eta; 


        //update the weights
        for (int l=0; l<mlp->Get_n_layers()-1; l++) { mlp->Get_layer(l) += - dW[l] / ((double)n_events_train); }

        //apply momentum
        for (auto& layer : dW) layer *= momentum; 

        //compute error with new weights
        double error(0.); 
        for (const TrainingData_t& data : training_data) error += rv_mag2( data.outputs - mlp->Eval(data.inputs) ); 
        error *= 1./((double)n_events_train * DoF_out); 
        
        error = sqrt(error);

        //print information about this epoch to stdout
        printf("\r -- epoch %i, error: %+.4e", i, error); cout << flush;  
        
        //update the graph & redraw
        if (i==0) { 
            for (int i=0; i<n_grad_iterations; i++) y_error[i] = log(error)/log(10.);
            graph = new TGraph(n_grad_iterations, x_epoch, y_error);
        } else    {
            graph->SetPointY(i, log(error)/log(10.)); 
        } 
        char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.2e;Epoch;log_{10}(RMS)", i+1, n_grad_iterations, error ); 
        graph->SetTitle(g_title);
        
        graph ->Draw(); 
        canvas->Modified(); 
        canvas->Update(); 
    }   

    /*
    printf("~~~~~~~~~~~~~~~~~ New mlp:"); 
    mlp->Print(); 
    printf("~~~~~~~~~~~~~~~~~ Old mlp:"); 
    mlp_target->Print(); 
    */

    //error between target / training MLPs 

    printf("\n\n~~~~~~~~~~~~~~~~~ Differences: (after)");
    for (int l=0; l<mlp->Get_n_layers()-1; l++) {

        RVec<double> diff = mlp->Get_layer(l) - mlp_target->Get_layer(l); 
        printf("\nLayer %i => %i", l, l+1); 
        int i_elem=0;
        double error=0;  
        for (int j=0; j<mlp->Get_layer_size(l+1); j++) { 

            error += pow( diff.at(i_elem), 2 ); 
            printf("\n  -  % .4e --- ", diff.at(i_elem++) );

            for (int k=0; k<mlp->Get_layer_size(l); k++) {
                error += pow( diff.at(i_elem), 2 ); 
                printf("% .4e ", diff.at(i_elem++));
            } 
        }
        printf( "   ---   rms of layer = %f\n ", sqrt(error/((double)diff.size())) ); 
        
    }
    cout << endl; 

    return 0; 
}