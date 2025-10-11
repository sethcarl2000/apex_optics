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
    for (const double& xx : v * v) ret += xx; 
    return ret; 
}

using WG = MultiLayerPerceptron::WeightGradient_t; 

//__________________________________________________________________________________________________________________
//Compute the gradient w/r/t each weight for the mlp (with respect to the given range of elements). 
void Compute_Jacobian_range(    RMatrix& J, 
                                const MultiLayerPerceptron* mlp, 
                                const std::vector<TrainingData_t>& data_vec, 
                                const size_t start, 
                                const size_t end,
                                const double dx =1e-10 )
{
    const char* const here = "Compute_Jacobian_range"; 
    
    const size_t n_weights = mlp->Get_n_weights(); 

    //check to make sure the Jacobian passed (J) is the right size
    if (J.GetNCols() != n_weights || J.GetNRows() != n_weights) {
        ostringstream oss; 
        oss << " in <" << here << ">: RMatrix 'J' passed is wrong size (" << J.GetNCols() << "x" << J.GetNRows() 
            << ") must be " << n_weights << "x" << n_weights << "."; 
        throw std::logic_error(oss.str()); 
        return;
    }

    //Check to make sure that the 'end' parameter given is not out-of-range. 
    if (end > data_vec.size()) {
        ostringstream oss; 
        oss << " in <" << here << ">: Argument 'size_t end' is invalid (" << end << "), input 'data_vec' size is:" << data_vec.size(); 
        throw std::logic_error(oss.str()); 
        return; 
    }

    //make a copy of the MLP which we can adjust the weights for 
    MultiLayerPerceptron mlp_cpy(*mlp); 

    for (size_t i=start; i<end; i++) { 
        
        const auto& data = data_vec[i]; 

        auto wg = mlp->Weight_gradient(data.inputs); 

        auto dZ = mlp->Eval(data.inputs) - data.outputs; 

        RVec<WG> dYi_dWu_dWv; dYi_dWu_dWv.reserve(n_weights); 
        //for each parameter, manually compute the Jacobian w/r/t the weights

        for (int l=0; l<mlp_cpy.Get_n_layers()-1; l++) {
            for (int j=0; j<mlp_cpy.Get_layer_size(l+1); j++) {
                for (int k=0; k<mlp_cpy.Get_layer_size(l)+1; k++) {
                    
                    mlp_cpy.Weight(l,j,k) += dx; 

                    WG wg_v = mlp_cpy.Weight_gradient(data.inputs); 
                    wg_v -= wg; 
                    wg_v *= 1./dx; 

                    dYi_dWu_dWv.push_back( std::move(wg_v) ); 

                    mlp_cpy.Weight(l,j,k) -= dx; 
                }
            }
        }

        int elem1 =0;  

        for (int l1=0; l1<mlp_cpy.Get_n_layers()-1; l1++) { 
            const size_t size1_next = mlp->Get_layer_size(l1+1); 
            const size_t size1_now  = mlp->Get_layer_size(l1) + 1; 
            for (int j1=0; j1<size1_next; j1++) {
                for (int k1=0; k1<size1_now; k1++) { //loop over all weights (row-wise)
                         
                    int elem2 =0;
                    
                    for (int l2=0; l2<mlp_cpy.Get_n_layers()-1; l2++) { 
                        const size_t size2_next = mlp->Get_layer_size(l2+1); 
                        const size_t size2_now  = mlp->Get_layer_size(l2) + 1; 
                        for (int j2=0; j2<size2_next; j2++) {
                            for (int k2=0; k2<size2_now; k2++) { //loop over all weights (column-wise)
                                
                                for (int i=0; i<mlp->Get_DoF_out(); i++) J.get(elem1, elem2) += 
                                    wg.at(i, l1, j1, k1) * wg.at(i, l2, j2, k2) + 
                                    dZ.at(i) * dYi_dWu_dWv.at(elem2).at(i, l1, j1, k1); 
                                
                                elem2++; 
                            }
                        }
                    }//loop over all weights (column-wise - elem2)

                    elem1++; 
                }
            }
        }//loop over all weights (row-wise - elem1)
    
    }//for (size_t i=start; i<end; i++)
}

//__________________________________________________________________________________________________________________
//Compute the gradient w/r/t each weight for the mlp (with respect to the given range of elements). 
void Compute_F_range(   ROOT::RVec<double>& dW, 
                        const MultiLayerPerceptron* mlp, 
                        const std::vector<TrainingData_t>& data_vec, 
                        size_t start, 
                        size_t end  )
{
    const char* const here = "Compute_F_range"; 

    //Check to make sure that the 'end' parameter given is not out-of-range. 
    if (end > data_vec.size()) {
        ostringstream oss; 
        oss << "Error in <" << here << ">: Argument 'size_t end' is invalid (" << end << "), input 'data_vec' size is:" << data_vec.size(); 
        throw std::logic_error(oss.str()); 
        return; 
    }

    
    for (size_t i=start; i<end; i++) { 
        
        const auto& data = data_vec[i]; 
        
        RVec<double> Z_err = mlp->Eval(data.inputs) - data.outputs; 

        //this is the gradient of each output coordinate (i), w/r/t each weight in the network. 
        auto weight_gradient = mlp->Weight_gradient(data.inputs); 

        //now, compute the gradient of the **loss function** w/r/t each weight: 
        int i_elem =0; 
        for (int l=0; l<mlp->Get_n_layers()-1; l++) {               // (l) - index of network layer
            for (int j=0; j<mlp->Get_layer_size(l+1); j++) {        // (j) - index of current layer output
                for (int k=0; k<mlp->Get_layer_size(l)+1; k++) {    // (k) - index of current layer input
                    
                    for (int i=0; i<mlp->Get_DoF_out(); i++) {      // (i) - index of output layer
                        dW.at( i_elem ) += Z_err[i] * weight_gradient.get(i,l,j,k);   
                    }
                    i_elem++; 
                }
            }
        }
    }//for (size_t i=start; i<end; i++)

}
//__________________________________________________________________________________________________________________

//__________________________________________________________________________________________________________________
// Compute the error w/r/t each element in the given range of the 'data_vec' for the current mlp.
double Evaluate_error_range(const std::vector<TrainingData_t>& data_vec, 
                            const MultiLayerPerceptron* mlp, 
                            size_t start, 
                            size_t end )
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
        error += rv_mag2( data.outputs - mlp->Eval(data.inputs) ); 
    }//for (size_t i=start; i<end; i++)

    return error; 
}

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

//____________________________________________________________________________________________________________________________________
int train_mlp_newtons_method(   const int n_grad_iterations = 10,
                                const char* path_infile = "",
                                const char* path_outfile = "",
                                const char* path_dbfile_starting_mlp = "",
                                RVec<int> mlp_structure={},
                                const char* tree_name="tracks_fp" )
{
    const char* const here = "train_new_mlp"; 
    
    //if there are more training events than this in the 'path_infile' root file, then only use this many. 
    const int max_events_train = 12e3; 

    //if we're training a pre-existing mlp, then we will NOT normalize the inputs (this was already done/undone in the first round of training!)
    const bool use_pre_existing_mlp = (string(path_dbfile_starting_mlp)!=""); 

    vector<string> branches_input   = {
        "x_sv", 
        "y_sv", 
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };

    vector<string> branches_output  = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 

    //this is the structure of the hidden layers. the eventual network will append an input layer and output layer on either side
    // of this set of hidden layers. 
    //
    // For example: for a network with the structure (4 inputs) => 6 => 6 => (3 outputs), this vector should be: 
    //  RVec<int> mlp_structure{6,6}; 
    //

    //add the input / output layers to this network. 
    mlp_structure.insert( mlp_structure.begin(), branches_input.size() ); 
    mlp_structure.push_back( branches_output.size() ); 

    //to create this data, we want to run in sequential mode. 
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
    vector<TrainingData_t> training_data; training_data.reserve(n_events_train); 

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
        
    vector<BranchLimits_t> lim_inputs, lim_outputs; 

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

        .Define("data", [&training_data, &lim_inputs, &lim_outputs, use_pre_existing_mlp](RVec<double>& inputs, RVec<double>& outputs)
        {   
            //normalize the inputs/outputs so that they all lie on the range x = [-1, +1]
            if (!use_pre_existing_mlp) {
                for (int i=0; i<inputs.size(); i++) {
                    double& val = inputs[i]; 
                    auto& limit = lim_inputs[i]; 

                    val = ( val - limit.mean )/( limit.max - limit.min ); 
                }
                for (int i=0; i<outputs.size(); i++) {
                    double& val = outputs[i]; 
                    auto& limit = lim_outputs[i]; 

                    val = ( val - limit.mean )/( limit.max - limit.min ); 
                }
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
    MultiLayerPerceptron* mlp;
    
    //check if we were handed a 'starting' polynomial, or should we start fresh with random weights? 
    if (string(path_dbfile_starting_mlp)!="") { 
        mlp = ApexOptics::Parse_mlp_from_file(path_dbfile_starting_mlp); 
        if (!mlp) {
            Error(here, "MLP was not parsed successfully from file: '%s'", path_dbfile_starting_mlp); 
            return -1; 
        }
    } else {
        //if not, start with a 'new' mlp, with randomly-initialized weights. 
        mlp = new MultiLayerPerceptron(mlp_structure); 
        
        //initialize with random, gaussian weights
        mlp->Add_gauss_noise(1.0); 
    }
    

    const int DoF_in  = mlp->Get_DoF_in(); 
    const int DoF_out = mlp->Get_DoF_out(); 
            
    ///////////////////////////////////////////////////////////////////////////////////////
    //  
    //  HYPERPARAMETERS - 
    //  
    //the extent to which 
    const double eta = 0.0001;

    const int update_period = 1; 
    //  
    //  
    //  
    ////////////////////////////////////////////////////////////////////////////////////////
    printf("Hyperparameters: ~~~\n"); 
    printf(" - Eta:               %.9f\n", eta);
    //mlp->Print(); 

    double x_epoch[n_grad_iterations];
    double y_error[n_grad_iterations];

    for (int i=0; i<n_grad_iterations; i++) { x_epoch[i]= i + 1; y_error[i]=0.; }

    TGraph* graph = nullptr;  
    
    char canv_title[200]; 
    sprintf(canv_title, "Eta: %.6f", eta); 
    auto canvas = new TCanvas("c", canv_title); 

        
    //number of threads 
    const size_t n_threads             = thread::hardware_concurrency(); 
    const size_t n_events_per_thread   = n_events_train / n_threads; 
    const size_t remainder             = n_events_train % n_threads; 
    
    //use all threads in parallel to evaluate the error of the passed mlp
    //______________________________________________________________________________________________________________
    auto Evaluate_error = [ &training_data, 
                            n_threads,
                            n_events_train,
                            n_events_per_thread,
                            remainder, 
                            DoF_out](const MultiLayerPerceptron* mlp)
    {   
        RVec<double> error_threads(n_threads, 0.); 
        
        size_t start = 0; 
        std::vector<thread> threads; threads.reserve(n_threads); 

        //launch all threads to evaluate the error
        for (size_t t=0; t<n_threads; t++) {
            //establish event range to run over
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0);
            
            threads.emplace_back([&training_data, &error_threads, start, end, mlp, t]
            {
                error_threads[t] = Evaluate_error_range(training_data, mlp, start, end); 
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
    double error = Evaluate_error(mlp); 
    double last_error = error; 

    RVec<RVec<double>> dW(mlp->Get_n_layers(), {}); 
    
    //zero out the weight-update vector
    for (int l=0; l<mlp->Get_n_layers()-1; l++)   
        dW[l] = RVec<double>( (mlp->Get_layer_size(l)+1) * mlp->Get_layer_size(l+1), 0. ); 
    

    printf("Running in parallel with %i threads...\n\n", (int)n_threads); 

    const int n_weights = mlp->Get_n_weights(); 


    printf("Pre-training error: % .4e\n", error); 


    //gradient descent trials. We will try to 'match' the inputs of the  
    for (int i=0; i<n_grad_iterations; i++) {

        //loop over all training data 

        //Create a dW_partial vector for each thread
        RVec<double> F_partial[n_threads]; 
        RMatrix J_partial[n_threads];  

        for (int t=0; t<n_threads; t++) {
            F_partial[t] = RVec<double>(n_weights, 0.); 
            J_partial[t] = RMatrix(n_weights, n_weights, 0.); 
        }

        //Create the partial J-matrix 
        
        //create & launch each thread (thank you to claude for teaching me how to use std::thread!)
        vector<thread> threads; threads.reserve(n_threads); 

        size_t start = 0; 
        for (size_t t=0; t<n_threads; t++) {
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0); 
            
            threads.emplace_back([&training_data, start, end, mlp, &F_partial, &J_partial, t]{

                Compute_F_range(F_partial[t],        mlp, training_data, start, end); 
                Compute_Jacobian_range(J_partial[t], mlp, training_data, start, end, 1e-10);  
            });
            start = end; 
        }

        //synchronize all the threads once they're done
        for (auto& thread : threads) thread.join();
        
        //combine all partial sums of 'J' and 'F' from each thread        
        RMatrix      &J = J_partial[0]; 
        RVec<double> &F = F_partial[0]; 

        for (int t=1; t<n_threads; t++) {
            J += J_partial[t]; 
            F += F_partial[t]; 
        }

        //now, apply the update. 
        J.Set_report_singular(false); 
        RVec<double> dX = J.Solve(F);

        dX *= eta; 
        
        //check for NaN in update. 
        bool good_update(true); 
        if (dX.size() != n_weights)  { good_update=false; } 
        else { 
            for (double& dx : dX) if (dx != dx) { good_update=false; break; } 
        }

        if (!good_update) {
            printf("Bad update; NaN update vec. breaking loop.\n"); 
            break; 
        }

        //apply the update.
        int i_elem =0;  
        for (int l=0; l<mlp->Get_n_layers()-1; l++) { for (double& weight : mlp->Get_layer(l)) weight += - dX[i_elem++]; }

        //compute error with new weights
        error = Evaluate_error(mlp);

        //print information about this epoch to stdout
        printf(" -- epoch %i, error: % .4e (%+.4e)      (progress %3.1f)\n", 
            i, 
            error, 
            error - last_error, 
            100.*((double)i+1)/((double)n_grad_iterations) ); cout << flush; 
        
        //update the graph & redraw
        
        if (i==0) {
            for (int i=0; i<n_grad_iterations; i++) y_error[i] = log(error)/log(10.);
            graph = new TGraph(n_grad_iterations, x_epoch, y_error);
        }    
    
        char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.4e;Epoch;log_{10}(RMS)", i+1, n_grad_iterations, error ); 
        graph->SetTitle(g_title);
        
        graph ->Draw(); 
        canvas->Modified(); 
        canvas->Update(); 
        
        graph->SetPointY(i, log(error)/log(10.));
        last_error = error; 
    
    }// for (int i=0; i<n_grad_iterations; i++) {
    
    graph->SetPointY(n_grad_iterations-1, log(error)/log(10.)); 
        
    char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.2e;Epoch;log_{10}(RMS)", n_grad_iterations, n_grad_iterations, error ); 
    graph->SetTitle(g_title);
    
    graph ->Draw(); 
    canvas->Modified(); 
    canvas->Update(); 

#if 0 
    //error between target / training MLPs. this is if we know what the 'toy' mlp data used to generate this data is. 
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
#endif 
    
    printf("\nCreating output file '%s'\n", path_outfile);

    //if we used a pre-existing mlp, then don't perform normalization. 
    if (use_pre_existing_mlp) {
        ApexOptics::Create_dbfile_from_mlp(path_outfile, mlp); 
        return 0; 
    }

    //now, create a new mlp, in which the inputs are (not!) normalized, as they have been for the training. 
    auto mlp_out = new MultiLayerPerceptron(mlp_structure); 

    //to begin with, copy all the weights just as they are in our training. 
    for (int l=0; l<mlp->Get_n_layers()-1; l++) mlp_out->Get_layer(l) = mlp->Get_layer(l); 


    //for the first and last layers, we need to un-do the normalization.
    //the 'normalization' is when we made is so that all of the input/output branches take on values in the range [-1,+1].   
    for (int j=0; j<mlp_out->Get_layer_size(1); j++) {
        for (int k=1; k<mlp_out->Get_layer_size(0)+1; k++) {
            
            auto& limit = lim_inputs.at(k-1); 

            mlp_out->Weight(0, j, 0) += - limit.mean * mlp->Weight(0, j, k) / (limit.max - limit.min); 

            mlp_out->Weight(0, j, k) = mlp->Weight(0, j, k) / (limit.max - limit.min); 
        }
    }
    
    //now, do the last output layer
    int last = mlp_out->Get_n_layers()-1;

    for (int j=0; j<mlp_out->Get_layer_size(last); j++) {

        auto& limit = lim_outputs.at(j); 

        mlp_out->Weight(last-1, j, 0) = mlp->Weight(last-1, j, 0) * (limit.max - limit.min)  +  limit.mean; 

        for (int k=1; k<mlp_out->Get_layer_size(last-1)+1; k++) {
            
            mlp_out->Weight(last-1, j, k) = mlp->Weight(last-1, j, k) * (limit.max - limit.min); 
        }
    }
    ApexOptics::Create_dbfile_from_mlp(path_outfile, mlp_out); 

    return 0; 
}