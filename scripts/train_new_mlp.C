#include "TROOT.h"
#include <iostream>
#include <chrono> 
#include <vector> 

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

int train_new_mlp(  const int n_events_train = 1e6, 
                    const int n_grad_iterations = 10)
{

    RVec<int> mlp_structure{3,3,3}; 

    //first, we create a 'dummy' mlp, and use it to create data: 
    auto mlp_target = new MultiLayerPerceptron(mlp_structure); 

    const int DoF_in  = mlp_target->Get_DoF_in(); 
    const int DoF_out = mlp_target->Get_DoF_out(); 

    //initialize with random, gaussian weights
    for (int l=0; l<mlp_target->Get_n_layers()-1; l++) {
        for (double& weight : mlp_target->Get_layer(l)) weight = gRandom->Gaus(); 
    }

    /*ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df(n_events_train);

    auto df_output = df

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

    df_output.Snapshot("training_data", "data/misc/mlp_train_data.root"); 
    
    */
    printf("making snapshot of %i training events...", n_events_train); cout << flush; 
    
    auto start  = chrono::high_resolution_clock::now(); 

    vector<TrainingData_t> training_data; training_data.reserve(n_events_train); 

    for (int i=0; i<n_events_train; i++) {

        TrainingData_t data{ .inputs={}, .outputs={} }; 

        data.inputs.reserve(DoF_in); 
        for (int j=0; j<DoF_in; j++) data.inputs.push_back( -1. + 2.*gRandom->Gaus() ); 

        data.outputs = mlp_target->Eval(data.inputs); 

        training_data.push_back(data); 
    }

    auto end    = chrono::high_resolution_clock::now(); 

    double time_elapsed = chrono::duration_cast<chrono::microseconds>( end - start ).count(); 
    
    printf("done.\ntime elapsed: %f seconds (%f us/event)\n", time_elapsed * 1e-6, time_elapsed/((double)n_events_train) ); 

    //now, we will attempt to 'reconstruct' this mlp with a new one... 
    vector<string> branches_in = {
        "x_inp",
        "y_inp",
        "z_inp"
    }; 

    vector<string> branches_out = {
        "x_out",
        "y_out",
        "z_out"
    };




    //we want this mlp to eventually match our target mlp 
    auto mlp = new MultiLayerPerceptron(mlp_structure); 

    //initialize with random, gaussian weights
    for (int l=0; l<mlp->Get_n_layers()-1; l++) {
        for (double& weight : mlp->Get_layer(l)) weight = gRandom->Gaus(); 
    }
            
    //the extent to which 
    double eta = 1.75e-1; 

    //the fraction of the 'existing' gradient which stays behind at the last step
    double momentum = 0.750; 

    printf("~~~~~~~~~~~~~~~~~ New mlp:"); 
    mlp->Print(); 

    double x_epoch[n_grad_iterations];
    double y_error[n_grad_iterations];
    for (int i=0; i<n_grad_iterations; i++) { x_epoch[i]=i; y_error[i]=0.; }

    TGraph* graph = nullptr;  
    
    char canv_title[200]; sprintf(canv_title, "Momentum: %f, Eta: %f", momentum, eta); 
    auto canvas = new TCanvas("c", canv_title); 

    RVec<RVec<double>> dW(mlp->Get_n_layers(), {}); 
    
        
    //for (int l=0; l<Get_n_layers()-1; l++) dW[l] = RVec<double>( (mlp->Get_layer_size(l)+1) * mlp->Get_layer_size(l+1) ); 
    //zero out the weight-update vector
    for (int l=0; l<mlp->Get_n_layers()-1; l++) dW[l] = RVec<double>( (mlp->Get_layer_size(l)+1) * mlp->Get_layer_size(l+1), 0. ); 
        
    //gradient descent trials. We will try to 'match' the inputs of the  
    for (int i=0; i<n_grad_iterations; i++) {

        //loop over all training data 
        for (const TrainingData_t& data : training_data) {

            RVec<double> Z_err = mlp->Eval(data.inputs) - data.outputs; 

            //this is the gradient of each output coordinate (i), w/r/t each weight in the network. 
            auto weight_gradient = mlp->Weight_gradient(data.inputs); 
            

            //now, compute the gradient of the **loss function** w/r/t each weight: 
            for (int l=0; l<mlp->Get_n_layers()-1; l++) {                   // (l) - index of network layer
                
                for (int j=0; j<mlp->Get_layer_size(l+1); j++) {        // (j) - index of current layer output
                    for (int k=0; k<mlp->Get_layer_size(l)+1; k++) {    // (k) - index of current layer input

                        for (int i=0; i<mlp->Get_DoF_out(); i++) {                  // (i) - index of output layer
                    
                            dW.at(l).at( j*(mlp->Get_layer_size(l)+1) + k ) += Z_err.at(i) * weight_gradient.at(i,l,j,k) * eta;   
                        }
                    }
                }
            }

        }

        //update the weights
        for (int l=0; l<mlp->Get_n_layers()-1; l++) { mlp->Get_layer(l) += - dW[l] / ((double)n_events_train); }

        //apply momentum
        for (auto& layer : dW) layer *= momentum; 

        //compute error with new weights
        double error(0.); 
        for (const TrainingData_t& data : training_data) error += rv_mag2( data.outputs - mlp->Eval(data.inputs) ); 
        error *= 1./((double)n_events_train); 
        
        //print information about this epoch to stdout
        printf("\r -- epoch %i, error: %+.4e", i, error); cout << flush;  
        
        //update the graph & redraw
        if (i==0) { 
            for (int i=0; i<n_grad_iterations; i++) y_error[i] = log(error);
            graph = new TGraph(n_grad_iterations, x_epoch, y_error);
        } else    {
            graph->SetPointY(i, log(error)); 
        } 
        char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), Error = %.2e;Epoch;log(Error)", i, n_grad_iterations, error ); 
        graph->SetTitle(g_title);
        
        graph->Draw(); 
        canvas->Modified(); 
        canvas->Update(); 
    }   

    printf("~~~~~~~~~~~~~~~~~ New mlp:"); 
    mlp->Print(); 
    printf("~~~~~~~~~~~~~~~~~ Old mlp:"); 
    mlp_target->Print(); 

    return 0; 
}