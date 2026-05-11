
#include "TROOT.h"
#include <chrono> 

using namespace std; 
using namespace ROOT::VecOps; 

//true when benchmarking poly, false when benchmarking MLP
#define BENCHMARK_POLY false

int benchmark_new_mlp(const int n_trials = 5,           //number of trials to run
                      const int n_iterations = 1e6)     //number of iterations per trial
{
    RVec<int> mlp_structure{5,8,8,8,4}; 

    MultiLayerPerceptron mlp(mlp_structure); 

#if BENCHMARK_POLY

    const int poly_order = 3; 
    
    NPoly pol_template(mlp.Get_DoF_in(), poly_order); 

    //fill the poly-array with polynomials (of the proper order), with random-gauss coefficients
    vector<NPoly> poly_vec; 
    for (int i=0; i<mlp.Get_DoF_out(); i++) { 
        NPoly poly(mlp.Get_DoF_in()); 
        for (int e=0; e<pol_template.Get_nElems(); e++) poly.Add_element( pol_template.Get_elemPowers(e), gRandom->Gaus() ); 
        poly_vec.push_back(poly); 
    }

    NPolyArray parr(poly_vec); 

    printf("Number of poly-array elements: %i\n", pol_template.Get_nElems() * mlp.Get_DoF_out()); 
    printf("poly-array DoF in/out: %i/%i, poly-array order: %i\n", parr.Get_DoF_in(), parr.Get_DoF_out(), poly_order); 
#else 
    //initialize network with random, gaussian weights
    for (int l=0; l<mlp.Get_n_layers()-1; l++) {
        for (double& weight : mlp.Get_layer(l)) weight = gRandom->Gaus(); 
    }
    cout << "MLP structure (inputs => hidden layers => outputs): " << mlp_structure << endl; 
#endif 

    RVec<double> X(mlp.Get_DoF_in(),    0.); 
    RVec<double> Z(mlp.Get_DoF_out(),   0.); 
    double z; 

    
    
    
    cout << "Starting evaluation... (" << n_iterations << " iterations/trial)\n" << flush; 

    vector<double> trial_times; //net time taken to run, in seconds

    for (int t=0; t<n_trials; t++) {
        
        //measure this, so we know how must of the elapsed time measured below is due to 'gRandom->Gaus()' 
        auto start_gaus = chrono::high_resolution_clock::now();  
        for (int i=0; i<n_iterations * (mlp.Get_DoF_in()); i++) { z = gRandom->Gaus(); }
        auto end_gaus = chrono::high_resolution_clock::now(); 
        double time_gaus = chrono::duration_cast<chrono::microseconds>( end_gaus - start_gaus ).count(); 

        printf(" -- Trial %2i...", t); cout << flush; 
        
        auto start = chrono::high_resolution_clock::now(); 
        
        for (int i=0; i<n_iterations; i++) {

            //initialize random input vector
            for (double& x : X) x = gRandom->Gaus(); 
#if BENCHMARK_POLY 
            Z = parr.Eval(X); 
#else   
            //straightforward evaluation
            //Z = mlp.Eval(X); 
            //computation of the 'weight gradient'
            auto wg = mlp.Weight_gradient(X); 
#endif  
        }

        auto end = chrono::high_resolution_clock::now(); 
        double trial_time = chrono::duration_cast<chrono::microseconds>( end - start ).count() - time_gaus; 
        
        printf("done.  %3.3f seconds elapsed (%.3f us/event)\n", trial_time/1e6, trial_time/((double)n_iterations)); 
        
        trial_times.push_back(trial_time/1e6); 
    }
    cout << "done with all trials." << endl; 

    double avg(0.); 
    for (double& x : trial_times) avg += x; 
    avg *= 1./((double)n_trials); 

    double stddev(0.); 
    for (double& x : trial_times) stddev += pow( x - avg, 2 );
    stddev = sqrt(stddev/n_trials);  

    printf(" --> Trial times: ~~~~~~~~~~\n --> Average: %f\n Std-dev: %f\n", avg, stddev ); 

    return 0; 
}