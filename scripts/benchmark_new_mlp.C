
#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

//true when benchmarking poly, false when benchmarking MLP
#define BENCHMARK_POLY false

int benchmark_new_mlp()
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

    const int n_trials = 1e7; 

    RVec<double> X(mlp.Get_DoF_in(),    0.); 
    RVec<double> Z(mlp.Get_DoF_out(),   0.); 
    double z; 

    
    //measure this, so we know how must of the elapsed time measured below is due to 'gRandom->Gaus()' 
    TStopwatch timer_gaus; 
    for (int i=0; i<n_trials * (mlp.Get_DoF_in()); i++) { z = gRandom->Gaus(); }
    double time_gaus = timer_gaus.RealTime(); 

    cout << "Starting evaluation..." << flush; 
    TStopwatch timer; 
    for (int i=0; i<n_trials; i++) {

        //initialize random input vector
        for (double& x : X) x = gRandom->Gaus(); 
#if BENCHMARK_POLY 
        Z = parr.Eval(X); 
#else 
        Z = mlp.Eval(X);
#endif  
    }
    double time_eval = timer.RealTime() - time_gaus; 
    cout << "done." << endl; 

    printf(" %f seconds elapsed to process %i events (%f us/event)\n", time_eval, n_trials, 1e6 * time_eval/((double)n_trials)); 

    return 0; 
}