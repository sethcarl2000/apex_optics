
#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

int benchmark_new_mlp()
{

    MultiLayerPerceptron mlp({5,5,5,5,4}); 

    //initialize network with random, gaussian weights
    for (int l=0; l<mlp.Get_n_layers()-1; l++) {
        for (double& weight : mlp.Get_layer(l)) weight = gRandom->Gaus(); 
    }

    const int n_trials = 1e7; 

    RVec<double> X(mlp.Get_DoF_in(), 0.); 
    RVec<double> Z(mlp.Get_DoF_out(), 0.); 

    NPoly poly(5); 
    NPoly pol_template(5,6); 
    
    
    for (int i=0; i<pol_template.Get_nElems(); i++) poly.Add_element( pol_template.Get_elemPowers(i), gRandom->Gaus() ); 
    cout << "Number of polynomial elements: " << poly.Get_nElems() << endl; 

    //measure this, so we know how must of the elapsed time measured below is due to 'gRandom->Gaus()' 
    TStopwatch timer_gaus; 
    for (int i=0; i<n_trials * (mlp.Get_DoF_in()); i++) { double x = gRandom->Gaus(); }
    double time_gaus = timer_gaus.RealTime(); 

    TStopwatch timer; 
    for (int i=0; i<n_trials; i++) {

        //initialize random input vector
        for (double& x : X) x = gRandom->Gaus(); 

        double z = poly.Eval(X); 
    }
    double time_eval = timer.RealTime() - time_gaus; 

    printf(" %f seconds elapsed to process %i events (%f us/event)\n", time_eval, n_trials, 1e6 * time_eval/((double)n_trials)); 

    return 0; 
}