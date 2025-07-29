
//////////////////////////////////////////////////////////////////////////
//
//  MultiLayerPerceptron
// 
//  My own implementation of the TMultiLayerPerceptron class, which I'm 
//  rebuilding primarily for two reasons: 
//  1. to rebuild with performance and thread-safety as priorites from the start
//  2. to enable efficient computation of the Jacobian and Hessian (for use with)
//     newton-iteration target-trajectory reconstruction. 
//
//////////////////////////////////////////////////////////////////////////

#include "MultiLayerPerceptron.h"
#include <cmath> 
#include <stdexcept>
#include <iostream> 
#include <stdio.h>
#include <limits> 

using namespace std; 
using namespace ROOT::VecOps; 

//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron::MultiLayerPerceptron(const ROOT::RVec<int>& _structure) 
    : fN_layers(_structure.size()),
    fLayer_size(_structure),
    fQuiet_nan{numeric_limits<double>::quiet_NaN()}
{   
    //if there is at least 2 layers (including input and output layesr, which there should be!)
    //initialize all weights for each layer to 0. 

    if (Get_n_layers() < 2) {
        throw std::invalid_argument("Network structure must have at least 2 layers!");
        return; 
    }

    fWeights.reserve(Get_n_layers()-1);

    for (int l=0; l<Get_n_layers()-1; l++) {

        //The connection to each layers is comprised of: 
        //      1. A matrix which describes the linear map between one layer's outputs, and the next's inputs. 
        //      2. A constant offset (weight) for each neuron. 
        // 
        //Here, were going to initailize the RMatrix which represents these linear maps. 
        //Each Input neron has one 'weight' attaching it to each output neuron, so if the current layer's size is N, and the
        // next layer's size is M, then we need N * M weights to attach all of the input neurons to all the output neurons, 
        // and then each output neruon also has a constant offset, so in total, the space 'between' each layer has 'N * M + M' 
        // free parameters.  
        fWeights.push_back(
            ROOT::RVec<double>( (fLayer_size[l] + 1) * fLayer_size[l+1], 0. )
        );  
    }   
    
}

//__________________________________________________________________________________________________________________________________
double& MultiLayerPerceptron::Weight(int l, int j, int k) 
{      
    //Check layer index
    if (l < 0 || l >= Get_n_layers()-1 ) {
        Error("Weight", "Invalid layer index (%i), must be [1,%i]", l, Get_n_layers()-2);
        return (fQuiet_nan=numeric_limits<double>::quiet_NaN());
    }

    //check 'row' (output) index
    if (j < 0 || j >= fLayer_size[l+1]) {
        Error("Weight", "Invalid 'j' (output) index (%i), must be [0,%i]", j, fLayer_size[l+1]-1);
        return (fQuiet_nan=numeric_limits<double>::quiet_NaN());
    }

    //check 'column' (input) index
    if (k < 0 || k >= fLayer_size[l]+1) {
        Error("Weight", "Invalid 'k' (input) index (%i), must be [0,%i]", k, fLayer_size[l]);
        return (fQuiet_nan=numeric_limits<double>::quiet_NaN());  
    }

    return fWeights[l][ j * (fLayer_size[l]+1) + k ]; 
}

//__________________________________________________________________________________________________________________________________
double MultiLayerPerceptron::Get_weight(int l, int j, int k) const
{      
    //Check layer index
    if (l < 0 || l >= Get_n_layers()-1 ) {
        Error("Get_weight", "Invalid layer index (%i), must be [1,%i]", l, Get_n_layers()-2);
        return numeric_limits<double>::quiet_NaN();
    }

    //check 'row' (output) index
    if (j < 0 || j >= fLayer_size[l+1]) {
        Error("Get_weight", "Invalid 'j' (output) index (%i), must be [0,%i]", j, fLayer_size[l+1]-1);
        return numeric_limits<double>::quiet_NaN();
    }

    //check 'column' (input) index
    if (k < 0 || k >= fLayer_size[l]+1) {
        Error("Get_weight", "Invalid 'k' (input) index (%i), must be [0,%i]", k, fLayer_size[l]);
        return numeric_limits<double>::quiet_NaN();
    }

    return fWeights[l][ j * (fLayer_size[l]+1) + k ]; 
}

//__________________________________________________________________________________________________________________________________
inline double MultiLayerPerceptron::Activation_fcn(double x) const
{
    //for right now, we're just going to use the cmath exp() function to make a sigmoid. 
    return x > 0. ? 1./( 1. + exp(-x)) : exp(x) / ( 1. + exp(x)); 
}

//__________________________________________________________________________________________________________________________________
void MultiLayerPerceptron::Print() const
{
    //
    printf("MultiLayerPerceptron.\n  -- structure: "); 
    printf(" (inputs: %i) => ", fLayer_size[0]);
    for (int i=1; i<Get_n_layers()-1; i++) printf("%i => ", fLayer_size[i]); 
    printf("(outputs: %i)\n", fLayer_size[Get_n_layers()-1]); 
    
    for (int l=0; l<Get_n_layers()-1; l++) {

        printf(" -- layer connection weights: l%i -> l%i\n", l, l+1);

        for (int j=0; j<fLayer_size[l+1]; j++) {
            printf("\n    ");
            for (int k=0; k<fLayer_size[l]; k++) printf("% .4e ", Get_weight(l, j, k) ); 
            
            printf(" -  % .4e", Get_weight(l, j, fLayer_size[l]) );     
        }
        printf("\n\n");
    }
    
    cout << endl;   

    return; 
}
//__________________________________________________________________________________________________________________________________
RVec<double> MultiLayerPerceptron::Eval(const RVec<double>& X) const 
{
    if ((int)X.size() != Get_DoF_in()) {
        Error("Eval", "Input vector wrong size (%i), expected (%i).", (int)X.size(), Get_DoF_in()); 
        return {}; 
    }
    
    int l=0; //the layer we're currently operating on (starting at 0) 

    //these vectors will propagate all values throughout the network.
    RVec<double> node_vals[Get_n_layers()]; 

    node_vals[0] = X; 
    for (int l=1; l<Get_n_layers(); l++) node_vals[l] = RVec<double>(fLayer_size[l], 0.); 

    //iterate through each layer, starting with the input layer
    for (const RVec<double>& weights : fWeights) {

        RVec<double>& input  = node_vals[l]; 
        RVec<double>& output = node_vals[l+1]; 

        int i_elem=0; 

        //iterate over all rows (elements of the output vector)
        for (int j=0; j<fLayer_size[l+1]; j++) {

            //iterate through all the columns (input vector elements + a constant )
            for (int k=0; k<fLayer_size[l]; k++) output[j] += weights[i_elem++] * input[k]; 
            
            //add the constant (which is the last element in each column of the 'weight matrix')
            output[j] += weights[i_elem++]; 
        }

        //do this so we don't apply the Activation function to the last layer
        if (l >= Get_n_layers()-1) return output; 

        //apply the activation function to each element of the output vector
        for (double& out_val : output) out_val = Activation_fcn(out_val); 

        //now, start again with the next row (or exit if we're done)
        l++; 
    }

    //if we got here, something has gone wrong...
    return {}; 
}
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
ClassImp(MultiLayerPerceptron); 