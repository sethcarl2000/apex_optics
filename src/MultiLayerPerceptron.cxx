
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
bool MultiLayerPerceptron::Check_index(int l, int j, int k) const 
{
    //Check layer index
    if ((l < 0 || l >= Get_n_layers()-1 )   || 
        (j < 0 || j >= fLayer_size[l+1])    ||
        (k < 0 || k >= fLayer_size[l]+1)) {
        
        Error("Check_index", "Invalid index: [l,j,k]=>[%i, %i, %i], can only be l=>[0,%i].",  l,j,k,  Get_n_layers()-1);

        Print(); 
        
        return false;
    }

    return true; 
}
RVec<double>& MultiLayerPerceptron::Get_layer(int l) 
{
    if (l < 0 || l >= Get_n_layers()-2 ) {
        throw logic_error("Tried to acess layer out-of-range; can only be l=[0, Get_n_layers()-2]"); 
        return fWeights[0];
    }
    return fWeights[l]; 
}
//__________________________________________________________________________________________________________________________________
double& MultiLayerPerceptron::Weight(int l, int j, int k) 
{      
    if (!Check_index(l,j,k)) return (fQuiet_nan=numeric_limits<double>::quiet_NaN()); 

    return fWeights[l][ j * (fLayer_size[l]+1) + k ]; 
}

//__________________________________________________________________________________________________________________________________
double MultiLayerPerceptron::Get_weight(int l, int j, int k) const
{      
    if (!Check_index(l,j,k)) return numeric_limits<double>::quiet_NaN(); 

    return fWeights[l][ j * (fLayer_size[l]+1) + k ]; 
}

//__________________________________________________________________________________________________________________________________
inline double MultiLayerPerceptron::Activation_fcn(double x) const
{
    //for right now, we're just going to use the cmath exp() function to make a sigmoid. 
    return ( x > 0. ? 1./( 1. + exp(-x)) : exp(x) / ( 1. + exp(x)) ) - 0.5;  
}

//__________________________________________________________________________________________________________________________________
inline double MultiLayerPerceptron::Activation_fcn_deriv(double x) const
{
    //Compute the derivative of the sigmoid 
    double S = Activation_fcn(x) + 0.5; 
    return S * ( 1. - S ); 
}


//__________________________________________________________________________________________________________________________________
RVec<double> MultiLayerPerceptron::Activation_fcn(const RVec<double>& X) const
{
    //for right now, we're just going to use the cmath exp() function to make a sigmoid. 
    return 1./( 1. + exp(X) ) - 0.5;  
}

//__________________________________________________________________________________________________________________________________
RVec<double> MultiLayerPerceptron::Activation_fcn_deriv(const RVec<double>& X) const
{
    //Compute the derivative of the sigmoid 
    auto Y = Activation_fcn(X); 
    return ( 0.5 + Y ) * ( 0.5 - Y ); 
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

            //add the constant (which is the last element in each column of the 'weight matrix')
            output[j] += weights[i_elem++]; 
            
            //iterate through all the columns (input vector elements + a constant )
            for (int k=0; k<fLayer_size[l]; k++) output[j] += weights[i_elem++] * input[k]; 
        }

        //do this so we don't apply the Activation function to the last layer
        if (l >= Get_n_layers()-2) return output; 

        //apply the activation function to each element of the output vector
        for (double& out_val : output) out_val = Activation_fcn(out_val); 

        //now, start again with the next row (or exit if we're done)
        l++; 
    }

    //if we got here, something has gone wrong...
    return {}; 
}
//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron::WeightGradient_t* MultiLayerPerceptron::Weight_gradient(const RVec<double>& X) const
{
    if ((int)X.size() != Get_DoF_in()) {
        Error("Weight_gradient", "Input vector wrong size (%i), expected (%i).", (int)X.size(), Get_DoF_in()); 
        return nullptr; 
    }

    //we want to allocate this on the heap to avoid copying it, but it would be more convenient to access it by reference 
    // in this function. once we're done, we'll return a ptr to it.
    WeightGradient_t *weight_gradient = new WeightGradient_t{ 
        .data       = RVec<double>(Get_n_layers()-1, {}), 
        .layer_size = fLayer_size 
    };  

    //structure is: 
    // outermost vector: each element is a coordinate of the output vector (i)
    // middle vector: each element is a distinct layer of weights (l)
    // inner vector:  each element is the gradient computed w/r/t a specific weight. so: 
    //
    //      dZ_i / dW^l_jk =  grad[l][ (i * layer_size[l+1] * (layer_size[l]+1)) + (j * (layer_size[l]+1)) + k ]
    //                                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        ^^^^^^^^^^^^^^^^^^   
    //                                      Total number of weights in layer 'l'        Total number of weights for output 'j' of layer 'l'
    // Where: 
    // Z_i is the output node the gradient is computed in respect to. 
    // dW^l_jk is the specific weight. 'l' is layer [0,N-1], 'k' is input variable index (x^l_k), and 'j' is output variable index (y^l_j)
    //   note each element with index k=0 is the constant-offset index. 
    //   so each layer is computed as: 
    // 
    //      y^l_j = w^l_jk * x^l_k + w^l_j0 
    //
    //   where the index 'k' is summed over on the RHS. 
    //
    RVec<RVec<double>>& grad = weight_gradient->data; 

    //the first step is actually quite similar to the feed-forward evaluation of the network. 
    //these vectors will propagate all values throughout the network.
    RVec<double> X_l[Get_n_layers()-1]; 
    RVec<double> Y_l[Get_n_layers()-1]; 

    for (int l=0; l<Get_n_layers()-1; l++) { 

        //initialize the gradient nested vector, and the 'layer buffers' 
        grad[l].reserve( Get_DoF_out() * fLayer_size[l+1] * (fLayer_size[l]+1) );
        //                               ^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^
        //              For layer 'l':   n. outputs         n. inputs + 1 constant      

        //initialize the layer buffers
        Y_l[l] = RVec<double>(fLayer_size[l+1], 0.); //this is the 'output' of each layer
        X_l[l].reserve(fLayer_size[l]);              //this is the 'input' to each layer
    }

    X_l[0] = X; 
    
    
    int l=0; 
    //iterate through each layer, starting with the input layer
    for (const RVec<double>& weights : fWeights) {

        RVec<double>& input  = X_l[l]; 
        RVec<double>& output = Y_l[l]; 

        int i_elem=0; 

        //iterate over all rows (elements of the output vector)
        for (int j=0; j<fLayer_size.at(l+1); j++) {

            //add the constant (which is the last element in each column of the 'weight matrix')
            output.at(j) += weights.at(i_elem++);

            //iterate through all the columns (input vector elements + a constant )
            for (int k=0; k<fLayer_size.at(l); k++) output[j] += weights.at(i_elem++) * input.at(k);  
        }

        if (l >= Get_n_layers()-2) break; 
        //apply the activation function to each element of the output vector
        X_l[l+1] = Activation_fcn(output); 
        
        //now, start again with the next row (or exit if we're done)
        l++; 
    }
    
    //now that we have cached all the layers, we're ready to start computing the gradient. 
    //we start with the last layer, and recursivley propagate all the way to the first layer. 
    
    //this 'A' matrix will be what we use to back-propagate thru all the layers, starting with the last. 
    RMatrix A = RMatrix::Square_identity(Get_DoF_out()); 

    int i_elem=0; 
    for (int l=Get_n_layers()-2; l>=0; l--) {

        //get the gradient for this layer
        auto& grad_l = grad[l];        
        
        for (int i=0; i<Get_DoF_out(); i++) {           //i -- index of the ouput we're computing the gradient w/r/t 
            for (int j=0; j<fLayer_size[l+1]; j++) {    //j -- index of the 'output' that this weight is associated with. 
                
                grad_l.push_back( A.at(i,j) );          // this is the 'weight' that is just a constant offset.      
                    
                for (int k=0; k<fLayer_size[l]; k++) {  //k -- index of the 'input' that this weight is associated with. 
                    
                    grad_l.push_back( A.at(i,j) * X_l[l].at(k) );
                } 
            }
        }

        //if we haven't reached the last layer yet, update the 'A' matrix
        if (l==0) break; 

        RVec<double> A_update_data; A_update_data.reserve(fLayer_size[l] * fLayer_size[l-1]);
        
        auto& weights = fWeights[l]; 

        for (int i=0; i<fLayer_size[l]; i++) {
            for (int j=0; j<fLayer_size[l-1]; j++) {
                A_update_data.push_back( 
                    weights.at( i*(fLayer_size[l]+1) + 1 + j ) * Activation_fcn_deriv(Y_l[l-1].at(j))
                ); 
            }
        }
        RMatrix A_update(fLayer_size[l], fLayer_size[l-1], A_update_data); 

        A = A * A_update; 
    }   
    
    return weight_gradient; 
}
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
ClassImp(MultiLayerPerceptron); 