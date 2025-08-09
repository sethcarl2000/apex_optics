
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
#include <sstream> 
#include <stdexcept> 

using namespace std; 
using namespace ROOT::VecOps; 

//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron::MultiLayerPerceptron(const ROOT::RVec<int>& _structure) 
    : fRd(std::random_device()),
    fGen(fRd()),
    fNormal_dist(std::normal_distribution<double>(0., 1.)),
    fN_layers(_structure.size()),
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
inline double MultiLayerPerceptron::Rand_gaus()
{
    return fNormal_dist(fGen); 
}
//__________________________________________________________________________________________________________________________________
void MultiLayerPerceptron::Add_gauss_noise(double stddev) 
{
    //add random gaussian noise to all weights
    for (auto& layer : fWeights) for (auto& weight : layer) weight += Rand_gaus() * stddev; 
    return; 
}
//__________________________________________________________________________________________________________________________________
bool MultiLayerPerceptron::Check_index(int l, int j, int k) const 
{
    //Check layer index
    if ((l < 0 || l >= Get_n_layers()-1 )   || 
        (j < 0 || j >= fLayer_size[l+1])    ||
        (k < 0 || k >= fLayer_size[l]+1)) {
        
        Error("Check_index", "Invalid index: [l,j,k]=>[%i, %i, %i], can only be l=>[0,%i].",  l,j,k,  Get_n_layers()-1);
        
        return false;
    }

    return true; 
}
//__________________________________________________________________________________________________________________________________
RVec<double>& MultiLayerPerceptron::Get_layer(int l) 
{
    if (l < 0 || l >= Get_n_layers()-1 ) {
        ostringstream oss; 
        oss << "in <MultiLayerPerceptron::Get_layer(int)>: Tried to get weights for layer " 
            << l << ", valid range is l=[0," << Get_n_layers()-2 << "]"; 
        throw logic_error(oss.str()); 
        return *(new RVec<double>{});
    }
    return fWeights[l]; 
}
//__________________________________________________________________________________________________________________________________
int MultiLayerPerceptron::Get_layer_size(int l) const
{
    if (l < 0 || l >= Get_n_layers() ) {
        ostringstream oss; 
        oss << "in <MultiLayerPerceptron::Get_layer_size(int)>: Tried to get weights for layer " 
            << l << ", valid range is l=[0," << Get_n_layers()-2 << "]"; 
        throw logic_error(oss.str());
        return -1;
    }
    return fLayer_size[l]; 
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
double& MultiLayerPerceptron::WeightGradient_t::at(int i, int l, int j, int k) 
{
    //Check layer index
    if (l < 0 || l >= layer_size.size()-1 ) {
        throw logic_error("WeightGradient_t::at() - invalid 'l' index passed; must be l=[0,layer_size.size()-2]"); 
        return *(new double(numeric_limits<double>::quiet_NaN()));
    }

    if (i < 0 || i >= DoF_out) {
        throw logic_error("WeightGradient_t::at() - invalid 'j' index passed; must be i=[0,DoF_out-1]"); 
        return *(new double(numeric_limits<double>::quiet_NaN()));
    }

    if (j < 0 || j >= layer_size[l+1]) {
        throw logic_error("WeightGradient_t::at() - invalid 'j' index passed; must be j=[0,layer_size[l+1]-1]"); 
        return *(new double(numeric_limits<double>::quiet_NaN()));
    }
    
    if (k < 0 || k >= layer_size[l]+1) {
        throw logic_error("WeightGradient_t::at() - invalid 'k' index passed; must be k=[0,layer_size[l]]"); 
        return *(new double(numeric_limits<double>::quiet_NaN()));
    }
    
    return data[l][ i * (layer_size[l+1] * (layer_size[l]+1))    +   j * (layer_size[l]+1)   +   k ]; 
    //                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        ^^^^^^^^^^^^^^^^^^^^^   
    //                  Total number of weights in layer 'l'         Total number of weights for output 'j' of layer 'l'
    
}
//__________________________________________________________________________________________________________________________________
double& MultiLayerPerceptron::WeightGradient_t::get(int i, int l, int j, int k) 
{
    //same as above, but does *not* perform bounds-checking 
    //   
    return data[l][ i * (layer_size[l+1] * (layer_size[l]+1))    +   j * (layer_size[l]+1)   +   k ]; 
    //                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        ^^^^^^^^^^^^^^^^^^^^^   
    //                  Total number of weights in layer 'l'         Total number of weights for output 'j' of layer 'l'
}
//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron::WeightGradient_t::WeightGradient_t(MultiLayerPerceptron::WeightGradient_t&& val) noexcept
{
    //move constructor
    data        = move(val.data); 
    layer_size  = move(val.layer_size);
    DoF_out = val.DoF_out;

    return; 
}
//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron::WeightGradient_t& MultiLayerPerceptron::WeightGradient_t::operator=(MultiLayerPerceptron::WeightGradient_t&& val) noexcept
{
    //move assignment operator
    if (this != &val) { //check to make sure we aren't trying to copy data to ourself
        data        = move(val.data); 
        layer_size  = move(val.layer_size);
        DoF_out = val.DoF_out; 
    }
    return *this;
}//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron::WeightGradient_t& MultiLayerPerceptron::WeightGradient_t::operator-=(const MultiLayerPerceptron::WeightGradient_t& rhs)
{
    //check if there's the same number of layers between the operands 
    const size_t n_layers = layer_size.size(); 
    if (n_layers == rhs.layer_size.size()) {
        
        int l=0; 
        for (; l<n_layers; l++) { 

            //check if layer sizes match
            if (layer_size[l] != rhs.layer_size[l]) break; 
            
            //if they do, subtract rhs from lhs.
            int i=0;  
            for (auto& weight : data[l]) weight -= rhs.data[l][i++]; 
        }
        //if we got to the last layer, it means all network layers were indeed the same size. 
        // if not, continue to the exception-throwing below (one or more of the layers in 'rhs' and 'lhs' were different sizes). 
        if (l==n_layers) return *this; 
    }
    
    //network structures were not congruent, throw an exception reporting this. 
    ostringstream oss; 
    oss << "in <MLP::WeightGradient_t::operator-=>: Layer size of operands does not match."
           "\n - LHS layer structure: ";
    for (const int& size : layer_size) oss << size << " "; 
    oss << "\n - RHS layer structure: ";
    for (const int& size : rhs.layer_size) oss << size << " "; 
    throw logic_error(oss.str()); 

    //return a reference to an empty weight gradient object 
    return *(new WeightGradient_t());  
}//__________________________________________________________________________________________________________________________________
inline double MultiLayerPerceptron::Activation_fcn(double x) const
{
    //for right now, we're just going to use the cmath exp() function to make a sigmoid. 
    return ( x > 0. ? 1./( 1. + exp(-x)) : exp(x) / ( 1. + exp(x)) ) - 0.5;  
}

//__________________________________________________________________________________________________________________________________
inline double MultiLayerPerceptron::Activation_fcn_deriv(double x) const
{
    //Compute the derivative of the sigmoid 
    double S = Activation_fcn(x); 
    return ( 0.5 - S ) * ( 0.5 + S );  
}


//__________________________________________________________________________________________________________________________________
RVec<double> MultiLayerPerceptron::Activation_fcn(const RVec<double>& X) const
{
    //for right now, we're just going to use the cmath exp() function to make a sigmoid. 
    return 1./( 1. + exp(-X) ) - 0.5;  
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
    printf("MultiLayerPerceptron.\n -- structure: "); 
    printf(" (inputs: %i) => ", fLayer_size[0]);
    for (int i=1; i<Get_n_layers()-1; i++) printf("%i => ", fLayer_size[i]); 
    printf("(outputs: %i)", fLayer_size[Get_n_layers()-1]); 
    
    int n_inputs = 0; for (int l=0; l<Get_n_layers()-1; l++) n_inputs += (fLayer_size[l]+1) * (fLayer_size[l+1]); 
    printf("; Total n. of weights: %i\n", n_inputs); 

    for (int l=0; l<Get_n_layers()-1; l++) {

        printf(" -- layer connection weights: l%i -> l%i\n", l, l+1);

        for (int j=0; j<fLayer_size[l+1]; j++) {
            printf("\n -  % .4e  --- ", Get_weight(l, j, 0) ); 
            for (int k=1; k<fLayer_size[l]+1; k++) printf("% .4e ", Get_weight(l, j, k) ); 
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
MultiLayerPerceptron::WeightGradient_t MultiLayerPerceptron::Weight_gradient(const RVec<double>& X) const
{
    if ((int)X.size() != Get_DoF_in()) {
        Error("Weight_gradient", "Input vector wrong size (%i), expected (%i).", (int)X.size(), Get_DoF_in()); 
        return WeightGradient_t(); 
    }

    //we want to allocate this on the heap to avoid copying it, but it would be more convenient to access it by reference 
    // in this function. once we're done, we'll return a ptr to it.
    WeightGradient_t weight_gradient;
    
    weight_gradient.data        = RVec<double>(Get_n_layers()-1, {}); 
    weight_gradient.layer_size  = fLayer_size; 
    weight_gradient.DoF_out     = Get_DoF_out(); 


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
    RVec<RVec<double>>& grad = weight_gradient.data; 

    //the first step is actually quite similar to the feed-forward evaluation of the network. 
    //these vectors will propagate all values throughout the network.
    RVec<double> X_l[Get_n_layers()-1]; 
    RVec<double> Y_l[Get_n_layers()-1]; 
    RVec<double> Y_l_deriv[Get_n_layers()-1]; 

    for (int l=0; l<Get_n_layers()-1; l++) { 

        //initialize the gradient nested vector, and the 'layer buffers' 
        grad[l].reserve( Get_DoF_out() * fLayer_size[l+1] * (fLayer_size[l]+1) );
        //                               ^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^
        //              For layer 'l':   n. outputs         n. inputs + 1 constant      

        //initialize the layer buffers
        Y_l[l] = RVec<double>(fLayer_size[l+1], 0.); //this is the 'output' of each layer
        X_l[l].reserve(fLayer_size[l]);              //this is the 'input' to each layer
        Y_l_deriv[l].reserve(fLayer_size[l+1]);      //this is the output of each layer, with the deriv. of the A-funct. applied
    }

    X_l[0] = X; 
    
    int l=0; 
    //iterate through each layer, starting with the input layer
    for (const RVec<double>& weights : fWeights) {

        RVec<double>& input  = X_l[l]; 
        RVec<double>& output = Y_l[l]; 

        int i_elem=0; 

        //iterate over all rows (elements of the output vector)
        for (int j=0; j<fLayer_size[l+1]; j++) {

            //add the constant (which is the last element in each column of the 'weight matrix')
            output[j] += weights[i_elem++];

            //iterate through all the columns (input vector elements + a constant )
            for (int k=0; k<fLayer_size[l]; k++) output[j] += weights[i_elem++] * input[k];  
        }

        Y_l_deriv[l] = Activation_fcn_deriv(output);

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
                
                grad_l.push_back( A.get(i,j) );          // this is the 'weight' that is just a constant offset.      
                    
                for (int k=0; k<fLayer_size[l]; k++) {  //k -- index of the 'input' that this weight is associated with. 
                    
                    grad_l.push_back( A.get(i,j) * X_l[l][k] );
                } 
            }
        }

        //if we haven't reached the last layer yet, update the 'A' matrix
        if (l==0) break; 

        RVec<double> A_update_data; A_update_data.reserve(fLayer_size[l+1] * fLayer_size[l]);
        
        auto& weights = fWeights[l]; 

        for (int j=0; j<fLayer_size[l+1]; j++) {
            for (int k=0; k<fLayer_size[l]; k++) {
                A_update_data.push_back( 
                    //weights.at( j*(fLayer_size[l]+1) + 1 + k ) * Activation_fcn_deriv( Y_l[l-1].at(k) )
                    weights[ j*(fLayer_size[l]+1) + 1 + k ] * Y_l_deriv[l-1][k] 
                );
            }
        }

        RMatrix A_l(fLayer_size[l+1], fLayer_size[l], A_update_data); 

        A = A * A_l; 
    }   
    
    return weight_gradient; 
}
//__________________________________________________________________________________________________________________________________
RMatrix MultiLayerPerceptron::Jacobian(const RVec<double>& X) const 
{
    const char* const here = "Jacobian"; 
    if (X.size() != (int)Get_DoF_in()) {
        Error(here, "Input is wrong DoF (%i), should be %i", (int)X.size(), Get_DoF_in()); 
        return RMatrix(0,0); 
    }

    //first, get the 'Y' vectors for each layer
    RVec<double> X_l[Get_n_layers()-1]; 
    RVec<double> Y_l[Get_n_layers()-1]; 
    RVec<double> Y_l_deriv[Get_n_layers()-1]; 

    for (int l=0; l<Get_n_layers()-1; l++) { 
        
        //initialize the layer buffers
        Y_l[l] = RVec<double>(fLayer_size[l+1], 0.); //this is the 'output' of each layer
        X_l[l].reserve(fLayer_size[l]);              //this is the 'input' to each layer
        Y_l_deriv[l].reserve(fLayer_size[l+1]);      //this is the output of each layer, with the deriv. of the A-funct. applied
    }

    X_l[0] = X; 
    
    int l=0; 
    //iterate through each layer, starting with the input layer
    for (const RVec<double>& weights : fWeights) {

        RVec<double>& input  = X_l[l]; 
        RVec<double>& output = Y_l[l]; 

        int i_elem=0; 

        //iterate over all rows (elements of the output vector)
        for (int j=0; j<fLayer_size[l+1]; j++) {

            //add the constant (which is the last element in each column of the 'weight matrix')
            output[j] += weights[i_elem++];

            //iterate through all the columns (input vector elements + a constant )
            for (int k=0; k<fLayer_size[l]; k++) output[j] += weights[i_elem++] * input[k];  
        }

        Y_l_deriv[l] = Activation_fcn_deriv(output);

        if (l >= Get_n_layers()-2) break; 
        //apply the activation function to each element of the output vector
        X_l[l+1] = Activation_fcn(output); 
        
        //now, start again with the next row (or exit if we're done)
        l++; 
    }

    //this 'A' matrix will be what we use to back-propagate thru all the layers, starting with the last. 
    RMatrix A = RMatrix::Square_identity(Get_DoF_out()); 

    int i_elem=0; 
    for (int l=Get_n_layers()-2; l>0; l--) {

        //if we haven't reached the last layer yet, update the 'A' matrix 
        RVec<double> A_update_data; A_update_data.reserve(fLayer_size[l+1] * fLayer_size[l]);
        
        auto& weights = fWeights[l]; 
        for (int j=0; j<fLayer_size[l+1]; j++) {
            for (int k=0; k<fLayer_size[l]; k++) {
                A_update_data.push_back( weights[ j*(fLayer_size[l]+1) + (k+1) ] * Y_l_deriv[l-1][k] );
            }
        }
        RMatrix Al(fLayer_size[l+1], fLayer_size[l], A_update_data);
        A = A * Al;
    }   

    RVec<double> A_update_data; A_update_data.reserve(fLayer_size[1] * fLayer_size[0]);
    
    auto& weights = fWeights[0]; 
    for (int j=0; j<fLayer_size[1]; j++) {
        for (int k=0; k<fLayer_size[0]; k++) {
            A_update_data.push_back( weights[ j*(fLayer_size[0]+1) + (k+1) ] );
        }
    }
    
    RMatrix Al(fLayer_size[1], fLayer_size[0], A_update_data); 
    return A * Al; 
}
//__________________________________________________________________________________________________________________________________
MultiLayerPerceptron* MultiLayerPerceptron::Concantenate(MultiLayerPerceptron *mlp1, MultiLayerPerceptron *mlp2) 
{
    const char* const here = "MultiLayerPerceptron::Concantenate"; 
    
    //check to see if passed ptrs are null
    if (!mlp1 || !mlp2) {
        fprintf(stderr, "Error in <%s>: One or both of input MLP ptrs are null\n"); 
        return nullptr;  
    }

    //check to see if mlp's are of compatable size
    if (mlp1->Get_DoF_out() != mlp2->Get_DoF_in()) { 
        fprintf(stderr, "Error in <%s>: Output DoF of first MLP (%i) does not match ouptut DoF of second MLP (%i)", 
            mlp1->Get_DoF_out(), 
            mlp2->Get_DoF_in()); 
        return nullptr;  
    }

    //create the new mlp structure, by sticking both structures together. 
    RVec<int> structure; 
    //we need to remove the last element from the first mlp, as its redundant

    for (int l=0; l<mlp1->Get_n_layers(); l++) structure.push_back(mlp1->Get_layer_size(l)); 
    for (int l=1; l<mlp2->Get_n_layers(); l++) structure.push_back(mlp2->Get_layer_size(l)); 
    
    MultiLayerPerceptron *mlp_out = new MultiLayerPerceptron(structure); 
    
    int l=0; 
    for (int i=0; i<mlp1->Get_n_layers()-1; i++) { mlp_out->Get_layer(l++) = mlp1->Get_layer(i); }
    for (int i=0; i<mlp2->Get_n_layers()-1; i++) { mlp_out->Get_layer(l++) = mlp2->Get_layer(i); }

    return mlp_out; 
}
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
ClassImp(MultiLayerPerceptron); 