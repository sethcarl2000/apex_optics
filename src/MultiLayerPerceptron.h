#ifndef MultiLayerPerceptron_h_
#define MultiLayerPerceptron_h_

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

#include "TROOT.h"
#include "TObject.h"
#include <ROOT/RVec.hxx>
#include "RMatrix.h"
#include <limits> 
#include <random> 

class MultiLayerPerceptron : public TObject {
public: 

    MultiLayerPerceptron(const ROOT::RVec<int>& _structure={}); 
    ~MultiLayerPerceptron() {}; 

    //get access to underlying weight by reference (so it can be altered). 
    // - performs checks on the validity of the input arguments. 
    double& Weight(int l, int j, int k); 

    //Add random gaussian noise to all weights with given stddev 
    void Add_gauss_noise(double stddev); 

    
    //const methods which will be guranteed thread-safe for quick evaluation at runtime

    //Get number of layers
    inline int Get_n_layers()   const { return fN_layers; }

    //Get number of input nodes
    inline int Get_DoF_in()     const { return fLayer_size[0]; }
    inline int Get_DoF_out()    const { return fLayer_size[fN_layers-1]; }

    //get value of a particular weight, as a const method. 
    double Get_weight(int l, int j, int k) const; 

    ROOT::RVec<double>& Get_layer(int l); 
    int Get_layer_size(int l) const; 

    //Evaluate the mlp for a given input value
    ROOT::RVec<double> Eval(const ROOT::RVec<double>& X) const; 

    struct WeightGradient_t{ 

        //default constructor
        WeightGradient_t() noexcept 
            : data({}), layer_size({}), DoF_out(0) {}; 
        //copy constructor
        WeightGradient_t(const WeightGradient_t &cpy) noexcept
            : data(cpy.data), layer_size(cpy.layer_size), DoF_out(cpy.DoF_out) {}; 
        //move constructor
        WeightGradient_t(WeightGradient_t&& val) noexcept; 

        //move assignment operator
        WeightGradient_t& operator=(WeightGradient_t&& val) noexcept; 


        ROOT::RVec<ROOT::RVec<double>> data; 
        ROOT::RVec<int> layer_size; 
        int DoF_out; 
        double& at(int i, int l, int j, int k); 
        double& get(int i, int l, int j, int k); 
    };

    //compute the gradient w/r/t each of the weights of the network. 
    MultiLayerPerceptron::WeightGradient_t Weight_gradient(const ROOT::RVec<double>& X) const; 

    //compute the jacobian matrix relating the outputs/inputs. 
    RMatrix Jacobian(const ROOT::RVec<double>& x) const; 

    //Print Network structure and all weights
    void Print() const; 

    bool Check_index(int l, int j, int k) const; 

private: 

    std::random_device fRd; 
    std::mt19937 fGen; 
    std::normal_distribution<double> fNormal_dist; 

    inline double Rand_gaus(); 

    //the nonlinear activation function 
    inline double Activation_fcn(double x) const; 
    inline double Activation_fcn_deriv(double x) const; 

    inline ROOT::RVec<double> Activation_fcn(const ROOT::RVec<double>& X) const; 
    inline ROOT::RVec<double> Activation_fcn_deriv(const ROOT::RVec<double>& X) const; 

    int fN_layers; 

    //each element in this vector represents a layer (including input, output, and all hidden layers!).
    //each element is the number of nodes in a given layer. 
    //Indices ranging from [0, fN_layers_total-1]. 
    // The [0]th layer is the number of input nodes, the [fN_layers-1]th layer is the number of output nodes. 
    ROOT::RVec<int> fLayer_size;  

    double fQuiet_nan; 

    //This will contain all weights for each layer. 
    ROOT::RVec<ROOT::RVec<double>> fWeights; 
    


    ClassDef(MultiLayerPerceptron,1); 
};

#endif
