#include "TROOT.h"

//test a dummy insatnce of the MultiLayerPerceptron class
int test_new_mlp() 
{
    MultiLayerPerceptron mlp({2,3,3}); 

    for (int l=0; l<mlp.Get_n_layers()-1; l++) {
        for (auto& weight : mlp.Get_layer(l)) weight = gRandom->Gaus(); 
    }    

    mlp.Print(); 
    
    auto wg = mlp.Weight_gradient({1,1});


    
    return 1; 
}