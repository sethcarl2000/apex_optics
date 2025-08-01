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

    int DoF_out = mlp.Get_DoF_out(); 

    for (int i=0; i<DoF_out; i++) {
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~output index :" << i << endl; 
        for (int l=0; l<mlp.Get_n_layers()-1; l++) {
            cout << "Layer " << l << endl; 
            auto& grad_i = wg->data.at(l); 

            int size_next = mlp.Get_layer_size(l+1); 
            int size_now  = mlp.Get_layer_size(l); 

            for (int j=0; j<size_next; j++) { 
                for (int k=0; k<size_now+1; k++) { 

                    printf("%+.3e ", wg->at(i, l, j, k)); 
                }
                cout << endl; 
            }
        } 
    }
    
    
    return 1; 
}