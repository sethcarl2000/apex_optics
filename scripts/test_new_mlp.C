#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

//test a dummy insatnce of the MultiLayerPerceptron class
int test_new_mlp() 
{
    
    RVec<double> X{ 
        gRandom->Gaus(), 
        gRandom->Gaus(), 
        gRandom->Gaus(), 
        gRandom->Gaus() 
    };

    cout << "Input X: " << X << endl; 

    MultiLayerPerceptron mlp({4,2,3,3}); 

    mlp.Get_layer(0) ={ 1,  2, 1, 1, 0, 
                        0,  0, 0, 1, 1 }; 

    mlp.Get_layer(1) ={ 1,  0, 1, 
                        2,  1, 0,  
                        1,  0, 1 }; 
    
    mlp.Get_layer(2) ={ 1,  0, 1, 1,
                        0,  0, 1, 0, 
                        1,  1, 1, 2 }; 
    
    
    //initialize each weight with a random gaussian
    /*for (int l=0; l<mlp.Get_n_layers()-1; l++) {
        for (auto& weight : mlp.Get_layer(l)) weight = gRandom->Gaus(); 
    } */   
    
    mlp.Print(); 
    
    //get a unique ptr to the 'WeightGradient_t' 
    auto wg = unique_ptr<MultiLayerPerceptron::WeightGradient_t>(mlp.Weight_gradient(X));
    
    //
    printf("\n\nGradient of selected elements (MLP Class output:)\n"); 
    printf(" dZ_0 / dW^2_02     =   %f\n", wg->at(0, 2, 0, 2));
    printf(" dZ_2 / dW^2_20     =   %f\n", wg->at(2, 2, 2, 0));
    printf(" dZ_1 / dW^1_21     =   %f\n", wg->at(1, 1, 2, 1)); 
    printf(" dZ_2 / dW^2_02     =   %f\n", wg->at(2, 2, 0, 2)); 

    auto aFcn    = [](double x){ return 1./(1. + exp(-x)) - 0.5; }; 
    auto aFcn_d  = [&aFcn](double x){ return (0.5 - aFcn(x))*(0.5 + aFcn(x)); }; 

    printf("  Result of manual calculation...\n"); 
    printf(" dZ_0 / dW^1_02     =   %f\n", aFcn(2. + aFcn(1 + 2*X[0] + X[1] + X[2])) );
    printf(" dZ_2 / dW^1_20     =   %f\n", 1. );
    printf(" dZ_1 / dW^0_21     =   %f\n", 0. ); 
    printf(" dZ_1 / dW^0_12     =   %f\n", aFcn_d(2. + X[0])*X[1] ); 
    cout << endl << endl; 
    
    RMatrix A_matrix_by_hand(3,3, { 0, 1, 1, 
                                    0, 1, 0, 
                                    1, 1, 2 }); 

    auto mat2 = RMatrix(3,3, {  aFcn_d(1 + X[1]), 0, 0, 
                                0, aFcn_d(2 + X[0]), 0,
                                0, 0, aFcn_d(1 + X[1])  } ); 
    
    
    cout << "A-matrix computed by hand:" << endl; 
    (A_matrix_by_hand * mat2).Print(); 


        

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
    
    //delete wg; 

    return 1; 
}