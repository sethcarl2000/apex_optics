#include "TROOT.h"
#include <cmath> 

using namespace std; 
using namespace ROOT::VecOps; 
using MLP = MultiLayerPerceptron; 

int test_mlp_derivatives()
{
    //test some manually-computed derivatives.
    RVec<double> X{ gRandom->Gaus(), gRandom->Gaus() };
    
    MLP mlp({2,2,2,2}); 

    mlp.Get_layer(0) = {
        0,  2, 0,
        1,  1, 1
    }; 

    mlp.Get_layer(1) = {
        0,  0, 1, 
        0,  1, 0
    }; 

    mlp.Get_layer(2) = {
        1,  1, 3, 
        0,  2, 0
    }; 

    auto AF     = [](double x)      { return 1./( 1. + exp(-x) ) - 0.5; }; 
    auto AF_d   = [&AF](double x)   { return (0.5 - AF(x))*(0.5 + AF(x)); }; 
    auto AF_dd  = [](double x)      { return -2. * pow( sinh(x), -3 ) * pow( cosh(x/2.), 4 ); }; 

    RMatrix J_manual(2,2, 0.); 

    double y1_x1 = X[0] + X[1] + 1; 
    y1_x1 = AF_d(AF(y1_x1)) * AF_d(y1_x1); 
    J_manual.at(0,0) = y1_x1 + 3 * AF_d(AF(2*X[0]))*AF_d(2*X[0])*2.; 
    J_manual.at(0,1) = y1_x1; 

    J_manual.at(1,0) = 2.*y1_x1; 
    J_manual.at(1,1) = 2.*y1_x1; 


    auto J = mlp.Jacobian(X); 

    printf("J (manual vs MLP-computed diff): "); 

    (J + (J_manual * -1.)).Print();

    MLP::HessianTensor_t H_manual; 
    H_manual.DoF_in  = 2; 
    H_manual.DoF_out = 2; 
    H_manual.data = RVec<double>(8, 0.); 
    
    double fdd_fd2 = AF_dd(AF(X[0] + X[1] + 1.)) * pow( AF_d(X[0] + X[1] + 1.), 2 ); 
    double fd_fdd  = AF_d(AF(X[0] + X[1] + 1.)) * AF_dd(X[0] + X[1] + 1.); 

    double fdd_fd2_2x = AF_dd(AF(2.*X[0])) * pow( AF_d(2.*X[0]), 2 ); 
    double fd_fdd_2x = AF_d(AF(2.*X[0])) * AF_dd(2.*X[0]); 
    

    H_manual.at(0, 0, 0) = fdd_fd2 + fd_fdd + 12.*fdd_fd2_2x + 12.*fd_fdd_2x; 
    H_manual.at(0, 0, 1) = fdd_fd2 + fd_fdd;  
    H_manual.at(0, 1, 0) = fdd_fd2 + fd_fdd;  
    H_manual.at(0, 1, 1) = fdd_fd2 + fd_fdd;  
    
    H_manual.at(1, 0, 0) = 2.*(fdd_fd2 + fd_fdd); 
    H_manual.at(1, 0, 1) = 2.*(fdd_fd2 + fd_fdd);  
    H_manual.at(1, 1, 0) = 2.*(fdd_fd2 + fd_fdd);  
    H_manual.at(1, 1, 1) = 2.*(fdd_fd2 + fd_fdd);  
    
    auto H = mlp.Hessian_tensor(X); 

    for (int i=0; i<2; i++) {
        printf("\nlayer %i \n", i); 
        for (int j=0; j<2; j++) {
            printf("\n - "); 
            for (int k=0; k<2; k++) printf("  %+.4e", H_manual.at(i,j,k) - H.at(i,j,k));  
        }
    }

    for (int i=0; i<2; i++) {
        printf("\nlayer %i \n", i); 
        for (int j=0; j<2; j++) {
            printf("\n - "); 
            for (int k=0; k<2; k++) printf("  %+.4e", H.at(i,j,k));  
        }
    }

    return 0; 
}