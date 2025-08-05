#include "TROOT.h"
#include "TF2.h"

using namespace std; 
using namespace ROOT::VecOps; 

MultiLayerPerceptron *mlp; 
int ind=0; 

double mlp_eval(double *x, double *par) { return mlp->Eval({x[0],x[1]}).at(par[0]); } 

int draw_dummy_mlp()
{   
    RVec<int> mlp_structure{2,5,5,4}; 

    //create the MLP, and give each weight a random, gaussian value
    mlp = new MultiLayerPerceptron(mlp_structure); 
    mlp->Add_gauss_noise(1.0); 

    //create this function which we will feed into the mlp 
    auto canv = new TCanvas("c", "canvas", 1200, 800); 
    canv->Divide(2,2, 0.01,0.01); 

    for (int i=0; i<4; i++) {

        canv->cd(i+1); 

        auto fcn = new TF2("mlp", mlp_eval, -2,2, -2,2, 1);
        
        fcn->SetParameter(0, i);  
        fcn->Draw("surf"); 
    }

    return 0; 
}