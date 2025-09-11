#include <NPolyArrayChain.h>
#include <TRandom3.h> 
#include <cstdio> 
#include <iostream> 
#include <NPoly.h> 
#include <NPolyArray.h> 
#include <vector> 
#include <ROOT/RVec.hxx>

using namespace std; 

int test_chain_jacobian()
{
    //make some random NPolyArrays 
    
    const vector<int> sizes{ 5, 5, 1, 5 };
    
    const int order = 2; 

    const double mag = 1.; 

    TRandom3 rand; 

    NPolyArrayChain chain; 

    NPolyArray symbolic = NPolyArrayChain::CreateIdentityNPolyArray(sizes[0], 1); 

    for (int i=0; i<sizes.size()-1; i++) {

        const int DoF_in  = sizes[i]; 
        const int DoF_out = sizes[i+1]; 

        vector<NPoly> polys{}; 

        for (int j=0; j<DoF_out; j++) {
        
            NPoly poly(DoF_in, order); 
            
            for (int k=0; k<poly.Get_nElems(); k++) poly.Get_elem(k)->coeff = rand.Gaus() * mag; 

            polys.push_back( poly ); 
        }

        chain.AppendArray( NPolyArray(polys) );

        symbolic = NPolyArray::Nest( NPolyArray(polys), symbolic ); 

        cout << "done with array " << i << endl; 
    }   

    ROOT::RVec<double> rand_input{}; 
    for (int i=0; i<sizes[0]; i++) rand_input.emplace_back( rand.Gaus() ); 

    RMatrix J_chain    = chain   .Jacobian( rand_input );
    RMatrix J_symbolic = symbolic.Jacobian( rand_input ); 

    RMatrix J_diff = J_chain + (J_symbolic * -1.); 

    J_diff.Print(); 

    return 0; 
}