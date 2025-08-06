#include "TROOT.h"

using namespace std; 
using namespace ROOT::VecOps; 

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
int create_toy_mlp_data(    const int n_events_train=1e5, 
                            const char* path_outfile="data/misc/mlp_toy_data.root",
                            const RVec<int>& hidden_layer_structure={},    
                            const char* tree_name="tracks_fp",
                            const char* path_toy_dbfile=""  )
{
    const char* const here = "create_toy_mlp_data"; 
 
    vector<string> branches_input   = {
        "x_sv", 
        "y_sv", 
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };
    const size_t DoF_in = branches_input.size(); 

    vector<string> branches_output  = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 
    const size_t DoF_out = branches_output.size(); 

    
    //add the input / output layers to this network. 
    RVec<int> mlp_structure = hidden_layer_structure; 
    mlp_structure.insert( mlp_structure.begin(), DoF_in ); 
    mlp_structure.push_back( DoF_out );
    
    auto mlp = new MultiLayerPerceptron(mlp_structure); 

    //randomize all weights in the network
    mlp->Add_gauss_noise(1.0); 

    cout << "Toy mlp: " << endl; 
    mlp->Print(); 

    //now, create the RDataFrame which 
    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df_create(n_events_train);

    //add all the inputs with the proper names
    vector<ROOT::RDF::RNode> nodes{ df_create };
    
    nodes.push_back( nodes.back().Define("X", [DoF_in] (){ RVec<double> v{}; v.reserve(DoF_in); return v; }, {}) ); 

    for (const auto& bname : branches_input) {

        auto new_node = nodes.back()
        
        //each input branch will be evenly-distributed over the interval x=[-1, +1]
        .Define(bname.data(), [](){ return -1. + 2.*gRandom->Rndm(); }, {}) 
        
        //now, add this input branch to the inptut vector
        .Redefine("X",  [](RVec<double>& X, double x){ X.push_back(x); return X; }, {"X", bname.data()}); 

        nodes.push_back(new_node); 
    }

    //now that we've added all inputs to the inptu vector "X", define the output vector "Z".
    nodes.push_back( nodes.back().Define("Z", [mlp](const RVec<double>& X){ return mlp->Eval(X); }, {"X"}) ); 

    //and define each output branch
    int i_output=0; 
    for (const auto& bname : branches_output) {

        nodes.push_back( nodes.back().Define(bname.data(), [i_output](const RVec<double>& Z){ return Z.at(i_output); }, {"Z"}) ); 
        i_output++; 
    }
    
    //now, were ready to make a 'snapshot' of all our branches. 
    //do this by putting all branch names into a single vector
    vector<string> all_branches_to_write = branches_input;
    copy( branches_output.begin(), branches_output.end(), back_inserter(all_branches_to_write) ); 


    printf("Info in <%s>: Creating snapshot of %i events...", here, n_events_train); cout << flush; 
    
    nodes.back().Snapshot(tree_name, path_outfile, all_branches_to_write); 
    
    cout << "done." << endl; 

    printf("Info in <%s>: Data written to '%s'\n", here, path_outfile); 


    return 0; 
}