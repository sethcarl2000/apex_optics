#include "TROOT.h"
#include <iostream>
#include <chrono> 
#include <vector> 
#include <thread>
#include <fstream>
#include <sstream>
#include <limits> 
#include <stdexcept> 
#include <ROOT/RResultPtr.hxx> 
#include <NPolyArray.h>
#include <NPolyArrayChain.h> 
#include <NPoly.h>
#include <ROOT/RVec.hxx> 
#include <RMatrix.h> 
#include <ROOT/RDataFrame.hxx> 
#include <ApexOptics.h> 
#include <TCanvas.h> 
#include <cstdlib> 
#include <TRandom3.h>

using namespace std; 
using namespace ROOT::VecOps; 

struct TrainingData_t {
    ROOT::RVec<double> inputs, outputs; 
}; 

double rv_mag2(const ROOT::RVec<double>& v) {
    double ret=0.; 
    for (const double& x : v) ret += x * x; 
    return ret; 
}

#if 0 
class NPolyArrayChainFactory {
private: 
    
    const bool is_RHRS{false}; 

    struct ChainLink {
        NPolyArrayChain::EArrayType type{NPolyArrayChain::kStatic};
        NPolyArray array{}; 
        int DoF_in{0}, DoF_out{0}, order{-1};
        std::string path_infile{""}, path_outfile{""}; 
        std::vector<std::string> outputs{}; 
    };

    std::vector<ChainLink> chain_links; 

public: 
    NPolyArrayChainFactory(const bool _is_RHRS=false) : is_RHRS{_is_RHRS}, arrays{} {};
    ~NPolyArrayChainFactory() {}; 

    //append an array from a file. 
    void AppendArray(const char* path_infile, const char* path_outfile,  
                const std::vector<std::string>& inputs,
                const std::vector<std::string>& outputs, const bool is_mutable=false) {
    
        ChainLink new_link{ 
            .type = (is_mutable ? NPolyArrayChain::kMutable : NPolyArrayChain::kStatic) 
        };  

        //try to parse this array. 
        try { 
            new_link.array = ApexOptics::Parse_NPolyArray_from_file(path_infile, outputs, inputs.size()); 
        
        } catch (const std::exception& e) {

            cerr << "Exception caught in <NPolyArrayChainFactory::AppendArray>, when trying to parse input file.\n what: " 
                 << e.what() << endl; 
            std::abort(); 
            return; 
        }

        if (!chain_links.empty()) {
            if (chain_links.back().array.Get_DoF_out() != new_link.array.Get_DoF_in()) return;
        }

        new_link.DoF_in  = inputs.size(); 
        new_link.DoF_out = outputs.size(); 
        new_link.path_infile  = string(path_infile); 
        new_link.path_outfile = string(path_outfile);  
        new_link.outputs = outputs; 

        chain_links.emplace_back( new_link ); 
    } 

    //append a buffer layer. 
    void AppendBuffer(int order) {
        
        //check if there has already been a buffer this is attached to. 
        ChainLink new_link{ 
            .type = NPolyArrayChain::kBuffer 
        }; 

        new_link.order = order; 

        //check if this buffer layer is being added on top of another (not allowed!)
        if (chain_links.back().type == NPolyArrayChain::kBuffer) {

            throw logic_error("in <NPolyArrayChainFactory::AppendBuffer>: Unable to 'stack' two buffer layers on top of each other.");
            return; 
        }
 
        chain_links.emplace_back( new_link ); 
    }; 

    //create an NPolyArrayChain, with the elements 'booked' with 'AppendBuffer' and 'AppendArray' 
    NPolyArrayChain CreatePolyArrayChain() {

        //create an empty NPolyArrayChain. 
        NPolyArrayChain chain; 

        //prepare to nest the polynomials, if applicable
        for (size_t i=0; i<chain_links.size()-1; i++) { 

            auto& link = chain_links[i]; 

            //if this is a buffer layer
            if (link.type == NPolyArrayChain::kBuffer) {
                
                const int DoF_next = chain_links[i+1].array.Get_DoF_in(); 
                
                chain.arrays.emplace_back( NPolyArrayChain::CreateIdentityNPolyArray(DoF_next, link.order), NPolyArrayChain::kBuffer ); 
            }

            //if this is a static or mutable layer
            if (link.type == NPolyArrayChain::kMutable || link.type == NPolyArrayChain::kStatic) {


            }
        }

        for (auto& link : chain_links) chain.arrays.emplace_back( link.array, link.type ); 
    }

    //save the resultant polynomails in their respective output files
    void SavePolyArrayChain() {
    
        if (chain_links.empty()) { return chain; }

        //if there's only one 'link' in the chain, then 
        if (chain_links.size()<2) {
            chain.arrays.push_back({ chain_links[0].array, chain_links[0].type }); 
            return chain; 
        }

        //prepare to nest the polynomials, if applicable
        for (size_t i=0; i<chain_links.size()-1; i++) { 
            
            //if this layer is a buffer, then nest it with the next
            if (chain_links[i].type == NPolyArrayChain::kBuffer) 
                chain_links[i+1].array = NPolyArray::Nest( chain_links[i+1].array, chain_links[i].array ); 
        
        }
    }
}; 
#endif 

enum EPolynomialType {
    kSingle=0,              //single polynomial mapping fp => sv
    kSingleWithBuffer,      //single polynomial mapping fp => sv, with 'buffer' polynomials on either side
    kSplit,                 //2 sepearate polyonimals mapping fp => q1, and q1 => sv. 
    kSplitWithBuffer        //2 sepearate polynomials mapping fp => q1 and q1 =>, with 3 buffer polynomials 
}; 


//__________________________________________________________________________________________________________________
//Compute the gradient w/r/t each weight for the mlp (with respect to the given range of elements). 
void Process_event_range(   ROOT::RVec<double>& dW, 
                            const NPolyArrayChain* chain, 
                            const std::vector<TrainingData_t>& data_vec,
                            size_t start, 
                            size_t end  )
{
    const char* const here = "Process_event_range"; 

    //Check to make sure that the 'end' parameter given is not out-of-range. 
    if (end > data_vec.size()) {
        ostringstream oss; 
        oss << "Error in <" << here << ">: Argument 'size_t end' is invalid (" << end << "), input 'data_vec' size is:" << data_vec.size(); 
        throw std::logic_error(oss.str()); 
        return; 
    }

    const int DoF_out = chain->arrays.back().first.Get_DoF_out(); 

    for (size_t evt=start; evt<end; evt++) {     

        const auto& data = data_vec[evt]; 

        //start with the output layer, and continue from there. 
        RVec<RVec<double>> X{ data.inputs };  X.reserve(chain->N_arrays()+1);
        int i=0; 
        for (const auto& array : chain->arrays) X.push_back( std::move(array.first.Eval(X.back())) ); 

        //compute the 'error' of the model
        auto dZ = data.outputs - X.back(); 

        RMatrix J = RMatrix::Identity(chain->Get_DoF_out()); 

        int i_elem=0; 

        //now, compute the jacboian for each layer, and the gradient w/r/t the coefficients of each 'output' layer
        for (int i = chain->N_arrays()-1; i >= 0; i--) {

            const auto& parr = chain->arrays[i].first;
            const bool is_modifiable = (chain->arrays[i].second != NPolyArrayChain::kStatic);  

            const int DoF_array = parr.Get_DoF_out(); 

            //if this is true, that means the elements of this array are 'modifiable', so we will compute the gradient w/r/t them. 
            if (is_modifiable) {

                for (int j=0; j<parr.Get_DoF_out(); j++) {

                    const NPoly* poly = parr.Get_poly(j); 
                    
                    for (double weight : poly->Eval_noCoeff(X[i])) {

                        for (int i_out=0; i_out<DoF_out; i_out++) dW[i_elem] += dZ[i_out] * J.get(i_out, j) * weight; 
                        i_elem++;  
                    }
                }//for (int j=0; j<parr.Get_DoF_out(); j++) 
            }//if (is_modifiable)

            if (i==0) break;  
            
            J *= parr.Jacobian(X[i]); 
        }
        
    }
    
}
//__________________________________________________________________________________________________________________

//__________________________________________________________________________________________________________________
// Compute the error w/r/t each element in the given range of the 'data_vec' for the current mlp.
double Evaluate_error_range(const std::vector<TrainingData_t>& data_vec, 
                            const NPolyArrayChain *chain,
                            const size_t start, 
                            const size_t end )
{
    const char* const here = "Evaluate_error_range"; 

    //Check to make sure that the 'end' parameter given is not out-of-range. 
    if (end > data_vec.size()) {
        ostringstream oss; 
        oss << "Error in <" << here << ">: Argument 'size_t end' is invalid (" << end << "), input 'data_vec' size is:" << data_vec.size(); 
        throw std::logic_error(oss.str()); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }
    double error=0.; 

    for (size_t i=start; i<end; i++) {
        auto& data = data_vec[i]; 
        error += rv_mag2( data.outputs - chain->Eval(data.inputs) ); 
    }

    return error; 
}


#define RUN_WITHOUT_MULTITHREADDING false

#define POLYNOMIAL_TYPE kSingleWithBuffer

//____________________________________________________________________________________________________________________________________
int train_polynomial_sandwich(  const int n_grad_iterations = 10,
                                const char* path_infile = "",
                                const char* path_dbfile = "",
                                const char* path_outfile = "",
                                const int input_layer_order = 1, 
                                const int output_layer_order = 1, 
                                const char* tree_name="tracks_fp" )
{
    const char* const here = "train_polynomial_sandwich"; 

    const bool is_RHRS = false; 
    
    //if there are more training events than this in the 'path_infile' root file, then only use this many. 
    const int max_events_train = 12.5e3; 

    //put in a number here which is orders of magnitude larger than the maximum reasonable error you expect. 
    //if the error of any particular iteration exceeds this, then quit. 
    const double max_error = 1e3; 


    vector<string> branches_input = {
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    };

    vector<string> branches_output  = {
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"    
    }; 

    //create the dataframe
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
    
    ROOT::RDataFrame df(tree_name, path_infile); 
    
    //check for existence of branches: 
    vector<string> column_names = df.GetColumnNames(); 
    vector<string> missing_columns{};  
    
    //make a vector of all the columns we need to find to proceed
    vector<string> all_columns_needed; 
    copy( branches_input.begin(),  branches_input.end(),  back_inserter(all_columns_needed) ); 
    copy( branches_output.begin(), branches_output.end(), back_inserter(all_columns_needed) ); 

    //loop through all branches needed, check if any are missing.
    for ( const string& column_needed : all_columns_needed ) { 
        
        bool is_found(false);    
        
        for ( const string& column_found : column_names ) if (column_found == column_needed) is_found = true; 
        
        if (!is_found) missing_columns.push_back(column_needed); 
    }

    //some missing columns were found. 
    if (missing_columns.size()>0) {
        fprintf(stderr, "Error in <%s>: From file '%s', missing columns: ", here, path_infile); 
        for (const string& col : missing_columns) {
            fprintf(stderr, "\n  '%s'", col.data()); 
        }
        cerr << "\n"; 
        return 1; 
    }

    //add noise to reduce the chances that a 'divide by zero' error is encountered; as some of the sieve-hole data 
    // has hole positions which are exactly 0. 
    const double noise_level = 1e-7; 
    TRandom3 rand; 



    //if more training events than this exist, cap them at 'max_events_train' 
    const int n_events_train = min<int>( *df.Count(), max_events_train ); 

    printf("Using %i training events\n", n_events_train); 

    //now that we know all the columns exist, we can continue; 
    
    vector<ROOT::RDF::RNode> input_nodes{ df 
        .Range(0, n_events_train)
        .Define("inputs",  [](){ return RVec<double>{}; }, {})
        .Define("outputs", [](){ return RVec<double>{}; }, {})
    }; 

    //add the input branches
    for (const auto& str : branches_input) {

        auto new_node = input_nodes.back() 

            .Redefine("inputs",  [&rand, noise_level](RVec<double>& V, double x) { V.push_back(x + rand.Gaus()*noise_level); return V; }, {"inputs", str.data()}); 

        input_nodes.push_back(new_node); 
    }
    //add the output branches
    for (const auto& str : branches_output) {
        
        auto new_node = input_nodes.back() 

            .Redefine("outputs", [&rand, noise_level](RVec<double>& V, double x) { V.push_back(x + rand.Gaus()*noise_level); return V; }, {"outputs", str.data()}); 

        input_nodes.push_back(new_node); 
    }
    
    cout << "Compiling vector of " << n_events_train << " training events..." << flush; 

    vector<TrainingData_t> training_data = *input_nodes.back() 

        .Define("data", [](RVec<double>& inputs, RVec<double>& outputs)
        {   
            return TrainingData_t{
                .inputs  = inputs,
                .outputs = outputs
            }; 
        }, {"inputs", "outputs"})
        
        .Take<TrainingData_t>("data"); //the 'Take' command makes RDataFrame compile a vector for us, of the given type. 
    
    cout << "done." << endl; 
    //___________________________________________________________________________________________________________________________
    NPolyArray poly_array = ApexOptics::Parse_NPolyArray_from_file(path_dbfile, branches_output, (int)branches_input.size()); 

    if (poly_array.Get_status() != NPolyArray::kGood) {
        Error(here, "NPolyArray from file '%s' did not parse successfully", path_dbfile); 
        return -1; 
    }

    

    //create the 'NPolyArrayChain' which will be the ~actual~ polynomial we want, as well as the 'bread' polynomials on either side
    NPolyArrayChain *sandwich = new NPolyArrayChain; 


    //add the first layer, before the NPolyArray (which is mutable)
    //sandwich->AppendArray( Create_Identity_NPolyArray(poly_array.Get_DoF_in(), input_layer_order), true ); 

    //This is where you specify the structure of the NPolyArrayChain that you want to train. 
    //You have the option of specifying the mutability of each array (whether or not its coefficients will be adjusted)
    // and you also have the option of adding 'mutable' buffer arrays at each step. 

    //here are some examples: 
    // 
    //  auto parr = ApexOptics::Parse_NPolyArray_from_file( path_dbfile, output_branches, input_branches.size() )
    //
    //  sandwich->InsertBufferArray( parr.Get_DoF_in() );
    //
    //  sandwich->AppendArray( parr, false );
    //
    //  sandwich->InsertBufferArray( parr.Get_DoF_out() ); 
    // 
    //
    //Line-by-line, this means:
    //
    //  auto parr = ApexOptics::Parse_NPolyArray_from_file( path_dbfile, output_branches, input_branches.size() )
    // --
    // -- Create a new NPolyArray object 'parr', which has its elements defined in the file at 'path_dbfile'. 
    // -- the outputs in that file are labeled by the elements of the vector<string> 'output_branches'.   
    // -- the 'input DoF' of the polynomial is given by the number of inputs: 'input_branches.size()'. 
    //
    //  sandwich->InsertBufferArray( parr.Get_DoF_in() ); 
    // --
    // -- Create a 'buffer' polynomial, which is linear-order, with the same number of inputs/outputs (which must match that
    // -- of the input DoF of our 'parr'. 
    //
    //  sandwich->AppendArray( parr, false ); 
    // --
    // -- Add our 'parr' which we parsed from the file to the chain of poylnomials. The 'false' argument means that we will not 
    // -- 'train' the elements polynomial (they will not be modified).
    //  
    //  sandwich->InsertBufferArray( parr.Get_DoF_out() ); 
    // --
    // -- Add another 'buffer' polynomial 'on top' of our array which we parsed from the file. Like the first buffer polynomial, 
    // -- this one will have its elements modified by the training process. 
    //

    const vector<string> branches_sv{
        "x_sv",
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };

    const vector<string> branches_q1{
        "x_q1",
        "y_q1",
        "dxdz_q1",
        "dydz_q1",
        "dpp_q1"
    }; 

    const vector<string> branches_fp{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp",
    }; 
    
    const char *path_dbfile_fp_sv, *path_dbfile_fp_q1, *path_dbfile_q1_sv; 
    const char *path_outfile_fp_sv, *path_outfile_fp_q1, *path_outfile_q1_sv; 
    NPolyArray parr_fp_sv, parr_fp_q1, parr_q1_sv; 

    switch (POLYNOMIAL_TYPE) {

        case kSingle : { 
            
            path_dbfile_fp_sv  = "data/csv/poly_prod_fp_sv_L_3ord.dat";
            path_outfile_fp_sv = "data/csv/poly_prod-trained_fp_sv_L_3ord.dat";
            parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file(path_dbfile_fp_sv, {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}, 4); 

            sandwich->AppendArray( parr_fp_sv, NPolyArrayChain::kMutable ); 
            break; 
        }

        case kSingleWithBuffer : {

            path_dbfile_fp_sv  = "data/csv/poly_WireAndFoil_noise-2.5e-3_fp_sv_L_4ord.dat";
            path_outfile_fp_sv = "data/csv/poly_WireAndFoil_noise-2.5e-3_trained_fp_sv_L_4ord.dat";
            parr_fp_sv = ApexOptics::Parse_NPolyArray_from_file(path_dbfile_fp_sv, {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}, 4); 

            sandwich->InsertBufferArray( parr_fp_sv.Get_DoF_in() ); 

            sandwich->AppendArray( parr_fp_sv, NPolyArrayChain::kStatic ); 

            sandwich->InsertBufferArray( parr_fp_sv.Get_DoF_out() ); 
            break; 
        }

        case kSplit : {

            path_dbfile_fp_q1  = "data/csv/poly_prod-buffer_fp_q1_L_3ord.dat";
            path_outfile_fp_q1 = "data/csv/poly_prod-trained_fp_q1_L_3ord.dat";
            parr_fp_q1 = ApexOptics::Parse_NPolyArray_from_file(path_dbfile_fp_q1, branches_q1, 4); 

            path_dbfile_q1_sv  = "data/csv/poly_prod-buffer_q1_sv_L_3ord.dat";
            path_outfile_q1_sv = "data/csv/poly_prod-trained_q1_sv_L_3ord.dat";
            parr_q1_sv = ApexOptics::Parse_NPolyArray_from_file(path_dbfile_q1_sv, branches_sv, 5); 
                                                        
            //add the output layer (which is mutable)
            sandwich->AppendArray( parr_fp_q1, NPolyArrayChain::kMutable ); 

            sandwich->AppendArray( parr_q1_sv, NPolyArrayChain::kMutable ); 
            break; 
        }

        case kSplitWithBuffer : {

            path_dbfile_fp_q1  = "data/csv/poly_prod_fp_q1_L_3ord.dat";
            path_outfile_fp_q1 = "data/csv/poly_prod-buffer_fp_q1_L_3ord.dat";
            parr_fp_q1 = ApexOptics::Parse_NPolyArray_from_file(path_dbfile_fp_q1, {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}, 4); 

            path_dbfile_q1_sv  = "data/csv/poly_prod_q1_sv_L_3ord.dat";
            path_outfile_q1_sv = "data/csv/poly_prod-buffer_q1_sv_L_3ord.dat";
            parr_q1_sv = ApexOptics::Parse_NPolyArray_from_file(path_dbfile_q1_sv, {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}, 5); 
                                                        
            //add the output layer (which is mutable)
            sandwich->InsertBufferArray( parr_fp_q1.Get_DoF_in() ); 
            
            sandwich->AppendArray( parr_fp_q1, NPolyArrayChain::kStatic ); 

            sandwich->InsertBufferArray( parr_q1_sv.Get_DoF_in() ); 

            sandwich->AppendArray( parr_q1_sv, NPolyArrayChain::kStatic ); 

            sandwich->InsertBufferArray( parr_q1_sv.Get_DoF_out() ); 
            break; 
        } 

    }

    
    const int DoF_in  = sandwich->Get_DoF_in(); 
    const int DoF_out = sandwich->Get_DoF_out(); 


    ///////////////////////////////////////////////////////////////////////////////////////
    //  
    //  HYPERPARAMETERS - 
    //  
    //the extent to which 
    double eta          = 1e-4;

    //the fraction of the 'existing' gradient which stays behind at the last step
    double momentum     = 1.0000 - 4e-3;
    
    //the number of epochs between graph updates: 
    const int update_period = 50; 
    //  
    //  
    //  
    ////////////////////////////////////////////////////////////////////////////////////////
    printf("Hyperparameters: ~~~\n"); 
    printf(" - Eta:               %.9f\n", eta);
    printf(" - Momentum:          %.9f\n", momentum);
    printf(" - (Drag-parameter):  %.3e\n\n", (1. - momentum)/eta); 
    //mlp->Print(); 

    double x_epoch[n_grad_iterations];
    double y_error[n_grad_iterations];

    for (int i=0; i<n_grad_iterations; i++) { x_epoch[i]= i + 1; y_error[i]=0.; }

    TGraph* graph = nullptr;  
    
    char canv_title[200]; 
    sprintf(canv_title, "Momentum: %.6f, Eta: %.6f", momentum, eta); 
    auto canvas = new TCanvas("c", canv_title); 

    

    //number of threads 
    const size_t n_threads             = ( RUN_WITHOUT_MULTITHREADDING ? 1 : std::thread::hardware_concurrency() ); 
    const size_t n_events_per_thread   = n_events_train / n_threads; 
    const size_t remainder             = n_events_train % n_threads; 
    
    //use all threads in parallel to evaluate the error of the passed mlp
    //______________________________________________________________________________________________________________
    auto Evaluate_error = [ &training_data, 
                            n_threads,
                            n_events_train,
                            n_events_per_thread,
                            remainder, 
                            DoF_out](const NPolyArrayChain* chain)
    {   
        RVec<double> error_threads(n_threads, 0.); 
        
        size_t start = 0; 
        std::vector<thread> threads; threads.reserve(n_threads); 

        //launch all threads to evaluate the error
        for (size_t t=0; t<n_threads; t++) {
            //establish event range to run over
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0);
            
            threads.emplace_back([&training_data, &error_threads, start, end, chain, t]
            {
                error_threads[t] = Evaluate_error_range(training_data, chain, start, end); 
            }); 
            start = end; 
        }

        //wait for the threads to finish
        for (auto &thread : threads) thread.join(); 

        //add the sum from all threads
        double error=0.; 
        for (auto &err : error_threads) error += err;

        //return the final error
        return sqrt( error / ((double)n_events_train * DoF_out)); 
    };
    //______________________________________________________________________________________________________________

    //the error for each generation. 
    double error = Evaluate_error(sandwich); 
    double last_error = error; 

    printf("starting error: %+.4e\n", error); 

    double best_error  = 1e30; 
    double start_error = error; 

    const size_t n_parameters = sandwich->Get_nElems_mutable(); 

    //initialize the update-vector
    RVec<double> best_parameters(n_parameters, 0.); 

    //__________________________________________________________________________________________________________________
    auto Save_parameters = [&sandwich,&best_error,&best_parameters](double error) 
    {
        if (error >= best_error) return false; 

        //this is the best error found so far. 
        best_error = error; 

        //save the current weights as the best yet found
        int i_elem=0;   
        for (int i=sandwich->N_arrays()-1; i>=0; i--) { 
            
            //check if this array is marked as mutable
            if (sandwich->arrays[i].second == NPolyArrayChain::kStatic) continue; 

            auto& parr = sandwich->arrays[i].first; 

            for (int j=0; j<parr.Get_DoF_out(); j++) { 
                
                auto poly = parr.Get_poly(j); 
                for (int k=0; k<poly->Get_nElems(); k++) best_parameters[i_elem++] = poly->Get_elem(k)->coeff;
            } 
        }//for (int i=sandwich->N_arrays()-1; i>=0; i--)

        return true; 
    }; 
    //__________________________________________________________________________________________________________________
    
    //__________________________________________________________________________________________________________________
    auto Set_parameters = [](NPolyArrayChain *chain, const RVec<double>& parameters) 
    {
        if (parameters.size() != (size_t)chain->Get_nElems_mutable()) {
            ostringstream oss; 
            oss << "in <Set_parameters>: RVec<double> arguemnt (2nd) is wrong size (" << parameters.size() << "), "
                   "must be size (" << chain->Get_nElems_mutable() << ")";  
            throw invalid_argument(oss.str()); 
            return; 
        }

        //save the current weights as the best yet found
        int i_elem=0;   
        for (int i=chain->N_arrays()-1; i>=0; i--) { 
            
            //check if this array is marked as mutable
            if (chain->arrays[i].second == NPolyArrayChain::kStatic) continue; 

            auto& parr = chain->arrays[i].first; 

            for (int j=0; j<parr.Get_DoF_out(); j++) { 
                
                auto poly = parr.Get_poly(j); 
                for (int k=0; k<poly->Get_nElems(); k++) poly->Get_elem(k)->coeff = parameters[i_elem++];
            } 
        }//for (int i=chain->N_arrays()-1; i>=0; i--)
    };
    //__________________________________________________________________________________________________________________
    
    
    //this will be the 'update vector' for each iteration
    RVec<double> dW( n_parameters, 0. ); 

    printf("Running in parallel with %i threads...\n\n", (int)n_threads); 
    
    
    //gradient descent trials. We will try to 'match' the inputs of the  
    for (int i=0; i<n_grad_iterations; i++) {

        //create a dW_partial vector for each thread
        RVec<double> dW_partial[n_threads]; 
        for (size_t t=0; t<n_threads; t++) dW_partial[t] = RVec<double>( n_parameters, 0.); 

        //create & launch each thread (thank you to claude for teaching me how to use std::thread!)
        vector<thread> threads; threads.reserve(n_threads); 

        size_t start = 0; 
        for (size_t t=0; t<n_threads; t++) {
            size_t end = start + n_events_per_thread + (t < remainder ? 1 : 0); 
            
            threads.emplace_back([&training_data, start, end, sandwich, &dW_partial, t]{
               
                Process_event_range(dW_partial[t], sandwich, training_data, start, end);  
            });
            start = end; 
        }

        //synchronize all the threads once they're done
        for (auto& thread : threads) thread.join();
        
        //apply momentum
        dW *= momentum; 

        //combine all partial results from each thread 
        for (size_t t=0; t<n_threads; t++) dW += dW_partial[t] * eta / ((double)n_events_train * DoF_out); 

        //add the parameter updates to the NPolyArrayChain
        //cout << dW << endl; 
    
        *sandwich += dW;  
        
        //compute error with new weights
        error = Evaluate_error(sandwich);

        //print information about this epoch to stdout
        printf("\r -- epoch %i, error: % .4e (%+.4e); current 'best' error saved: % .4e  -- (progress %3.1f)", 
            i, 
            error, 
            error - last_error, 
            best_error,
            100.*((double)i+1)/((double)n_grad_iterations) ); cout << flush; 
        
        if (error > max_error) {
            printf("\nMaximum error exceeded (%.4e > %.4e). iteration loop terminated.\n", error, max_error); 
            break; 
        }

        if (error != error) {
            printf("\nNan Error reached. iteration loop terminated.\n"); 
            break; 
        }

        //update the graph & redraw
        if (i % update_period == 0) {

            Save_parameters(error); 
            
            if (i==0) {
                for (int i=0; i<n_grad_iterations; i++) y_error[i] = log(error)/log(10.);
                graph = new TGraph(n_grad_iterations, x_epoch, y_error);
            }    
        
            char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.4e;Epoch;log_{10}(RMS)", i+1, n_grad_iterations, error ); 
            graph->SetTitle(g_title);
            
            graph ->Draw(); 
            canvas->Modified(); 
            canvas->Update(); 
        }
        graph->SetPointY(i, log(error)/log(10.));
        last_error = error; 
    
    }// for (int i=0; i<n_grad_iterations; i++) {
    
    graph->SetPointY(n_grad_iterations-1, log(error)/log(10.)); 
        
    char g_title[200]; sprintf(g_title, "Epoch (%4i/%i), RMS = %.4e;Epoch;log_{10}(RMS)", n_grad_iterations, n_grad_iterations, error ); 
    graph->SetTitle(g_title);
    
    graph ->Draw(); 
    canvas->Modified(); 
    canvas->Update(); 

    if (best_error > start_error) {
        Error(here, "Best error found in interation loop is worse than starting error. no new MLP saved.");
        return 1; 
    }

    //perform one last check to see if the final gradient iteration was the best error yet found
    Save_parameters( Evaluate_error(sandwich) ); 

    //set the parameters of each polynomial to be those of the best-error
    Set_parameters( sandwich, best_parameters ); 

    cout << endl; 
    printf("Best error recorded: %.5e\n", best_error); 

    map<string, NPoly*> polymap;  


    switch (POLYNOMIAL_TYPE) { 

        case kSingle : {

            map<string,NPoly*> polymap; 
            int i=0; for (const auto& str : branches_sv) polymap[str] = parr_fp_sv.Get_poly(i++); 
            
            ApexOptics::Create_dbfile_from_polymap(is_RHRS, string(path_outfile_fp_sv), polymap ); 
            break;
        }

        case kSingleWithBuffer : {

            parr_fp_sv = NPolyArray::Nest( parr_fp_sv, sandwich->arrays[0].first ); 
            parr_fp_sv = NPolyArray::Nest( sandwich->arrays[2].first, parr_fp_sv ); 

            map<string,NPoly*> polymap; 
            int i=0; for (const auto& str : branches_sv) polymap[str] = parr_fp_sv.Get_poly(i++); 
            
            ApexOptics::Create_dbfile_from_polymap(is_RHRS, string(path_outfile_fp_sv), polymap ); 
            break;
        }

        case kSplit : {

            map<string,NPoly*> polymap; 
            int i=0; for (const auto& str : branches_q1) polymap[str] = parr_fp_q1.Get_poly(i++); 

            ApexOptics::Create_dbfile_from_polymap(is_RHRS, string(path_outfile_fp_q1), polymap); 

            polymap.clear(); 
            i=0; for (const auto& str : branches_sv) polymap[str] = parr_q1_sv.Get_poly(i++); 

            ApexOptics::Create_dbfile_from_polymap(is_RHRS, string(path_outfile_q1_sv), polymap); 
            break;
        }

        case kSplitWithBuffer : {

            parr_fp_q1 = NPolyArray::Nest( parr_fp_q1, sandwich->arrays[0].first ); 
            
            map<string,NPoly*> polymap; 
            int i=0; for (const auto& str : branches_q1) polymap[str] = parr_fp_q1.Get_poly(i++); 

            ApexOptics::Create_dbfile_from_polymap(is_RHRS, string(path_outfile_fp_q1), polymap); 

            parr_q1_sv = NPolyArray::Nest( parr_q1_sv, sandwich->arrays[2].first );
            parr_q1_sv = NPolyArray::Nest( sandwich->arrays[4].first, parr_q1_sv );
            
            polymap.clear(); 
            i=0; for (const auto& str : branches_sv) polymap[str] = parr_q1_sv.Get_poly(i++); 

            ApexOptics::Create_dbfile_from_polymap(is_RHRS, string(path_outfile_q1_sv), polymap); 
            break; 
        }
    }

    return 0; 

}