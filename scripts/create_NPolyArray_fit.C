//stdlib headers
#include <memory>
#include <string> 
#include <sstream> 
#include <map>
#include <fstream>
#include <iomanip>
#include <thread> 
#include <cmath> 
#include <memory> 
#include <vector> 
#include <utility> 
//ROOT headers 
#include <TFile.h>
#include <TVector3.h>
#include <TParameter.h>
#include <TSystem.h> 
#include <TDatime.h> 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
//Apex-optics headers
#include <ApexOptics.h> 
#include "include/RDFNodeAccumulator.h"

//should we use multithreadding when element-pruning? 
#define MULTITHREAD_PRUNING true 

using namespace std;

//sieve coords  {"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"}
//q1 coords     {"x_q1","y_q1","dxdz_q1","dydz_q1","dpp_q1"}
//q1-fwd coords {"fwd_x_q1","fwd_y_q1","fwd_dxdz_q1","fwd_dydz_q1","fwd_dpp_q1"} 
//fp coords     {"x_fp","y_fp","dxdz_fp","dydz_fp"}

//Taking the RDataFrame, this computes the maximum value each element contributes over the whole dataset. 
vector<pair<const NPoly::NPolyElem*,double>> Get_maximum_element_values(ROOT::RDF::RNode df, const NPoly& poly, const vector<string>& inputs); 

//creates db '.dat' files for polynomials which are meant to map from focal-plane coordinates to sieve coordinates. 
//if you want to avoid creating (or overwriting) a db-file, just enter "" as the db-file name, and only the histograms will be drawn instead. 
int create_NPolyArray_fit(  const int poly_order     =4,
                            const char* path_infile  ="data/mc/mc_R_production.root",
                            const char* stem_outfile ="data/poly/mc-production_sv_fp_R_4ord.dat",  
                            const vector<string> inputs  ={"x_sv","y_sv","dxdz_sv","dydz_sv","dpp_sv"},
                            const vector<string> outputs ={"x_fp","y_fp","dxdz_fp","dydz_fp"},
                            const double pruning_constant=-1., 
                            const char* tree_name    ="tracks_fp") 
{
    if (inputs.empty() || outputs.empty()) {
        Error(__func__, "inputs and/or outputs are empty");
        return -1; 
    }
    
    //check if we can load the apex optics lib
    if (gSystem->Load("libApexOptics") < 0) {
        Error(__func__, "libApexOptics could not be loaded."); 
        return 1; 
    }

    auto infile = new TFile(path_infile, "READ");

    //check if the infile could be opened
    if (!infile || infile->IsZombie()) {
        Error(__func__, "root file '%s' could not be opened.", path_infile); 
        return 1; 
    }

    //check if we can find the 'is_RHRS' parameter. Fatal error if not! 
    TParameter<bool>* param_is_RHRS = (TParameter<bool>*)infile->Get("is_RHRS"); 
    if (!param_is_RHRS) {
        Error(__func__, "Could not find TParameter<bool> 'is_RHRS' in file '%s'.", path_infile); 
        return 1; 
    }
    const bool is_RHRS = param_is_RHRS->GetVal(); 

    infile->Close(); 
    delete infile; 

    ROOT::EnableImplicitMT(); 
    Info(__func__, "Multi-threadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize()); 

    //Now, we are ready to process the tree using RDataFrame

    auto df_ptr = unique_ptr<ROOT::RDataFrame>(nullptr); 
    try { df_ptr = unique_ptr<ROOT::RDataFrame>(new ROOT::RDataFrame(tree_name, path_infile)); }
    catch (const std::invalid_argument& e) {
        Error(__func__, "Trying to create RDataFrame threw std::invalid_argument exception.\n what(): %s", e.what()); 
        return -1; 
    } 
    ROOT::RDataFrame& df = *df_ptr; 
    
    const size_t DoF_in  = inputs.size();  
    const size_t DoF_out = outputs.size();  

    cout << "Creating polynomials for inputs => outputs..." << flush; 
    
    //for each of the branches in the 'branches_sv' created above, make a polynomial which takes all the branches
    // of the 'branches_fp' vec 
    map<string,NPoly*> polymap = ApexOptics::Create_NPoly_fit(df, poly_order, inputs, outputs);

    cout << "done." << endl;    

    //check if we're going to do pruning
    if (pruning_constant > 0.) {

        printf(
            "Performing polynimal pruning. minimum element value over training dataset: %.4e\n", 
            pruning_constant
        ); 
        
        for (auto poly_it = polymap.begin(); poly_it != polymap.end(); poly_it++) {
            
            printf(
                "Pruning polynomial '%s'...", 
                poly_it->first.data()
            );
            cout << flush; 

            auto poly = poly_it->second; 
            
            //create a polynomial with the correct, 'pruned' DoF. 
            NPoly *poly_pruned = new NPoly(poly->Get_nDoF()); 

            auto&& max_elem_values = Get_maximum_element_values(df, *poly, inputs);
            
            for (const auto& elem_and_value : max_elem_values) {

                //if the maximum value this element takes on over the whole training dataset 
                // is above the user's specified miminum value, then we can keep it. 
                if (elem_and_value.second > pruning_constant) {
                    poly_pruned->Add_element(*elem_and_value.first); 
                }
            }

            printf(
                "done.\n"
                "number of elemnents in polynomial (before pruning / after pruning):  %i / %i\n", 
                poly->Get_nElems(), poly_pruned->Get_nElems()
            ); cout << flush; 

            //now, replace the originial polynomial with the pruned one. 
            poly_it->second = poly_pruned; 
        }
    }

    string path_outfile(stem_outfile);

    //create the output file ___________________________________________
    if (path_outfile != "") {
        cout << "Creating dbfiles for polynomials..." << flush; 
        
        // first, we're going to write some of our own data to the file.
        // the 'trunc' flag means that we're going to wipe the contents of the file and start over.  
        fstream outfile(path_outfile, ios::out | ios::trunc); 

        if (!outfile.is_open()) {
            Error(__func__, "Unable to open file '%s'.", path_outfile.c_str()); 
            return -1; 
        }

        //now that we can open the file, let's add some metadata to it. 
        ostringstream metadata; 
        
        //add basic metadata 
        metadata <<   
            "# Polynomial metadata: \n"
            "# -- Order: " << poly_order << "\n"
            "# -- Trained on data: '" << path_infile << "'\n";
          
        //add the branches; 
        metadata << "# -- Input branches: { "; 
        for (const auto& br : inputs) metadata << br << " "; 
        metadata << "}\n"; 

        //add a timestamp
        TDatime dtime; 
        metadata << 
            "# -- Time created: " << dtime.AsString() << "\n"; 

        if (pruning_constant > 0.) {
            metadata << "# -- Polynomial-pruning element min value: " << pruning_constant << "\n"; 
        }

        cout << metadata.str() << endl; 

        outfile << metadata.str() << endl; 
        outfile.close(); 

        ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap, true); 

        cout << "done." << endl; 

    } else { 

        cout << "Skipping db-file creation." << endl; 
    }
    //____________________________________________________________________________

    //delete our poly models
    for (auto it = polymap.begin(); it != polymap.end(); it++ ) delete it->second;
    
    return 0;
}


//________________________________________________________________________________________________________________________________
vector<pair<const NPoly::NPolyElem*,double>> Get_maximum_element_values(ROOT::RDF::RNode df, const NPoly& poly, const vector<string>& inputs)
{
    //first, let's get a vector (of events) for training event
    using namespace ROOT::VecOps; 

    RDFNodeAccumulator rna(df); 

    //create the (empty) inputs vector
    rna.Define("X_inputs", [](){ return RVec<double>{}; }, {}); 

    //now, fill the input vector with all our inputs
    for (const auto& input : inputs) {
        rna.Overwrite("X_inputs", [](RVec<double>& X, double inp){ 
            X.push_back(inp); return X; 
        }, {"X_inputs", input.data()}); 
    }    

    //now, tell the RDataFrame to give us a vector of all these events (which are each vectors themeselves)
    vector<RVec<double>> events = *rna.Get().Take<RVec<double>>("X_inputs"); 

    //and, finally, we can loop over all events and get the max values for each one. 
    const int n_threads = MULTITHREAD_PRUNING ? thread::hardware_concurrency() : 1;

    //this vector will hold the maximum values for each element. 
    vector<pair<const NPoly::NPolyElem*, double>> max_vals[n_threads];
    
    //initialize the vectors 
    const int n_elems = poly.Get_nElems(); 
    for (int t=0; t<n_threads; t++) {
        max_vals[t].reserve(n_elems); 
        for (int i=0; i<n_elems; i++) max_vals[t].push_back({poly.Get_elem(i), 0.}); 
    } 

    vector<thread> threads; 
    
    size_t start=0; 
    const size_t n_events = events.size(); 
    const size_t n_events_per_thread = n_events/n_threads; 
    
    for (int t=0; t<n_threads; t++) {

        //decide how many events to process with this thread 
        size_t end = start + n_events_per_thread; 

        //account for the fact that the number of threads may not evenly divide the number of events
        end += t < (n_events % n_threads) ? 1 : 0;  

        //launch all the threads 
        threads.emplace_back([t,start,end,&max_vals,&events,&poly,n_elems]{
            for (size_t i=start; i<end; i++) {
                
                //get one particular traning-data 'event'
                const auto& event = events[i]; 

                //evaluate each polynomial's ouptut for this event
                const auto&& elem_vals = poly.Eval_noCoeff(event); 

                //find the maximum value each element produces for this event. 
                double val;
                for (int e=0; e<n_elems; e++) {
                    val = fabs(elem_vals[e]); 
                    if (val > max_vals[t][e].second) max_vals[t][e].second = val;    
                }
            }
        }); 
    }
    
    //now, join all the threads
    for (auto& t : threads) t.join(); 

    //find the maximum value found for each element of all threads 
    auto& max_vals_found = max_vals[0]; 
    
    for (size_t t=1; t<n_threads; t++) 
    {
        for (int e=0; e<n_elems; e++) 
            max_vals_found[e].second = max<double>( max_vals_found[e].second, max_vals[t][e].second ); 
    }
    
    //now, return these max values found
    return max_vals_found; 
} 
