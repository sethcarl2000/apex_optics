#include "TROOT.h"
#include <ROOT/RDataFrame.hxx> 
#include <vector> 
#include <iostream> 
#include <cstdio> 
#include <cmath> 

using namespace std; 
using namespace ROOT::VecOps; 

//contains data we need to take from the input file
struct HoleData_t {
    
    //sieve-coordinates (except momentum! we can't get this from v-wire data). 
    double x_sv, y_sv, dxdz_sv, dydz_sv; 

    //coefficients of the polynomial used to reconstruct these focal-plane coordinates (see usage below). 
    RVec<double> a_y_fp, a_dxdz_fp, a_dydz_fp; 

    double x_fp_min{-0.6}, x_fp_max{0.6}; 
}; 

int convert_holefit_data(size_t n_events_to_generate, const bool is_RHRS, const char* path_infile, const char* path_outfile) 
{
    const char* const here = "convert_holefit_data"; 

    //converts data produced by the file 'make_data_from_holefit.C' macro to a format which is friendly to the polynomial-fitting
    //macros which are called 'fitpoints_[...].C'. and ApexOptics::Create_NPoly_fit

    //for this first one, run in sequential mode. 
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df_get_vector("hole_data", path_infile); 

    cout << "Parsing holes from input file..." << flush; 

    vector<HoleData_t> hole_data = *df_get_vector 

        .Define("holedata_vec", []( double x_sv, 
                                    double y_sv, 
                                    double dxdz_sv, 
                                    double dydz_sv, 
                                    RVec<double> a_y_fp, 
                                    RVec<double> a_dxdz_fp,
                                    RVec<double> a_dydz_fp,// ) 
                                    double x_fp_min, 
                                    double x_fp_max )
        {
            return HoleData_t{
                .x_sv       = x_sv, 
                .y_sv       = y_sv, 
                .dxdz_sv    = dxdz_sv, 
                .dydz_sv    = dydz_sv, 
                .a_y_fp     = a_y_fp,
                .a_dxdz_fp  = a_dxdz_fp,
                .a_dydz_fp  = a_dydz_fp,
                .x_fp_min   = x_fp_min,
                .x_fp_max   = x_fp_max 
            }; 
        }, {"x_sv", "y_sv", "dxdz_sv", "dydz_sv", "a_y_fp", "a_dxdz_fp", "a_dydz_fp", "x_fp_min", "x_fp_max"})
        
        //this 'Take' command tells the RDataFrame to put each event in the column 'holedata_vec' into a vector, 
        // and hand that vector back to us. 
        .Take<HoleData_t>("holedata_vec"); 
    
    printf("done. %zi holes parsed from input file.", hole_data.size()); cout << endl;  

    
    //create a random number generator which we will use to pick random x_fp in the range given below
    TRandom3 rand; 

    //now, we proceed with actually making our polynomial-fit data. 
    const double x_fp_min = -0.600; 
    const double x_fp_max = +0.600; 
   
    //decide how many events to generate. 
    const size_t n_events_per_hole = n_events_to_generate / hole_data.size() + 1; 
    
    n_events_to_generate = n_events_per_hole * hole_data.size(); 
    //we will generate the nearest multiple of the number of holes rounded up. 
   
        
    //make sure we run sequentially (no parallel processing)
    ROOT::DisableImplicitMT(); 

    size_t n_events_this_hole=0; 
    auto hole_data_iterator     = hole_data.begin(); 
    auto beyond_last_element    = hole_data.end(); 

    ROOT::RDataFrame df(n_events_to_generate); 
    
    //now, we actually generate all the events we want. 
    cout << "Generating " << n_events_to_generate << " events..." << flush; 
        
    auto df_output = df 

        .Define("hole_data", [&hole_data_iterator, beyond_last_element, &n_events_this_hole, n_events_per_hole]()
        {
            //check if we've gone out-of-range
            if (hole_data_iterator == beyond_last_element) {
                throw logic_error(  "in <RDF::Define('holedata')>: hole_data_iterator tried to access data beyond "
                                    "the last element." ); 
                return HoleData_t{};
            }

            //Get the 'hole-data' struct. 
            HoleData_t hole_data = *hole_data_iterator;

            //if we've made enough events for this hole, then move on to the next one. 
            if (++n_events_this_hole >= n_events_per_hole) {
                hole_data_iterator++;
                n_events_this_hole = 0;
            }

            return hole_data; 
        }, {})

        //define the sieve-coordinates from the hole data. 
        .Define("x_sv",    [](const HoleData_t& hole_data){ return hole_data.x_sv; },       {"hole_data"})
        .Define("y_sv",    [](const HoleData_t& hole_data){ return hole_data.y_sv; },       {"hole_data"})
        .Define("dxdz_sv", [](const HoleData_t& hole_data){ return hole_data.dxdz_sv; },    {"hole_data"})
        .Define("dydz_sv", [](const HoleData_t& hole_data){ return hole_data.dydz_sv; },    {"hole_data"})

        //generate x_fp uniformly in the range defined above. 
        .Define("x_fp", [x_fp_min, x_fp_max, &rand](const HoleData_t& hole_data)
        {
            return hole_data.x_fp_min + (hole_data.x_fp_max - hole_data.x_fp_min) * rand.Rndm(); 
        }, {"hole_data"})

        //this last one is a dummy model (which is not too far from accurate), because most of the 'polynomial-fit' 
        // scripts are expecting this branch, but we don't have this data from vertical-wire runs. 
        .Define("dpp_sv",  [x_fp_max](double x_fp){ return x_fp * ( 0.04 / x_fp_max ); },   {"x_fp"})

        //these coordinates are defined by polyomials, with x_fp as an input. 
        .Define("y_fp", [](const HoleData_t& hole_data, double x_fp){

            double ret =0.; 
            const RVec<double>& coeffs = hole_data.a_y_fp; 
            for (size_t i=0; i<coeffs.size(); i++) ret += coeffs[i] * pow( x_fp, i ); 
            return ret; 

        }, {"hole_data", "x_fp"})

        .Define("dxdz_fp", [](const HoleData_t& hole_data, double x_fp){

            double ret =0.; 
            const RVec<double>& coeffs = hole_data.a_dxdz_fp; 
            for (size_t i=0; i<coeffs.size(); i++) ret += coeffs[i] * pow( x_fp, i ); 
            return ret; 

        }, {"hole_data", "x_fp"})

        .Define("dydz_fp", [](const HoleData_t& hole_data, double x_fp){

            double ret =0.; 
            const RVec<double>& coeffs = hole_data.a_dydz_fp; 
            for (size_t i=0; i<coeffs.size(); i++) ret += coeffs[i] * pow( x_fp, i ); 
            return ret; 

        }, {"hole_data", "x_fp"})

        .Snapshot("tracks_fp", path_outfile, {
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "dpp_sv",
            
            "x_fp",
            "y_fp",
            "dxdz_fp",
            "dydz_fp"
        }); 

    cout << "done." << endl; 


    //add 'is_RHRS' parameter
    cout << "Adding parameters..." << flush; 
    auto file = new TFile(path_outfile, "UPDATE"); 

    auto param_is_RHRS = new TParameter<bool>("is_RHRS", is_RHRS); 
    
    auto param_x_fp_min = new TParameter<double>("x_fp_min", x_fp_min); 
    auto param_x_fp_max = new TParameter<double>("x_fp_max", x_fp_max); 

    auto param_n_sieve_holes = new TParameter<int>("n_sieve_holes", (int)hole_data.size()); 
    
    param_is_RHRS       ->Write();
    param_x_fp_min      ->Write(); 
    param_x_fp_max      ->Write(); 
    param_n_sieve_holes ->Write(); 

    file->Close(); 
    delete file; 

    cout << "done." << endl; 
    return 0; 
}   
