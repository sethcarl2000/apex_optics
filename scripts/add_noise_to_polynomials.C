#include "TROOT.h"
#include "TRandom.h"
#include <vector> 
#include <string> 
#include <map> 
#include <NPoly.h>
#include <ApexOptics.h>
#include <NPolyArrayChain.h> 
#include <NPolyArray.h> 
#include <stdexcept> 
#include <TRandom3.h> 

using namespace std; 

int add_noise_to_polynomials(   const bool is_RHRS, 
                                const char* path_infile, 
                                const char* path_outfile, 
                                const double noise_mag=1e-3, 
                                const int noise_order=1     )  {

    const char* const here = "add_noise_to_polynomial"; 
    
    const vector<string> inputs{
        "x_fp",
        "y_fp",
        "dxdz_fp",
        "dydz_fp"
    }; 
    const size_t DoF_in = inputs.size(); 

    vector<string> outputs{
        "x_sv", 
        "y_sv",
        "dxdz_sv",
        "dydz_sv",
        "dpp_sv"
    };
    const size_t DoF_out = outputs.size(); 

    //try to parse the NPolyArray
    NPolyArray parr; 
    try {

        parr = ApexOptics::Parse_NPolyArray_from_file(path_infile, outputs, (int)inputs.size()); 
    
    } catch (const std::exception& e) {

        Error(here, "NPolyArray from file '%s' did not parse successfully.\n what(): %s", path_infile, e.what()); 
        return -1; 
    } 

    //now, add random noise. 
    TRandom3 rand; 
    auto Create_random_array = [noise_mag, noise_order, &rand](const int DoF)
    {   
        vector<NPoly> polys; 
        for (int i=0; i<DoF; i++) {

            NPoly poly(DoF, noise_order); 

            //set all element coefficients to random, gauss noise  
            for (int j=0; j<poly.Get_nElems(); j++) poly.Get_elem(j)->coeff = rand.Gaus() * noise_mag; 
            
            //Now, set all 'linear' elments to 1. in the limit noise_mag -> 0, this means that the polynomial
            ROOT::RVec<int> powers(DoF, 0); 
            powers[i] = 1; 

            poly.Find_element(powers)->coeff += 1.; 

            polys.push_back(poly); 
        }

        return NPolyArray(polys); 
    };

    //now, sandwich the input polynomial with a random NPolyArray on either side. 
    auto noisy_array1 = NPolyArray::Nest( parr, Create_random_array(DoF_in) ); 

    auto noisy_array = NPolyArray::Nest( Create_random_array(DoF_out), noisy_array1 ); 

    //now, save it in a file 
    map<string, NPoly*> polymap; 

    int i_pol=0; for (const auto& str : outputs) polymap[str] = noisy_array.Get_poly(i_pol++); 

    ApexOptics::Create_dbfile_from_polymap(is_RHRS, path_outfile, polymap);

    return 0; 
}