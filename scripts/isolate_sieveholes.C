////////////////////////////////////////////////////////////////////////////////
//
//  This is a macro which creates an application which allows the creation of training data 
//  by selecting & isolating sieve-hole data 
//
////////////////////////////////////////////////////////////////////////////////
#include <TColor.h> 
#include <TGClient.h> 
#include <isolate_sieveholes/PickSieveHoleApp.h> 

//_____________________________________________________________________________________________________________________________________
int isolate_sieveholes( const bool is_RHRS, 
                        const char* path_infile,
                        const char* target_name = "V2",
                        const char* coord_x ="dxdz_sv",
                        const char* coord_y ="dydz_sv",
                        const char* drawing_option ="col", 
                        const unsigned int palette =kSunset)
{
    new PickSieveHoleApp(gClient->GetRoot(), 
                         1400, 
                         700, 
                         is_RHRS,
                         path_infile, 
                         target_name, 
                         coord_x, 
                         coord_y, 
                         1, 
                         drawing_option,
                         palette);

    return 0; 
}
