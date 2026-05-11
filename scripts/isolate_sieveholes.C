////////////////////////////////////////////////////////////////////////////////
//
//  This is a macro which creates an application which allows the creation of training data 
//  by selecting & isolating sieve-hole data 
//
////////////////////////////////////////////////////////////////////////////////
#include <TColor.h> 
#include <TGClient.h> 
#include <isolate_sieveholes/PickSieveHoleApp.h> 
#include <ModularOpticsModel.h> 
#include "opticsModel/forward.h"

//_____________________________________________________________________________________________________________________________________
int isolate_sieveholes( const bool is_RHRS, 
                        const char* path_infile,
                        const char* target_name = "V2",
                        const char* path_poly="data/poly/fits_30Dec/V1-v05_fp_sv_R_4ord.dat",
                        const char* coord_x ="dxdz_sv",
                        const char* coord_y ="dydz_sv",
                        const char* drawing_option ="col", 
                        const unsigned int palette =kSunset)
{
    auto app = new PickSieveHoleApp(gClient->GetRoot(), 
                         1400, 
                         700, 
                         is_RHRS,
                         path_infile, 
                         target_name, 
                         coord_x, 
                         coord_y, 
                         drawing_option,
                         palette);

    ForwardModel model(path_poly); 

    app->SetOpticsModel(&model); 

    app->SetDrawRange_y(-0.0375, +0.0475); 

    app->LaunchApplication(); 

    return 0; 
}
