#include <ApexOptics.h>
#include <TVector3.h> 
#include <cstdio> 
#include <iostream> 

void examine_coord_diff()
{  
    bool is_RHRS=false; 

    TVector3 x0 = ApexOptics::HCS_to_SCS(is_RHRS, TVector3(0., 0., 0.)); 

    TVector3 ux = ApexOptics::HCS_to_SCS(is_RHRS, TVector3(1., 0., 0.)); 
    TVector3 uy = ApexOptics::HCS_to_SCS(is_RHRS, TVector3(0., 1., 0.)); 
    TVector3 uz = ApexOptics::HCS_to_SCS(is_RHRS, TVector3(0., 0., 1.)); 

    //transformation matrix: 
    printf("trans. matrix:\n");
    printf("%+.8f & %+.8f & %+.8f \\\\\n", (ux - x0)[0], (uy - x0)[0], (uz - x0)[0]);
    printf("%+.8f & %+.8f & %+.8f \\\\\n", (ux - x0)[1], (uy - x0)[1], (uz - x0)[1]); 
    printf("%+.8f & %+.8f & %+.8f \n", (ux - x0)[2], (uy - x0)[2], (uz - x0)[2]);  

    printf("const. offset vector:\n"); 
    printf("%+.8f \\\\ \n%+.8f \\\\ \n%+.8f\n", 1e3*x0[0], 1e3*x0[1], 1e3*x0[2]); 
}