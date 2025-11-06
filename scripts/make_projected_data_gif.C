#include <TCanvas.h> 
#include <TPad.h> 
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <ApexOptics.h> 
#include <TRandom3.h> 
#include <TVector3.h> 
#include <Interactive3dHist.hxx>

using namespace std; 
using ApexOptics::Trajectory_t; 

int make_projected_data_gif(const char* path_infile)
{
    //make a gif of 'projecting' 
    ROOT::RDataFrame df("tracks_fp", path_infile);

    const double z_min = -750.e-3; 
    const double z_max = +750.e-3; 

    int stat_multiplier = 25; 

    //this is kinda funky, but its a way we can make use of the interactive hist method
    vector<Trajectory_t> vec_events = *df
        
        .Define("Xsv", [](double x, double y, double dxdz, double dydz) 
        {
            return Trajectory_t{x, y, dxdz, dydz}; 
        }, {"x_sv","y_sv","dxdz_sv","dydz_sv"})

        .Take<Trajectory_t>("Xsv"); 

    ROOT::RDataFrame df_proj((size_t)vec_events.size()*stat_multiplier); 

    TRandom3 rand; 

    size_t i_event=0; 
    const size_t n_events   = vec_events.size()*stat_multiplier; 
    
    auto df_draw = df_proj

        .Define("projected_point", [&]()
        {
            double z = z_min + (z_max - z_min)*rand.Rndm(); 
            
            Trajectory_t traj = vec_events.at(i_event++/stat_multiplier); 

            return TVector3(
                traj.x + traj.dxdz*z, 
                traj.y + traj.dydz*z, 
                z
            ); 
        }, {})

        .Define("x", [](TVector3 v){ return v.x(); }, {"projected_point"})
        .Define("y", [](TVector3 v){ return v.y(); }, {"projected_point"})
        .Define("z", [](TVector3 v){ return v.z(); }, {"projected_point"});

    new Interactive3dHist(df_draw, 
        {"x", 200, -0.06,0.06},
        {"y", 200, -0.06,0.06},
        {"z", 200, z_min,z_max}
    );

    return 0;  
}