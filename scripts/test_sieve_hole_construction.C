
#include <map>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm> 
#include <random>
#include "TROOT.h"

using namespace std; 

namespace ApexTargetGeometry {
  
  inline double Get_sieve_angle(bool _is_RHRS) {
    return ( _is_RHRS ? -5.372 : 5.366 ) * TMath::Pi() / 180.;
  }
  
  inline double Get_HRS_angle(bool _is_RHRS) {
    return ( _is_RHRS ? -12.50 : 12.50 ) * TMath::Pi() / 180.;
  }

  //units in mm. 
  inline TVector3 Get_APEX_Target_center() {
    return TVector3( 0., 0., -1053.7952 );
  }
    
  //units in mm. These are in Target coordinates (TCS), obtained by rotating hall
  // coordinates first by '-sieve_angle' about the y-axis, then by pi/2 about the
  // z-axis. 
  inline TVector3 Get_sieve_pos(bool _is_RHRS) {
    return TVector3(  _is_RHRS ?  -1.101 :  -1.301,
                      _is_RHRS ?  -3.885 :   6.672,
                      _is_RHRS ? 794.609 : 795.766 ); 
  } 
  
}; 


struct SieveHole {
    
    int row,col; 
    double x,y,radius_front,radius_back;   
    
    //defining (overloading) this operator lets us use the std::find() function on a vector<SieveHole> 
    bool operator==(const SieveHole& rhs) { return ((row==rhs.row) && (col==rhs.col)); }

    static std::vector<SieveHole> Create_sieve_holes(bool is_RHRS=false) {
        return {}; 
    } 

}; 

struct Track_t {
    TVector3 pos,mom; 
};

int test_sieve_hole_construction(const bool is_RHRS=true, const int vwire_num=2, const int n_events=1e5)
{   
    const char* const here = "test_sieve_hole_construction"; 

    //create a random number generator
    std::random_device rd;
    auto twister   = mt19937(rd());  
    
    auto real_dist = std::uniform_real_distribution<double>(0., 1.); 
    auto Get_rnd_range = [&real_dist,&twister](double min, double max)
    { 
        return min + (max-min)*real_dist(twister); 
    }; 

    vector<SieveHole> big_holes; 

    //create RHRS sieve-holes
    vector<SieveHole> sieve_holes; 

    const double dx = 5.842; //all units in mm

    //the sign is different here, because the L & R-sieves are mirror-images of
    // of one another
    const double dy = is_RHRS ? 4.826 : -4.826;  //
    
    //the first row is 8 rows above (-x) the center hole
    const double x0 = -dx * 8.; 
    //the first column is 7-spaces +y from the center
    const double y0 =  dy * 7.; 

    //sieve thickness
    const double sieve_thickness = 12.7; 

    //two possible hole radii
    const double holeR_small = 0.6985; 
    const double holeR_big   = 1.3462; 

    const double holeR_widened = 0.9525; //radius of the holes whose exits are widened. 

    const int nRows=17; 
    for (int row=0; row<nRows; row++) { 
        //17 rows in total, numbered 0-16.
        // row 0 is the highest in HCS, so lowest-X in TCS.
        // Recall that in TCS, the central 'big hole' is the origin of x & y.
        
        //even rows are shifted in +y half a col-spacing
        bool evenRow = (row % 2==0); 
        int nCols;
        if (evenRow) { nCols = 15; }
        else  {
            if (row==1 || row==15) {nCols = 12;} //each odd row has a 'gap' missing, except these two
            else                   {nCols = 11;}
        } 
        
        for (int col=0; col<nCols; col++) { 
        
            //check wheter this hole is a 'big' hole
            bool bigHole ((row==8 && col==7) || (row==12 && col==3)); 

            //get this hole's x-position
            double x = x0 + ((double)row)*dx;
            double y = y0 - ((double)col)*dy;

            if (!evenRow) { 
                y += -dy/2.; //odd-holes are shifted in y a bit
                //skip over the gap which happens in some rows, but not these. 
                if ( (row!=1 && row!=15) && col>5 ) y += -dy; 
            }
            
            //we're ready to define our new hole. 
            SieveHole new_hole{
                .row    = row,
                .col    = col,
                .x      = x,
                .y      = y,
                .radius_front   = (bigHole ? holeR_big : holeR_small), 
                .radius_back    = (bigHole ? holeR_big : holeR_small)
            }; 

            //check to see if this is one of the holes where the exit-hole is wider than the entrance-hole. 
            //this is true of the top-3 and bottom-3 rows, but on those rows, it is NOT so for the lateral-most 3 holes. 
            if (row <= 2 || row >= 14) {
                if (col <= 12) {
                    new_hole.radius_back = holeR_widened; 
                }
            }

            sieve_holes.push_back(new_hole); 
            
            if (bigHole) { big_holes.push_back(new_hole); }

        }//for (int col=0; col<nCols; col++) 
    }//for (int row=0; row<nRows; row++) 
    //_______________________________________________________________________________________________________________

    //test the generation of particles at random positions
    vector<TVector3> vwire_positions{
        TVector3( -3.225,    0.,    -196.214  ),
        TVector3( -0.725,    0.,       3.786  ),
        TVector3(  1.725,    0.,     203.786  )
    }; 

    if (vwire_num < 1 &&  vwire_num > 3) {
        Error(here, "Invalid V-wire number: %i, can only be 1->3", vwire_num); 
        return -1; 
    }
    
    TVector3 vertex = vwire_positions.at(vwire_num-1); 
    
    //change the wire pos to sieve coordinates
    vertex.RotateY( ApexTargetGeometry::Get_sieve_angle(is_RHRS) ); 
    vertex.RotateZ( TMath::Pi()/2. ); 

    //add the offsets of the sieve-plane center
    vertex += -1.* ApexTargetGeometry::Get_sieve_pos(is_RHRS); 
 

    const vector<SieveHole> *holes = &sieve_holes; 

    
    //area ratio of 'big hole' to 'small holes'. 
    //so... the large holes are **approximately** 4 times the area of the other holes. so I will give each 'big' hole 3+ 'ballots' in the raffle of which holes to populate. 
    //the gist of this subroutine is that we want to 'throw' tracks at each hole, in an amount which is proportional to the opening-face
    //of the hole. to do this efficiently, we choose all holes randomly, by having an 'index' generated over a random distribution. 
    // however, we must account for the fact that there are two 'big' holes. To account for this, we give each of the big holes an 
    // extra 'weight', which amounts to an extra chance at selection, each time a random hole is picked. 
    double big_hole_extra_area = (holeR_big*holeR_big)/(holeR_small*holeR_small) - 1.; 

    auto Random_hole = [holes, &big_holes, &Get_rnd_range, big_hole_extra_area]() 
    {
        double index = Get_rnd_range( -2.*big_hole_extra_area, (double)holes->size()-1 ); 
        
        if (index < -big_hole_extra_area) return big_holes.at(0); 
        if (index < 0. )                  return big_holes.at(1); 

        return holes->at( (int)index ); 
    }; 


    ROOT::EnableImplicitMT(); 
    printf("Multi-threadding enabled. Thread pool size: %i\n", ROOT::GetThreadPoolSize() ); 

    printf("Generating %.2e tracks...", (double)n_events); cout << endl;


    ROOT::RDataFrame df(n_events);

    auto df_output = df 

        .Define("Xsv",     [Random_hole, Get_rnd_range, vertex, sieve_thickness]()
        {   
            TVector3 pos_sieve, momentum; 

            //keep looping until we find a 'good' target trajectory
            bool found=false; 

            while (!found) {
                
                SieveHole hole = Random_hole(); 
                
                int tries=0; 

                while (1) { //choose a random spot to 'aim' at
                    
                    pos_sieve = TVector3(
                        Get_rnd_range( -hole.radius_front, hole.radius_front ),
                        Get_rnd_range( -hole.radius_front, hole.radius_front ),
                        0. 
                    ); 

                    //check if this 'spot' is within the radius of the hole 
                    if ( pos_sieve.Mag() < hole.radius_front ) break; 
                }

                double x_exit = pos_sieve.x();
                double y_exit = pos_sieve.y(); 

                pos_sieve[0] += hole.x; 
                pos_sieve[1] += hole.y; 

                momentum = (pos_sieve - vertex).Unit(); 
                
                //now, see if this track can 'thread the needle' from the entrance hole, out of the exit hole. 
                x_exit += (momentum.x()/momentum.z()) * sieve_thickness; 
                y_exit += (momentum.y()/momentum.z()) * sieve_thickness; 
                
                //If the hole hits the wall of the sieve hole before it exists, then throw this one out. 
                if ( sqrt(x_exit*x_exit + y_exit*y_exit) < hole.radius_back ) break; 

                //while (tries++ < 25)
            
            }//while (!found)
                
            

            Track_t track_sv{
                .pos = pos_sieve, 
                .mom = momentum
            }; 

            return track_sv; 
        }) 

        .Define("x_sv", [](Track_t &trk){ return trk.pos.x(); }, {"Xsv"})
        .Define("y_sv", [](Track_t &trk){ return trk.pos.y(); }, {"Xsv"});

    
    string h_title = (is_RHRS ? "Sieve Hole drawing - RHRS" : "Sieve Hole drawing - LHRS"); 
    
    char buff[50]; sprintf(buff, " (V%i)", vwire_num);
    h_title += string(buff); 

    h_title += ";x_sv;y_sv"; 
    

    auto hist = df_output.Histo2D({"h", h_title.data(), 300, -50, 50, 300, -50, 50}, "x_sv", "y_sv"); 
    
    gStyle->SetPalette(kSunset); 
    
    hist->DrawCopy("col2"); 
    
    cout << "done." << endl; 

    //uncomment this to draw all sieve holes. 
    /*new TCanvas; 
    gStyle->SetOptStat(0); 

    string h_title = (is_RHRS ? "Sieve Hole drawing - RHRS;x_sv;y_sv" : "Sieve Hole drawing - LHRS;x_sv;y_sv"); 

    auto hist = new TH2D("h", h_title.data(), 200, -60e-3, 60e-3, 200, -60e-3, 60e-3); 
    hist->Draw(); 

    for (const SieveHole& hole : sieve_holes) {

        auto circ_f = new TEllipse( hole.x, hole.y, hole.radius_front, hole.radius_front ); 
        circ_f->Draw(); 

        auto circ_b = new TEllipse( hole.x, hole.y, hole.radius_back, hole.radius_back ); 
        circ_b->Draw();  
    }*/ 

    return 0; 
}