
#include <map>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm> 
#include "TROOT.h"

using namespace std; 


struct SieveHole {
    
    int row,col; 
    double x,y,radius_front,radius_back;   
    
    //defining (overloading) this operator lets us use the std::find() function on a vector<SieveHole> 
    bool operator==(const SieveHole& rhs) { return ((row==rhs.row) && (col==rhs.col)); }

    static std::vector<SieveHole> Create_sieve_holes(bool is_RHRS=false) {
        return {}; 
    } 

}; 

int test_sieve_hole_construction(bool is_RHRS=true)
{

    //create RHRS sieve-holes
    vector<SieveHole> sieve_holes; 

    const double dx = 5.842e-3; //all units in meters

    //the sign is different here, because the L & R-sieves are mirror-images of
    // of one another
    const double dy = is_RHRS ? 4.826e-3 : -4.826e-3;  //
    
    //the first row is 8 rows above (-x) the center hole
    const double x0 = -dx * 8.; 
    //the first column is 7-spaces +y from the center
    const double y0 =  dy * 7.; 

    //two possible hole radii
    const double holeR_small = 0.6985e-3; 
    const double holeR_big   = 1.3462e-3; 

    const double holeR_widened = 0.9525e-3; //radius of the holes whose exits are widened. 

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

        }//for (int col=0; col<nCols; col++) 
    }//for (int row=0; row<nRows; row++) 
    //_______________________________________________________________________________________________________________

    new TCanvas; 
    gStyle->SetOptStat(0); 

    string h_title = (is_RHRS ? "Sieve Hole drawing - RHRS;x_sv;y_sv" : "Sieve Hole drawing - LHRS;x_sv;y_sv"); 

    auto hist = new TH2D("h", h_title.data(), 200, -60e-3, 60e-3, 200, -60e-3, 60e-3); 
    hist->Draw(); 

    for (const SieveHole& hole : sieve_holes) {

        auto circ_f = new TEllipse( hole.x, hole.y, hole.radius_front, hole.radius_front ); 
        circ_f->Draw(); 

        auto circ_b = new TEllipse( hole.x, hole.y, hole.radius_back, hole.radius_back ); 
        circ_b->Draw();  

    }

    return 0; 
}