#include "TROOT.h"
#include <vector>
#include <regex> 
#include <utility> 
#include <sstream> 
#include <string> 
#include <iostream> 
#include <algorithm> 
#include <limits> 
#include <cstdio> 
#include <fstream> 
#include <iostream> 

using namespace std; 
using namespace ROOT::VecOps; 

struct SieveHole {

    int row,col; 
    double x,y,radius_front,radius_back; 
    bool is_big; 

    //defining (overloading) this operator lets us use the std::find() function on a vector<SieveHole> 
    bool operator==(const SieveHole& rhs) { return ((row==rhs.row) && (col==rhs.col)); }
}; 

//constructs a container of SieveHole structs with accurate positions, row/column indices, and sizes. 
// all units in mm for sieve hole positions and sizes. 
vector<SieveHole> Construct_sieve_holes(bool is_RHRS) 
{
    vector<SieveHole> sieve_holes; 

    const double dx = 5.842e-3; //all units in meters

    //the sign is different here, because the L & R-sieves are mirror-images of
    // of one another
    const double dy = is_RHRS ? 4.826e-3 : -4.826e-3;
    
    //the first row is 8 rows above (-x) the center hole
    const double x0 = -dx * 8.; 
    //the first column is 7-spaces +y from the center
    const double y0 =  dy * 7.; 

    //sieve thickness
    const double sieve_thickness = 12.7e-3; 

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
                .radius_back    = (bigHole ? holeR_big : holeR_small),
                .is_big = bigHole
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
    
    return sieve_holes; 
}


int make_data_from_holefit(const bool is_RHRS, const char* path_infile, const char* path_outfile)
{   
    const char* const here = "make_data_from_holefit"; 

    //each fp-coordinate polynomal will have polynomial_order+1 coefficients. 
    const int polynomial_order=3; 

    //I will pass a nan-int as row/col for invalid hole names. 
    const int nan_int = numeric_limits<int>::quiet_NaN(); 

    //parse all holes in file
    vector<SieveHole> holes = Construct_sieve_holes(is_RHRS); 

    auto Get_sieveHole = [&holes, nan_int](string hole_name) 
    {
        int row, col; 

        //thanks to claude for showing me how to use the std::regex library
        regex pattern("r(\\d+)c(\\d+)"); 
        smatch matches; 

        //get the row and column of the sieve hole. 
        if (regex_match(hole_name, matches, pattern)) {
            row = std::stoi(matches[1].str()); 
            col = std::stoi(matches[2].str()); 
        } else {
            return SieveHole{.row=-1, .col=-1}; 
        }

        //this function returns 'true' if a sieve hole has the correct row & column we just parsed from the .dat file
        auto check_row_col_match = [row,col](const SieveHole& hole) {
            return (hole.row==row && hole.col==col); 
        }; 

        //check our collection of sieve holes to see if we can find the one with this name
        auto it = std::find_if( holes.begin(), holes.end(), check_row_col_match ); 
        if (it == holes.end()) {
            Error("Get_sieveHole", "Following hole parsed from file is not a valid sieve-hole: %s", hole_name.c_str()); 
            return SieveHole{.row=-1, .col=-1}; 
        }

        return *it; 
    }; 

    //now, we can parse the sieve-holes from the file. 
    
    //first row of the file is the react-vertex. 
    ifstream infile(path_infile); 
    if (!infile.is_open()) {
        Error(here, "File '%s' could not be opened.", path_infile); 
        return -1; 
    }


    string line, token; 
    istringstream iss; 

    //line 1: should be the token 'react-vertex', followed by 3 doubles. 
    if (!getline(infile, line)) {
        Error(here, "Unexpected end-of-file reached while parsing line 'react-vertex'");
        return -1; 
    }
    iss = istringstream(line); 
    TVector3 position_vtx; 
    iss >> token >> position_vtx[0] >> position_vtx[1] >> position_vtx[2]; 
    if (token != "react-vertex") {
        Error(here, "token 'react-vertex' missing from file: %s", path_infile); 
        return -1; 
    }

    //line 2: should be the token 'rast-span', followed by 3 doubles. 
    if (!getline(infile, line)) {
        Error(here, "Unexpected end-of-file reached while parsing line 'rast-span'");
        return -1; 
    }
    iss = istringstream(line); 
    TVector3 rast_span; 
    iss >> token >> rast_span[0] >> rast_span[1] >> rast_span[2]; 
    if (token != "rast-span") {
        Error(here, "token 'rast-span' missing from file: %s", path_infile); 
        return -1; 
    }


    //stores the sieve-hole, with fpcoord
    struct SieveHole_with_fpcoord_t {
        SieveHole hole;
        RVec<double> a_y_fp, a_dxdz_fp, a_dydz_fp;  
    };
    vector<SieveHole_with_fpcoord_t> holes_with_fpcoord{}; 

    //now, all remaining lines have a specific format, which will be detailed below. 
    int line_num=2; 
    while (getline(infile, line)) {
        
        line_num++; 

        iss = istringstream(line); 

        string hole_name; 
        iss >> hole_name; 

        //if we've reached the last line, quit. 
        if (hole_name=="eof") break; 
        
        //get the sieve-hole associated with this line. 
        SieveHole hole = Get_sieveHole(hole_name); 

        if (hole.row==-1) {
            //this happens if 'hole_name' is not a valid sieve-hole name. 
            Error(here, "invalid hole-name encountered: '%s'.\n In file '%s', line %i", hole_name.c_str(), path_infile, line_num);
            return -1;  
        }

        //so, a vaild sieve-hole has been found. let's proceed: 
        double x_sv, y_sv; 

        //the next two tokens are the x & y coordinates of the sieve-hole. 
        iss >> x_sv; 
        iss >> y_sv; 

        //lets parse the polynomial for each fp-coordinate.
        SieveHole_with_fpcoord_t hole_data{
            .hole       = hole, 
            .a_y_fp     = RVec<double>(polynomial_order+1, 0.),
            .a_dxdz_fp  = RVec<double>(polynomial_order+1, 0.),
            .a_dydz_fp  = RVec<double>(polynomial_order+1, 0.)
        }; 

        for (double& a : hole_data.a_y_fp)    iss >> a; 
        for (double& a : hole_data.a_dxdz_fp) iss >> a; 
        for (double& a : hole_data.a_dydz_fp) iss >> a; 
        
        holes_with_fpcoord.push_back(hole_data); 
    }

    printf("%zi holes parsed from file '%s'.\n.", holes_with_fpcoord.size(), path_infile); 
    infile.close(); 

    //now, we are ready to make our own data. 

    //this ensures that the RDataFrame will operate in sequential mode, so it will NOT try to process the loop in parallel. 
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
    
    //create our Dataframe. 
    ROOT::RDataFrame df(holes_with_fpcoord.size()); 

    int i_elem=0; 


    //There's a subtle point here:
    // The react-vertex is given to us in a coordinate system
    // whose axes are correctly orentied along the sieve-coordinate system axes, but its origin is given relative to the hall-coordinate
    // system origin. Displacements in Sieve-coordinates are given relative to the sieve's central hole, so we must subtract the 
    // centeral sieve-hole's coordinates, to properly convert the coordinate below to sieve-coordinates. 
    NPoly dummy(4,2); 

    position_vtx += -ApexOptics::Get_sieve_pos(is_RHRS); 

    auto df_out = df

        //position vertex, in scs. 
        .Define("position_vtx_scs", [&position_vtx](){ return position_vtx; }, {})

        .Define("raster_rms_scs", [&rast_span](){ return rast_span; }, {})

        .Define("hole_data", [&holes_with_fpcoord, &i_elem]()
        {
            SieveHole_with_fpcoord_t hole_data = holes_with_fpcoord[i_elem++]; 
            return hole_data; 
        }, {})

        .Define("hole_row", [](const SieveHole_with_fpcoord_t& data){ return data.hole.row; }, {"hole_data"})
        .Define("hole_col", [](const SieveHole_with_fpcoord_t& data){ return data.hole.col; }, {"hole_data"})

        .Define("x_sv",     [](const SieveHole_with_fpcoord_t& data){ return data.hole.x; },   {"hole_data"})
        
        .Define("y_sv",     [](const SieveHole_with_fpcoord_t& data){ return data.hole.y; },   {"hole_data"})
        
        .Define("dxdz_sv",  [](const SieveHole_with_fpcoord_t& data, TVector3 vtx)
        { 
            //the '0. - vtx.z()' here is the 'dz', as z=0 in sieve-coordinates is the face of the sieve. 
            return ( data.hole.x - vtx.x() )/( 0. - vtx.z() ); 
        }, {"hole_data", "position_vtx_scs"})

        .Define("dydz_sv",  [](const SieveHole_with_fpcoord_t& data, TVector3 vtx)
        { 
            //the '0. - vtx.z()' here is the 'dz', as z=0 in sieve-coordinates is the face of the sieve. 
            return ( data.hole.y - vtx.y() )/( 0. - vtx.z() ); 
        }, {"hole_data", "position_vtx_scs"})
        
        .Define("a_y_fp",    [](const SieveHole_with_fpcoord_t& data){ return data.a_y_fp; },    {"hole_data"})
        .Define("a_dxdz_fp", [](const SieveHole_with_fpcoord_t& data){ return data.a_dxdz_fp; }, {"hole_data"})
        .Define("a_dydz_fp", [](const SieveHole_with_fpcoord_t& data){ return data.a_dydz_fp; }, {"hole_data"})

        .Snapshot("hole_data", path_outfile, {
            "position_vtx_scs",
            "raster_rms_scs", 
            "hole_row", 
            "hole_col",
            "x_sv",
            "y_sv",
            "dxdz_sv",
            "dydz_sv",
            "a_y_fp",
            "a_dxdz_fp",
            "a_dydz_fp"
        }); 

    return 0; 
}