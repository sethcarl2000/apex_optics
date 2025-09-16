#ifndef SieveHole_h_
#define SieveHole_h_

///////////////////////////////////////////////
//
//  This defines a few structs: 
//
//  SieveHole
//  SieveHoleData
//  FPCoordPolynomial
//
///////////////////////////////////////////////

//this struct contains info about sieve holes
struct SieveHole {

    int row,col; 
    double x,y,radius_front,radius_back; 
    bool is_big; 

    //defining (overloading) this operator lets us use the std::find() function on a vector<SieveHole> 
    bool operator==(const SieveHole& rhs) const { return ((row==rhs.row) && (col==rhs.col)); }

    //defining this operator lets us use so-called 'ordered sets', like std::map, which can employ clever algorithms 
    // to search for a given element very qucikly
    bool operator<(const SieveHole& rhs) const { 
        if (row < rhs.row) return true;
        return col < rhs.col; 
    }; 
}; 

#endif 