#include "PolynomialCut.h"
#include "PolynomialCutApp.h"
#include "TGClient.h"
#include <sstream> 
#include <string> 
#include <stdexcept> 
#include <iostream> 
#include <fstream> 
#include <cstdio> 
#include <utility> 

using namespace std; 


//_____________________________________________________________________________________________________________________
bool PolynomialCut::IsInside(double x, double y) const
{
    //check to make sure this PolynomialCut object has a valid status before proceeding. 
    if (GetStatus() != kGood) {

        const char* status;
        switch (GetStatus()) {
            case kError : status="kError"; break; 
            case kNot_init : status="kNot_init"; break; 
        } 
        ostringstream oss; 
        oss << "in <PolynomialCut::IsInside>: Status of PolynomialCut object is '" << status << "'." << endl; 
        throw logic_error(oss.str()); 
        return false; 
    }
    //now, we will raycast to each of the segments in the 'segment' vector. 
    //we will cast a ray upward, as that is the easiest to compute. 
    int n_segments_match=0; 

    for (const Segment_t& seg : fSegments) {

        //first, check x 
        double t = ( x - seg.x0 ) / seg.sx; 

        if (t < 0. || t >= 1.) continue; //does not match this segment (x out of range)

        //now, check y 
        double y_segment = seg.y0 + seg.sy * t; 

        if (y > y_segment) continue; //does not match this segment (y out of range) 
        
        n_segments_match++; 
    }   
    //if the number of intersections with all segments is even, then the point is outside. otherwise, it is inside. 
    if (n_segments_match % 2 == 0) { return false; } else { return true; }
}
//_____________________________________________________________________________________________________________________
void PolynomialCut::AddVertices(const vector<Vertex_t>& vertices)
{
    std::copy( vertices.begin(), vertices.end(), std::back_inserter(fVertices) ); 
    InitSegments(); 
}
//_____________________________________________________________________________________________________________________
void PolynomialCut::InitSegments() 
{
    //this returns 'true' if any segments collide with one another, and false if not.
    const char* const here = "PolynomialCut::InitSegments"; 

    if (fVertices.size() < 3) {
        return; 
    }
    fSegments.clear(); fSegments.reserve(fVertices.size()); 

    for (size_t i=0; i<fVertices.size()-1; i++) fSegments.push_back(MakeSegment(fVertices[i], fVertices[i+1])); 
    
    //add one last segment to complete the loop, by making a segment connecting the last -> first segments 
    fSegments.push_back(MakeSegment(fVertices.back(), fVertices.front())); 

    if (fVertices.size() <= 3 && GetStatus() == PolynomialCut::kNot_init) {
        fStatus = PolynomialCut::kGood; 
        return; 
    }

    //now, check each vertex to see if they collide with any other. 
    //first check the segment between the last vertex and the new vertex 
    for (size_t i=0; i<fSegments.size(); i++) {
        Segment_t seg_old = fSegments[i]; 

        for (size_t j=i+1; j<fSegments.size(); j++) {
            Segment_t seg_new = fSegments[j]; 

            //if this new segment intersects with any of the pre-existing segments, then throw an exception. 
            if (DoSegmentsIntersect(seg_old, seg_new)) {

                ostringstream oss; 

                oss <<  "in <" << here << ">: Vertices "
                        "(vert-index:["<< j <<","<< (j+1==fSegments.size() ? 0 : j+1) <<"])"
                        " collide with vertices"
                        "(vert-index:["<< i <<","<< i+1 <<"])";  

                fStatus = kError; 
                throw InvalidVertexException(oss.str()); 
                return; 
            }   
        }
    }

    if (GetStatus() == PolynomialCut::kNot_init) fStatus = PolynomialCut::kGood; 
}
//_____________________________________________________________________________________________________________________
bool PolynomialCut::DoSegmentsIntersect(const Segment_t& s1, const Segment_t& s2) const 
{
    //To solve this linear system: 
    //
    //      x1 + s1x * t1 = x2 + s2x * t2
    //      y1 + s1y * t1 = y2 + s2y * t2 

    //which can be witten as: 
    // 
    //      [ s1x  -s2x ] [ t1 ]  =  [ x2 - x1 ]
    //      [ s1y  -s2y ] [ t2 ]  =  [ y2 - y1 ]

    //or, 
    //      [ a  b ] [ t1 ] = [ x ]
    //      [ c  d ] [ t2 ] = [ y ]
    
    //               [ t1 ] = [ a  b ]^-1 [ x ]
    //               [ t2 ] = [ c  d ]    [ y ]

    //compute the point at which these two segemnts intersect 
    double a =  s1.sx; 
    double b = -s2.sx; 
    double c =  s1.sy; 
    double d = -s2.sy; 

    double det = a*d - b*c; 

    //if the determinant is zero (or NaN), then the segments are parallel 
    if ( det != det || fabs(det) < 1e-20 ) return false; 

    double x = s2.x0 - s1.x0; 
    double y = s2.y0 - s1.y0; 
    
    //compute 't'
    double t1 = (  d*x - b*y )/det; 
    double t2 = ( -c*x + a*y )/det; 

    //now, the way the segments are designed, the segments are only 'valid' with t = [0,1). 
    if (t1 < 0. || t1 >= 1. ||
        t2 < 0. || t2 >= 1.) return false; 

    return true; 
}
//_____________________________________________________________________________________________________________________
void PolynomialCut::Create_dbfile(const char* path_outfile) const
{
    const char* const here = "PolynomialCut::Create_dbfile"; 

    //create a new output file. if one exists, then delete all contents and overwrite it. 
    fstream file(path_outfile, ios::out | ios::trunc); 

    //check if file could be opened
    if (!file.is_open()) { 
        ostringstream oss; 
        oss << "in <" << here << ">: unable to open output file"; 
        throw DBFileException(oss.str().c_str(), path_outfile); 
        return;  
    }

    if (GetStatus() != kGood) {
        ostringstream oss; 
        string status="undefined status"; 
        switch (GetStatus()) {
            case kError    : status="kError"; break; 
            case kNot_init : status="kNot_init"; break; 
        }
        oss << "in <" << here << ">: Tried to write PolynomialCut object with status '" << status << "'"; 
        throw logic_error(oss.str()); 
    }

    //add the preamble to the dbfile
    file << "# PolynomialCut vertex file. the data format is space-delimited, and is as follows:\n"
            "#\n"
            "# [vertex-index] [vertex-x] [vertex-y]\n"
            "#\n"
            "# note that any lines startig with '#' will be ignored, so you can add comments to indicate\n"
            "# the circumstances under which this data was created, and how it should be used\n"
            "# (like which branches 'x' and 'y' correspond to.)\n"; 

    int i=0; 
    for (const Vertex_t& vertex : fVertices) {

        char line[200]; sprintf(line, "%-3i %+.9e %+.9e\n", i++, vertex.x, vertex.y); 
        file << line; 
    }

    file.close(); 
}
//_____________________________________________________________________________________________________________________
void PolynomialCut::Parse_dbfile(const char* path_infile)
{
    const char* const here = "PolynomialCut::Parse_dbfile"; 

    //delete all pre-existing vertices
    fVertices.clear(); 

    ifstream infile(path_infile); 
    if (!infile.is_open()) {
        ostringstream oss; 
        oss << "in <" << here << ">: unable to open input file"; 
        throw DBFileException(oss.str().c_str(), path_infile); 
        return; 
    }

    //now, we can parse the file
    string token, line; 
    istringstream iss; 

    vector<Vertex_t> vertices; 

    while (getline(infile, line)) {

        iss = istringstream(line); 
        
        //check if this line is a comment line. if so, skip. 
        iss >> token; 
        if (!token.empty() && token == "#") continue; 

        Vertex_t vertex; 
        iss >> vertex.x >> vertex.y; 

        vertices.push_back(vertex); 
    }


    //we'll do a try-catch block here, as an exception will be thrown in this step if any of the segments overlap 
        try {   
            
            AddVertices(vertices); 

        } catch (const InvalidVertexException& e) {

            ostringstream oss; oss << "in <" << here << ">: vertices in input file are invalid; one or more segements intersect"; 
            fStatus = kError; 
            throw DBFileException(oss.str().c_str(), path_infile); 
            return; 
        }
}
//_____________________________________________________________________________________________________________________
void PolynomialCut::InteractiveApp(TH2* hist, const char* drawing_option, unsigned int palette)
{
    //launch the interactive app
    new PolynomialCutApp(gClient->GetRoot(), 1200, 800, hist, drawing_option, palette);
}
//_____________________________________________________________________________________________________________________

ClassImp(PolynomialCut); 