#ifndef PolynomialCut_h_
#define PolynomialCut_h_
///////////////////////////////////////////////////////////////////
//
//  PolynomialCut
//
//  This is a class for performing an arbitrary polynomial cut on a TH2 histogram. 
//  
///////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include <vector> 
#include <exception> 

class PolynomialCut {
public: 

    struct Vertex_t  { double x, y; }; 
    struct Segment_t { double x0, y0, sx, sy; };
    
    PolynomialCut() : fStatus{PolynomialCut::kNot_init}, fVertices{}, fSegments{} {}; 
    ~PolynomialCut() {}; 

    //Check if a point is inside the region enclosed by this polynomial 
    bool IsInside(double x, double y) const; 
    
    //add a vertex to a list of vertices 
    void AddVertices(const std::vector<Vertex_t>& vertices); 

    std::vector<Vertex_t> GetVertices() const { return fVertices; }

    //create a dbfile at the given path, with all vertices. check implementation in PolynomialCut.cxx to see output file format. 
    void Create_dbfile(const char* path_outfile) const; 

    void Parse_dbfile(const char* path_infile); 


    enum EStatus { kError=-1, kGood=0, kNot_init=1 };
    EStatus GetStatus() const { return fStatus; }

    //this is a custom exception class, which will be thrown if a new vertex added is invalid. 
    class InvalidVertexException : public std::exception {
    private: 
        std::string fMessage;
        ClassDef(PolynomialCut::InvalidVertexException,1); 
    public:   
        explicit InvalidVertexException(std::string message) : fMessage{message} {}; 

        const char* what() const noexcept override { return fMessage.c_str(); }; 
    }; 

    //exception thrown when there is some sort of fatal error with the input/output file parsing/writing.  
    class DBFileException : public std::exception {
    private: 
        std::string fError, fPath;  
        ClassDef(PolynomialCut::DBFileException,1);  
    public: 
        explicit DBFileException(const char* error, const char* path) 
            : fError(error), fPath(path) {};

        const char* what() const noexcept override { 
            std::string *report = new std::string("\n -- error:           " + fError + 
                                                  "\n -- with file path: '" + fPath + "'\n"); 
            return report->c_str(); 
        }
    }; 

private: 

    //create a segment from two vertices
    Segment_t MakeSegment(const Vertex_t& v1, const Vertex_t& v2) const {
        return Segment_t{ 
            .x0 = v1.x, 
            .y0 = v1.y,
            .sx = v2.x - v1.x, 
            .sy = v2.y - v1.y
        }; 
    };

    //Check if two segements intersect 
    bool DoSegmentsIntersect(const Segment_t& seg1, const Segment_t& seg2) const; 

    //Create segments from the list of all vertices 
    void InitSegments(); 

    EStatus fStatus;
    std::vector<Vertex_t>  fVertices; 
    std::vector<Segment_t> fSegments; 

ClassDef(PolynomialCut,1);
};

#endif 