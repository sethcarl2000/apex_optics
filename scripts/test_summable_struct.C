
//this is meant to test if I can use the ROOT::RDataFrame::Sum() method on a column which holds user-defined class. 
#include "TROOT.h"
#include <stdexcept>

using namespace std; 
using namespace ROOT::VecOps; 

struct SummableStruct_t{ 
    RVec<double> data; 

    void operator+=(const SummableStruct_t& rhs) {
        if (rhs.data.size() == data.size()) {
            data += rhs.data; 
        } else {
            throw std::logic_error("<SummableStruct_t::operator+=>: Number of elements in 'data' vectors does not match between opearands."); 
        }   
    }
}; 

int test_summable_struct()
{ 
    ROOT::RDataFrame df(1e5); 

    auto df_output = df 

        .Define("struct",   []()
        { 
            return SummableStruct_t{ .data={ 
                gRandom->Gaus(),
                gRandom->Gaus(),
                gRandom->Gaus()
            }}; 
        }, {})

        .Define("str_mag",  [](const SummableStruct_t& str)
        {
            return sqrt( str.data[0]*str.data[0] + str.data[1]*str.data[1] + str.data[2]*str.data[2] );    
        }, {"struct"});


    auto hist = df_output.Histo1D<double>({"h", "mag", 200, 0., 5.}, "str_mag"); 

    
    auto df_sum = df_output.Sum<SummableStruct_t>("str_mag"); 

    cout << "average mag: " << (df_sum)->data << endl; 
        

    hist->DrawCopy();


    return 0; 
}