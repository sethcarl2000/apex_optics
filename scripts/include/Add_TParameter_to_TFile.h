#ifndef Add_TParameter_to_TFile_H
#define Add_TParameter_to_TFile_H

#include <TParameter.h> 

template<typename T> void Add_TParameter_to_TFile(const char* name, T val)
{
    auto param = new TParameter<T>(name, val); 
    param->Write(); 
}


#endif 