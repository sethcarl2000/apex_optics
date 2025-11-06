#ifndef Get_TParameter_from_TFile_H
#define Get_TParameter_from_TFile_H

#include <TParameter.h>
#include <TFile.h> 
#include <stdexcept>
#include <optional>  

//This function fetches the TParameter with name 'param_name' from the TFile under rel. path 'path_file'. 
//_______________________________________________________________________________________________________
template<typename T> std::optional<T> Get_TParameter_from_TFile(const char* path_file, const char* param_name)
{
    const char* const here = "Get_TParameter_from_TFile"
    //try to open the file in read-only mode    
    TFile* file; 
    TParameter<T> *param_ptr; 

    T param; 

    try{ 
        file = new TFile(path_file, "READ");
        if (!file) throw std::exception(Form("Unable to open TFile under path '%s'", path_file)); 

        param_ptr = (TParameter<T>*)file->Get(param_name); 
        if (!param_ptr) throw std::exception(Form("Unable to find parameter '%s' in file '%s'", param_name, path_file)); 
            
        param = param_ptr->GetVal();
        
        file->Close(); 
        delete file; 

    } catch (const std::exception& e) {

        Error(here, "Something went wrong trying to fetch parameter '%s' from file '%s'.\n what(): %s", param_name, path_file, e.what()); 
        return std::nullopt; 
    }

    return std::optional<T>{param}; 
}
//_______________________________________________________________________________________________________

#endif