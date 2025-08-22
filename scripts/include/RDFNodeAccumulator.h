#ifndef RDFNodeAccumulator_h_
#define RDFNodeAccumulator_h_

#include <vector>
#include <string> 
#include <stdexcept> 
#include <sstream>
#include <ROOT/RDataFrame.hxx>

//A small helper class which is meant to avoid some of the awkward syntax assocaited with RDataFrame creation. 
class RDFNodeAccumulator {
public: 
    enum EStatus {
        kGood    =0,
        kWarning =1, 
        kError =-1
    }; 

private: 
    EStatus fStatus;
    std::vector<ROOT::RDF::RNode> fNodes;     
    std::ostringstream fErrorMsg; 

public: 
    RDFNodeAccumulator(ROOT::RDF::RNode start) : fNodes{start} {}; 
    ~RDFNodeAccumulator() {}; 

    template<typename F> void Define(const char* new_branch, F expression, const vector<std::string>& inputs) {
        
        if (GetStatus() != kGood) {
            //fErrorMsg << "\n - <Define>: branch '" << new_branch << "' cannot be added, status is not 'kGood'";
            return;  
        }

        //check if this branch is already defined. if so, error! 
        if (IsBranchDefined(new_branch)) {
            //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
            fErrorMsg << "\n - <Define>: branch '" << new_branch << "' is already defined."
                         "\n             use 'DefineIfMissing' to skip definition if already defined, or 'Overwrite' to overwrite preexisting definition.";
            fStatus = kError; 
            return;  
        }

        try {        
            fNodes.emplace_back( fNodes.back().Define(new_branch, expression, inputs) ); 
      
        } catch (const std::exception& e) {
            //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
            fErrorMsg << "\n - <Define>: exception caught defining branch'" << new_branch << "'. what(): " << e.what(); 
            fStatus = kError; 
            return; 
        }
    }; 

    //overwrite a previously defined branch, or define it for the first time if it doesn't exist 
    template<typename F> void Overwrite(const char* new_branch, F expression, const vector<std::string>& inputs) {
        
        if (GetStatus() != kGood) {
            fErrorMsg << "\n - <Overwrite>: branch '" << new_branch << "' cannot be added, status is not ok.";
            return;
        }

        try {
            if (IsBranchDefined(new_branch)) { 
                fNodes.emplace_back( fNodes.back().Redefine(new_branch, expression, inputs) ); 
            } else {
                fNodes.emplace_back( fNodes.back().Define(new_branch, expression, inputs) ); 
            }
            
        } catch (const std::exception& e) {
            //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
            fErrorMsg << "\n - <Overwrite>: exception caught defining branch'" << new_branch << "'. "
                         "\n                what(): " << e.what(); 
            fStatus = kError; 
            return; 
        }
    }; 

    //same as above, but only define if this column is not yet defined in the dataframe
    template<typename F> void DefineIfMissing(const char* new_branch, F expression, const std::vector<std::string>& inputs) {

        if (!IsBranchDefined(new_branch)) Define(new_branch, expression, inputs); 
    }
    

    //check if a branch is defined in the dataframe as it currently exists
    bool IsBranchDefined(std::string branch) {
        for (const std::string& column : fNodes.back().GetColumnNames()) { if (branch == column) return true; }
        return false; 
    }
    bool IsBranchDefined(const char* branch) { 
        std::string b_str(branch); 
        return IsBranchDefined(b_str);  
    }

    //returns false if any of the branch-names in 'branches' are missing. otherwise, returns true. 
    bool HasBranches(const std::vector<std::string>& branches) {
        for (const auto& str : branches) { if (!IsBranchDefined(str)) return false; }
        return true; 
    }
    //bool IsBranchDefined(const char* branch) const { string str(branch); return IsBranchDefined(str); }

    //these do exactly the same as the functions above with the same names, but allow the input branch name to be a
    // std::string
    template<typename F> void Define(std::string new_branch, F expression, const std::vector<std::string>& inputs) {
        Define(new_branch.c_str(), expression, inputs); 
    }
    
    template<typename F> void DefineIfMissing(std::string new_branch, F expression, const std::vector<std::string>& inputs) {
        DefineIfMissing(new_branch.c_str(), expression, inputs); 
    }

    template<typename F> void Overwrite(std::string new_branch, F expression, const std::vector<std::string>& inputs) {
        Overwrite(new_branch.c_str(), expression, inputs); 
    }

    ROOT::RDF::RNode& Get() { return fNodes.back(); }

    EStatus GetStatus() const { return fStatus; }
    std::string GetErrorMsg() const { return fErrorMsg.str(); }

    //assignment operator 
    ROOT::RDF::RNode& operator=(ROOT::RDF::RNode node) {
        fNodes.emplace_back(node); 
        return fNodes.back(); 
    }
    ClassDef(RDFNodeAccumulator,1); 
};

#endif 