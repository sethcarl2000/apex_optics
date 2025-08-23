#ifndef RDFNodeAccumulator_h_
#define RDFNodeAccumulator_h_

#include <vector>
#include <string> 
#include <stdexcept> 
#include <sstream>
#include <cstdlib> 
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
    std::vector<ROOT::RDF::RNode> fNodes;     
    std::ostringstream fErrorMsg; 
    EStatus fStatus{kGood};
    bool fAbortOnError{true}; 

public: 
    RDFNodeAccumulator(ROOT::RDF::RNode start); 
    ~RDFNodeAccumulator(); 

    //define a new branch, with lambda func. 'expression', and inputs 'inputs' 
    template<typename F> void Define(const char* new_branch, F expression, const vector<std::string>& inputs); 

    //overwrite a previously defined branch, or define it for the first time if it doesn't exist 
    template<typename F> void Overwrite(const char* new_branch, F expression, const vector<std::string>& inputs); 

    //same as above, but only define if this column is not yet defined in the dataframe
    template<typename F> void DefineIfMissing(const char* new_branch, F expression, const std::vector<std::string>& inputs);

    //check if a branch is defined in the dataframe as it currently exists
    bool IsBranchDefined(std::string branch); 


    //bool IsBranchDefined(const char* branch) const { string str(branch); return IsBranchDefined(str); }

    //these do exactly the same as the functions above with the same names, but allow the input branch name to be a
    // std::string
    template<typename F> void Define(std::string new_branch, F expression, const std::vector<std::string>& inputs);
    
    template<typename F> void DefineIfMissing(std::string new_branch, F expression, const std::vector<std::string>& inputs);

    template<typename F> void Overwrite(std::string new_branch, F expression, const std::vector<std::string>& inputs);

    bool IsBranchDefined(const char* branch); 

    ROOT::RDF::RNode& Get() { return fNodes.back(); }

    EStatus GetStatus() const { return fStatus; }
    std::string GetErrorMsg() const { return fErrorMsg.str(); }

    //called when we want to print the stored error message and abort. 
    void PrintMsgAndAbort() {
        cerr << "<RDFNodeAccumulator>: Aborted execution. Error message: " << GetErrorMsg() << endl; 
        std::abort(); 
    }

    //set whether an abort should be called when an error is encountered. 
    void SetAbortOnError(bool _val) { fAbortOnError=_val; }

    //assignment operator 
    ROOT::RDF::RNode& operator=(ROOT::RDF::RNode node) {
        fNodes.emplace_back(node); 
        return fNodes.back(); 
    }
    ClassDef(RDFNodeAccumulator,1); 
};

//now, we get actual class definitions. 

//__________________________________________________________________________________________________________________________________
RDFNodeAccumulator::RDFNodeAccumulator(ROOT::RDF::RNode start) : fNodes{ start } {/* noop */}

//__________________________________________________________________________________________________________________________________
RDFNodeAccumulator::~RDFNodeAccumulator() {/* noop */}

//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Define(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
{
    //check if the status is ok. if not, then abort. 
    if (GetStatus() != kGood) {
        //fErrorMsg << "\n - <Define>: branch '" << new_branch << "' cannot be added, status is not 'kGood'";
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }

    //check if this branch is already defined. if so, error! 
    if (IsBranchDefined(new_branch)) {
        //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
        fErrorMsg << "\n - <Define>: branch '" << new_branch << "' is already defined."
                     "\n             use 'DefineIfMissing' to skip definition if already defined, or 'Overwrite' to overwrite preexisting definition.";
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }

    try {        
        fNodes.emplace_back( fNodes.back().Define(new_branch, expression, inputs) ); 
    
    } catch (const std::exception& e) {
        //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
        fErrorMsg << "\n - <Define>: exception caught defining branch '" << new_branch << "'. what(): " << e.what(); 
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }
}

//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Overwrite(const char* new_branch, F expression, const vector<std::string>& inputs) {
        
    if (GetStatus() != kGood) {
        fErrorMsg << "\n - <Overwrite>: branch '" << new_branch << "' cannot be added, status is not ok.";
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }

    try {
        if (IsBranchDefined(new_branch)) { 
            fNodes.emplace_back( fNodes.back().Redefine(new_branch, expression, inputs) ); 
        } else {
            fNodes.emplace_back( fNodes.back().Define(new_branch, expression, inputs) ); 
        }
        
    } catch (const std::exception& e) {
        //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
        fErrorMsg << "\n - <Overwrite>: exception caught defining branch '" << new_branch << "'. "
                     "\n                what(): " << e.what(); 
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }
}
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::DefineIfMissing(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
{
    if (!IsBranchDefined(new_branch)) Define(new_branch, expression, inputs); 
}
//__________________________________________________________________________________________________________________________________
bool RDFNodeAccumulator::IsBranchDefined(std::string branch) 
{
    for (const std::string& column : fNodes.back().GetColumnNames()) { if (branch == column) return true; }
    return false; 
}
//__________________________________________________________________________________________________________________________________

//these functions are the same as above, but with a different signature

//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Define(std::string new_branch, F expression, const std::vector<std::string>& inputs)
{
    Define(new_branch.c_str(), expression, inputs); 
} 
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Overwrite(std::string new_branch, F expression, const vector<std::string>& inputs) 
{
    Overwrite(new_branch.c_str(), expression, inputs); 
}
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::DefineIfMissing(std::string new_branch, F expression, const std::vector<std::string>& inputs) 
{
    DefineIfMissing(new_branch.c_str(), expression, inputs); 
}
//__________________________________________________________________________________________________________________________________
bool RDFNodeAccumulator::IsBranchDefined(const char* branch)
{
    std::string branch_str(branch); 
    return IsBranchDefined(branch_str); 
}
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________

ClassImp(RDFNodeAccumulator); 

#endif 