/// @file Template definiton for polymorphic class 'ModularOpticsModel'   
#ifndef ModularOpticsModel_H
#define ModularOpticsModel_H

#include <ArmMode.h>
#include <ROOT/RDataFrame.hxx>
#include <functional>

/// @brief Parent class for modular optics models 
class ModularOpticsModel {
protected: 

    ROOT::RDF::ColumnNames_t fInputs, fOutputs; 

public: 
    ModularOpticsModel() {}; 
    ~ModularOpticsModel() {}; 

    /// @brief defines output branches for optics model(s)
    /// @param node_in RDF node which has all requisite input branches
    /// @return RDF node with desired target variables defined
    virtual ROOT::RDF::RNode DefineOutputs(ROOT::RDF::RNode node_in) const = 0;

    ROOT::RDF::ColumnNames_t GetInputs()  const { return fInputs; }
    ROOT::RDF::ColumnNames_t GetOutputs() const { return fOutputs; }
};


#endif