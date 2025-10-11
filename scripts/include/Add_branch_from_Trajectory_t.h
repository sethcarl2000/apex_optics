#ifndef Add_branch_from_Trajectory_t_H
#define Add_branch_from_Trajectory_t_H

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <map> 
#include <string>
#include <vector> 
#include <ApexOptics.h> 

//______________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Trajectory_t struct.
ROOT::RDF::RNode Add_branch_from_Trajectory_t(
    ROOT::RDF::RNode df, 
    const char* branch_in, 
    std::map<std::string, double ApexOptics::Trajectory_t::*> branches
)
{
    const int n_nodes = branches.size() + 1; 
    ROOT::RVec<ROOT::RDF::RNode> nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double ApexOptics::Trajectory_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Trajectory_t::coord' 
        nodes.push_back( nodes.back()
            
            .Define(branch_name, [coord](const ApexOptics::Trajectory_t& track) { return track.*coord; }, {branch_in})
        ); 
    }
    return nodes.back(); 
}

#endif 