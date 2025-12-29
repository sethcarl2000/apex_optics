#ifndef Add_branch_from_Trajectory_t_H
#define Add_branch_from_Trajectory_t_H

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <map> 
#include <string>
#include <sstream> 
#include <vector> 
#include <stdexcept> 
#include <ApexOptics.h> 

//______________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Trajectory_t struct.
ROOT::RDF::RNode Add_branch_from_Trajectory_t(
    ROOT::RDF::RNode df, 
    const char* branch_in, 
    std::map<std::string, double ApexOptics::Trajectory_t::*> branches
)
{
    using ApexOptics::Trajectory_t; 

    const int n_nodes = branches.size() + 1; 
    ROOT::RVec<ROOT::RDF::RNode> nodes{ df }; 

    int i_branch=0; 
    for (auto it = branches.begin(); it != branches.end(); it++) {
        
        //name of this output branch
        const char* branch_name = it->first.data(); 

        double Trajectory_t::*coord = it->second; 

        //define a new branch with name 'branch_name' which corresponds to 'Trajectory_t::coord' 
        nodes.push_back( nodes.back()
            
            .Define(branch_name, [coord](const Trajectory_t& track) { return track.*coord; }, {branch_in})
        ); 
    }
    return nodes.back(); 
}

//______________________________________________________________________________________________________________________
//this is a helper function which automates the creation of branches, which are just members of the Trajectory_t struct.
// this adds the branches in-order. so if you provide the parameters {"x_sv","y_sv","dxdz_sv"}, then they will be assigned to
// the branches Trajectory_t::x, Trajectory_t::y, Trajectory_t::dxdz.
ROOT::RDF::RNode Add_branches_from_Trajectory_t(
    ROOT::RDF::RNode df, 
    const char* branch_in, 
    std::vector<std::string> branches
)
{
    using namespace std;  
    using ApexOptics::Trajectory_t; 

    //number of branches we're going to define 
    const int n_branches = branches.size(); 
 
    //check if there's more than 5 branches to define (not allowed!)
    if (n_branches > 5) {
        ostringstream oss; 

        oss << "in <" << __func__ << ">: Attempted to define "<<n_branches<<" branches from a Trajectory_t struct, "
            "but there can be no more than 5; it only has 5 members!";
        throw invalid_argument(oss.str()); 
        return df; 
    }

    ROOT::RVec<ROOT::RDF::RNode> nodes{ df }; 

    const vector<double Trajectory_t::*> Trajectory_t_branches{
        &Trajectory_t::x,
        &Trajectory_t::y,
        &Trajectory_t::dxdz,
        &Trajectory_t::dydz,
        &Trajectory_t::dpp
    }; 


    for (int i=0; i<n_branches; i++) {
        
        //name of this output branch
        const char* branch_name = branches[i].data(); 

        double Trajectory_t::*coord = Trajectory_t_branches[i]; 

        //define a new branch with name 'branch_name' which corresponds to 'Trajectory_t::coord' 
        nodes.push_back( nodes.back()
            .Define(branch_name, [coord](const Trajectory_t& track) { return track.*coord; }, {branch_in})
        );
    }
    return nodes.back(); 
}

#endif 