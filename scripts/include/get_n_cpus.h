#ifndef get_n_cpus_h
#define get_n_cpus_h

#include <string> 
#include <cstdlib>
#include <thread> 

//tries to see if we are in a slurm job. If we are, then report the # of alloted cpus (if we are not, report 'hardware concurrency')
unsigned int get_n_cpus()
{
    const char* n_cpus_env = std::getenv("SLURM_CPUS_PER_TASK");

    if (!n_cpus_env) {
        //we do NOT appear to be in a slurm job...
        return std::thread::hardware_concurrency(); 
    } else {
        //we ARE in a slurm job. report our allotment (hardware_concurrency will tell us about CPUs that are not ours to be using!)
        return std::stoi( std::string{n_cpus_env} );
    }
}


#endif