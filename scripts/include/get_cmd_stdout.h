#ifndef get_cmd_stdout_H
#define get_cmd_stdout_H

#include <fstream>
#include <TSystem.h> 

std::ifstream process_cmd(std::string cmd, std::string& cmd_stdout, std::string& cmd_stderr) 
{
    cmd = cmd + " > cmd.log"; 
    gSystem->Exec(cmd.data()); 

    return std::ifstream("cmd.log"); 
}

#endif 