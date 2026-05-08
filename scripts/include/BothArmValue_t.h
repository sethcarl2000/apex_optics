#ifndef BothArmValue_t_H
#define BothArmValue_t_H

#include <ArmMode.h> 

//very simple holder for quantities that have one copy for each arm. 
template<typename T> struct BothArmValue_t {
    T R; 
    T L; 

    T& operator() (ArmMode::Bit flag) { 
        if (flag & ArmMode::kRHRS) { return R; } else { return L; }
    }
};

#endif 