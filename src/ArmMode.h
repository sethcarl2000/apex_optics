#ifndef ArmMode_H
#define ArmMode_H

#include <stdint.h> 

/// @brief Bit flag to indicate which arm to use. 
namespace ArmMode {
    enum Bit : int8_t { 
        kNone=0, 
        kRHRS=1<<0,
        kLHRS=1<<1, 
        kBoth=(kRHRS | kLHRS)
    };
}

#endif 