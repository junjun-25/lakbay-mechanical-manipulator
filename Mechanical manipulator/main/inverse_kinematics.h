#ifndef INVERSE_KINEMATICS_H 
#define INVERSE_KINEMATICS_H 

#include "robot_config.h"
#include <math.h>
#include <stdbool.h>

extern volatile bool within_limits;

// Utility
inline float clampf(float v, float lo, float hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}
void IK(float X, float Y, float Z, float qx, float qy, float qz, float qw, float *q_prev, float *q_final);
void forward_kinematics(float q[DOF], float *T06);



#endif 