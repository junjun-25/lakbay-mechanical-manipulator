#ifndef ROBOT_CONFIG_H
#define ROBOT_CONFIG_H

#include <math.h>
// ====================================================
// === ROBOT CONFIGURATION PARAMETERS
// ====================================================
#define DOF 6

// ==================================================
// ===== LINK LENGTH ======
#define a1 100.0f
#define a2 500.0f
#define a3 500.0f
#define a4 100.0f
#define a5 100.0f
#define a6 100.0f

#define d4  (a3 + a4)
#define d6  (a4 + a5)


// =========================================
// === Numerical IK Refinement (LM) =======
// =========================================
#define MAX_ITER 30
#define POS_TOL  0.5f  
#define ORI_TOL  0.01f
#define ALPHA    1.0f 


// ============================================
// ============= JOINT LIMITS =================
// ============================================

static const float q_min[DOF] = {
    -165.0f * ((float)M_PI / 180.0f),
    -70.0f  * ((float)M_PI / 180.0f),
    -160.0f * ((float)M_PI / 180.0f),
    -180.0f * ((float)M_PI / 180.0f),
    -100.0f * ((float)M_PI / 180.0f),
    -180.0f * ((float)M_PI / 180.0f)
};

static const float q_max[DOF] = {
     165.0f * ((float)M_PI / 180.0f),
     170.0f * ((float)M_PI / 180.0f),
     160.0f * ((float)M_PI / 180.0f),
     180.0f * ((float)M_PI / 180.0f),
     100.0f * ((float)M_PI / 180.0f),
     180.0f * ((float)M_PI / 180.0f)
};






#endif // ROBOT_CONFIG_H
