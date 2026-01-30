#include <stdio.h>
#include <math.h>
#include <esp_dsp.h>
#include "inverse_kinematics.h"

volatile bool within_limits = true;

#define MAX_JUMP_RAD (10.0f * M_PI / 180.0f)  

void trans_mat(float theta, float alpha, float a, float d, float *A) {
    float ct = cosf(theta);
    float st = sinf(theta);
    float ca = cosf(alpha);
    float sa = sinf(alpha);

    A[0]  = ct;      A[1]  = -st*ca;  A[2]  = st*sa;  A[3]  = a*ct;
    A[4]  = st;      A[5]  = ct*ca;   A[6]  = -ct*sa; A[7]  = a*st;
    A[8]  = 0.0f;    A[9]  = sa;      A[10] = ca;     A[11] = d;
    A[12] = 0.0f;    A[13] = 0.0f;    A[14] = 0.0f;   A[15] = 1.0f;
}

void forward_kinematics(float q[DOF], float *T06)
{
    const float dh_table[DOF * 4] = {
        q[0],              (float)M_PI / 2, 0.0f, a1,
        q[1],              0.0f,             a2,   0.0f,
        q[2] + (float)M_PI / 2, (float)M_PI / 2, 0.0f, 0.0f,
        q[3],             -(float)M_PI / 2,  0.0f, d4,
        q[4],              (float)M_PI / 2,  0.0f, 0.0f,
        q[5],              0.0f,              0.0f, d6
    };

    /* Identity matrix */
    float T[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

    float Ti[16];
    float Ttemp[16];

    for (int i = 0; i < DOF; i++) {
        trans_mat(
            dh_table[i * 4 + 0],
            dh_table[i * 4 + 1],
            dh_table[i * 4 + 2],
            dh_table[i * 4 + 3],
            Ti
        );

        /* Ttemp = T * Ti */
        dspm_mult_4x4x4_f32(T, Ti, Ttemp);

        /* Copy back */
        for (int k = 0; k < 16; k++) {
            T[k] = Ttemp[k];
        }
    }

    /* Output final transform */
    for (int i = 0; i < 16; i++) {
        T06[i] = T[i];
    }
}


void R03(float q1, float q2, float q3, float *R) {
    float c1 = cosf(q1);
    float s1 = sinf(q1);

    float q23 = q2 + q3;
    float c23 = cosf(q23);
    float s23 = sinf(q23);

    R[0] = -s23 * c1; R[1] =  s1; R[2] =  c23 * c1;
    R[3] = -s23 * s1; R[4] = -c1; R[5] =  c23 * s1;
    R[6] =  c23; R[7] =  0.0f; R[8] =  s23;
}

void inverse_position_a(float wx, float wy, float wz, float *q) {
    float t1_a = atan2f(wy, wx);

    float r1 = fmaxf(sqrtf(fmaxf(0.0f, wx*wx + wy*wy)), 1e-12f);
    float r2 = wz - a1;
    float r3 = fmaxf(sqrtf(fmaxf(0.0f, r1*r1 + r2*r2)), 1e-12f);
    float phi1 = atan2f(r2, r1);

    float phi2cos = clampf(((d4*d4) - (a2*a2) - (r3*r3)) / (-2.0f * a2 * r3), -1.0f, 1.0f);
    float phi2 = atan2f(sqrtf(fmaxf(0.0f, 1.0f - phi2cos*phi2cos)), phi2cos);
    float t2_a = phi1 + phi2;

    float phi3cos = clampf(((r3*r3) - (a2*a2) - (d4*d4)) / (-2.0f * a2 * d4), -1.0f, 1.0f);
    float phi3 = atan2f(sqrtf(fmaxf(0.0f, 1.0f - phi3cos*phi3cos)), phi3cos);
    float t3_a = phi3 - (float)M_PI;

    q[0] = t1_a;
    q[1] = t2_a;
    q[2] = t3_a;
}

void inverse_position_b(float wx, float wy, float wz, float *q){
    float t1_b1 = atan2f(wy, wx) - (float)M_PI;
    float t1_b2 = atan2f(wy, wx) + (float)M_PI;

    float r1 = fmaxf(sqrtf(fmaxf(0.0f, wx*wx + wy*wy)), 1e-12f);
    float r4 = wz - a1;
    float r5 = fmaxf(sqrtf(fmaxf(0.0f, r1*r1 + r4*r4)), 1e-12);
    float phi1 = atan2f(r1, r4);

    float phi2cos = clampf(((d4*d4) - (a2*a2) - (r5*r5)) / (-2.0f * a2 * r5), -1.0f, 1.0f);
    float phi2 = atan2f(sqrtf(fmaxf(0.0f, 1.0f - phi2cos*phi2cos)), phi2cos);
    float phi3 = (float)M_PI/2.0f - (phi1 + phi2);
    float t2_b = (float)M_PI - phi3;

    float phi4cos = clampf(((r5*r5) - (a2*a2) - (d4*d4)) / (-2.0f * a2 * d4), -1.0f, 1.0f);
    float phi4 = atan2f(sqrtf(fmaxf(0.0f, 1.0f - phi4cos*phi4cos)), phi4cos);
    float t3_b = phi4 - (float)M_PI;

    q[0] = t1_b1;
    q[1] = t1_b2;
    q[2] = t2_b;
    q[3] = t3_b;
}

void inverse_orientation(float q1, float q2, float q3, float *R06, float *q_wrist){
    float R0_3[9]; 
    float R36[9]; 
    
    R03(q1, q2, q3, R0_3);

    float R0_3inv[9] = {R0_3[0],R0_3[3],R0_3[6],
                         R0_3[1],R0_3[4],R0_3[7],
                         R0_3[2],R0_3[5],R0_3[8]};

    dspm_mult_3x3x3_f32(R0_3inv, R06, R36);

    float s = sqrtf(R36[2]*R36[2] + R36[5]*R36[5]);
    float t4_up, t5_up, t6_up;
    float t4_down, t5_down, t6_down;

    if (s < 1e-6f) {
        // Singular case: t5 is 0 or pi
        t5_up = (R36[8] > 0) ? 0.0f : M_PI;
        t4_up = 0.0f;     
        t6_up = atan2f(R36[3], R36[0]); 
        t5_down = (R36[8] > 0) ? -0.0f : -M_PI;
        t4_down = 0.0f;
        t6_down = atan2f(R36[3], R36[0]);
    } else {
        t5_up = atan2f(s, R36[8]);
        t4_up = atan2f(R36[5], R36[2]);
        t6_up = atan2f(R36[7], -R36[6]);

        t5_down = atan2f(-s, R36[8]);
        t4_down = atan2f(-R36[5], -R36[2]);
        t6_down = atan2f(-R36[7], R36[6]);
    }

    q_wrist[0] = t4_up;
    q_wrist[1] = t5_up;
    q_wrist[2] = t6_up;

    q_wrist[3] = t4_down;
    q_wrist[4] = t5_down;
    q_wrist[5] = t6_down;
}

void ik_solutions(float *T06, float *solutions){
    float R06[9];

    R06[0] = T06[0]; R06[1] = T06[1]; R06[2] = T06[2];
    R06[3] = T06[4]; R06[4] = T06[5];  R06[5] = T06[6]; 
    R06[6] = T06[8];  R06[7] = T06[9];  R06[8] = T06[10];

    float X = T06[3];
    float Y = T06[7];
    float Z = T06[11];

    float wx = X - d6 * R06[2];
    float wy = Y - d6 * R06[5];
    float wz = Z - d6 * R06[8];

    float arm_a[3];
    float arm_b[4];

    inverse_position_a(wx, wy, wz, arm_a);
    inverse_position_b(wx, wy, wz, arm_b);

    float wrist_a[6];
    float wrist_b[6];
    float wrist_c[6];

    inverse_orientation(arm_a[0], arm_a[1], arm_a[2],R06, wrist_a);
    inverse_orientation(arm_b[0], arm_b[2], arm_b[3],R06, wrist_b);
    inverse_orientation(arm_b[1], arm_b[2], arm_b[3],R06, wrist_c);
    
    solutions[0]  = arm_a[0]; solutions[1]  = arm_a[1]; solutions[2]  = arm_a[2];
    solutions[3]  = wrist_a[0]; solutions[4]  = wrist_a[1]; solutions[5]  = wrist_a[2];

    solutions[6]  = arm_a[0]; solutions[7]  = arm_a[1]; solutions[8]  = arm_a[2];
    solutions[9]  = wrist_a[3]; solutions[10] = wrist_a[4]; solutions[11] = wrist_a[5];

    solutions[12] = arm_b[0]; solutions[13] = arm_b[2]; solutions[14] = arm_b[3];
    solutions[15] = wrist_b[0]; solutions[16] = wrist_b[1]; solutions[17] = wrist_b[2];

    solutions[18] = arm_b[0]; solutions[19] = arm_b[2]; solutions[20] = arm_b[3];
    solutions[21] = wrist_b[3]; solutions[22] = wrist_b[4]; solutions[23] = wrist_b[5];

    solutions[24] = arm_b[1]; solutions[25] = arm_b[2]; solutions[26] = arm_b[3];
    solutions[27] = wrist_c[0]; solutions[28] = wrist_c[1]; solutions[29] = wrist_c[2];

    solutions[30] = arm_b[1]; solutions[31] = arm_b[2]; solutions[32] = arm_b[3];
    solutions[33] = wrist_c[3]; solutions[34] = wrist_c[4]; solutions[35] = wrist_c[5];

}

void filter_solutions(float *solutions, const float *q_min, const float *q_max, int *valid_count, float *valid_solutions) {
    *valid_count = 0;

    for (int i = 0; i < 6; i++){
        int valid = 1;
        for (int j = 0; j < DOF; j++){
            float q = solutions[i * DOF + j];
            if (q < q_min[j] || q > q_max[j]) {
                valid = 0;
                break;
            }
        }
        if (valid) {
            for (int j = 0; j < DOF; j++){
                valid_solutions[*valid_count * DOF + j] = solutions[i * DOF + j];
            }
            (*valid_count)++;
        }
    }
}

void final(float *valid_solutions, int valid_count, float *q_prev, float *q_final) {
    if (valid_count == 0) {
        for (int j = 0; j < DOF; j++) {
            q_final[j] = clampf(q_prev[j], q_min[j], q_max[j]);
        }
        within_limits = false;
        return;
    }

    int best_idx = 0;
    float best_dist = 1e30f;

    // Find closest solution to q_prev
    for (int i = 0; i < valid_count; i++) {
        float dist = 0.0f;
        for (int j = 0; j < DOF; j++) {
            float diff = valid_solutions[i * DOF + j] - q_prev[j];
            dist += diff * diff;
        }

        if (dist < best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }

    for (int j = 0; j < DOF; j++) {
        float delta = fabsf(valid_solutions[best_idx * DOF + j] - q_prev[j]);
        if (delta > MAX_JUMP_RAD) {
            // Too large jump → keep previous joint values
            for (int k = 0; k < DOF; k++)
                q_final[k] = q_prev[k];
            within_limits = false;
            return;
        }
    }


    for (int j = 0; j < DOF; j++)
        q_final[j] = valid_solutions[best_idx * DOF + j];

    within_limits = true; 
}

static inline void normalize_quat(float *qx, float *qy, float *qz, float *qw)
{
    float n = (*qx)*(*qx) + (*qy)*(*qy) + (*qz)*(*qz) + (*qw)*(*qw);

    // avoid divide-by-zero
    if (n < 1e-12f) {
        // fallback: identity rotation
        *qx = 0.0f;
        *qy = 0.0f;
        *qz = 0.0f;
        *qw = 1.0f;
        return;
    }

    float inv = 1.0f / sqrtf(n);
    *qx *= inv;
    *qy *= inv;
    *qz *= inv;
    *qw *= inv;
}



void IK(float X, float Y, float Z, float qx, float qy, float qz, float qw, float *q_prev ,float *q_final){

    float R[9];
    float R06[9];
    float T[16];
    float solutions[36];
    float valid_solutions[36];
    int valid_count = 0;

    normalize_quat(&qx, &qy, &qz, &qw);

    // --- Quaternion → Rotation Matrix (3x3) --
    R[0] = 1.0f - 2.0f*(qy*qy + qz*qz);
    R[1] = 2.0f*(qx*qy - qw*qz);
    R[2] = 2.0f*(qx*qz + qw*qy);

    R[3] = 2.0f*(qx*qy + qw*qz);
    R[4] = 1.0f - 2.0f*(qx*qx + qz*qz);
    R[5] = 2.0f*(qy*qz - qw*qx);

    R[6] = 2.0f*(qx*qz - qw*qy);
    R[7] = 2.0f*(qy*qz + qw*qx);
    R[8] = 1.0f - 2.0f*(qx*qx + qy*qy);

    // urdf frame offset
    static const float Ru[9] = {0,0,1, 0,-1,0, 1,0,0};

    // Ru is orthogonal so Ru = RuT
    // R06 = R* RuT
    dspm_mult_3x3x3_f32(R, Ru, R06);



    T[0] = R06[0]; T[1] = R06[1]; T[2] = R06[2]; T[3] = X;
    T[4] = R06[3]; T[5] = R06[4]; T[6] = R06[5]; T[7] = Y;
    T[8] = R06[6]; T[9] = R06[7]; T[10] = R06[8]; T[11] = Z;
    T[12] = 0.0f; T[13] = 0.0f; T[14] = 0.0f; T[15] = 1.0f;


    ik_solutions(T, solutions);
    filter_solutions(solutions, q_min, q_max, &valid_count, valid_solutions);


    final(valid_solutions, valid_count, q_prev, q_final);
}
