#include <stdio.h>
#include <math.h>
#include <freertos/FreeRTOS.h>  
#include <freertos/task.h>
#include "esp_dsp.h"
#include "matrix.h"
#include <stdlib.h>


#define DOF 6
#define a1 100.0f
#define a2 500.0f
#define a3 500.0f
#define a4 100.0f
#define a5 100.0f
#define a6 100.0f

#define d4 (a3 + a4)
#define d6 (a4 + a5)

float dt = 0.1f;


float xdot[DOF] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

float xddot[DOF] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
//float xcmd[DOF] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

static float T0[16], T01[16], T02[16], T03[16], T04[16], T05[16], T06[16];
static float J[36];
static float JT[36];


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

void forward_kinematics(float q[DOF], float *T0, float *T01, float *T02, float *T03, float *T04, float *T05, float *T06) {
    const float dh_table[DOF*4] = {
        q[0], (float)M_PI/2, 0.0f, a1,
        q[1], 0.0f,   a2,   0.0f,
        q[2] +(float)M_PI/2, (float)M_PI/2, 0.0f, 0.0f,
        q[3], -(float)M_PI/2, 0.0f, d4,
        q[4], (float)M_PI/2, 0.0f, 0.0f,
        q[5], 0.0f,   0.0f, d6
    };

    T0[0] = 1.0f; T0[1] = 0.0f; T0[2] = 0.0f; T0[3] = 0.0f;
    T0[4] = 0.0f; T0[5] = 1.0f; T0[6] = 0.0f; T0[7] = 0.0f; 
    T0[8] = 0.0f; T0[9] = 0.0f; T0[10] = 1.0f; T0[11] = 0.0f; 
    T0[12] = 0.0f; T0[13] = 0.0f; T0[14] = 0.0f; T0[15] = 1.0f; 

    float T12[16];
    float T23[16];
    float T34[16];
    float T45[16];
    float T56[16];

    trans_mat(dh_table[0], dh_table[1], dh_table[2], dh_table[3], T01);
    trans_mat(dh_table[4], dh_table[5], dh_table[6], dh_table[7], T12);
    trans_mat(dh_table[8], dh_table[9], dh_table[10], dh_table[11], T23);
    trans_mat(dh_table[12], dh_table[13], dh_table[14], dh_table[15], T34);
    trans_mat(dh_table[16], dh_table[17], dh_table[18], dh_table[19], T45);
    trans_mat(dh_table[20], dh_table[21], dh_table[22], dh_table[23], T56);

    dspm_mult_4x4x4_f32(T0, T01, T01);   
    dspm_mult_4x4x4_f32(T01, T12, T02);  
    dspm_mult_4x4x4_f32(T02, T23, T03);
    dspm_mult_4x4x4_f32(T03, T34, T04);
    dspm_mult_4x4x4_f32(T04, T45, T05);
    dspm_mult_4x4x4_f32(T05, T56, T06);
}

void jacobian_matrix(float *T0, float *T01, float *T02, float *T03, float *T04, float *T05, float *T06, float *J)
{
    // z-axis of each frame (rotation axis)
    float z0[3] = {T0[2], T0[6], T0[10]};
    float z1[3] = {T01[2], T01[6], T01[10]};
    float z2[3] = {T02[2], T02[6], T02[10]};
    float z3[3] = {T03[2], T03[6], T03[10]};
    float z4[3] = {T04[2], T04[6], T04[10]};
    float z5[3] = {T05[2], T05[6], T05[10]};

    // position of each frame
    float dn0[3] = {T0[3], T0[7], T0[11]};
    float dn1[3] = {T01[3], T01[7], T01[11]};
    float dn2[3] = {T02[3], T02[7], T02[11]};
    float dn3[3] = {T03[3], T03[7], T03[11]};
    float dn4[3] = {T04[3], T04[7], T04[11]};
    float dn5[3] = {T05[3], T05[7], T05[11]};
    float dn6[3] = {T06[3], T06[7], T06[11]}; 

    float *z[6] = {z0, z1, z2, z3, z4, z5};
    float *d[7] = {dn0, dn1, dn2, dn3, dn4, dn5, dn6};

    for(int i=0; i<6; i++)
    {
        float lv[3];  // linear velocity component: cross(z_i, (d6 - d_i))
        lv[0] = z[i][1]*(dn6[2]-d[i][2]) - z[i][2]*(dn6[1]-d[i][1]);
        lv[1] = z[i][2]*(dn6[0]-d[i][0]) - z[i][0]*(dn6[2]-d[i][2]);
        lv[2] = z[i][0]*(dn6[1]-d[i][1]) - z[i][1]*(dn6[0]-d[i][0]);

        // set linear part (row-major)
        J[0*6 + i] = lv[0];
        J[1*6 + i] = lv[1];
        J[2*6 + i] = lv[2];

        // set angular part (row-major)
        J[3*6 + i] = z[i][0];
        J[4*6 + i] = z[i][1];
        J[5*6 + i] = z[i][2];
    }

}


void print_matrix4x4(float *T) {
    for(int row = 0; row < 4; row++) {
        for(int col = 0; col < 4; col++) {
            printf("%8.3f ", T[row*4 + col]); // row-major
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix6x6(float *J) {
    for(int row = 0; row < 6; row++) {
        for(int col = 0; col < 6; col++) {
            printf("%8.6f ", J[row*6 + col]);
        }
        printf("\n");
    }
    printf("\n");
}

void ee_vel(float q_prev[DOF], float qdot[DOF], float xdot[DOF])
{    
    forward_kinematics(q_prev, T0, T01, T02, T03, T04, T05, T06);
    jacobian_matrix(T0, T01, T02, T03, T04, T05, T06, J);
    dspm_mult_f32(J, qdot, xdot, DOF, DOF, 1);
}

void tau2F(float q_prev[DOF], float tau[DOF], float F[DOF]){
    
    forward_kinematics(q_prev, T0, T01, T02, T03, T04, T05, T06);
    jacobian_matrix(T0, T01, T02, T03, T04, T05, T06, J);

    for(int row = 0; row < 6; row++){
        for(int col = 0; col < 6; col++){
            JT[col * 6 + row] = J[row * 6 + col];
        }
    }

    float lamda = 0.05*0.05;

    float I[36] = {lamda, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                   0.0f, lamda, 0.0f, 0.0f, 0.0f, 0.0f,
                   0.0f, 0.0f, lamda, 0.0f, 0.0f, 0.0f,
                   0.0f, 0.0f, 0.0f, lamda, 0.0f, 0.0f,
                   0.0f, 0.0f, 0.0f, 0.0f, lamda, 0.0f,
                   0.0f, 0.0f, 0.0f, 0.0f, 0.0f, lamda};

    float J1[36];

    dspm_mult_f32(JT, J, J1, 6, 6, 6);

    float J2[36];
    for (int i = 0; i < 36; i++) {
        J2[i] = J1[i] + I[i];
    }

    Matrix Aug = NewMatrix(6, 12);
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            Aug->index[i][j] = J2[i * 6 + j];
            Aug->index[i][6 + j] = (i == j) ? 1.0f : 0.0f;
        }
    }
    // Perform Gauss–Jordan elimination to invert Jacobian
    ReducedRowEchelonForm(Aug);

    // extract inverse of Jacobian from the right half
    float J3[36];
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            J3[i * 6 + j] = Aug->index[i][6 + j];
        }
    }

    //Damped least squared psuedo inverse ng Jacobian
    float J_damped[36];
    dspm_mult_f32(J, J3, J_damped, 6, 6, 6);

    dspm_mult_f32(J_damped, tau, F, 6, 6, 1);
    FreeMatrix(&Aug);

}


void lp_filter(){

}

void admittance_control(float *Md, float *Dd, float *Kd, 
                        float *F_filtered, float *xd, 
                        float *xc)
{
    static float xddot_filtered[DOF] = {0};  // persistent filter state
    float xddot[DOF];
    float x[DOF] = {0};
    float xcmd[DOF] = {0};
    float xdot_filtered[DOF] = {0};

    float alpha = 0.9f;   // LPF factor (0.0–1.0), tune this

    for (int i = 0; i < DOF; i++) {
        // Admittance control law (raw acceleration)
        xddot[i] = (1.0f / Md[i]) * (F_filtered[i] 
                    - (Dd[i] * (xdot[i])) 
                    - (Kd[i] * (xc[i] - xd[i])));

        if (!isfinite(xddot[i])) {
            xddot[i] = 0.0f;
        }

        // Low-pass filter on acceleration
        xddot_filtered[i] = alpha * xddot_filtered[i] + (1.0f - alpha) * xddot[i];


        // Integrate filtered acceleration → velocity → position
        xdot[i] += xddot_filtered[i] * dt;

    // Low-pass filter on velocity
        xdot_filtered[i] = alpha * xdot_filtered[i] + (1.0f - alpha) * xdot[i];

        xcmd[i] += xdot_filtered[i] * dt;

        x[i] = xd[i] + xcmd[i];
    }

    printf("xddot (raw)\n");
    for(int i = 0; i < DOF; i++){
        printf("%f \n", xddot[i]);
    }
    printf("xddot_filtered\n");
    for(int i = 0; i < DOF; i++){
        printf("%f \n", xddot_filtered[i]);
    }

    printf("xdot\n");
    for(int i = 0; i < DOF; i++){
        printf("%f \n", xdot[i]);
    }

    printf("xcmd\n");
    for(int i = 0; i < DOF; i++){
        printf("%f \n", x[i]);
    }
}


void app_main(void)
{
int count = 0;
    while(count < 1000){
    float q_prev[DOF] = {0.0f, 0.0f, -0.0f, 0.0f, 0.0f, 0.0f};
    float F[DOF];
    float tau[DOF] = {0.0f, 0.0f, -0.0f, 0.0f, 0.0f, 0.0f};
    float F_filter[DOF] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    float Md[6];
    float Dd[6];
    float Kd[6];

    for (int i = 0; i < 6; i++) {
        Md[i] = 1.0f; 
        Dd[i] = 10.0f;  
        Kd[i] = 30.0f;  
    } // for tuning


    tau2F(q_prev, tau, F);


    printf("Forces: \n");
    for(int i = 0; i < 6; i++){
        printf("%f N \n", F[i]);
    }


    float xd[DOF] = {1300.0f, 0.0, 100.0f, 0.0, 0.0, 0.0};
    float xc[DOF] = {1300.0f, 0.0, 100.0f, 0.0, 0.0, 0.0};

    admittance_control(Md, Dd, Kd, F_filter, xd, xc);


    vTaskDelay(10/portTICK_PERIOD_MS);

    count++;
    }
}
