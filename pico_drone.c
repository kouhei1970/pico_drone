#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pico/stdlib.h"
#include "hardware/i2c.h"


u_int8_t euler2dcm(float *euler,float *dcm)
{
    float phi=euler[0];
    float theta=euler[1];
    float psi=euler[2];
    
    float e11= cos(theta)*cos(psi);
    float e12= cos(theta)*sin(psi);
    float e13=-sin(theta);
    
    float e21= sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
    float e22= sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi);
    float e23= sin(phi)*cos(theta);
    
    float e31= cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);
    float e32= cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);
    float e33= cos(phi)*cos(theta);
    
    dcm[0]=e11;
    dcm[1]=e12;
    dcm[2]=e13;

    dcm[3]=e21;
    dcm[4]=e22;
    dcm[5]=e23;

    dcm[6]=e31;
    dcm[7]=e32;
    dcm[8]=e33;

    return true;
}


u_int8_t dcm2euler(float *dcm, float *euler)
{
    float phi=atan2(dcm[1*3+2], dcm[2*3+2]);
    float theta=atan2(-dcm[0*3+2], sqrt(dcm[1*3+2]*dcm[1*3+2] + dcm[2*3+2]*dcm[2*3+2]));
    float psi=atan2(dcm[0*3+1], dcm[0*3+0]);

    euler[0]=phi;
    euler[1]=theta;
    euler[2]=psi;

    return true;
}


u_int8_t quat2dcm(float *q, float *dcm)
{
    float q1=q[0];
    float q2=q[1];
    float q3=q[2];
    float q4=q[3];
    
    float e11=   q1*q1 - q2*q2 - q3*q3 + q4*q4;
    float e12= 2 * (q1*q2 + q3*q4);
    float e13= 2 * (q1*q3 - q2*q4);
    
    float e21= 2 * (q1*q2 - q3*q4);
    float e22= - q1*q1 + q2*q2 - q3*q3 + q4*q4;
    float e23= 2 * (q2*q3 + q1*q4);
    
    float e31= 2 * (q1*q3 + q2*q4);
    float e32= 2 * (q2*q3 - q1*q4);
    float e33= - q1*q1 - q2*q2 + q3*q3 + q4*q4;

    dcm[0]=e11;
    dcm[1]=e12;
    dcm[2]=e13;

    dcm[3]=e21;
    dcm[4]=e22;
    dcm[5]=e23;

    dcm[6]=e31;
    dcm[7]=e32;
    dcm[8]=e33;

    return true;
}

u_int8_t dcm2quat(float *dcm, float *q)
{
    float e11=dcm[0*3+0];
    float e12=dcm[0*3+1];
    float e13=dcm[0*3+2];
    float e21=dcm[1*3+0];
    float e22=dcm[1*3+1];
    float e23=dcm[1*3+2];
    float e31=dcm[2*3+0];
    float e32=dcm[2*3+1];
    float e33=dcm[2*3+2];
    
    float q1=0.5*sqrt(1 + e11 - e22 - e33);
    float q2=0.5*sqrt(1 - e11 + e22 - e33);
    float q3=0.5*sqrt(1 - e11 - e22 + e33);
    float q4=0.5*sqrt(1 + e11 + e22 + e33);
    
    //q=[q1, q2, q3, q4]
    //idx=q.index(max(q))
    
    u_int8_t idx;

    //一番大きいqを探す
    idx=3;
    if (q1>q2){
        if (q1>q3){
            if (q1>q4) idx=0;
        }
    }
    if (q2>q3){
        if (q2>q4) idx=1;        
    }
    if (q3>q4)idx=2;


    if (idx==0){
        q[1]=( e12+e21)/4/q1;
        q[2]=( e13+e31)/4/q1;
        q[3]=( e23-e32)/4/q1;
    }
    else if (idx==1){
        q[0]=( e12+e21)/4/q2;
        q[2]=( e23+e32)/4/q2;
        q[3]=(-e13+e31)/4/q2;
    }
    else if (idx==2){
        q[0]=( e13+e31)/4/q3;
        q[1]=( e23+e32)/4/q3;
        q[3]=( e12-e21)/4/q3;
    }
    else if (idx==3){
        q[0]=( e23-e32)/4/q4;
        q[1]=(-e13+e31)/4/q4;
        q[2]=( e12+e21)/4/q4;
    }

    return true;
}

int main() {
    stdio_init_all();
    printf("Pico_drone ver 0\n");
    return true;
}


/*
def euler2dcm(euler):
    phi=euler[0]
    theta=euler[1]
    psi=euler[2]
    
    e11= np.cos(theta)*np.cos(psi)
    e12= np.cos(theta)*np.sin(psi)
    e13=-np.sin(theta)
    
    e21= np.sin(phi)*np.sin(theta)*np.cos(psi) - np.cos(phi)*np.sin(psi)
    e22= np.sin(phi)*np.sin(theta)*np.sin(psi) + np.cos(phi)*np.cos(psi)
    e23= np.sin(phi)*np.cos(theta)
    
    e31= np.cos(phi)*np.sin(theta)*np.cos(psi) + np.sin(phi)*np.sin(psi)
    e32= np.cos(phi)*np.sin(theta)*np.sin(psi) - np.sin(phi)*np.cos(psi)
    e33= np.cos(phi)*np.cos(theta)
    
    return((e11,e12,e13),(e21,e22,e23),(e31,e32,e33))

def dcm2euler(e):
    phi=np.arctan2(e[1][2], e[2][2])
    theta=np.arctan2(-e[0][2], np.sqrt(e[1][2]**2 + e[2][2]**2))
    psi=np.arctan2(e[0][1], e[0][0])
    
    return (phi, theta, psi)

def quat2dcm(q):
    q1=q[0]
    q2=q[1]
    q3=q[2]
    q4=q[3]
    
    e11=   q1**2 - q2**2 - q3**2 + q4**2
    e12= 2 * (q1*q2 + q3*q4)
    e13= 2 * (q1*q3 - q2*q4)
    
    e21= 2 * (q1*q2 - q3*q4)
    e22= - q1**2 + q2**2 - q3**2 + q4**2
    e23= 2 * (q2*q3 + q1*q4)
    
    e31= 2 * (q1*q3 + q2*q4)
    e32= 2 * (q2*q3 - q1*q4)
    e33= - q1**2 - q2**2 + q3**2 + q4**2
    
    return((e11,e12,e13),(e21,e22,e23),(e31,e32,e33))

def dcm2quat(e):
    e11=e[0][0]
    e12=e[0][1]
    e13=e[0][2]
    e21=e[1][0]
    e22=e[1][1]
    e23=e[1][2]
    e31=e[2][0]
    e32=e[2][1]
    e33=e[2][2]
    
    q1=0.5*np.sqrt(1 + e11 - e22 - e33)
    q2=0.5*np.sqrt(1 - e11 + e22 - e33)
    q3=0.5*np.sqrt(1 - e11 - e22 + e33)
    q4=0.5*np.sqrt(1 + e11 + e22 + e33)
    
    q=[q1, q2, q3, q4]
    idx=q.index(max(q))
    
    if idx==0:
        q[1]=( e12+e21)/4/q1
        q[2]=( e13+e31)/4/q1
        q[3]=( e23-e32)/4/q1
    elif idx==1:
        q[0]=( e12+e21)/4/q2
        q[2]=( e23+e32)/4/q2
        q[3]=(-e13+e31)/4/q2
    elif idx==2:
        q[0]=( e13+e31)/4/q3
        q[1]=( e23+e32)/4/q3
        q[3]=( e12-e21)/4/q3
    elif idx==3:
        q[0]=( e23-e32)/4/q4
        q[1]=(-e13+e31)/4/q4
        q[2]=( e12+e21)/4/q4

    return q
    */