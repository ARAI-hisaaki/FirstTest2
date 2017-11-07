//Class_Inclide
#include "mbed.h"

#define MV_C3       4
#define MV_R3       4
#define MV_C6       7
#define MV_R6       7
#define MV_C9       10
#define MV_R9       10
#define MV_C12      13
#define MV_R12      13
#define MV_C15      16
#define MV_R15      16  
#ifndef CHANGEME_H_
#define CHANGEME_H_
 

class CalculeMatrix{
    private:
        void pivot(int ,double [MV_C15][30],int );
        void sweep(int ,double [MV_C15][30],int );
    
    public:
        void Addition_Mat(double ,double [MV_C15][MV_R15],double ,double [MV_C15][MV_R15],double [MV_C15][MV_R15]);
        void Addition_Mat3x3(double ,double [MV_C3][MV_R3],double ,double [MV_C3][MV_R3],double [MV_C3][MV_R3]);
        void Multi_Mat(double a [MV_C15][MV_R15],double b [MV_C15][MV_R15],double c [MV_C15][MV_R15]);
        void Fusion_3x3Mat(double [MV_C15][MV_R15],int ,int ,double [MV_C3][MV_R3]);
        void Trans_Mat(double [MV_C15][MV_R15],double [MV_C15][MV_R15]);
        void Trans_Mat3x3(double [MV_C3][MV_R3],double [MV_C3][MV_R3]);
        void Inv_Mat(double [MV_C15][MV_R15],double [MV_C15][MV_R15]);
    
};

#endif
