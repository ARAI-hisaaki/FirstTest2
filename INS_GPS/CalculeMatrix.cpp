//Class_Inclide
#include "mbed.h"
#include "CalculeMatrix.h"

void CalculeMatrix::Addition_Mat(double A,double a [MV_C15][MV_R15],double B,double b [MV_C15][MV_R15],double c [MV_C15][MV_R15]){

        int columnA,columnB,rowA,rowB;//列（縦）の数a,b   行（横）の数a,b

        columnA = (int)a [1][0];
        rowA = (int)a [0][1];
        columnB = (int)b [1][0];
        rowB = (int)b [0][1];

        if(columnA != columnB ||rowA != rowB){
            //pc.printf("\r\n  Addition_Mat  a [1][0] = %d    a [0][1] = %d    b [1][0] = %d    b [0][1] = %d  \r\n",columnA,rowA,columnB,rowB);
        }
        for(int i=1;i<=columnA;i++){
            for(int j=1;j<=rowA;j++){
                c [i][j] = A*a[i][j] + B*b [i][j];
            }
        }
        c [1][0] = columnA;
        c [0][1] = rowA;
}

void CalculeMatrix::Addition_Mat3x3(double A,double a [MV_C3][MV_R3],double B,double b [MV_C3][MV_R3],double c [MV_C3][MV_R3]){

        int columnA,columnB,rowA,rowB;//列（縦）の数a,b   行（横）の数a,b

        columnA = (int)a [1][0];
        rowA = (int)a [0][1];
        columnB = (int)b [1][0];
        rowB = (int)b [0][1];

        if(columnA != columnB ||rowA != rowB){
            //pc.printf("\r\n  Addition_Mat  a [1][0] = %d    a [0][1] = %d    b [1][0] = %d    b [0][1] = %d  \r\n",columnA,rowA,columnB,rowB);
        }
        for(int i=1;i<=columnA;i++){
            for(int j=1;j<=rowA;j++){
                c [i][j] = A*a[i][j] + B*b [i][j];
            }
        }
        c [1][0] = columnA;
        c [0][1] = rowA;
}

void CalculeMatrix::Multi_Mat(double a [MV_C15][MV_R15],double b [MV_C15][MV_R15],double c [MV_C15][MV_R15]){
    
    int columnA,columnB,rowA,rowB;//列（縦）の数a,b   行（横）の数a,b
    double Me_ij=0.0;

    columnA = (int)a [1][0];
    rowA = (int)a [0][1];
    columnB = (int)b [1][0];
    rowB = (int)b [0][1];

    if(rowA != columnB){

    }

    for(int i=1;i<=columnA;i++){
        for(int j=1;j<=rowB;j++){

            for(int k=1;k<=rowA;k++){
                Me_ij = Me_ij + a [i][k]*b[k][j];
            }

            c [i][j] = Me_ij;
            Me_ij = 0.0;
        }
    }

    c [0][0] = a [0][0];
    c [1][0] = columnA;
    c [0][1] = rowB;
}

void CalculeMatrix::Fusion_3x3Mat(double a [MV_C15][MV_R15],int Column,int Row,double b [MV_C3][MV_R3]){
    
    /*
    int columnA,columnB,rowA,rowB;//列（縦）の数a,b   行（横）の数a,b
    columnA = (int)a [1][0];
    rowA = (int)a [0][1];
    columnB = (int)b [1][0];
    rowB = (int)b [0][1];

    if(columnB != 3 || rowB != 3){
        //System.out.printf("\n  Fusion_3X3Mat   columnB = %d   rowB = %d  \n",columnB,rowB);
    }
    */

    for(int i=1;i<=3;i++){
        for(int j=1;j<=3;j++){
            a [3*(Column-1) + i][3*(Row-1) + j] = b [i][j];
        }
    }
}

void CalculeMatrix::Trans_Mat(double a [MV_C15][MV_R15],double b [MV_C15][MV_R15]){

    int columnA,rowA;

    columnA = (int)a [1][0];
    rowA = (int)a [0][1];

    for(int i=1;i<=columnA;i++){
        for(int j=1;j<=rowA;j++){
            b [j][i] = a [i][j];
        }
    }

    b [0][0] = a [0][0];
    b [0][1] = a [1][0];
    b [1][0] = a [0][1];
}

void CalculeMatrix::Trans_Mat3x3(double a [MV_C3][MV_R3],double b [MV_C3][MV_R3]){

    int columnA,rowA;

    columnA = (int)a [1][0];
    rowA = (int)a [0][1];

    for(int i=1;i<=columnA;i++){
        for(int j=1;j<=rowA;j++){
            b [j][i] = a [i][j];
        }
    }

    b [0][0] = a [0][0];
    b [0][1] = a [1][0];
    b [1][0] = a [0][1];
}


void CalculeMatrix::Inv_Mat(double a [MV_C15][MV_R15],double b [MV_C15][MV_R15]){

    int N = ((int)a [1][0]);
    double A [15][30];

    //A = new double [N][2*N];
    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            A [i][j] = a [i+1][j+1];

            if(i == j){
                A [i][j+N] = 1.0;
            }else{
                A [i][j+N] = 0.0;
            }

        }
    }

    for(int k=0; k<N; k++){
        pivot(k,A,N);
        sweep(k,A,N);
    }

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            b [i+1][j+1] = A [i][j+N];
        }
    }

    b [0][0] = a [0][0];
    b [1][0] = a [1][0];
    b [0][1] = a [0][1];
}

void CalculeMatrix::pivot(int k,double a[MV_C15][30],int N){
    double max,copy;
    int ip=k;
        
    max=abs(a[k][k]);
        
    for(int i=k+1; i<N; i++){
        if(max<abs(a[i][k])){
            ip=i;
            max=abs(a[i][k]);
        }
    }

    if(ip!=k){
        for(int j=0; j<2*N; j++){
            copy    =a[ip][j];
            a[ip][j]=a[k][j];
            a[k][j] =copy;
        }
    }
}

void CalculeMatrix::sweep(int k,double a[MV_C15][30],int N){
    double piv,mmm;

    piv=a[k][k];

    for(int j=0; j<2*N; j++){
            a[k][j]=a[k][j]/piv;
    }
    
    for(int i=0; i<N; i++){
            mmm=a[i][k];

            if(i!=k){
                for(int j=k; j<2*N; j++){
                    a[i][j]=a[i][j]-mmm*a[k][j];
                }
            }
    }   
}


