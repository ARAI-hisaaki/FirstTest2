#include "mbed.h"
#include "NormalDistribution.h"

Serial pcDis(USBTX, USBRX);
Timer Distribution;

double normal_distribution::CalculeDistribution(double a,double b){
    Distribution.start();
    //aが平均，bが分散の正規分布を算出．

    double Max_Fx,Fx,ran,ran2;
    double c = 0;
    
    ran = rand()/32767.0;
    ran2 = rand()/32767.0;
    //System.out.printf("%f      ",rand);

    Max_Fx = 1 / sqrt(2*Pi*b);
    Fx = ran*Max_Fx;

    if(ran2 >= 0.5){
        c = a + sqrt( -2*b * log( Fx*sqrt(2*Pi*b) ) );
    }else if(ran2 < 0.5){
        c = a - sqrt( -2*b * log( Fx*sqrt(2*Pi*b) ) );
    }
    
    Distribution.reset();

    return c;
}