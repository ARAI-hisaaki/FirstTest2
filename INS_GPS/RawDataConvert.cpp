//Class_Inclide
#include "mbed.h"
#include "RawDataConvert.h"

void DataConvert::Convers_Accel(double temp,double A_Raw [3],double A_Bias [3],double A [3]){
    /*
     * 三軸の加速度の大きさを算出
     */
    double Ax_Lsb,Ay_Lsb,Az_Lsb;
    double Ax_Zero,Ay_Zero,Az_Zero;

    Ax_Lsb = 1.0/9.80665;Ay_Lsb = 1.0/9.80665;Az_Lsb = 1.0/9.80665;
    Ax_Zero = 0.0;Ay_Zero = 0.0;Az_Zero = 0.0;

    /*
    Ax_Lsb = 1.0;Ay_Lsb = 1.0;Az_Lsb = 1.0;
    Ax_Zero = 0.0;Ay_Zero = 0.0;Az_Zero = 0.0;
    */


    A [0] = (A_Raw [0] - Ax_Zero)/Ax_Lsb + A_Bias [0];
    A [1] = (A_Raw [1] - Ay_Zero)/Ay_Lsb + A_Bias [1];
    A [2] = (A_Raw [2] - Az_Zero)/Az_Lsb + A_Bias [2];
}

void DataConvert::Convers_Gyro(double temp,double G_Raw [3],double G_Bias [3],double G [3]){
    /*
     * 三軸の角速度の大きさを算出
     */
    double Gx_Lsb,Gy_Lsb,Gz_Lsb;
    double Gx_Zero,Gy_Zero,Gz_Zero;


    Gx_Lsb = 1.0/( Pi/180.0 );Gy_Lsb = 1.0/( Pi/180.0 );Gz_Lsb = 1.0/( Pi/180.0 );
    Gx_Zero = 0.0;Gy_Zero = 0.0;Gz_Zero = 0.0;

    /*
    Gx_Lsb = 1.0;Gy_Lsb = 1.0;Gz_Lsb = 1.0;
    Gx_Zero = 0.0;Gy_Zero = 0.0;Gz_Zero = 0.0;
    */


    G [0] = (G_Raw [0] - Gx_Zero)/Gx_Lsb + G_Bias [0];
    G [1] = (G_Raw [1] - Gy_Zero)/Gy_Lsb + G_Bias [1];
    G [2] = (G_Raw [2] - Gz_Zero)/Gz_Lsb + G_Bias [2];
}

void DataConvert::Comvers_GPS(double X_NMEA [3],double X [3]){
    /*
     * NMEAデータをWGS84の座標データに変換
     * Xの緯度は地心緯度ではなく測地緯度
     * 高度は気圧のデータ
     */

    int XN_0,XE_0;
    double XN_1,XE_1;
    double XN,XE,XH;

    XN_0 = (int)(X_NMEA [0]/100.0);
    XN_1 = X_NMEA [0] - ((double)XN_0)*100.0;
    XN_1 = XN_1/60.0;
    XN = (double)XN_0 + XN_1;
    XN = XN/180.0*Pi;

    XE_0 = (int)(X_NMEA [1]/100.0);
    XE_1 = X_NMEA [1] - ((double)XE_0)*100.0;
    XE_1 = XE_1/60.0;
    XE = (double)XE_0 + XE_1;
    XE = XE/180.0*Pi;

    XH = X_NMEA [2];

    X [0] = Comvers_GeodeticLatitude(XN);
    X [1] = XE;
    X [2] = XH;

    //System.out.printf("%f %f  %f  \n",X [0],X [1],X [2]);

}

double DataConvert::Comvers_GeodeticLatitude(double Geocentric_Latitude){
    //地心緯度を測地緯度に変換

    //double b;//WGS84による，短半径
    //double f;//偏平率
    double Geodetic_Latitude;

    //b = r_e*(1-f);

    Geodetic_Latitude = atan(1/(1-f)/(1-f)*tan(Geocentric_Latitude));

    return Geodetic_Latitude;
}

double DataConvert::Comvers_GeocentricLatitude(double Geodetic_Latitude){
    //測地緯度を地心緯度に変換
    double Geocentric_Latitude;
    
    Geocentric_Latitude = atan((1-f)*(1-f)*tan(Geodetic_Latitude));
    
    return Geocentric_Latitude;

}