#define r_e     6378137.0                           //地球の長半径
#define f       0.00335281                          //偏平率
#define We      0.000072921151467                   //地球の自転速度
#define e3      1000.0
#define e6      1000000.0
#define Pi      3.14159265358979323846264338
#define e       2.71828182845904523536028747135

class DataConvert{
    private:
        
    
    public:
        void Convers_Accel(double ,double [3],double [3],double [3]);
        void Convers_Gyro(double ,double [3],double [3],double [3]);
        void Convers_Mag(double ,double [3],double [3],double [3]);
        void Comvers_GPS(double [3],double [3]);
        double Comvers_GeodeticLatitude(double );
        double Comvers_GeocentricLatitude(double );

};