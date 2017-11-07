#define r_e     6378137.0
#define f       0.00335281
#define We      0.000072921151467
#define e3      1000.0
#define e6      1000000.0
#define Pi      3.14159265358979323846264338
#define e       2.71828182845904523536028747135

class Initialization{
    private:
        void VectorRectification(double ,double ,double [3],double [3],double [3]);
        double Signum(double);
        double Cul_Sqrt(double,double);
    
    public:
        void ComputeDCM(double [3],double [3][3],double [3],double [3],double [3]);
        void ComputeDCM_FromAccelAzimuthal(double [3],double ,double [3],double [3][3]);
        void ComputeDCM_FromQuaternion(double [4],double [3][3],double [3][3]);
        void ComputeQuaternion_FromDCM(double [3][3],double [4]);
        void ComputeEulerAngle_FromDCM(double [3][3],double [3]);
        void ComputeGravity(double [3],double [3]);
        
};