#define r_e     6378137.0
#define f       0.00335281
#define We      0.000072921151467
#define e3      1000.0
#define e6      1000000.0
#define Pi      3.14159265358979323846264338
#define e       2.71828182845904523536028747135

class KalmanFilter{

    private:
        void ErrorDynamics_Element(double [16][16],double ,double [3],double [3],double [3],double [3][3],double ,double );
        void ErrorDynamics_Transition(double [16][16],double T,double [4][4],double [4][4],double [4][4],double [4][4],double [4][4],double [4][4],double [4][4],double [4][4],double [4][4],double [4][4],double ,double );
        void Noise_CovarianceMatrix(double[16][16],double ,double [3][3],double ,double ,double ,double );
        void MeasurementVector(double [16][16],double [3],double [3],double [3][3],double [3],double [3],double [3],double [3]);
        void MeasurementNoise_CovarianceMatrix(double [16][16],double [16],double [16],double [16]);
        void MeasurementDesign_Matrix(double [16][16],double [3][3],double [3]);
    
    public:
        void PredictStep(double [16][16],double [16][16],double ,double [3],double [3],double [3],double [3][3],double ,double ,double ,double );
        void MeasurementUpdate(double [16][16],double [16][16],double [3],double [3],double [3][3],double [3],double [3],double [3],double [3],double [3],double [3],double [3]);
        void State_MeasurementUpdate(double [16][16],double [3],double [3],double [3][3],double [3][3],double [4],double [3],double [3]);

};