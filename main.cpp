// Arai Hisaaki(2015 CORE B3)
// This program is a GPS/INS navigation program (alpha version) executed on the STM32-F401RE.

#include "mbed.h"
//sensor
#include "MPU9250.h"
#include "MS5607I2C.h"              //MS5607
#include "SDFileSystem.h"           //SD
//INS_GPS        
#include "CalculeMatrix.h"
#include "Initialization.h"
#include "KalmanFilter.h"
#include "RawDataConvert.h"
#include "Strapdown_INS.h"
#include "CalculeMatrix.h"

#define e3      1000.0
#define e6      1000000.0
#define Pi      3.14159265358979323846264338
#define e       2.71828182845904523536028747135

#define IMU_ratio       0.015625    //IMU [s]   
#define IMU_Size        64          //Hz
#define Tem_Size        1           //IMU_Size*IMU_ratio*Kal_M_ratio
#define Kal_M_ratio     1           //Measurement_update_time   [s]
#define i_nmea_s        200

//Pin Setting 
Serial pc(SERIAL_TX, SERIAL_RX);                        //PC
SPI spi(D11, D12, D13);                                 //MPU9250
mpu9250_spi imu(spi,D10);                               //MPU9250
SDFileSystem    sd(PC_12, PC_11, PC_10, PD_2, "sd");    //SDcard
MS5607I2C       ms5607(D14, D15, false);                //MS5607
Serial GPS(D8,D2);                                      //GPS

//Function Setting 
void Ini();
void Sensor_Read();
void KalP();
void KalM();
void Receive_GPS();
void Receive_Command();
void Command_Library();
Strapdown_INS       INS;
KalmanFilter        Kalman;
DataConvert         RDC;
Initialization      Initia;
CalculeMatrix       CM;
Ticker              ADtimer;
Timer               Aime;
Timer           COM_Time;
int MPU_Num,IMU_Num,Tem_Num,Kal_Num;
int KalP_Num0,KalP_Num1;
int KalMU_flag,DelIMU;

//MPU9250
double Temp [Tem_Size],Temp_ver;             //Sencer_Temp
double Ax [IMU_Size],Ay [IMU_Size],Az [IMU_Size];         //Accel
double Gx [IMU_Size],Gy [IMU_Size],Gz [IMU_Size];         //Gyro
double Al [IMU_Size];
double Time_Sencer;      //Sencer_Time

//GPS
float UTC_Time,Latitude,Longitude,MSL_Alti,Geoid_Sepa,Speed,Course,UTC_Date,Magnetic,HDOP;
char C_Status,C_Latitude,C_MSL,C_Geoid,C_Longitude,C_Magnetic,C_Mode;
int PFI,Satellite;
char NMEA [i_nmea_s];
int ReceiveON = 0;
int i_nmea = 0;

//INS
double delta_t;                                             //Time Step
double Ab [3],Wb [3],Mb [3],Sencer_temp;                    //IMU data
double Xn [3],Vn [3],qbn [4],Cbn [3][3],Cnb [3][3];         //INS quantity
double Xn1 [3],Vn1 [3],qbn1 [4],Cbn1 [3][3],Cnb1 [3][3];    //INS quantity(観測更新に使うINS状態量を保存)

//Kalman Filter
double Xn_GPS [3],Vn_GPS [3],Mn_GPS [3];                    //GPS and Barometric Altitude Sencer Data,Magnet Data
double qx_GPS [3],qv_GPS [3],qm_GPS [3];                    //Dispersion of GPS, Barometric Altitude and Magnet Sencer
double Tba,Tbg,qa,qg;                                       //IMU Sencer time constant
double A_Bias [3],G_Bias [3];                               //IMU(Accel,Gyro)  bias
double Xhat [16][16],P [16][16];                            //KalmaFilter quantity,KalmaFilter Covariance Matrix
double P0 [16][16],dP = 0.0;


//Command Valiable
char COM [10];
int RCOM = 0,icom;

int main() {
    pc.baud(115200);
    GPS.baud(115200/2);
    GPS.attach(Receive_GPS, Serial::RxIrq);
    //pc.baud(921600);
    //MPU9250_Setting
    pc.printf("%d   \r\n",IMU_Num);
    if(imu.init(1,BITS_DLPF_CFG_98HZ)){  //INIT the mpu9250
        pc.printf("\r\nCouldn't initialize MPU9250 via SPI!");
    }    
    pc.printf("\r\nWHOAMI=0x%2x\r\n",imu.whoami()); //output the I2C address to know if SPI is working, it should be 104   
    pc.printf("Gyro_scale=%u\r\n",imu.set_gyro_scale(BITS_FS_2000DPS));    //Set full scale range for gyros 
    pc.printf("Acc_scale=%u\r\n",imu.set_acc_scale(BITS_FS_16G));          //Set full scale range for accs
    pc.printf("AK8963 WHIAM=0x%2x\r\n",imu.AK8963_whoami());
    /*
    imu.init(1,BITS_DLPF_CFG_98HZ);
    imu.whoami();
    imu.set_gyro_scale(BITS_FS_2000DPS);
    imu.set_acc_scale(BITS_FS_16G);
    imu.AK8963_whoami();
    imu.AK8963_calib_Magnetometer();
    */
    pc.printf("%d   \n",IMU_Num);
    MPU_Num = 0;
    IMU_Num = 0;
    Tem_Num = 0;
    Kal_Num = 0;
    KalMU_flag = 0;
    delta_t = IMU_ratio;
    wait(2);
    Ini();
    //Timer Setting
    Aime.start();
    ADtimer.attach(&Sensor_Read,IMU_ratio);//
    
    while(1){
        
        KalP();
        
        if(KalMU_flag == 1){    
            KalP();    
            KalM();    
            for(int i=1;i<16;i++){
                for(int j=1;j<16;j++){
                    dP += P [i][j] - P0 [i][j];
                    P0 [i][j] = P [i][j];
                }
            }
            /*
            pc.printf("Kal_NUM = %d dP10^6 = %f Bias %f %f %f %f %f %f\r\n",
            Kal_Num,
            dP*e6,
            A_Bias [0],
            A_Bias [1],
            A_Bias [2],
            G_Bias [0],
            G_Bias [1],
            G_Bias [2]
            );
            */
            //pc.printf("%f   %f\r\n",Ab [0],A_Bias [0]);
            Kal_Num ++;
            KalMU_flag = 0;
            dP = 0.0;
            
        }
        
        
    }
    
}

void Ini(){
    double D [16][16],U_T [16][16];
    double Sigma_p,Sigma_v,Sigma_e,Sigma_a,Sigma_g;
    
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            Cbn [i][j] = 0.0;
            Cnb [i][j] = 0.0;
        }
    }
    
    for(int i=0;i<16;i++){
        for(int j=0;j<16;j++){
            D [i][j] = 0.0;
            U_T [i][j] = 0.0;
            Xhat [i][j] = 0.0;
            P [i][j] = 0.0;
        }
    }
    
    Xhat [1][0] = 15.0;Xhat [0][1] = 1.0;
    P [1][0] = 15.0;P [0][1] = 15.0;
    D [1][0] = 15.0;D [0][1] = 15.0;
    U_T [1][0] = 15.0;U_T [0][1] = 15.0;
    
    Xn_GPS [0] = Latitude;
    Xn_GPS [1] = Longitude;
    Xn_GPS [2] = ms5607.getAltitude();;
    RDC.Comvers_GPS(Xn_GPS,Xn_GPS);
    
    /*
    Xn [0] = Pi/180*38;
    Xn [1] = Pi/180*135;
    Xn [2] = 30;
    */
    
    Xn [0] = Xn_GPS [0];
    Xn [1] = Xn_GPS [1];
    Xn [2] = Xn_GPS [2];
    
    Vn [0] = 0.0;
    Vn [1] = 0.0;
    Vn [2] = 0.0;

    qbn [0] = 0.29569;
    qbn [1] = 0.653963;
    qbn [2] = - 0.288898;
    qbn [3] = 0.651530;

    Initia.ComputeDCM_FromQuaternion(qbn,Cbn,Cnb);

    A_Bias [0] = 0.0;
    A_Bias [1] = 0.0;
    A_Bias [2] = 0.0;
    
    G_Bias [0] = 0.0;
    G_Bias [1] = 0.0;
    G_Bias [2] = 0.0;
    
    Mb [0] = -1.0;
    Mb [0] = 0.0;
    Mb [0] = 0.0;
    
    double CEP_Xn,CEP_Vn;
    
    CEP_Xn = 0.001;
    qx_GPS [0] = (2*Pi*CEP_Xn/40000/e3)*(2*Pi*CEP_Xn/40000/e3);
    qx_GPS [1] = (2*Pi*CEP_Xn/40000/e3)*(2*Pi*CEP_Xn/40000/e3);
    qx_GPS [2] = (CEP_Xn/1.1774)*(CEP_Xn/1.1774);

    CEP_Vn = 0.001;
    qv_GPS [0] = (CEP_Vn/1.1774)*(CEP_Vn/1.1774);
    qv_GPS [1] = (CEP_Vn/1.1774)*(CEP_Vn/1.1774);
    qv_GPS [2] = (CEP_Vn/1.1774)*(CEP_Vn/1.1774);

    qm_GPS [0] = 7.943/e3;
    qm_GPS [1] = 7.943/e3;
    qm_GPS [2] = 7.943/e3;
    
    Sigma_p = 1;//CEP [m]
    Sigma_v = 0.5;//CEP [m/s]
    Sigma_e = 0.01;//

    Sigma_a = 0.603/e6;//CEP:0.1[m/s^2]
    Sigma_g = 2.5/e6;//CEP:0.5[度/s]
    
    for(int i=1;i<=15;i++){
        int j = i;
        U_T [i][j] = 1.0;
    }

    D [1][1] = Sigma_p*Sigma_p/r_e;
    D [2][2] = Sigma_p*Sigma_p/r_e/cos(Xn [0]);
    D [3][3] = Sigma_p*Sigma_p;
    D [4][4] = Sigma_v*Sigma_v;
    D [5][5] = Sigma_v*Sigma_v;
    D [6][6] = Sigma_v*Sigma_v;
    D [7][7] = Sigma_e*Sigma_e;
    D [8][8] = Sigma_e*Sigma_e;
    D [9][9] = Sigma_e*Sigma_e;
    D [10][10] = Sigma_a*Sigma_a;
    D [11][11] = Sigma_a*Sigma_a;
    D [12][12] = Sigma_a*Sigma_a;
    D [13][13] = Sigma_g*Sigma_g;
    D [14][14] = Sigma_g*Sigma_g;
    D [15][15] = Sigma_g*Sigma_g;
    
    CM.Multi_Mat(D, U_T, P);
}

void Sensor_Read(){
    double XnL [3],VnL [3],qbnL [4],CbnL [3][3],CnbL [3][3]; 
    
    imu.read_acc();
    imu.read_rot();
        
    Ax [IMU_Num] = imu.accelerometer_data[0];
    Ay [IMU_Num] = imu.accelerometer_data[1];
    Az [IMU_Num] = imu.accelerometer_data[2];
    Gx [IMU_Num] = imu.gyroscope_data[0];
    Gy [IMU_Num] = imu.gyroscope_data[1];
    Gz [IMU_Num] = imu.gyroscope_data[2];
    Al [IMU_Num] = ms5607.getAltitude();
    //pc.printf("%f   \r\n",Az [IMU_Num]);
    //センサ座標系から機体座標系に変換
    double A_Raw [3],G_Raw [3];
    A_Raw [0] = - Az [IMU_Num];
    A_Raw [1] = - Ax [IMU_Num];
    A_Raw [2] = - Ay [IMU_Num];
    G_Raw [0] = - Gz [IMU_Num];
    G_Raw [1] = - Gx [IMU_Num];
    G_Raw [2] = - Gy [IMU_Num];
    
    RDC.Convers_Accel(0.0,A_Raw,A_Bias,Ab);
    RDC.Convers_Gyro(0.0,G_Raw,G_Bias,Wb);
    
    INS.INS_TimeUpdate(delta_t,Xn,Vn,Ab,Wb,qbn,Cbn,Cnb,XnL,VnL,qbnL,CbnL,CnbL);
    INS.INS_DataUpdate(Xn,Vn,qbn,Cbn,Cnb,XnL,VnL,qbnL,CbnL,CnbL);
    
    MPU_Num ++;
    IMU_Num ++;
    
    if(IMU_Num == IMU_Size){
        IMU_Num = 0;
        Temp [Tem_Num] = imu.Temperature;
        //pc.printf("%d   %f  %f\r\n",Tem_Num,Xn [2],A_Raw [0]);
        Tem_Num ++;
        
        if(Tem_Num == Tem_Size){
            Tem_Num = 0;
            KalMU_flag = 1;
            DelIMU = IMU_Num;
            
            //観測更新に使うINS状態量を保存
            INS.INS_DataUpdate(Xn1,Vn1,qbn1,Cbn1,Cnb1,Xn,Vn,qbn,Cbn,Cnb);  
            
            for(int i=0;i<Tem_Size;i++){
                Temp_ver +=  Temp [i];
            }
            Temp_ver = Temp_ver/Tem_Size;
        }
    }
    
}

void KalP(){
    double KalP_TimeStep;
    
    KalP_Num1 = IMU_Num;
    if(KalP_Num0 != KalP_Num1){
            
        if(KalP_Num1 > KalP_Num0){
            KalP_TimeStep = delta_t*(double)(KalP_Num1 - KalP_Num0);
        }else{
            KalP_TimeStep = delta_t*(double)(KalP_Num1 + IMU_Size - KalP_Num0);
        }
            
        Kalman.PredictStep(Xhat, P, KalP_TimeStep, Xn, Vn, Ab, Cbn, Tba, Tbg, qa, qg);
        KalP_Num0 = KalP_Num1;
    }
}

void KalM(){
    double XnL [3],VnL [3],qbnL [4],CbnL [3][3],CnbL [3][3]; 
    double Ab1 [3],Wb1 [3];

    Xn_GPS [0] = Latitude;
    Xn_GPS [1] = Longitude;
    Xn_GPS [2] = Al [DelIMU];
    
    RDC.Comvers_GPS(Xn_GPS,Xn_GPS);

    Vn_GPS [0] = 0.514444*Speed*cos(Course/180*Pi);
    Vn_GPS [1] = 0.514444*Speed*sin(Course/180*Pi);
    Vn_GPS [2] = - (Al [DelIMU] - Al [DelIMU - 8]) * 8.0;
    
    /*
    Mn_GPS [0] = 0;
    Mn_GPS [1] = 0;
    Mn_GPS [2] = 1;
    */

    Kalman.MeasurementUpdate(Xhat, P, Xn1, Vn1, Cbn1, Mb, Xn_GPS, Vn_GPS, Mn_GPS, qx_GPS, qv_GPS, qm_GPS);
    Kalman.State_MeasurementUpdate(Xhat, Xn1, Vn1, Cbn1, Cnb1, qbn1, A_Bias, G_Bias);
    
    for(int DelStep = DelIMU;DelStep < IMU_Num;DelStep ++){
        double A_Raw [3],G_Raw [3];
        A_Raw [0] = - Az [DelStep];
        A_Raw [1] = - Ax [DelStep];
        A_Raw [2] = - Ay [DelStep];
        G_Raw [0] = - Gz [DelStep];
        G_Raw [1] = - Gx [DelStep];
        G_Raw [2] = - Gy [DelStep];
        
        RDC.Convers_Accel(0.0,A_Raw,A_Bias,Ab1);
        RDC.Convers_Gyro(0.0,G_Raw,G_Bias,Wb1);
    
        INS.INS_TimeUpdate(delta_t,Xn1,Vn1,Ab1,Wb1,qbn1,Cbn1,Cnb1,XnL,VnL,qbnL,CbnL,CnbL);
        INS.INS_DataUpdate(Xn1,Vn1,qbn1,Cbn1,Cnb1,XnL,VnL,qbnL,CbnL,CnbL);
    }
    INS.INS_DataUpdate(Xn,Vn,qbn,Cbn,Cnb,Xn1,Vn1,qbn1,Cbn1,Cnb1);
    
}
void Receive_GPS(){
    char NMEA_s = GPS.getc();
    //pc.printf("%c",NMEA_s);
    
    if(ReceiveON == 0 && NMEA_s == '$'){
        NMEA [i_nmea] = NMEA_s;
        ReceiveON = 1;
        i_nmea ++;
        
    }else if(ReceiveON == 1){
        NMEA [i_nmea] = NMEA_s;
        i_nmea ++;
        //pc.printf("%c",NMEA_s);
        if(NMEA_s == ','){
            
            if(strcmp(NMEA,"$GPGGA,") == 0){
                ReceiveON = 3;
                pc.printf("X,");
            }else if(strcmp(NMEA,"$GPRMC,") == 0){
                ReceiveON = 2;
                pc.printf("V,");
            }else{
                ReceiveON = 0;
                for(int i=0;i < i_nmea;i++){
                    NMEA [i] = NULL;
                }
                i_nmea = 0;
            }
        }
    }else if(ReceiveON == 2){
        NMEA [i_nmea] = NMEA_s;
        i_nmea ++;
        int NMEA_SIZE = i_nmea;
        if(NMEA_s == '$'){
            float RMCXX;
            if(sscanf(NMEA,"$GPRMC,%f,%c,%f,%c,%f,%c,%f,%f,%f",&UTC_Time,&Latitude,&C_Latitude,&Longitude,&C_Longitude,&RMCXX,&Speed,&Course) >= 1){
                                        //pc.printf("AAAAAAA ");
                                        pc.printf("%f,%f,\r\n",Speed,Course);
                                        //pc.printf("19 ");
            }
            ReceiveON = 0;
            for(int i=0;i<NMEA_SIZE;i++){
                NMEA [i] = NULL;
            }
            i_nmea = 0;

        }
        
        
    }else if(ReceiveON == 3){
        NMEA [i_nmea] = NMEA_s;
        i_nmea ++;
        int NMEA_SIZE = i_nmea;
        if(NMEA_s == '$'){
            //pc.printf("BBBBBBB ");
            if(sscanf(NMEA,"$GPGGA,%f,%f,%c,%f,%c,%d,%d,%f,%f,%c,%f,%c",
            &UTC_Time,&Latitude,&C_Latitude,&Longitude,&C_Longitude,&PFI,&Satellite,&HDOP,&MSL_Alti,&C_MSL,&Geoid_Sepa,&C_Geoid) >= 1){
                                        
                                        pc.printf("%f,%f,%f,%f,%f,",UTC_Time,Latitude,Longitude,MSL_Alti,Geoid_Sepa);
                                        //pc.printf("19 ");
            }
            ReceiveON = 0;
            for(int i=0;i<NMEA_SIZE;i++){
                NMEA [i] = NULL;
            }
            i_nmea = 0;
 
        }
        
    }else{
        ReceiveON = 0;
        i_nmea = 0;
    } 

}

void Receive_Command(){
    int InvalidCOM = 1;
    
    if(RCOM == 0){
        COM [0] = pc.getc();
        if(COM [0] == '$'){
            RCOM = 1;
            icom = 1;
            COM_Time.start();
            pc.printf("%c",COM [0]);
            InvalidCOM = 0;
        }else{
            pc.printf("%c",COM [0]);
            pc.printf("\r\nInvalid  %d  %f\r\n",icom,COM_Time.read());
        }
    }else if(RCOM == 1 && COM_Time.read() < 10.0){
        COM [icom] = pc.getc();
        if(COM [icom] == '#'){
            pc.printf("%c\r\n",COM [icom]);
            InvalidCOM = 0;
            RCOM = 0;
            icom = 0;
            Command_Library();
            COM_Time.stop();
            COM_Time.reset();
        }else if(icom == 9){
            pc.printf("%c",COM [9]);
            pc.printf("\r\nInvalid  %d  %f\r\n",icom,COM_Time.read());
        }else{
            pc.printf("%c",COM [icom]);
            InvalidCOM = 0;
            icom ++;
        }
    }else if(RCOM == 1 && COM_Time.read() >= 10.0){
        pc.printf("%c",pc.getc());
        pc.printf("\r\nInvalid Time Out\r\n");
    }else{
        pc.printf("%c",pc.getc());
        pc.printf("\r\nInvalid  %d  %f\r\n",icom,COM_Time.read());
    }
    
    if(InvalidCOM == 1){
        RCOM = 0;
        COM_Time.stop();
        COM_Time.reset();
        //char COM [10] = {};
        for(int i=0;i<9;i++){
            COM [i] = NULL;
        }
    }
    
}


void Command_Library(){
    int EffectiveCOM = 0;

    //Command = "$CFM#" : Change Flight Mode
    if(strcmp(COM,"$CFM#") == 0){
        pc.printf("Change Flight Mode\r\n");
        //R_MODE = 10;
        EffectiveCOM = 1;
    }
    
    
    if(EffectiveCOM == 0){
        pc.printf("\r\nInvalid  %d  %f\r\n",icom,COM_Time.read());
    }else if(EffectiveCOM == 1){
        for(int i=0;i<9;i++){
            COM [i] = NULL;
        }
    }
}