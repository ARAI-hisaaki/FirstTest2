//Class_Inclide
#include "mbed.h"
#include "CalculeMatrix.h"
#include "Strapdown_INS.h"
#include "KalmanFilter.h"
#include "Initialization.h"

//Function
CalculeMatrix Cal_Mat;
Initialization Initi_kal;

void KalmanFilter::PredictStep(double Xhat [16][16],double P [16][16],double T,double Xn [3],double Vn [3],
double Ab [3],double Cbn [3][3],double Tba,double Tbg,double qa,double qg){
    
    /*
     * 時刻t+k/2のINSの状態を代入し，
     * 時刻ｔで事後推定値だったXhatに，時刻t+kの事前推定値を代入.
     * 時刻ｔで事後誤差共分散行列だったPに，時刻t+kの事前誤差共分散行列を代入.
     *
     * k=1の場合，つまりINS時間更新と同周期で予測ステップを繰り返す場合は，時刻tのINS状態量を使う
     *
     * INS状態量の誤差をカルマンフィルタの状態量にしているため，INS状態量の観測更新後のXhatの成分はすべて0
     * よって，本来はXhatの時間更新の計算はいらない
     */

    double PrioriXhat [16][16];
    double Phi [16][16],Phi_T [16][16],Q [16][16];
    double PrioriP_A [16][16],PrioriP_A2 [16][16];
    
    for(int i=0;i<16;i++){
        for(int j=0;j<16;j++){
            PrioriXhat [i][j] = 0.0;
            Phi [i][j] = 0.0;
            Phi_T [i][j] = 0.0;
            Q [i][j] = 0.0;
            PrioriP_A [i][j] = 0.0;
            PrioriP_A2 [i][j] = 0.0;
        }
    }
    
    PrioriXhat [1][0] = 15;PrioriXhat [0][1] = 1;
    Phi [1][0] = 15;Phi [0][1] = 15;
    Phi_T [1][0] = 15;Phi_T [0][1] = 15;
    Q [1][0] = 15;Q [0][1] = 15;
    PrioriP_A [1][0] = 15;PrioriP_A [0][1] = 15;
    PrioriP_A2 [1][0] = 15;PrioriP_A2 [0][1] = 15;
    
    
    //計算に必要な各行列の算出

    ErrorDynamics_Element(Phi,T,Xn,Vn,Ab,Cbn,Tba,Tbg);

    Noise_CovarianceMatrix(Q,T,Cbn,Tba,Tbg,qa,qg);
    
    //事前推定値PrioriXhatの算出，およびXhatへの代入
    //Cal_Mat.Multi_Mat(Phi, Xhat, PrioriXhat);
    //Cal_Mat.Addition_Mat(0, PrioriXhat, 1, PrioriXhat, Xhat);

    //事前誤差共分散行列Pの算出
    Cal_Mat.Multi_Mat(Phi, P, PrioriP_A);
    Cal_Mat.Trans_Mat(Phi, Phi_T);
    Cal_Mat.Multi_Mat(PrioriP_A, Phi_T, PrioriP_A2);
    Cal_Mat.Addition_Mat(1, PrioriP_A2, 1, Q, P);
    
    
}

void KalmanFilter::MeasurementUpdate(double Xhat [16][16],double P [16][16],double Xn [3],double  Vn [3],double Cbn [3][3],
double Mb [3],double Xn_GPS [3],double Vn_GPS [3],double Mn_GPS [3],double qx_GPS [3],double qv_GPS [3],double qm_GPS [3]){

    /*
     * 座標の観測値Xn_GPSは，緯度経度はGPS，高度は気圧
     * 速度の観測値Vn_GPSは，緯度経度方向はGPS，高度方向は気圧
     * 磁力の観測値Mn_GPSは，国土地理院の地磁気地図より求めたデータを用いる
     *
     * よって，これらの分散，qx_GPS，qv_GPS，qm_GPSは三軸ごとに異なる．
     */

    double y [16][16],G [16][16],C [16][16],C_T [16][16],R [16][16];
    double PosterioriXhat [16][16];
    double PosterioriP [16][16];

    double PC_T [16][16],G_A[16][16],G_A2 [16][16],G_A2inv [16][16];
    double CXhat [16][16],PosterioriXhat_A [16][16],PosterioriXhat_A2 [16][16];
    double GC [16][16],PosterioriP_A [16][16];

    double II [16][16];
    
    for(int i=0;i<16;i++){
        for(int j=0;j<16;j++){
            y [i][j] = 0.0;
            G [i][j] = 0.0;
            C [i][j] = 0.0;
            C_T [i][j] = 0.0;
            R [i][j] = 0.0;
            PosterioriXhat [i][j] = 0.0;
            PosterioriP [i][j] = 0.0;
            PC_T [i][j] = 0.0;
            G_A [i][j] = 0.0;
            G_A2 [i][j] = 0.0;
            G_A2inv [i][j] = 0.0;
            CXhat [i][j] = 0.0;
            PosterioriXhat_A [i][j] = 0.0;
            PosterioriXhat_A2 [i][j] = 0.0;
            GC [i][j] = 0.0;
            PosterioriP_A [i][j] = 0.0;
            II [i][j] = 0.0;
        }
    }
    

    //地磁気有り
    /*
    y [1][0] = 9;y [0][1] = 1;
    G [1][0] = 15;G [0][1] = 9;
    C [1][0] = 9;C [0][1] = 15;
    C_T [1][0] = 15;C_T [0][1] = 9;
    R [1][0] = 9;R [0][1] = 9;
    */


    //地磁気なし
    y [1][0] = 6;y [0][1] = 1;
    G [1][0] = 15;G [0][1] = 6;
    C [1][0] = 6;C [0][1] = 15;
    C_T [1][0] = 15;C_T [0][1] = 6;
    R [1][0] = 6;R [0][1] = 6;



    PosterioriXhat [1][0] = 15;PosterioriXhat [0][1] = 1;
    PosterioriP [1][0] = 15;PosterioriP [0][1] = 15;


    //計算に必要な各行列の算出

    II [1][0] = 15;II [0][1] = 15;

    for(int i=1;i<=15;i++){
        II [i][i] = 1;
    }

    MeasurementVector(y,Xn,Vn,Cbn,Mb,Xn_GPS,Vn_GPS,Mn_GPS);

    MeasurementDesign_Matrix(C,Cbn,Mb);

    Cal_Mat.Trans_Mat(C, C_T);

    MeasurementNoise_CovarianceMatrix(R,qx_GPS,qv_GPS,qm_GPS);


    //カルマンゲインGの算出
    Cal_Mat.Multi_Mat(P, C_T, PC_T);
    Cal_Mat.Multi_Mat(C, PC_T, G_A);
    Cal_Mat.Addition_Mat(1, G_A, 1, R, G_A2);
    Cal_Mat.Inv_Mat(G_A2, G_A2inv);
    Cal_Mat.Multi_Mat(PC_T, G_A2inv, G);

    //事後推定値PosterioriXhatの算出，およびXhatへの代入
    Cal_Mat.Multi_Mat(C, Xhat, CXhat);
    Cal_Mat.Addition_Mat(1, y, -1, CXhat, PosterioriXhat_A);
    Cal_Mat.Multi_Mat(G, PosterioriXhat_A, PosterioriXhat_A2);
    Cal_Mat.Addition_Mat(1, Xhat, 1, PosterioriXhat_A2, PosterioriXhat);
    Cal_Mat.Addition_Mat(0, PosterioriXhat, 1, PosterioriXhat, Xhat);

    //事後誤差共分散行列PosterioriPの算出，およびPへの代入
    Cal_Mat.Multi_Mat(G, C, GC);
    Cal_Mat.Addition_Mat(1, II, -1, GC, PosterioriP_A);
    Cal_Mat.Multi_Mat(PosterioriP_A, P, PosterioriP);
    Cal_Mat.Addition_Mat(0, II, 1, PosterioriP, P);

}

void KalmanFilter::State_MeasurementUpdate(double Xhat [16][16],double Xn [3],double Vn [3],double Cbn [3][3],double Cnb [3][3],double qbn [4],double A_Bias [3],double G_Bias [3]){
    /*
     * 推定されたINS状態量の誤差を修正し，それに伴い，カルマンフィルタの状態量を0にする
     */

    double Cqbn [16][16],Cqbn_inv [16][16],Delta_qbn [16][16],Delta_e [16][16];
    double q1_tiX2,q2_tiX2,q3_tiX2,q4_tiX2;
    double qbn_nolm;

    for(int i=0;i<16;i++){
        for(int j=0;j<16;j++){
            Cqbn [i][j] = 0.0;
            Cqbn_inv [i][j] = 0.0;
            Delta_qbn [i][j] = 0.0;
            Delta_e [i][j] = 0.0;
        }
    }
    
    Xn [0] = Xn [0] + Xhat [1][1];
    Xn [1] = Xn [1] + Xhat [2][1];
    Xn [2] = Xn [2] + Xhat [3][1];

    Vn [0] = Vn [0] + Xhat [4][1];
    Vn [1] = Vn [1] + Xhat [5][1];
    Vn [2] = Vn [2] + Xhat [6][1];

    q1_tiX2 = 2*qbn [0];
    q2_tiX2 = 2*qbn [1];
    q3_tiX2 = 2*qbn [2];
    q4_tiX2 = 2*qbn [3];

    Cqbn [1][0] = 4;    Cqbn [0][1] = 4;

    Cqbn [1][1] = q4_tiX2;      Cqbn [1][2] = - q3_tiX2;    Cqbn [1][3] = q2_tiX2;      Cqbn [1][4] = - q1_tiX2;
    Cqbn [2][1] = q3_tiX2;      Cqbn [2][2] = q4_tiX2;      Cqbn [2][3] = - q1_tiX2;    Cqbn [2][4] = - q2_tiX2;
    Cqbn [3][1] = - q2_tiX2;    Cqbn [3][2] = q1_tiX2;      Cqbn [3][3] = q4_tiX2;      Cqbn [3][4] = - q3_tiX2;
    Cqbn [4][1] = q1_tiX2;      Cqbn [4][2] = q2_tiX2;      Cqbn [4][3] = q3_tiX2;      Cqbn [4][4] = q4_tiX2;

    Cal_Mat.Inv_Mat(Cqbn, Cqbn_inv);

    Delta_e [1][0] = 4;Delta_e [0][1] = 1;

    Delta_e [1][1] = Xhat [7][1];
    Delta_e [2][1] = Xhat [8][1];
    Delta_e [3][1] = Xhat [9][1];
    Delta_e [4][1] = 0;

    Cal_Mat.Multi_Mat(Cqbn_inv,Delta_e,Delta_qbn);

    qbn [0] =qbn [0] +  Delta_qbn [1][1];
    qbn [1] =qbn [1] +  Delta_qbn [2][1];
    qbn [2] =qbn [2] +  Delta_qbn [3][1];
    qbn [3] =qbn [3] +  Delta_qbn [4][1];

    qbn_nolm = sqrt(qbn [0]*qbn [0] + qbn [1]*qbn [1] + qbn [2]*qbn [2] + qbn [3]*qbn [3]);

    qbn [0] = qbn [0]/qbn_nolm;
    qbn [1] = qbn [1]/qbn_nolm;
    qbn [2] = qbn [2]/qbn_nolm;
    qbn [3] = qbn [3]/qbn_nolm;

    Initi_kal.ComputeDCM_FromQuaternion(qbn, Cbn, Cnb);

    A_Bias [0] = A_Bias [0] + Xhat [10][1];
    A_Bias [1] = A_Bias [1] + Xhat [11][1];
    A_Bias [2] = A_Bias [2] + Xhat [12][1];

    G_Bias [0] = G_Bias [0] + Xhat [13][1];
    G_Bias [1] = G_Bias [1] + Xhat [14][1];
    G_Bias [2] = G_Bias [2] + Xhat [15][1];

    Xhat [1][1] = 0.0;
    Xhat [2][1] = 0.0;
    Xhat [3][1] = 0.0;
    Xhat [4][1] = 0.0;
    Xhat [5][1] = 0.0;
    Xhat [6][1] = 0.0;
    Xhat [7][1] = 0.0;
    Xhat [8][1] = 0.0;
    Xhat [9][1] = 0.0;
    Xhat [10][1] = 0.0;
    Xhat [11][1] = 0.0;
    Xhat [12][1] = 0.0;
    Xhat [13][1] = 0.0;
    Xhat [14][1] = 0.0;
    Xhat [15][1] = 0.0;

}

void KalmanFilter::ErrorDynamics_Element(double Phi [16][16],double T,double Xn [3],double Vn [3],double Ab [3],double Cbn [3][3],double Tba,double Tbg){

    double Frr [4][4],Frv [4][4];//d_Xn
    double Fvr [4][4],Fvv [4][4],Fve [4][4],Fvb [4][4];//d_Vn
    double Fer [4][4],Fev [4][4],Fee [4][4],Feb [4][4];//d_Ebn  (DCM誤差の表現法)
    double Fba [4][4];//d_Bacc
    double Fbg [4][4];//d_Bgyro

    double L,h;//Xn
    double vN,vE,vD;//Vn
    double aN,aE,aD;//Ab
    double Wnx,Wny,Wnz;//航法座標系の慣性空間に対する回転（高速で緯度経度が変化すると発生）

    double r,G,m1,mu;
    double ee,r_m,r_p;

    double Cos_L,Sin_L,Tan_L;
    
    L = Xn [0];h = Xn [2];
    vN = Vn [0];vE = Vn [1];vD = Vn [2];
    aN = Cbn [0][0]*Ab [0] + Cbn [0][1]*Ab [1] + Cbn [0][2]*Ab [2];
    aE = Cbn [1][0]*Ab [0] + Cbn [1][1]*Ab [1] + Cbn [1][2]*Ab [2];
    aD = Cbn [2][0]*Ab [0] + Cbn [2][1]*Ab [1] + Cbn [2][2]*Ab [2];

    Cos_L = cos(L);Sin_L = sin(L);Tan_L = tan(L);

    ee = f*(2.0 - f);
    r_m = r_e*(1.0 - ee)/pow((1.0 - ee*Sin_L*Sin_L),1.5);//大圏半径
    r_p = r_e/sqrt(1.0 - ee*Sin_L*Sin_L);//子午半径

    Wnx = vE/(r_p+h) + 2*We*Cos_L;
    Wny = - vN/(r_m + h);
    Wnz = - vE/(r_p + h)*Tan_L -2*We*Sin_L;

    r = r_e + h;
    G = 6.672598/e6/e3/100.0;
    m1 = 5.9736*e6*e6*e6*e6;
    mu = G*m1;


    //以下，誤差に関するシステムダイナミクス行列の各成分の行列

    Frr [1][0] = 3;Frr [0][1] = 3;

    Frr [1][1] = 0;                             Frr [1][2] = 0; Frr [1][3] = - vN/( r*r );
    Frr [2][1] = vE*Sin_L/( r*Cos_L*Cos_L );    Frr [2][2] = 0; Frr [2][3] = - vE/( r*r*Cos_L );
    Frr [3][1] = 0;                             Frr [3][2] = 0; Frr [3][3] = 0;


    Frv [1][0] = 3;Frv [0][1] = 3;

    Frv [1][1] = 1.0/r; Frv [1][2] = 0;                 Frv [1][3] = 0;
    Frv [2][1] = 0;     Frv [2][2] = 1.0/( r*Cos_L );   Frv [2][3] = 0;
    Frv [3][1] = 0;     Frv [3][2] = 0;                 Frv [3][3] = -1.0;


    Fvr [1][0] = 3;Fvr [0][1] = 3;

    Fvr [1][1] = - vE*( vE/(r*Cos_L*Cos_L) + 2.0*We*Cos_L ) - r*We*We*cos(2.0*L);  Fvr [1][2] = 0; Fvr [1][3] = (vE*Tan_L*Tan_L - vN*vD)/( r*r ) - We*We*sin(2.0*L)/2.0;
    Fvr [2][1] = vN*( vE/(r*Cos_L*Cos_L) + 2.0*We*Cos_L ) - 2.0*vD*We*Sin_L;       Fvr [2][2] = 0; Fvr [2][3] = - (vN*vE*Tan_L + vE*vD)/( r*r );
    Fvr [3][1] = 2.0*vE*We*Sin_L + r*We*We*sin(2.0*L);                             Fvr [3][2] = 0; Fvr [3][3] = (vN*vN + vE*vE)/( r*r ) - We*We*Cos_L*Cos_L - 2.0*mu/( r*r*r );


    Fvv [1][0] = 3;Fvv [0][1] = 3;

    Fvv [1][1] = vD/r;                      Fvv [1][2] = - 2.0*( vE*Tan_L/r + We*Sin_L );   Fvv [1][3] = vN/r;
    Fvv [2][1] = vE*Tan_L/r + 2.0*We*Sin_L; Fvv [2][2] = ( vD+vN*Tan_L )/r;                 Fvv [2][3] = vE/r + 2.0*We*Cos_L;
    Fvv [3][1] = - 2.0*vN/r;                Fvv [3][2] = - 2.0*( vE/r + We*Cos_L );         Fvv [3][3] = 0;


    Fve [1][0] = 3;Fve [0][1] = 3;

    Fve [1][1] = 0;     Fve [1][2] = aD;    Fve [1][3] = - aE;
    Fve [2][1] = -aD;   Fve [2][2] = 0;     Fve [2][3] = aN;
    Fve [3][1] = aE;    Fve [3][2] = - aN;  Fve [3][3] = 0;


    Fvb [1][0] = 3;Fvb [0][1] = 3;

    Fvb [1][1] = Cbn [0][0];    Fvb [1][2] = Cbn [0][1];    Fvb [1][3] = Cbn [0][2];
    Fvb [2][1] = Cbn [1][0];    Fvb [2][2] = Cbn [1][1];    Fvb [2][3] = Cbn [1][2];
    Fvb [3][1] = Cbn [2][0];    Fvb [3][2] = Cbn [2][1];    Fvb [3][3] = Cbn [2][2];


    Fer [1][0] = 3;Fer [0][1] = 3;

    Fer [1][1] = We*Sin_L;                              Fer [1][2] = 0; Fer [1][3] = vE/( r*r );
    Fer [2][1] = 0;                                     Fer [2][2] = 0; Fer [2][3] = - vN/( r*r );
    Fer [3][1] = ( vE + vE*Tan_L*Tan_L )/r + We*Cos_L;  Fer [3][2] = 0; Fer [3][3] = - vE*Tan_L*Tan_L/( r*r );


    Fev [1][0] = 3;Fev [0][1] = 3;

    Fev [1][1] = 0;     Fev [1][2] = - 1/r;     Fev [1][3] = 0;
    Fev [2][1] = 1/r;   Fev [2][2] = 0;         Fev [2][3] = 0;
    Fev [3][1] = 0;     Fev [3][2] = Tan_L/r;   Fev [3][3] = 0;


    Fee [1][0] = 3;Fee [0][1] = 3;

    Fee [1][1] = 0;     Fee [1][2] = Wnz;   Fee [1][3] = - Wny;
    Fee [2][1] = - Wnz; Fee [2][2] = 0;     Fee [2][3] = Wnx;
    Fee [3][1] = Wny;   Fee [3][2] = - Wnx; Fee [3][3] = 0;


    Feb [1][0] = 3;Feb [0][1] = 3;

    Feb [1][1] = Cbn [0][0];    Feb [1][2] = Cbn [0][1];    Feb [1][3] = Cbn [0][2];
    Feb [2][1] = Cbn [1][0];    Feb [2][2] = Cbn [1][1];    Feb [2][3] = Cbn [1][2];
    Feb [3][1] = Cbn [2][0];    Feb [3][2] = Cbn [2][1];    Feb [3][3] = Cbn [2][2];


    Fba [1][0] = 3;Fba [0][1] = 3;

    Fba [1][1] = - 1.0/Tba; Fba [1][2] = 0;         Fba [1][3] = 0;
    Fba [2][1] = 0;         Fba [2][2] = - 1.0/Tba; Fba [2][3] = 0;
    Fba [3][1] = 0;         Fba [3][2] = 0;         Fba [3][3] = - 1.0/Tba;


    Fbg [1][0] = 3;Fbg [0][1] = 3;

    Fbg [1][1] = - 1.0/Tbg; Fbg [1][2] = 0;         Fbg [1][3] = 0;
    Fbg [2][1] = 0;         Fbg [2][2] = - 1.0/Tbg; Fbg [2][3] = 0;
    Fbg [3][1] = 0;         Fbg [3][2] = 0;         Fbg [3][3] = - 1.0/Tbg;

    ErrorDynamics_Transition(Phi,T,Frr,Frv,Fvr,Fvv,Fve,Fvb,Fer,Fev,Fee,Feb,Tba,Tbg);
}

void KalmanFilter::ErrorDynamics_Transition(double Phi [16][16],double T,double Frr [4][4],double Frv [4][4],double Fvr [4][4],double Fvv [4][4],
            double Fve [4][4],double Fvb [4][4],double Fer [4][4],double Fev [4][4],double Fee [4][4],double Feb [4][4],double Tba,double Tbg){


    /*
     * システムダイナミクスから状態遷移行列を算出
     */

    double Frr1 [4][4],Frv1 [4][4];//d_Xn
    double Fvr1 [4][4],Fvv1 [4][4],Fve1 [4][4],Fvb1 [4][4];//d_Vn
    double Fer1 [4][4],Fev1 [4][4],Fee1 [4][4],Feb1 [4][4];//d_Ebn  (DCM誤差の表現法)
    double Fba1 [4][4];//d_Bacc
    double Fbg1 [4][4];//d_Bgyro

    double Exp_Tba,Exp_Tbg;
    double I  [4][4];

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            Frr1 [i][j] = 0.0;
            Frv1 [i][j] = 0.0;
            Fvr1 [i][j] = 0.0;
            Fvv1 [i][j] = 0.0;
            Fve1 [i][j] = 0.0;
            Fvb1 [i][j] = 0.0;
            Fer1 [i][j] = 0.0;
            Fev1 [i][j] = 0.0;
            Fee1 [i][j] = 0.0;
            Feb1 [i][j] = 0.0;
            Fba1 [i][j] = 0.0;
            Fbg1 [i][j] = 0.0;
            I [i][j] = 0.0;
        }
    }

    Phi [1][0] = 15.0;Phi [0][1] = 15.0;

    I [1][0] = 3.0;I [0][1] = 3.0;
    I [1][1] = 1.0;I [2][2] = 1.0;I [3][3] = 1.0;

    Exp_Tba = pow(e, - T/Tba);
    Exp_Tbg = pow(e, - T/Tbg);


    //システムダイナミクスの各要素から，遷移行列Phiの各要素を作成

    Cal_Mat.Addition_Mat3x3(1, I, T, Frr, Frr1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Frv, Frv1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Fvr, Fvr1);
    Cal_Mat.Addition_Mat3x3(1, I, T, Fvv, Fvv1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Fve, Fve1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Fvb, Fvb1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Fer, Fer1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Fev, Fev1);
    Cal_Mat.Addition_Mat3x3(1, I, T, Fee, Fee1);
    Cal_Mat.Addition_Mat3x3(0, I, T, Feb, Feb1);


    Fba1 [1][0] = 3;Fba1 [0][1] = 3;

    Fba1 [1][1] = Exp_Tba;  Fba1 [1][2] = 0;        Fba1 [1][3] = 0;
    Fba1 [2][1] = 0;        Fba1 [2][2] = Exp_Tba;  Fba1 [2][3] = 0;
    Fba1 [3][1] = 0;        Fba1 [3][2] = 0;        Fba1 [3][3] = Exp_Tba;


    Fbg1 [1][0] = 3;Fbg1 [0][1] = 3;

    Fbg1 [1][1] = Exp_Tbg;  Fbg1 [1][2] = 0;        Fbg1 [1][3] = 0;
    Fbg1 [2][1] = 0;        Fbg1 [2][2] = Exp_Tbg;  Fbg1 [2][3] = 0;
    Fbg1 [3][1] = 0;        Fbg1 [3][2] = 0;        Fbg1 [3][3] = Exp_Tbg;

    //遷移行列Phiの各要素を，Phiに代入

    Cal_Mat.Fusion_3x3Mat(Phi, 1, 1, Frr1);
    Cal_Mat.Fusion_3x3Mat(Phi, 1, 2, Frv1);
    Cal_Mat.Fusion_3x3Mat(Phi, 2, 1, Fvr1);
    Cal_Mat.Fusion_3x3Mat(Phi, 2, 2, Fvv1);
    Cal_Mat.Fusion_3x3Mat(Phi, 2, 3, Fve1);
    Cal_Mat.Fusion_3x3Mat(Phi, 2, 4, Fvb1);
    Cal_Mat.Fusion_3x3Mat(Phi, 3, 1, Fer1);
    Cal_Mat.Fusion_3x3Mat(Phi, 3, 2, Fev1);
    Cal_Mat.Fusion_3x3Mat(Phi, 3, 3, Fee1);
    Cal_Mat.Fusion_3x3Mat(Phi, 3, 5, Feb1);
    Cal_Mat.Fusion_3x3Mat(Phi, 4, 4, Fba1);
    Cal_Mat.Fusion_3x3Mat(Phi, 5, 5, Fbg1);
}

void KalmanFilter::Noise_CovarianceMatrix(double Q [16][16],double T,double Cbn [3][3],double Tba,double Tbg,double qa,double qg){
        
    //Qは時刻tから時刻t+kまでのプロセスノイズの共分散行列

    double Qvv [4][4],Qvb [4][4];
    double Qee [4][4],Qeb [4][4];
    double Qvb_T [4][4],Qba [4][4];
    double Qeb_T [4][4],Qbg [4][4];

    double Qvv_A,Qvb_A;
    double Qee_A,Qeb_A;
    double Qba_A;
    double Qbg_A;

    double I  [4][4];
    double Cbn_X [4][4];

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            Qvv [i][j] = 0.0;
            Qvb [i][j] = 0.0;
            Qee [i][j] = 0.0;
            Qeb [i][j] = 0.0;
            Qvb_T [i][j] = 0.0;
            Qba [i][j] = 0.0;
            Qeb_T [i][j] = 0.0;
            Qbg [i][j] = 0.0;
            I [i][j] = 0.0;
            Cbn_X [i][j] = 0.0;
        }
    }

    Q [1][0] = 15;Q [0][1] = 15;

    I [1][0] = 3.0;I [0][1] = 3.0;
    I [1][1] = 1.0;I [2][2] = 1.0;I [3][3] = 1.0;

    Cbn_X [1][0] = 3.0;Cbn_X [0][1] = 3.0;
        
    Cbn_X [1][1] = Cbn [0][0];Cbn_X [1][2] = Cbn [0][1];Cbn_X [1][3] = Cbn [0][2];
    Cbn_X [2][1] = Cbn [1][0];Cbn_X [2][2] = Cbn [1][1];Cbn_X [2][3] = Cbn [1][2];
    Cbn_X [3][1] = Cbn [2][0];Cbn_X [3][2] = Cbn [2][1];Cbn_X [3][3] = Cbn [2][2];

    //プロセスノイズの共分散行列Qの各要素の作成

    Qvv_A = T*T*T*qa/3.0;
    Qvb_A = T*T*( 1.0/2.0 - T/(Tba*3.0) + T*T/(6.0*Tba*Tba) )*qa;
    Qee_A = T*T*T*qg/3.0;
    Qeb_A = T*T*( 1.0/2.0 - T/(Tbg*3.0) + T*T/(6.0*Tbg*Tbg) )*qg;
    Qba_A = T*( 1.0 - T/Tba + 2.0*T*T/(3.0*Tba*Tba) )*qa;
    Qbg_A = T*( 1.0 - T/Tbg + 2.0*T*T/(3.0*Tbg*Tbg) )*qg;

    Cal_Mat.Addition_Mat3x3(0, I, Qvv_A, I, Qvv);
    Cal_Mat.Addition_Mat3x3(0, I,Qvb_A, Cbn_X, Qvb);
    Cal_Mat.Addition_Mat3x3(0, I, Qee_A, I, Qee);
    Cal_Mat.Addition_Mat3x3(0, I, Qeb_A, Cbn_X, Qeb);
    Cal_Mat.Addition_Mat3x3(0, I, Qba_A, I, Qba);
    Cal_Mat.Addition_Mat3x3(0, I, Qbg_A, I, Qbg);
    Cal_Mat.Trans_Mat3x3(Qvb, Qvb_T);
    Cal_Mat.Trans_Mat3x3(Qeb, Qeb_T);


    //各要素をQに代入

    Cal_Mat.Fusion_3x3Mat(Q, 2, 2, Qvv);
    Cal_Mat.Fusion_3x3Mat(Q, 2, 4, Qvb);
    Cal_Mat.Fusion_3x3Mat(Q, 3, 3, Qee);
    Cal_Mat.Fusion_3x3Mat(Q, 3, 5, Qeb);
    Cal_Mat.Fusion_3x3Mat(Q, 4, 2, Qvb_T);
    Cal_Mat.Fusion_3x3Mat(Q, 4, 4, Qba);
    Cal_Mat.Fusion_3x3Mat(Q, 5, 3, Qeb_T);
    Cal_Mat.Fusion_3x3Mat(Q, 5, 5, Qbg);
}

void KalmanFilter::MeasurementVector(double y [16][16],double Xn [3],double Vn [3],double Cbn [3][3],double Mb [3],
            double Xn_GPS [3],double Vn_GPS [3],double Mn_GPS [3]){

    /*
     * 観測ベクトルの算出
     *
     * Mn_GPSは実際は国土地理院の地磁気地図より求めたデータを用いる
     */

    y [1][1] = Xn_GPS [0] - Xn [0];
    y [2][1] = Xn_GPS [1] - Xn [1];
    y [3][1] = Xn_GPS [2] - Xn [2];

    y [4][1] = Vn_GPS [0] - Vn [0];
    y [5][1] = Vn_GPS [1] - Vn [1];
    y [6][1] = Vn_GPS [2] - Vn [2];

    y [7][1] = Mn_GPS [0] - (Cbn [0][0]*Mb [0] + Cbn [0][1]*Mb [1] + Cbn [0][2]*Mb [2]);
    y [8][1] = Mn_GPS [1] - (Cbn [1][0]*Mb [0] + Cbn [1][1]*Mb [1] + Cbn [1][2]*Mb [2]);
    y [9][1] = Mn_GPS [2] - (Cbn [2][0]*Mb [0] + Cbn [2][1]*Mb [1] + Cbn [2][2]*Mb [2]);

}

void KalmanFilter::MeasurementNoise_CovarianceMatrix(double R [16][16],double qx_GPS [16],double qv_GPS [16],double qm_GPS [16]){

    R [1][1] = qx_GPS [0];
    R [2][2] = qx_GPS [1];
    R [3][3] = qx_GPS [2];

    R [4][4] = qv_GPS [0];
    R [5][5] = qv_GPS [1];
    R [6][6] = qv_GPS [2];

    R [7][7] = qm_GPS [0];
    R [8][8] = qm_GPS [1];
    R [9][9] = qm_GPS [2];
}

void KalmanFilter::MeasurementDesign_Matrix(double C [16][16],double Cbn [3][3],double Mb [3]){

    double x,y,z;
    double A [4][4];

    double I  [4][4];

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            A [i][j] = 0.0;
            I [i][j] = 0.0;
        }
    }
    
    I [1][0] = 3.0;I [0][1] = 3.0;
    I [1][1] = 1.0;I [2][2] = 1.0;I [3][3] = 1.0;

    x = Cbn [0][0]*Mb [0] + Cbn [0][1]*Mb [1] + Cbn [0][2]*Mb [2];
    y = Cbn [1][0]*Mb [0] + Cbn [1][1]*Mb [1] + Cbn [1][2]*Mb [2];
    z = Cbn [2][0]*Mb [0] + Cbn [2][1]*Mb [1] + Cbn [2][2]*Mb [2];

    A [1][0] =  3;A [0][1] = 3;

    A [1][1] = 0;               A [1][2] =  - z*( - 1.0 );  A [1][3] = y*( - 1.0 );
    A [2][1] = z*( - 1.0 );     A [2][2] =  0;              A [2][3] = - x*( - 1.0 );
    A [3][1] = - y*( - 1.0 );   A [3][2] =  x*( - 1.0 );    A [3][3] = 0;


    Cal_Mat.Fusion_3x3Mat(C, 1, 1, I);
    Cal_Mat.Fusion_3x3Mat(C, 2, 2, I);
    Cal_Mat.Fusion_3x3Mat(C, 3, 3, A);
}
