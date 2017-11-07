//Class_Inclide
#include "mbed.h"
#include "Strapdown_INS.h"


double ee,r_m,r_p,L,h,vN,vE;
double Sin_L,Cos_L,Tan_L;

void Strapdown_INS::INS_TimeUpdate(double delta_t,double Xn [3],double Vn [3],double Ab [3],double Wb [3],double qbn [4],double Cbn [3][3],double Cnb [3][3],double Xn1 [3],double Vn1 [3],double qbn1 [4],double Cbn1 [3][3],double Cnb1 [3][3]){
    
    double Gn [3];//Gn0 = gN,Gn1 = gE,Gn2 = gD;
    
    L = Xn [0];//測地緯度，GPSデータは変換の必要がある．
    h = Xn [2];//測地高度，
    vN = Vn [0];
    vE = Vn [1];
    
    Sin_L = sin(L);
    Cos_L = cos(L);
    Tan_L = tan(L);

    ee = f*(2 - f);
    r_m = r_e*(1.0 - ee)/pow((1.0 - ee*Sin_L*Sin_L),1.5);//子午半径
    r_p = r_e/sqrt(1.0 - ee*Sin_L*Sin_L);//大圏半径
    
    Quaternion_TimeUpdate(delta_t,Xn,Vn,Wb,qbn,qbn1);
    DCM_TimeUpdate(qbn1,Cbn1,Cnb1);
    Gravity(Gn,Xn);
    Velocity_TimeUpdate(delta_t,Xn,Vn,Ab,Wb,qbn,Cbn,Gn,Vn1);
    Position_TimeUpdate(delta_t,Xn,Vn,Xn1);    
}

void Strapdown_INS::Quaternion_TimeUpdate(double delta_t,double Xn [3],double Vn [3],double Wb [3],double qbn [4],double qbn1 [4]){
    
    double d_phi_n [3],Wn [3],d_theta_b [3],phi_b_nolm;
    double kappa,lambda;
    double delta [3],epsilon[3];
    double Phi [4][4]; 
    
    //Cooning_Correction
    d_theta_b [0] = Wb [0]*delta_t;
    d_theta_b [1] = Wb [1]*delta_t;
    d_theta_b [2] = Wb [2]*delta_t;

    //クォータニオン時間更新のための状態遷移行列Phiの算出
    phi_b_nolm = sqrt(d_theta_b [0]*d_theta_b [0] + d_theta_b [1]*d_theta_b [1] + d_theta_b [2]*d_theta_b [2]);
    kappa = 1 - phi_b_nolm*phi_b_nolm/24.0;
    lambda = - phi_b_nolm*phi_b_nolm/4.0 + phi_b_nolm*phi_b_nolm*phi_b_nolm*phi_b_nolm/192.0;

    Wn [0] = vE/(r_p+h) + We*Cos_L;
    Wn [1] = - vN/(r_m + h);
    Wn [2] = - vE/(r_p + h)*Tan_L - We*Sin_L;

    d_phi_n [0] = Wn [0]*delta_t;
    d_phi_n [1] = Wn [1]*delta_t;
    d_phi_n [2] = Wn [2]*delta_t;

    delta [0] = kappa*d_theta_b [0] - d_phi_n [0];
    delta [1] = kappa*d_theta_b [1] - d_phi_n [1];
    delta [2] = kappa*d_theta_b [2] - d_phi_n [2];

    epsilon [0] = kappa*d_theta_b [0] + d_phi_n [0];
    epsilon [1] = kappa*d_theta_b [1] + d_phi_n [1];
    epsilon [2] = kappa*d_theta_b [2] + d_phi_n [2];

    Phi [0][0] = (2.0 + lambda)/2.0;    Phi [0][1] = (epsilon [2])/2.0;     Phi [0][2] = ( - epsilon [1])/2.0;  Phi [0][3] = (epsilon [0])/2.0;
    Phi [1][0] = ( - epsilon [2])/2.0;  Phi [1][1] = (2.0 + lambda)/2.0;    Phi [1][2] = (epsilon [0])/2.0;     Phi [1][3] = (epsilon [1])/2.0;
    Phi [2][0] = (epsilon [1])/2.0;     Phi [2][1] = ( - epsilon [0])/2.0;  Phi [2][2] = (2.0 + lambda)/2.0;    Phi [2][3] = (epsilon [2])/2.0;
    Phi [3][0] = ( - delta [0])/2.0;    Phi [3][1] = ( - delta [1])/2.0;    Phi [3][2] = ( - delta [2])/2.0;    Phi [3][3] = (2.0 + lambda)/2.0;

    //クォータニオン時間更新
    qbn1 [0] = Phi [0][0]*qbn [0] + Phi [0][1]*qbn [1] + Phi [0][2]*qbn [2] + Phi [0][3]*qbn [3];
    qbn1 [1] = Phi [1][0]*qbn [0] + Phi [1][1]*qbn [1] + Phi [1][2]*qbn [2] + Phi [1][3]*qbn [3];
    qbn1 [2] = Phi [2][0]*qbn [0] + Phi [2][1]*qbn [1] + Phi [2][2]*qbn [2] + Phi [2][3]*qbn [3];
    qbn1 [3] = Phi [3][0]*qbn [0] + Phi [3][1]*qbn [1] + Phi [3][2]*qbn [2] + Phi [3][3]*qbn [3];   
}


void Strapdown_INS::DCM_TimeUpdate(double qbn1 [4],double Cbn1 [3][3],double Cnb1 [3][3]){
    
    double q1,q2,q3,q4;

    q1 = qbn1 [0];
    q2 = qbn1 [1];
    q3 = qbn1 [2];
    q4 = qbn1 [3];

    Cbn1 [0][0] = q1*q1 - q2*q2 - q3*q3 + q4*q4;    Cbn1 [0][1] = 2*(q1*q2 - q3*q4);                Cbn1 [0][2] = 2*(q3*q1 + q2*q4);
    Cbn1 [1][0] = 2*(q1*q2 + q3*q4);                Cbn1 [1][1] = q2*q2 - q1*q1 - q3*q3 + q4*q4;    Cbn1 [1][2] = 2*(q2*q3 - q1*q4);
    Cbn1 [2][0] = 2*(q3*q1 - q2*q4);                Cbn1 [2][1] = 2*(q2*q3 + q1*q4);                Cbn1 [2][2] = q3*q3 - q1*q1 - q2*q2 + q4*q4;

    Cnb1 [0][0] = q1*q1 - q2*q2 - q3*q3 + q4*q4;    Cnb1 [0][1] = 2*(q1*q2 + q3*q4);                Cnb1 [0][2] = 2*(q3*q1 - q2*q4);
    Cnb1 [1][0] = 2*(q1*q2 - q3*q4);                Cnb1 [1][1] = q2*q2 - q1*q1 - q3*q3 + q4*q4;    Cnb1 [1][2] = 2*(q2*q3 + q1*q4);
    Cnb1 [2][0] = 2*(q3*q1 + q2*q4);                Cnb1 [2][1] = 2*(q2*q3 - q1*q4);                Cnb1 [2][2] = q3*q3 - q1*q1 - q2*q2 + q4*q4;
}


void Strapdown_INS::Gravity(double Gn [3],double Xn [3]){

    double g0,g2,g4,gn,gn1,gn2,gn3;

    //重力計算に使われる係数
    g0 = 9.7803267715;
    g2 = 0.0052790414;
    g4 = 0.0000232718;
    gn = 1.63/e6/100.0;
    gn1 = 3.1571/e6/10.0;
    gn2 = 2.1027/e6/e3;
    gn3 = 7.3749/e6/e6/100.0;
    
    Gn [0] = - gn*h*Sin_L*Cos_L;
    Gn [1] = 0;
    Gn [2] = g0*(1.0 + g2*Sin_L*Sin_L + g4*pow(Sin_L,4.0) )*(1.0 - (gn1 - gn2*Sin_L*Sin_L)*h + gn3*h*h);
}


void Strapdown_INS::Velocity_TimeUpdate(double delta_t,double Xn [3],double Vn [3],double Ab [3],double Wb [3],double qbn [4],double Cbn [3][3],double Gn [3],double Vn1 [3]){

    double L_dot,l_dot;
    double We_n [3],delta_Vb [3];
    double Vn1_A [3],Vn1_B [3],Vn1_C [3];

    L_dot = vN/(r_m + h);
    l_dot = vE/(r_p + h)*Cos_L;

    We_n [0] = (l_dot + 2.0*We)*Cos_L;
    We_n [1] = - L_dot;
    We_n [2] = - (l_dot + 2.0*We)*Sin_L;

    delta_Vb [0] = Ab [0]*delta_t;
    delta_Vb [1] = Ab [1]*delta_t;
    delta_Vb [2] = Ab [2]*delta_t;

    Vn1_A [0] = Cbn [0][0]*delta_Vb [0] + Cbn [0][1]*delta_Vb [1] + Cbn [0][2]*delta_Vb [2];
    Vn1_A [1] = Cbn [1][0]*delta_Vb [0] + Cbn [1][1]*delta_Vb [1] + Cbn [1][2]*delta_Vb [2];
    Vn1_A [2] = Cbn [2][0]*delta_Vb [0] + Cbn [2][1]*delta_Vb [1] + Cbn [2][2]*delta_Vb [2];

    Vn1_B [0] = (We_n [1]*Vn [2] - We_n [2]*Vn [1])*delta_t;
    Vn1_B [1] = (We_n [2]*Vn [0] - We_n [0]*Vn [2])*delta_t;
    Vn1_B [2] = (We_n [0]*Vn [1] - We_n [1]*Vn [0])*delta_t;

    Vn1_C [0] = Gn [0]*delta_t;
    Vn1_C [1] = Gn [1]*delta_t;
    Vn1_C [2] = Gn [2]*delta_t;

    Vn1 [0] = Vn1_A [0] - Vn1_B [0] + Vn1_C [0] + Vn [0];
    Vn1 [1] = Vn1_A [1] - Vn1_B [1] + Vn1_C [1] + Vn [1];
    Vn1 [2] = Vn1_A [2] - Vn1_B [2] + Vn1_C [2] + Vn [2];
}


void Strapdown_INS::Position_TimeUpdate(double delta_t,double Xn [3],double Vn [3],double Xn1 [3]){

    Xn1 [0] = Xn [0] + delta_t*Vn [0]/(r_m + h);
    Xn1 [1] = Xn [1] + delta_t*Vn [1]/(r_p + h)/Cos_L;
    Xn1 [2] = Xn [2] - delta_t*Vn [2];
}


void Strapdown_INS::INS_DataUpdate(double Xn [3],double Vn [3],double qbn [4],double Cbn [3][3],double Cnb [3][3],double Xn1 [3],double Vn1 [3],double qbn1 [4],double Cbn1 [3][3],double Cnb1 [3][3]){

    Xn [0] = Xn1 [0];
    Xn [1] = Xn1 [1];
    Xn [2] = Xn1 [2];

    Vn [0] = Vn1 [0];
    Vn [1] = Vn1 [1];
    Vn [2] = Vn1 [2];

    qbn [0] = qbn1 [0];
    qbn [1] = qbn1 [1];
    qbn [2] = qbn1 [2];
    qbn [3] = qbn1 [3];

    Cbn [0][0] = Cbn1 [0][0];Cbn [0][1] = Cbn1 [0][1];Cbn [0][2] = Cbn1 [0][2];
    Cbn [1][0] = Cbn1 [1][0];Cbn [1][1] = Cbn1 [1][1];Cbn [1][2] = Cbn1 [1][2];
    Cbn [2][0] = Cbn1 [2][0];Cbn [2][1] = Cbn1 [2][1];Cbn [2][2] = Cbn1 [2][2];

    Cnb [0][0] = Cnb1 [0][0];Cnb [0][1] = Cnb1 [0][1];Cnb [0][2] = Cnb1 [0][2];
    Cnb [1][0] = Cnb1 [1][0];Cnb [1][1] = Cnb1 [1][1];Cnb [1][2] = Cnb1 [1][2];
    Cnb [2][0] = Cnb1 [2][0];Cnb [2][1] = Cnb1 [2][1];Cnb [2][2] = Cnb1 [2][2];
}

