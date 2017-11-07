//Class_Inclide
#include "mbed.h"
#include "Initialization.h"
#include "CalculeMatrix.h"



CalculeMatrix Cal_Mat_Initia;

void Initialization::ComputeDCM(double Xn [3],double Cbn [3][3],double Ab_Vector [3],double Mb_Vector [3],double Mn_Vector [3]){

    /*
     * 二つの方向ベクトルと，それらの外積で求めたベクトルの，計三つのベクトルの
     * 機体座標系での表現と，航法座標系での表現から
     * DCMを計算
     *
     * しかし，実際はセンサのバイアス誤差等の影響で，
     * 内積が機体座標系と航法座標系で一致しないので，
     * 地磁気センサに問題があるとして，それを修正したうえで行う
     * 修正するベクトル量の算出は，エグイ
     */

    double AbMb_Inner,GnMn_Inner,Delta_MbVector [3],Mb_Vector_Hat [3],MbVectorHat_Nolm;
    double Gn [3],Gn_Vector [3],Gn_nolm;
    double A [16][16],B [16][16],C [16][16],InvA [16][16];//Cbn A = B;
    double AbMb [3],GnMn [3];

    for(int i=0;i<16;i++){
        for(int j=0;j<16;j++){
            A [i][j] = 0.0;
            B [i][j] = 0.0;
            C [i][j] = 0.0;
            InvA [i][j] = 0.0;
        }
    }

    ComputeGravity(Gn,Xn);

    Gn_nolm = sqrt(Gn [0]*Gn [0] + Gn [1]*Gn [1] + Gn [2]*Gn [2]);
    Gn_Vector [0] = Gn [0]/Gn_nolm;
    Gn_Vector [1] = Gn [1]/Gn_nolm;
    Gn_Vector [2] = Gn [2]/Gn_nolm;

    //内積計算
    AbMb_Inner = Ab_Vector [0]*Mb_Vector [0] + Ab_Vector [1]*Mb_Vector [1] + Ab_Vector [2]*Mb_Vector [2];
    GnMn_Inner = Gn_Vector [0]*Mn_Vector [0] + Gn_Vector [1]*Mn_Vector [1] + Gn_Vector [2]*Mn_Vector [2];

    //地磁気ベクトルの修正量算出
    VectorRectification(AbMb_Inner,GnMn_Inner,Ab_Vector,Mb_Vector,Delta_MbVector);

    //地磁気ベクトル修正
    Mb_Vector_Hat [0] = Mb_Vector [0] - Delta_MbVector [0];
    Mb_Vector_Hat [1] = Mb_Vector [1] - Delta_MbVector [1];
    Mb_Vector_Hat [2] = Mb_Vector [2] - Delta_MbVector [2];

    MbVectorHat_Nolm = sqrt( Mb_Vector_Hat [0]*Mb_Vector_Hat [0] + Mb_Vector_Hat [1]*Mb_Vector_Hat [1] + Mb_Vector_Hat [2]*Mb_Vector_Hat [2] );

    Mb_Vector_Hat [0] = Mb_Vector_Hat [0]/MbVectorHat_Nolm;
    Mb_Vector_Hat [1] = Mb_Vector_Hat [1]/MbVectorHat_Nolm;
    Mb_Vector_Hat [2] = Mb_Vector_Hat [2]/MbVectorHat_Nolm;

    //System.out.printf("Nolm = %f    \n",Mb_Vector_Hat [0]*Mb_Vector_Hat [0] + Mb_Vector_Hat [1]*Mb_Vector_Hat [1] + Mb_Vector_Hat [2]*Mb_Vector_Hat [2]);
    //System.out.printf("Inner    AbMb_Hat = %f   GnMn = %f   \n",(Ab_Vector [0]*Mb_Vector_Hat [0] + Ab_Vector [1]*Mb_Vector_Hat [1] + Ab_Vector [2]*Mb_Vector_Hat [2])*e6,GnMn_Inner*e6);

    //外積算出
    AbMb [0] = Ab_Vector [1]*Mb_Vector_Hat [2] - Ab_Vector [2]*Mb_Vector_Hat [1];
    AbMb [1] = Ab_Vector [2]*Mb_Vector_Hat [0] - Ab_Vector [0]*Mb_Vector_Hat [2];
    AbMb [2] = Ab_Vector [0]*Mb_Vector_Hat [1] - Ab_Vector [1]*Mb_Vector_Hat [0];

    GnMn [0] = Gn_Vector [1]*Mn_Vector [2] - Gn_Vector [2]*Mn_Vector [1];
    GnMn [1] = Gn_Vector [2]*Mn_Vector [0] - Gn_Vector [0]*Mn_Vector [2];
    GnMn [2] = Gn_Vector [0]*Mn_Vector [1] - Gn_Vector [1]*Mn_Vector [0];

    //行列生成
    A [1][0] = 3;A [0][1] = 3;

    A [1][1] = Ab_Vector [0];   A [1][2] = Mb_Vector_Hat [0];   A [1][3] = AbMb [0];
    A [2][1] = Ab_Vector [1];   A [2][2] = Mb_Vector_Hat [1];   A [2][3] = AbMb [1];
    A [3][1] = Ab_Vector [2];   A [3][2] = Mb_Vector_Hat [2];   A [3][3] = AbMb [2];


    B [1][0] = 3;B [0][1] = 3;

    B [1][1] = Gn_Vector [0];   B [1][2] = Mn_Vector [0];   B [1][3] = GnMn [0];
    B [2][1] = Gn_Vector [1];   B [2][2] = Mn_Vector [1];   B [2][3] = GnMn [1];
    B [3][1] = Gn_Vector [2];   B [3][2] = Mn_Vector [2];   B [3][3] = GnMn [2];

    //逆行列計算
    Cal_Mat_Initia.Inv_Mat(A,InvA);

    Cal_Mat_Initia.Multi_Mat(B,InvA,C);

    Cbn [0][0] = C [1][1];  Cbn [0][1] = C [1][2];  Cbn [0][2] = C [1][3];
    Cbn [1][0] = C [2][1];  Cbn [1][1] = C [2][2];  Cbn [1][2] = C [2][3];
    Cbn [2][0] = C [3][1];  Cbn [2][1] = C [3][2];  Cbn [2][2] = C [3][3];

}
void Initialization::ComputeDCM_FromAccelAzimuthal(double Ab [3],double Theta,double Xn [3],double Cbn [3][3]){
    /*
    * Thetaは方位角，単位は[deg]．（北が0°，東が90°）
    * まず，重力ベクトルに垂直で，機軸ベクトル1と基底ベクトル1の方向が同じ
    * g座標系での表現をし，それを航法座標系に変換する
    */
    double Gn [3],Gn_Vector [3],Gn_nolm,Cos_TheG;
    double Theta_r,Thetag_r;                    //方位角のrad単位，g座標系での方位角[rad]
    double Cos_T,Sin_T,Cos_Tg,Sin_Tg,Cos_Gg,Sin_Gg;
    double Ab_Vector [3],Ab_nolm;                //観測した重力ベクトルの機体座標系での単位ベクトル
    double Vb1 [3],Vn1 [3];                       //機軸方向ベクトルの機体座標系での表現と，航法座標系での表現
    double Xg1 [3],Xg2 [3],Xg3 [3];        //機体座標系から見た，g座標系の基底ベクトル
    double Xg1_nolm,Xg2_nolm,Xg3_nolm;
    double Cgb [3][3],Cbg [3][3],Cgn [3][3];                          //機体座標系bからg座標系への変換，gからnへの変換
    double A [16][16],B [16][16],C [16][16];//Cbn A = B;

    Theta_r = Theta*Pi/180;
    ComputeGravity(Gn, Xn);

    //航法座標系における単位重力ベクトル
    Gn_nolm = sqrt(Gn [0]*Gn [0] + Gn [1]*Gn [1] + Gn [2]*Gn [2]);
    Gn_Vector [0] = Gn [0]/Gn_nolm;
    Gn_Vector [1] = Gn [1]/Gn_nolm;
    Gn_Vector [2] = Gn [2]/Gn_nolm;

    //重力に垂直な平面と，航法座標系の地平面に垂直な平面の成すCos角の算出
    Cos_TheG = Gn [2]/Gn_nolm;

    //g座標系での方位角[rad]の算出
    Thetag_r = atan(Cos_TheG*tan(Theta_r));

    //機体座標家で観測されたの重力ベクトル
    Ab_nolm = sqrt(Ab [0]*Ab [0] + Ab [1]*Ab [1] + Ab [2]*Ab [2]);
    Ab_Vector [0] = Ab [0]/Ab_nolm;
    Ab_Vector [1] = Ab [1]/Ab_nolm;
    Ab_Vector [2] = Ab [2]/Ab_nolm;

    Vb1 [0] = 1.0;
    Vb1 [1] = 0;
    Vb1 [2] = 0;

    //g座標系での算出
    //Xg3算出
    Xg3 [0] = Ab_Vector [0];
    Xg3 [1] = Ab_Vector [1];
    Xg3 [2] = Ab_Vector [2];

    //Xg2算出
    Xg2 [0] = Ab_Vector [1]*Vb1 [2] - Ab_Vector [2]*Vb1 [1];
    Xg2 [1] = Ab_Vector [2]*Vb1 [0] - Ab_Vector [0]*Vb1 [2];
    Xg2 [2] = Ab_Vector [0]*Vb1 [1] - Ab_Vector [1]*Vb1 [0];
    Xg2_nolm = sqrt(Xg2 [0]*Xg2 [0] + Xg2 [1]*Xg2 [1] + Xg2 [2]*Xg2 [2]);

    Xg2 [0] = Xg2 [0]/Xg2_nolm;
    Xg2 [1] = Xg2 [1]/Xg2_nolm;
    Xg2 [2] = Xg2 [2]/Xg2_nolm;

    //Xg1算出
    Xg1 [0] = Xg2 [1]*Xg3 [2] - Xg2 [2]*Xg3 [1];
    Xg1 [1] = Xg2 [2]*Xg3 [0] - Xg2 [0]*Xg3 [2];
    Xg1 [2] = Xg2 [0]*Xg3 [1] - Xg2 [1]*Xg3 [0];
    Xg1_nolm = sqrt(Xg1 [0]*Xg1 [0] + Xg1 [1]*Xg1 [1] + Xg1 [2]*Xg1 [2]);

    Xg1 [0] = Xg1 [0]/Xg1_nolm;
    Xg1 [1] = Xg1 [1]/Xg1_nolm;
    Xg1 [2] = Xg1 [2]/Xg1_nolm;

    //機体座標系bからg座標系への方向余弦行列
    Cgb [0][0] = Xg1 [0];Cgb [0][1] = Xg2 [0];Cgb [0][2] = Xg3 [0];
    Cgb [1][0] = Xg1 [1];Cgb [1][1] = Xg2 [1];Cgb [1][2] = Xg3 [1];
    Cgb [2][0] = Xg1 [2];Cgb [2][1] = Xg2 [2];Cgb [2][2] = Xg3 [2];

    //重力のベクトルは地面に垂直として計算
    //3-2-1系でg座標系から航法座標系への方向余弦行列
    Cos_T = cos( - Theta_r);
    Sin_T = sin( - Theta_r);

    Cgn [0][0] = Cos_T;     Cgn [0][1] = Sin_T; Cgn [0][2] = 0.0;
    Cgn [1][0] = - Sin_T;   Cgn [1][1] = Cos_T; Cgn [1][2] = 0.0;
    Cgn [2][0] = 0.0;       Cgn [2][1] = 0.0;   Cgn [2][2] = 1.0;

    //行列生成
    A [1][0] = 3;A [0][1] = 3;
    
    A [1][1] = Cgb [0][0];  A [1][2] = Cgb [0][1];  A [1][3] = Cgb [0][2];
    A [2][1] = Cgb [1][0];  A [2][2] = Cgb [1][1];  A [2][3] = Cgb [1][2];
    A [3][1] = Cgb [2][0];  A [3][2] = Cgb [2][1];  A [3][3] = Cgb [2][2];
    
    B [1][0] = 3;B [0][1] = 3;
    
    B [1][1] = Cgn [0][0];  B [1][2] = Cgn [0][1];  B [1][3] = Cgn [0][2];
    B [2][1] = Cgn [1][0];  B [2][2] = Cgn [1][1];  B [2][3] = Cgn [1][2];
    B [3][1] = Cgn [2][0];  B [3][2] = Cgn [2][1];  B [3][3] = Cgn [2][2];

    Cal_Mat_Initia.Inv_Mat(A, A);                //CgbをCbgに変換

    //逆行列計算
    Cal_Mat_Initia.Multi_Mat(B,A,C);

    Cbn [0][0] = C [1][1];  Cbn [0][1] = C [1][2];  Cbn [0][2] = C [1][3];
    Cbn [1][0] = C [2][1];  Cbn [1][1] = C [2][2];  Cbn [1][2] = C [2][3];
    Cbn [2][0] = C [3][1];  Cbn [2][1] = C [3][2];  Cbn [2][2] = C [3][3];
}

void Initialization::ComputeDCM_FromQuaternion(double qbn [4],double Cbn [3][3],double Cnb [3][3]){
    
    double q1,q2,q3,q4;

    q1 = qbn [0];
    q2 = qbn [1];
    q3 = qbn [2];
    q4 = qbn [3];

    Cbn [0][0] = q1*q1 - q2*q2 - q3*q3 + q4*q4;    Cbn [0][1] = 2*(q1*q2 - q3*q4);                Cbn [0][2] = 2*(q3*q1 + q2*q4);
    Cbn [1][0] = 2*(q1*q2 + q3*q4);                Cbn [1][1] = q2*q2 - q1*q1 - q3*q3 + q4*q4;    Cbn [1][2] = 2*(q2*q3 - q1*q4);
    Cbn [2][0] = 2*(q3*q1 - q2*q4);                Cbn [2][1] = 2*(q2*q3 + q1*q4);                Cbn [2][2] = q3*q3 - q1*q1 - q2*q2 + q4*q4;

    Cnb [0][0] = q1*q1 - q2*q2 - q3*q3 + q4*q4;    Cnb [0][1] = 2*(q1*q2 + q3*q4);                Cnb [0][2] = 2*(q3*q1 - q2*q4);
    Cnb [1][0] = 2*(q1*q2 - q3*q4);                Cnb [1][1] = q2*q2 - q1*q1 - q3*q3 + q4*q4;    Cnb [1][2] = 2*(q2*q3 + q1*q4);
    Cnb [2][0] = 2*(q3*q1 + q2*q4);                Cnb [2][1] = 2*(q2*q3 - q1*q4);                Cnb [2][2] = q3*q3 - q1*q1 - q2*q2 + q4*q4;

}

void Initialization::VectorRectification(double AbMb_Inner,double GnMn_Inner,double Ab_Vector [3],double Mb_Vector [3],double Delta_MbVector [3]){
    double X,Y,Z,W;
    double Mx,My,Mz;

    double AbMb_CrossProduct [3],AbMb_dm0 [3],AbMb_dm1 [3];
    double Alpha,Beta,Gamma;
    double A,B,C,D;
    double I,J,K;

    double d_mx0, d_my0, d_mz0;
    double d_mx1, d_my1, d_mz1;

    double Check_Inner0,Check_Inner1;

    X = Ab_Vector [0];
    Y = Ab_Vector [1];
    Z = Ab_Vector [2];
    W = AbMb_Inner - GnMn_Inner;

    Mx = Mb_Vector [0];
    My = Mb_Vector [1];
    Mz = Mb_Vector [2];

    AbMb_CrossProduct [0] = Ab_Vector [1]*Mb_Vector [2] - Ab_Vector [2]*Mb_Vector [1];
    AbMb_CrossProduct [1] = Ab_Vector [2]*Mb_Vector [0] - Ab_Vector [0]*Mb_Vector [2];
    AbMb_CrossProduct [2] = Ab_Vector [0]*Mb_Vector [1] - Ab_Vector [1]*Mb_Vector [0];

    Alpha = AbMb_CrossProduct [0];
    Beta = AbMb_CrossProduct [1];
    Gamma = AbMb_CrossProduct [2];

    A = Beta/(Z*Beta -Y*Gamma);
    B = (X*Beta - Y*Alpha)/(Z*Beta -Y*Gamma);
    C = Gamma/(Z*Beta -Y*Gamma);
    D = (X*Gamma - Z*Alpha)/(Z*Beta -Y*Gamma);

    I = 1.0 + B*B + D*D;
    J = W*(A*B + C*D) + Mx -B*My + D*Mz;
    K = W*W*(A*A + C*C) - 2.0*W*(A*My -C*Mz);

    //d_mx0 = ( J + Math.sqrt(J*J - I*K) )/I;
    //d_mx1 = ( J - Math.sqrt(J*J - I*K) )/I;

    d_mx0 = ( - K)/( - J + sqrt(J*J - I*K) );
    d_mx1 = ( + K)/( + J + sqrt(J*J - I*K) );

    d_my0 = + A*W -B*d_mx0;
    d_my1 = + A*W -B*d_mx1;
    d_mz0 = - C*W -D*d_mx0;
    d_mz1 = - C*W -D*d_mx1;

    AbMb_dm0 [0] = Ab_Vector [1]*(Mb_Vector [2] - d_mz0) - Ab_Vector [2]*(Mb_Vector [1] - d_my0);
    AbMb_dm0 [1] = Ab_Vector [2]*(Mb_Vector [0] - d_mx0) - Ab_Vector [0]*(Mb_Vector [2] - d_mz0);
    AbMb_dm0 [2] = Ab_Vector [0]*(Mb_Vector [1] - d_my0) - Ab_Vector [1]*(Mb_Vector [0] - d_mx0);

    AbMb_dm1 [0] = Ab_Vector [1]*(Mb_Vector [2] - d_mz1) - Ab_Vector [2]*(Mb_Vector [1] - d_my1);
    AbMb_dm1 [1] = Ab_Vector [2]*(Mb_Vector [0] - d_mx1) - Ab_Vector [0]*(Mb_Vector [2] - d_mz1);
    AbMb_dm1 [2] = Ab_Vector [0]*(Mb_Vector [1] - d_my1) - Ab_Vector [1]*(Mb_Vector [0] - d_mx1);

    Check_Inner0 = AbMb_dm0 [0]*AbMb_CrossProduct [0] + AbMb_dm0 [1]*AbMb_CrossProduct [1] + AbMb_dm0 [2]*AbMb_CrossProduct [2];
    Check_Inner1 = AbMb_dm1 [0]*AbMb_CrossProduct [0] + AbMb_dm1 [1]*AbMb_CrossProduct [1] + AbMb_dm1 [2]*AbMb_CrossProduct [2];

    if( Signum( Check_Inner0 ) > 0 ){
        Delta_MbVector [0] = d_mx0;
        Delta_MbVector [1] = d_my0;
        Delta_MbVector [2] = d_mz0;
        //System.out.printf(" 0   Check_Inner0 = %f   \n",Check_Inner0);
    }else if( Signum( Check_Inner1 ) > 0 ){
        Delta_MbVector [0] = d_mx1;
        Delta_MbVector [1] = d_my1;
        Delta_MbVector [2] = d_mz1;
        //System.out.printf(" 1   Check_Inner1 = %f   \n",Check_Inner1);
    }else{
        Delta_MbVector [0] = 0;
        Delta_MbVector [1] = 0;
        Delta_MbVector [2] = 0;
        //System.out.printf(" 2   Check_Inner0 = %f   Check_Inner1 = %f   \n",Check_Inner0,Check_Inner1);
    }
    //System.out.printf("%f   %f  %f  \n",Delta_MbVector [0]*e6,Delta_MbVector [1]*e6,Delta_MbVector [2]*e6);
}

void Initialization::ComputeQuaternion_FromDCM(double Cbn [3][3],double q [4]){
    /*
     * DCMからクォータニオンを計算
     */
    double q1,q2,q3,q4,q1_Abs,q2_Abs,q3_Abs;
    double qAlpha,qBeta,qGamma;
    double trCbn;

    int Alpha,Beta,Gamma;

    /*
    double q_nolm;
    q_nolm = q [0]*q [0] + q [1]*q [1] + q [2]*q [2] + q [3]*q [3];
    */

    Alpha = 5;
    Beta = 5;
    Gamma = 5;
    qAlpha = 0;
    qBeta = 0;
    qGamma = 0;

    trCbn = Cbn [0][0] + Cbn [1][1] + Cbn [2][2];
    
    double B = -1.0/e6/e6;
    
    q1 = Cul_Sqrt( Cbn [0][0]/2.0 + (1.0 - trCbn)/4.0 , B);
    q2 = Cul_Sqrt( Cbn [1][1]/2.0 + (1.0 - trCbn)/4.0 , B);
    q3 = Cul_Sqrt( Cbn [2][2]/2.0 + (1.0 - trCbn)/4.0 , B);
    q4 = Cul_Sqrt((1.0 + trCbn), B)/2.0;

    q1_Abs = abs(q1);
    q2_Abs = abs(q2);
    q3_Abs = abs(q3);


    //System.out.printf("Cbn        %f  %f  %f  %f  \n",Cbn [0][0],Cbn [1][1],Cbn [2][2],q_nolm);
    //System.out.printf("q        %f  %f  %f  %f  \n",q1,q2,q3,q4);
    //System.out.printf("q_Abs        %f  %f  %f  %f  \n",q1_Abs,q2_Abs,q3_Abs,q4);


    if(q1_Abs>=q2_Abs && q1_Abs>=q3_Abs){
        qAlpha = q1_Abs;
        qBeta = q2_Abs;
        qGamma = q3_Abs;

        Alpha = 1;
        Beta = 2;
        Gamma = 3;
    }else if(q2_Abs>=q3_Abs && q2_Abs>=q1_Abs){
        qAlpha = q2_Abs;
        qBeta = q3_Abs;
        qGamma = q1_Abs;

        Alpha = 2;
        Beta = 3;
        Gamma = 1;
    }else if(q3_Abs>=q1_Abs && q3_Abs>=q2_Abs){
        qAlpha = q3_Abs;
        qBeta = q1_Abs;
        qGamma = q2_Abs;

        Alpha = 3;
        Beta = 1;
        Gamma = 2;
    }else{
        //System.out.printf("Error    ComputeQuaternion_FromDCM   \n");
    }

    //System.out.printf("qAlpha     %f  %f  %f  %f  \n",qAlpha,qBeta,qGamma,Cbn [Gamma -1][Beta -1] - Cbn [Beta -1][Gamma -1]);

    if(q4 != 0){
        qAlpha = Signum(Cbn [Gamma -1][Beta -1] - Cbn [Beta -1][Gamma -1])*qAlpha;
        qBeta = Signum( qAlpha*(Cbn [Beta -1][Alpha -1] + Cbn [Alpha -1][Beta -1]) )*qBeta;
        qGamma = Signum( qAlpha*(Cbn [Gamma -1][Alpha -1] + Cbn [Alpha -1][Gamma -1]) )*qGamma;
    }else{
        qBeta = Signum( qAlpha*(Cbn [Beta -1][Alpha -1] + Cbn [Alpha -1][Beta -1]) )*qBeta;
        qGamma = Signum( qAlpha*(Cbn [Gamma -1][Alpha -1] + Cbn [Alpha -1][Gamma -1]) )*qGamma;
    }

    q [Alpha - 1] = qAlpha;
    q [Beta - 1] = qBeta;
    q [Gamma - 1] = qGamma;
    q [3] = q4;
}

void Initialization::ComputeEulerAngle_FromDCM(double Cbn [3][3],double eAngle [3]){

    double Theta0,Theta1,Theta2;
    double CosT2;

    if(abs(Cbn [0][0]) < 1.0/e6 && abs(Cbn [1][0]) < 1.0/e6){
        eAngle [0] = 0.0;
        eAngle [1] = 90.0;
        eAngle [2] = 0.0;
    }else{

        Theta1 = asin(Cbn [2][0]);
        CosT2 = cos(Theta1);

        Theta0 = atan2( - Cbn [2][1]/CosT2,Cbn [2][2]/CosT2);
        Theta2 = atan2( - Cbn [1][0]/CosT2,Cbn [0][0]/CosT2);

        eAngle [0] = Theta0*180.0/Pi;
        eAngle [1] = Theta1*180.0/Pi;
        eAngle [2] = Theta2*180.0/Pi;
    }

}

void Initialization::ComputeGravity(double Gn [3],double Xn [3]){

    double g0,g2,g4;
    double gn,gn1,gn2,gn3;
    double h,L;
    
    double Sin_L,Cos_L;
    
    //重力計算に使われる係数
    g0 = 9.7803267715;
    g2 = 0.0052790414;
    g4 = 0.0000232718;
    gn = 1.63/e6/100.0;
    gn1 = 3.1571/e6/10.0;
    gn2 = 2.1027/e6/e3;
    gn3 = 7.3749/e6/e6/100.0;
    
    
    h = Xn [2];
    L = Xn [0];
    
    Sin_L = sin(L);
    Cos_L = cos(L);
    
    Gn [0] = - gn*h*Sin_L*Cos_L;
    Gn [1] = 0;
    Gn [2] = g0*(1.0 + g2*Sin_L*Sin_L + g4*pow(Sin_L,4.0) )*(1.0 - (gn1 - gn2*Sin_L*Sin_L)*h + gn3*h*h);
}

double Initialization::Cul_Sqrt(double A,double B){
    //数値計算上，値が0になるはずの時，わずかに負の側になってしまう場合がある．それを避けてルートをとる関数
    if(A > B && A < 0.0){
        return 0.0;
    }else if(A >= 0.0){
        return sqrt(A);
    }else{
        //System.out.printf("Error    CulSqrt \n");
        return 0.0;
    }
}

double Initialization::Signum(double A){
    double sgn;
    
    if(A > 0.0){
        sgn = 1.0;
    }else if(A < 0.0){
        sgn = - 1.0;
    }else{
        sgn = 0.0;
    }
    
    return sgn;
}
