/// \file ODS.h
/// This file contains some basic algorithms for Orbit Dynamics and Safety
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

#include <cstdlib>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "ODSLib/ODSConstant.h"
#include "ODSLib/ODSRightFun.h"

/// All the functions are in the namespace ODS
/// 
namespace ODS {

////////////////////////////////////////////////////////////////////////////////
/// Coordinate Transformation functions
/// 坐标变换相关函数
////////////////////////////////////////////////////////////////////////////////

/// Coordinate rotation matrix
/// 坐标旋转矩阵
template <typename T>
bool RotationAxis(const int axis, const T alpha, Eigen::Matrix<T, 3, 3> & Mtx);

/// Compute transformation matrix from VVLH to ICS.
/// 计算从VVLH坐标系到输入直角坐标系（地心惯性系或地固系）的转换矩阵
template <typename T>
bool VVLH2ICSMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx);
// bool VVLHToICSMtx(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, Eigen::Matrix3d &mtx);

/// Compute transformation matrix from ICS to VVLH
/// 计算从直角坐标系到VVLH坐标系的转换矩阵
template <typename T>
bool ICS2VVLHMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx);
// bool ICSToVVLHMtx(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel, Eigen::Matrix3d &mtx);

/// Transformate the coordinates from ICS to VVLH
/// 从直角坐标系转换到VVLH坐标系
template <typename T>
bool ICS2VVLH(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> & Chaser, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & RelState);
// bool ICSToVVLH(const Eigen::VectorXd &Target, const Eigen::VectorXd &Chaser, Eigen::VectorXd &RelState);

/// Transformate Cartesian coordinates to classical elements
/// 计算从直角坐标到轨道根数
template <typename T>
bool Cart2Elem(const Eigen::Matrix<T, Eigen::Dynamic, 1> & cart, 
    Eigen::Array<T, Eigen::Dynamic, 1> & elem);

/// Transformate Cartesian coordinates to classical elements
/// 计算从轨道根数到直角坐标
template <typename T>
bool Elem2Cart(const Eigen::Array<T, Eigen::Dynamic, 1> & elem, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & cart);

////////////////////////////////////////////////////////////////////////////////
/// Math functions
/// 数学相关函数
////////////////////////////////////////////////////////////////////////////////

/// Compute square
/// 求平方函数
template <typename T>
T Sqr(T x) { return x * x; }

/// Get the root of function for distance_point_ellipse by dichotomy
/// 点到椭圆距离涉及非线性方程的二分法求解
template <typename T>
T GetRoot(T r0, T z0, T z1, T g);

/// Get the root of function for distance_point_ellipsoid by dichotomy
/// 点到椭球距离涉及非线性方程的二分法求解
template <typename T>
T GetRoot(T r0, T r1, T z0, T z1, T z2, T g);
// double GetRoot(double r0, double r1, double z0, double z1, double z2, double g);

/// Compute distance from a point to an ellipse
/// 点到椭圆距离求解
template <typename T>
T DistancePointEllipse(T e0, T e1, T y0, T y1, T &x0, T &x1);

/// Compute distance from a point to an ellipsoid
/// 点到椭球距离求解
template <typename T>
T DistancePointEllipsoid(T e0, T e1, T e2, T y0, T y1, T y2, T &x0, T &x1, T &x2);

//********************************************************************
/// 计算点和椭球之间的最近距离
/// @author	孙振江
/// @date	2017.1.2
/// Input
/// @param	pointPos		点的坐标，在椭球主轴坐标系中表示
/// @param	ellipAx1		椭球的主轴长度x
/// @param	ellipAx2		椭球的主轴长度y
/// @param	ellipAx3		椭球的主轴长度z
/// Output
/// @param	distanceVec		最小距离矢量
/// @param	posAtEllip		椭球面上最小距离对应的点，在椭球主轴坐标系中表示
/// @return	true			计算成功
/// 		false			点不在椭球外部
//********************************************************************
bool DistancePointEllipsoid(const Eigen::Vector3d &pointPos, double ellipAx1, double ellipAx2, double ellipAx3, Eigen::Vector3d &distanceVec, Eigen::Vector3d &posAtEllip);

//********************************************************************
/// 计算椭球和椭球之间的最近距离
/// @author	Wang Hua(更新仿真库还不是很熟练，暂时添加一个函数：牛智勇)
/// @date	2012-04-13/2017.1.3
/// Input
/// @param	ellipAx11			第一个椭球的主轴长度x
/// @param	ellipAx12			第一个椭球的主轴长度y
/// @param	ellipAx13			第一个椭球的主轴长度z
/// @param	ellipAx21			第二个椭球的主轴长度x
/// @param	ellipAx22			第二个椭球的主轴长度y
/// @param	ellipAx23			第二个椭球的主轴长度z
/// @param	mtx1To2				从第一个椭球的主轴坐标系转到第二个椭球的主轴坐标系的矩阵
/// @param	ellip2AtEllip1Pos	第二个椭球中心在第一个椭球主轴坐标系中的位置，在第一个椭球主轴坐标系中表示
/// Output
/// @param	distanceVec12		最小距离矢量，从第一个椭球到第二个椭球，在第一个椭球主轴坐标系中表示
/// @param	posAtEllip1			第一个椭球面上最小距离对应的点，在第一个椭球主轴坐标系中表示
/// @param	posAtEllip2			第二个椭球面上最小距离对应的点，在第二个椭球主轴坐标系中表示
/// @return	true				计算成功
///			false				椭球接触或相交，或者没有在规定次数内收敛
//********************************************************************
bool DistanceEllipsoid(double ellipAx11, double ellipAx12, double ellipAx13, double ellipAx21, double ellipAx22, double ellipAx23,
                       const Eigen::MatrixXd &mtx1To2, const Eigen::Vector3d &ellip2AtEllip1Pos, Eigen::Vector3d &distanceVec12, Eigen::Vector3d &posAtEllip1, Eigen::Vector3d &posAtEllip2);

//********************************************************************
/// 对特征值和特征向量按照降序排列
/// @Author	孙振江
/// @Date	2017.1.3
/// @Input
/// @In/Out
/// @Param	eigenvalue
/// @Param	eigenvector
/// @Output
/// @Return
//********************************************************************
void EigenDesSort(Eigen::VectorXd &eigenvalue, Eigen::MatrixXd &eigenvector);

//********************************************************************
/// ODE积分一步，变步长四阶龙格库塔
/// @Author	孙振江
/// @Date	2017.1.11
/// @Input
/// @Param	f       计算右函数类
/// @Param	t       初始时刻
/// @Param	h       积分步长
/// @Param	eps     积分精度
/// @In/Out
/// @Param	X       初始和终端状态
/// @Output
/// @Return
//********************************************************************
template <typename T>
void RK4Step(const CRightFun<T> &f, double t, double h, double eps, Eigen::Matrix<T, Eigen::Dynamic, 1> &y)
{
    int n = y.size();
    int m, i, j, k;
    double hh, p, dt, x, tt, q, a[4];
    Eigen::Matrix<T, Eigen::Dynamic, 1> g, b, c, d, e;
    g.resize(n);
    b.resize(n);
    c.resize(n);
    d.resize(n);
    e.resize(n);
    hh = h;
    m = 1;
    p = 1.0 + eps;
    x = t;
    c = y;
    while (p >= eps)
    {
        a[0] = hh / 2.0;
        a[1] = a[0];
        a[2] = hh;
        a[3] = hh;
        g = y;
        y = c;
        dt = h / m;
        t = x;
        for (j = 0; j <= m - 1; j++)
        {
            f(t, y, d);
            b = y;
            e = y;
            for (k = 0; k <= 2; k++)
            {
                y = e + a[k] * d;
                b = b + a[k + 1] * d / 3.0;
                tt = t + a[k];
                f(tt, y, d);
            }
            y = b + hh * d / 6.0;
            t = t + dt;
        }
        p = 0.0;
        for (i = 0; i <= n - 1; i++)
        {
            q = abs(cons(y[i]) - cons(g[i]));
            if (q > p)
                p = q;
        }
        hh = hh / 2.0;
        m = m + m;
    }
    return;
}

//********************************************************************
/// 勒让德-高斯(Legendre-Gauss)数值积分函数
/// 《C常用算法程序集》
/// @Author	孙振江
/// @Date	2017.1.5
/// @Input
/// @Param	f       被积函数类
/// @Param	a       定积分下限
/// @Param	b       定积分上限
/// @Param	eps     定积分精度
/// @In/Out
/// @Output
/// @Return	积分结果
//********************************************************************
double IntLRGS(const CFun11 &f, double a, double b, double eps = 1.0e-10);

//********************************************************************
/// 三重积分：勒让德-高斯(Legendre-Gauss)方法
/// 《C常用算法程序集》
/// @Author	孙振江
/// @Date	2017.1.5
/// @Input
/// @Param	func       被积函数类
/// @Param	funcy1     y轴下限计算函数
/// @Param	funcy2     y轴上限计算函数
/// @Param	funcz1     z轴下限计算函数
/// @Param	funcz2     z轴上限计算函数
/// @Param	x1         x轴积分下限
/// @Param	x2         x轴积分上限
/// @Param	js         各轴分段数量
/// @In/Out
/// @Output
/// @Return	积分结果
//********************************************************************
double IntLRGS3D(const CFun31 &func, const CFun11 &funcy1, const CFun11 &funcy2, const CFun21 &funcz1, const CFun21 &funcz2,
                 double x1, double x2, int js[3]);


//********************************************************************
/// RK78 ODE积分函数（还有问题，待修改）
/// @Author	孙振江
/// @Date	2017.1.4
/// @Input
/// @Param	f       计算右函数类
/// @Param	T0      初始时刻
/// @Param	T1      终端时刻
/// @In/Out
/// @Param	X       初始和终端状态
/// @Output
/// @Return
//********************************************************************
// void RK78(const CRightFun<double> &f, double T0, double T1, Eigen::VectorXd &X);
// void ODS::RK78(const CRightFun<double>& f, double T0, double T1, VectorXd& X)
// {
//     VectorXd Y0 = X;
//     int N = Y0.size();
// 
//     double ERREST;
//     double H0 = 0.001; double HS = 0.1; double H1 = 100.0;
//     double EPS = 1.e-12; double BS = 20 * EPS;
// 
//     int I, J, K;
// 
//     // AlgebraicMatrix<double> Z(N,16);
// 
//     std::vector<VectorXd> Z;
//     for (I = 0; I < N; I++)
//     {
//         VectorXd tmp(16);
//         tmp.fill(0);
//         Z.push_back(tmp);
//     }
// 
//     VectorXd Y1(N);
//     VectorXd Y1cons(N);
// 
//     double VIHMAX = 0.0, T, H;
//     double RFNORM, HH0, HH1;
// 
//     double HSQR = 1.0 / 9.0;
//     double A[13], C[13], D[13];
//     double B[13][12];
// 
// 
// 
//     A[0] = 0.0; A[1] = 1.0 / 18.0; A[2] = 1.0 / 12.0; A[3] = 1.0 / 8.0; A[4] = 5.0 / 16.0; A[5] = 3.0 / 8.0;
//     A[6] = 59.0 / 400.0; A[7] = 93.0 / 200.0; A[8] = 5490023248.0 / 9719169821.0; A[9] = 13.0 / 20.0; A[10] = 1201146811.0 / 1299019798.0; A[11] = 1.0;
//     A[12] = 1.0;
// 
//     B[0][0] = 0.0; B[0][1] = 0.0; B[0][2] = 0.0; B[0][3] = 0.0; B[0][4] = 0.0;
//     B[0][5] = 0.0; B[0][6] = 0.0; B[0][7] = 0.0; B[0][8] = 0.0; B[0][9] = 0.0;
//     B[0][10] = 0.0; B[0][11] = 0.0;
// 
//     B[1][0] = 1.0 / 18.0; B[1][1] = 0.0; B[1][2] = 0.0; B[1][3] = 0.0; B[1][4] = 0.0;
//     B[1][5] = 0.0; B[1][6] = 0.0; B[1][7] = 0.0; B[1][8] = 0.0; B[1][9] = 0.0;
//     B[1][10] = 0.0; B[1][11] = 0.0;
// 
//     B[2][0] = 1.0 / 48.0; B[2][1] = 1.0 / 16.0; B[2][2] = 0.0; B[2][3] = 0.0; B[2][4] = 0.0;
//     B[2][5] = 0.0; B[2][6] = 0.0; B[2][7] = 0.0; B[2][8] = 0.0; B[2][9] = 0.0;
//     B[2][10] = 0.0; B[2][11] = 0.0;
// 
//     B[3][0] = 1.0 / 32.0; B[3][1] = 0.0; B[3][2] = 3.0 / 32.0; B[3][3] = 0.0; B[3][4] = 0.0;
//     B[3][5] = 0.0; B[3][6] = 0.0; B[3][7] = 0.0; B[3][8] = 0.0; B[3][9] = 0.0;
//     B[3][10] = 0.0; B[3][11] = 0.0;
// 
//     B[4][0] = 5.0 / 16.0; B[4][1] = 0.0; B[4][2] = -75.0 / 64.0; B[4][3] = 75.0 / 64.0; B[4][4] = 0.0;
//     B[4][5] = 0.0; B[4][6] = 0.0; B[4][7] = 0.0; B[4][8] = 0.0; B[4][9] = 0.0;
//     B[4][10] = 0.0; B[4][11] = 0.0;
// 
//     B[5][0] = 3.0 / 80.0; B[5][1] = 0.0; B[5][2] = 0.0; B[5][3] = 3.0 / 16.0; B[5][4] = 3.0 / 20.0;
//     B[5][5] = 0.0; B[5][6] = 0.0; B[5][7] = 0.0; B[5][8] = 0.0; B[5][9] = 0.0;
//     B[5][10] = 0.0; B[5][11] = 0.0;
// 
//     B[6][0] = 29443841.0 / 614563906.0; B[6][1] = 0.0; B[6][2] = 0.0; B[6][3] = 77736538.0 / 692538347.0; B[6][4] = -28693883.0 / 1125000000.0;
//     B[6][5] = 23124283.0 / 1800000000.0; B[6][6] = 0.0; B[6][7] = 0.0; B[6][8] = 0.0; B[6][9] = 0.0;
//     B[6][10] = 0.0; B[6][11] = 0.0;
// 
//     B[7][0] = 16016141.0 / 946692911.0; B[7][1] = 0.0; B[7][2] = 0.0; B[7][3] = 61564180.0 / 158732637.0; B[7][4] = 22789713.0 / 633445777.0;
//     B[7][5] = 545815736.0 / 2771057229.0; B[7][6] = -180193667.0 / 1043307555.0; B[7][7] = 0.0; B[7][8] = 0.0; B[7][9] = 0.0;
//     B[7][10] = 0.0; B[7][11] = 0.0;
// 
//     B[8][0] = 39632708.0 / 573591083.0; B[8][1] = 0.0; B[8][2] = 0.0; B[8][3] = -433636366.0 / 683701615.0; B[8][4] = -421739975.0 / 2616292301.0;
//     B[8][5] = 100302831.0 / 723423059.0; B[8][6] = 790204164.0 / 839813087.0; B[8][7] = 800635310.0 / 3783071287.0; B[8][8] = 0.0; B[8][9] = 0.0;
//     B[8][10] = 0.0; B[8][11] = 0.0;
// 
//     B[9][0] = 246121993.0 / 1340847787.0; B[9][1] = 0.0; B[9][2] = 0.0; B[9][3] = -37695042795.0 / 15268766246.0; B[9][4] = -309121744.0 / 1061227803.0;
//     B[9][5] = -12992083.0 / 490766935.0; B[9][6] = 6005943493.0 / 2108947869.0; B[9][7] = 393006217.0 / 1396673457.0; B[9][8] = 123872331.0 / 1001029789.0; B[9][9] = 0.0;
//     B[9][10] = 0.0; B[9][11] = 0.0;
// 
//     B[10][0] = -1028468189.0 / 846180014.0; B[10][1] = 0.0; B[10][2] = 0.0; B[10][3] = 8478235783.0 / 508512852.0; B[10][4] = 1311729495.0 / 1432422823.0;
//     B[10][5] = -10304129995.0 / 1701304382.0; B[10][6] = -48777925059.0 / 3047939560.0; B[10][7] = 15336726248.0 / 1032824649.0; B[10][8] = -45442868181.0 / 3398467696.0; B[10][9] = 3065993473.0 / 597172653.0;
//     B[10][10] = 0.0; B[10][11] = 0.0;
// 
//     B[11][0] = 185892177.0 / 718116043.0; B[11][1] = 0.0; B[11][2] = 0.0; B[11][3] = -3185094517.0 / 667107341.0; B[11][4] = -477755414.0 / 1098053517.0;
//     B[11][5] = -703635378.0 / 230739211.0; B[11][6] = 5731566787.0 / 1027545527.0; B[11][7] = 5232866602.0 / 850066563.0; B[11][8] = -4093664535.0 / 808688257.0; B[11][9] = 3962137247.0 / 1805957418.0;
//     B[11][10] = 65686358.0 / 487910083.0; B[11][11] = 0.0;
// 
//     B[12][0] = 403863854.0 / 491063109.0; B[12][1] = 0.0; B[12][2] = 0.0; B[12][3] = -5068492393.0 / 434740067.0; B[12][4] = -411421997.0 / 543043805.0;
//     B[12][5] = 652783627.0 / 914296604.0; B[12][6] = 11173962825.0 / 925320556.0; B[12][7] = -13158990841.0 / 6184727034.0; B[12][8] = 3936647629.0 / 1978049680.0; B[12][9] = -160528059.0 / 685178525.0;
//     B[12][10] = 248638103.0 / 1413531060.0; B[12][11] = 0.0;
// 
//     C[0] = 14005451.0 / 335480064.0; C[1] = 0.0; C[2] = 0.0; C[3] = 0.0; C[4] = 0.0; C[5] = -59238493.0 / 1068277825.0;
//     C[6] = 181606767.0 / 758867731.0; C[7] = 561292985.0 / 797845732.0; C[8] = -1041891430.0 / 1371343529.0; C[9] = 760417239.0 / 1151165299.0; C[10] = 118820643.0 / 751138087.0; C[11] = -528747749.0 / 2220607170.0;
//     C[12] = 1.0 / 4.0;
// 
//     D[0] = 13451932.0 / 455176623.0; D[1] = 0.0; D[2] = 0.0; D[3] = 0.0; D[4] = 0.0; D[5] = -808719846.0 / 976000145.0;
//     D[6] = 1757004468.0 / 5645159321.0; D[7] = 656045339.0 / 265891186.0; D[8] = -3867574721.0 / 1518517206.0; D[9] = 465885868.0 / 322736535.0; D[10] = 53011238.0 / 667516719.0; D[11] = 2.0 / 45.0;
//     D[12] = 0.0;
// 
//     for (I = 0; I < N; I++)
//     {
//         Z[I][0] = Y0[I];
//         Z[I][1] = 0.0;
//     }
// 
//     H = abs(HS); HH0 = abs(H0); HH1 = abs(H1);
//     T = T0; RFNORM = 0.0; ERREST = 0.0;
// 
//     while (T != T1) {
// 
//         // compute new stepsize
//         if (RFNORM != 0) { H = H*min(4.0, exp(HSQR*log(EPS / RFNORM))); }
//         if (abs(H)>abs(HH1)) { H = HH1; }
//         else if (abs(H)<abs(HH0)*0.99) {
//             H = HH0;
//             cout << "--- WARNING, MINIMUM STEPSIZE REACHED IN RK" << endl;
//         }
// 
//         if ((T + H - T1)*H>0) { H = T1 - T; }
// 
//         for (J = 0; J<13; J++) {
// 
//             for (I = 0; I<N; I++) {
// 
//                 Y0[I] = 0.0; // EVALUATE RHS AT 13 POINTS
// 
//                 for (K = 0; K<J; K++) { Y0[I] = Y0[I] + Z[I][K + 3] * B[J][K]; }
// 
//                 Y0[I] = H*Y0[I] + Z[I][0];
//             }
// 
//             //Y1 = f(Y0, T + H*A[J]);
//             f(T + H*A[J], Y0, Y1);
// 
//             for (I = 0; I<N; I++) { Z[I][J + 3] = Y1[I]; }
//         }
// 
//         for (I = 0; I<N; I++) {
// 
//             Z[I][1] = 0.0; Z[I][2] = 0.0; // EXECUTE 7TH,8TH ORDER STEPS
// 
//             for (J = 0; J<13; J++) {
//                 Z[I][1] = Z[I][1] + Z[I][J + 3] * D[J];
//                 Z[I][2] = Z[I][2] + Z[I][J + 3] * C[J];
//             }
// 
//             Y1[I] = (Z[I][2] - Z[I][1])*H;
//             Z[I][2] = Z[I][2] * H + Z[I][0];
//         }
// 
// 
//         for (I = 0; I<N; I++) { Y1cons[I] = Y1[I]; }
// 
//         //RFNORM = normtmp(N, Y1cons); // ESTIMATE ERROR AND DECIDE ABOUT BACKSTEP
//         for (I = 0; I<N; I++) {
//             if (Y1cons[I]<0) { Y1cons[I] = -Y1cons[I]; }
//             RFNORM = max(RFNORM, X[I]);
//         }
// 
//         if ((RFNORM>BS) && (abs(H / H0)>1.2)) {
//             H = H / 3.0;
//             RFNORM = 0;
//         }
//         else {
//             for (I = 0; I<N; I++) { Z[I][0] = Z[I][2]; }
//             T = T + H;
//             VIHMAX = max(VIHMAX, H);
//             ERREST = ERREST + RFNORM;
//         }
//     }
// 
//     for (I = 0; I<N; I++) { Y1[I] = Z[I][0]; }
// 
//     X = Y1;
// 
// }

}
