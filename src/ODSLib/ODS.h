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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 数学相关函数
//********************************************************************
/// 求平方函数
/// @Author	孙振江
/// @Date	2016.12.25
/// @Input
/// @Param	x
/// @Output
/// @Return	x*x
//********************************************************************
template <typename T>
T Sqr(T x) { return x * x; }

//********************************************************************
/// 点到椭圆距离涉及非线性方程的二分法求解
/// @Author	孙振江
/// @Date	2016.12.25
/// @Input
/// @Param	r0
/// @Param	z0
/// @Param	z1
/// @Param	g
/// @Output
/// @Param
/// @Return	s
//********************************************************************
double GetRoot(double r0, double z0, double z1, double g);

//********************************************************************
/// 点到椭球距离涉及非线性方程的二分法求解
/// @Author	孙振江
/// @Date	2016.12.25
/// @Input
/// @Param	r0, r1
/// @Param	z0, z1, z2
/// @Param	g
/// @Output
/// @Param
/// @Return	s
//********************************************************************
double GetRoot(double r0, double r1, double z0, double z1, double z2, double g);

//********************************************************************
/// 点到椭圆距离求解
/// Distance from a Point to an Ellipse, an Ellipsoid, or a Hyperellipsoid
/// Geometric Tools, LLC, https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
/// @Author	孙振江
/// @Date   2016.12.25
/// @Input
/// @Param  e0          椭圆第一个主轴，e0 >= e1
/// @Param  e1          椭圆第二个主轴
/// @Param  y0          点的第一个坐标，y0 >= 0
/// @Param  y1          点的第二个坐标，y1 >= 0
/// @Output
/// @Param	x0          椭圆上最近点第一坐标
/// @Param	x1          椭圆上最近点第二坐标
/// @Return	distance    最近距离
//********************************************************************
double DistancePointEllipse(double e0, double e1, double y0, double y1, double &x0, double &x1);

//********************************************************************
/// 点到椭球距离求解
/// Distance from a Point to an Ellipse, an Ellipsoid, or a Hyperellipsoid
/// Geometric Tools, LLC, https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
/// @Author	孙振江
/// @Date   2016.12.25
/// @Input
/// @Param  e0          椭球第一个主轴，e0 >= e1 >= e2 >0
/// @Param  e1          椭球第二个主轴
/// @Param  e2          椭球第三个主轴
/// @Param  y0          点的第一个坐标，y0 >= 0
/// @Param  y1          点的第二个坐标，y1 >= 0
/// @Param  y2          点的第三个坐标，y2 >= 0
/// @Output
/// @Param	x0          椭圆上最近点第一坐标
/// @Param	x1          椭圆上最近点第二坐标
/// @Param	x2          椭圆上最近点第三坐标
/// @Return	distance    最近距离
//********************************************************************
double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2, double &x0, double &x1, double &x2);

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
void RK78(const CRightFun<double> &f, double T0, double T1, Eigen::VectorXd &X);

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
}
