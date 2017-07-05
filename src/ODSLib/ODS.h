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

/// Compute closet distance from a point to an ellipsoid
/// 计算点和椭球之间的最近距离
template <typename T>
bool DistancePointEllipsoid(const Eigen::Matrix<T, 3, 1> &pointPos, 
                            T ellipAx1, T ellipAx2, T ellipAx3, 
                            Eigen::Matrix<T, 3, 1> &distanceVec, 
                            Eigen::Matrix<T, 3, 1> &posAtEllip);

/// Compute closet distance between two ellipsoids
/// 计算椭球和椭球之间的最近距离
template <typename T>
bool DistanceEllipsoid(T ellipAx11, T ellipAx12, T ellipAx13, 
                       T ellipAx21, T ellipAx22, T ellipAx23,
                       const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx1To2, 
                       const Eigen::Matrix<T, 3, 1> &ellip2AtEllip1Pos, 
                       Eigen::Matrix<T, 3, 1> &distanceVec12, 
                       Eigen::Matrix<T, 3, 1> &posAtEllip1, 
                       Eigen::Matrix<T, 3, 1> &posAtEllip2);

/// Sort eigen values and vectors in descending order
/// 对特征值和特征向量按照降序排列
template <typename T>
void EigenDesSort(Eigen::Matrix<T, Eigen::Dynamic, 1> &eigenvalue, 
                  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &eigenvector);

/// Integrate one step for ODE by 4-order adaptive Runge-Kutta
/// ODE积分一步，变步长四阶龙格库塔
template <typename T>
void RK4Step(const CRightFun<T> &f, double t, double h, double eps, 
             Eigen::Matrix<T, Eigen::Dynamic, 1> &y);

/// RK78 ODE积分函数（还有问题，待修改）
// void RK78(const CRightFun<double> &f, double T0, double T1, Eigen::VectorXd &X);

/// Numerical integration of Legendre-Gauss algorithm
/// 勒让德-高斯(Legendre-Gauss)数值积分函数
template <typename T>
T IntLRGS(const CFun11 &f, T a, T b, double eps = 1.0e-10);

/// Triple integration of Legendre-Gauss algorithm
/// 三重积分：勒让德-高斯(Legendre-Gauss)方法
template <typename T>
T IntLRGS3D(const CFun31 &func, const CFun11 &funcy1, const CFun11 &funcy2, 
            const CFun21 &funcz1, const CFun21 &funcz2, T x1, T x2, int js[3]);

}
