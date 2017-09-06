/// \file ODS.h
/// This file contains some basic algorithms for Orbit Dynamics and Safety
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

#include <cstdlib>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "ODSLib/ODS_Constant.h"
#include "ODSLib/ODS_RightFun.h"

/// All the functions are in the namespace ODS
/// 
namespace ODS {

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

}

#include "ODSLib/ODS_Math_t.h"