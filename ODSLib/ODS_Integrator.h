/// \file ODS.h
/// This file contains some basic algorithms for Orbit Dynamics and Safety
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

#include <cstdlib>
#include <cmath>
#include <vector>

#include <eigen3/Eigen/Dense>

#include "ODSLib/ODS_Constant.h"
#include "ODSLib/ODS_RightFun.h"

/// All the functions are in the namespace ODS
/// 
namespace ODS {

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

#include "ODSLib/ODS_Integrator_t.h"