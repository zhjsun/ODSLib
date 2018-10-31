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

////////////////////////////////////////////////////////////////////////////////
/// Coordinate Transformation functions
/// 坐标变换相关函数
////////////////////////////////////////////////////////////////////////////////

/// Coordinate rotation matrix
/// 坐标旋转矩阵
template <typename T>
bool RotationAxis(const int axis, const T alpha, Eigen::Matrix<T, 3, 3> & Mtx);


////////////////////////////////////////////////////////////////////////////////

/// Compute transformation matrix from VVLH to ICS.
/// 计算从VVLH坐标系到输入直角坐标系（地心惯性系或地固系）的转换矩阵
template <typename T>
bool VVLH2ICSMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx);

/// Compute transformation matrix from ICS to VVLH
/// 计算从直角坐标系到VVLH坐标系的转换矩阵
template <typename T>
bool ICS2VVLHMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx);

/// Transformate the coordinates from ICS to VVLH
/// 从直角坐标系转换到VVLH坐标系
template <typename T>
bool ICS2VVLH(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> & Chaser, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & RelState);

/// Transformate the coordinates from VVLH to ICS
/// 从VVLH坐标系转换到直角坐标系
template <typename T>
bool VVLH2ICS(const Eigen::Matrix<T, Eigen::Dynamic, 1> &Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &RelState, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> &Chaser );

/// Compute transformation matrix from LVLH to ICS.
/// 计算从LVLH坐标系到输入直角坐标系（地心惯性系或地固系）的转换矩阵
template <typename T>
bool LVLH2ICSMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx);

/// Compute transformation matrix from ICS to LVLH
/// 计算从直角坐标系到VVLH坐标系的转换矩阵
template <typename T>
bool ICS2LVLHMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx);

/// Transformate the coordinates from ICS to LVLH
/// 从直角坐标系转换到LVLH坐标系
template <typename T>
bool ICS2LVLH(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> & Chaser, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & RelState);

/// Transformate the coordinates from LVLH to ICS
/// 从LVLH坐标系转换到直角坐标系
template <typename T>
bool LVLH2ICS(const Eigen::Matrix<T, Eigen::Dynamic, 1> &Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &RelState, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> &Chaser );


////////////////////////////////////////////////////////////////////////////////

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

}

#include "ODSLib/ODS_Transformation_t.h"