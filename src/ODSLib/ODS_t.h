/// \file ODS_t.h
/// This file contains the implementations for the methods in "ODS.h"
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

#include "ODSLib/ODS.h"

namespace ODS {

/// Coordinate rotation matrix
/// 坐标旋转矩阵
/// 
/// \Param[in]	axis    the rotation axis [1, 2, 3]
/// \Param[in]  alpha   the rotation angle (rad)
/// \Param[out] Mtx     the rotation matrix
/// \Return             true = normal; false = error
template <typename T>
bool RotationAxis(const int axis, const T alpha, Eigen::Matrix<T, 3, 3> & Mtx)
{
    assert(axis >= 1 && axis <= 3);
    switch (axis)
    {
    case 1:
        Mtx << 1, 0, 0, 0, cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha);
        break;
    case 2:
        Mtx << cos(alpha), 0, -sin(alpha), 0, 1, 0, sin(alpha), 0, cos(alpha);
        break;
    default:
        Mtx << cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha), 0, 0, 0, 1;
        break;
    }

    return true;
}

/// \brief Return transformation matrix from VVLH to ICS.
/// 计算从VVLH坐标系到输入直角坐标系（地心惯性系或地固系）的转换矩阵
/// 
/// VVLH(Vehicle Velocity Local Horizontal): 
/// the Z axis is along the negative position vector
/// the Y axis is along the negative orbit normal
/// the X axis is toward velocity 
/// 
/// \Param[in]  pos	    the position of vehicle in ICS (m)
/// \Param[in]  vel	    the velocity of vehicle in ICS (m/s)
/// \Param[out] mtx	    the transformation matrix from VVLH to ICS
/// \Return	            true = normal; false = error
//********************************************************************
template <typename T>
bool VVLH2ICSMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx) {

    T r, s, n;
    Eigen::Matrix<T, 3, 1> ss, nn;
    r = pos.norm();
    nn = pos.cross(vel);
    n = nn.norm();
    ss = pos.cross(nn);
    s = ss.norm();

    Eigen::Matrix<T, 3, 3> mtxTmp;
    mtxTmp.row(0) = -ss / s;
    mtxTmp.row(1) = -nn / n;
    mtxTmp.row(2) = -pos / r;

    mtx = mtxTmp.inverse();

    return true;
}

/// \brief Return transformation matrix from ICS to VVLH
/// 计算从直角坐标系到VVLH坐标系的转换矩阵
/// 
/// The definition of VVLH coordinate system can be found in "VVLHToICSMtx"
/// 
/// \Param[in]  pos     the position of vehicle 地心惯性系或地固系中的飞行器位置[m]
/// \Param[in]  vel     velocity of vehicle 地心惯性系或地固系中的飞行器速度[m/s]
/// \Param[out] mtx     ICS到VVLH的转移矩阵
/// \Return	            true = normal; false = error
template <typename T>
bool ICS2VVLHMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx) {

    T r, s, n;
    Eigen::Matrix<T, 3, 1> ss, nn;
    r = pos.norm();
    nn = pos.cross(vel);
    n = nn.norm();
    ss = pos.cross(nn);
    s = ss.norm();

    mtx.row(0) = -ss / s;
    mtx.row(1) = -nn / n;
    mtx.row(2) = -pos / r;

    return true;
}

}