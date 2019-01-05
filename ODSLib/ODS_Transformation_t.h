/// \file ODS_t.h
/// This file contains the implementations for the methods in "ODS.h"
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

#include <iostream>

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

/// \brief Transformate the coordinates from ICS to VVLH
/// 从直角坐标系转换到VVLH坐标系
/// 
/// The definition of VVLH coordinate system can be found in "VVLHToICSMtx"
/// 
/// \Param[in]  Target      the coordinates of target in ICS (m, m/s)
/// \Param[in]  Chaser      the coordinates of chaser in ICS (m, m/s)
/// \Param[out] RelState    the coordinates of chaser in target's VVLH (m, m/s)
/// \Return	                true = normal; false = error
template <typename T>
bool ICS2VVLH(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> & Chaser, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & RelState) {

    Eigen::Matrix<T, 3, 3> Mo, Mw;                                // 坐标变换矩阵
    Eigen::Matrix<T, 3, 1> RTar, VTar, RCha, VCha, RRel, VRel, w;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(6, 6);
    Eigen::Matrix<T, Eigen::Dynamic, 1> state(6);

    RTar = Target.head(3);
    VTar = Target.tail(3);
    RCha = Chaser.head(3);
    VCha = Chaser.tail(3);

    ICS2VVLHMtx(RTar, VTar, Mo);
    w = -Mo.row(1)*VTar.norm() / RTar.norm();
    Mw << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;

    M.topLeftCorner(3, 3) = Mo;
    M.topRightCorner(3, 3).fill(0);
    M.bottomLeftCorner(3, 3) = -Mo*Mw;
    M.bottomRightCorner(3, 3) = Mo;
    state = Chaser - Target;
    RelState = M*state;

    return true;
}

/// \brief Transformate the coordinates from VVLH to ICS
/// 从VVLH坐标系转换到直角坐标系
/// 
/// The definition of VVLH coordinate system can be found in "VVLHToICSMtx"
/// 
/// \Param[in]  Target      the coordinates of target in ICS (m, m/s)
/// \Param[in]  RelState    the coordinates of chaser in target's VVLH (m, m/s)
/// \Param[out] Chaser      the coordinates of chaser in ICS (m, m/s)
/// \Return	                true = normal; false = error
template <typename T>
bool VVLH2ICS(const Eigen::Matrix<T, Eigen::Dynamic, 1> &Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &RelState, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> &Chaser) {

    Eigen::Matrix<T, 3, 3> Mo, Mw;                                // 坐标变换矩阵
    Eigen::Matrix<T, 3, 1> RTar, VTar, RCha, VCha, RRel, VRel, w;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(6, 6), M_inv(6, 6);
    Eigen::Matrix<T, Eigen::Dynamic, 1> state(6);

    RTar = Target.head(3);
    VTar = Target.tail(3);
    RCha = Chaser.head(3);
    VCha = Chaser.tail(3);

    ICS2VVLHMtx(RTar, VTar, Mo);
    w = -Mo.row(1)*VTar.norm() / RTar.norm();
    Mw << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;

    M.topLeftCorner(3, 3) = Mo;
    M.topRightCorner(3, 3).fill(0);
    M.bottomLeftCorner(3, 3) = -Mo*Mw;
    M.bottomRightCorner(3, 3) = Mo;
    M_inv = M.inverse();

    Chaser = M_inv*RelState + Target;

    return true;
}

/// \brief Return transformation matrix from LVLH to ICS.
/// 计算从LVLH坐标系到输入直角坐标系（地心惯性系或地固系）的转换矩阵
/// 
/// LVLH(Local Vertical Local Horizontal): 
/// the X axis is along the radiation direction toward velocity 
/// the Y axis is co-velocity direction
/// the Z axis is along the orbital normal direction
/// 
/// \Param[in]  pos	    the position of vehicle in ICS (m)
/// \Param[in]  vel	    the velocity of vehicle in ICS (m/s)
/// \Param[out] mtx	    the transformation matrix from VVLH to ICS
/// \Return	            true = normal; false = error
template <typename T>
bool LVLH2ICSMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx) {

    T r, s, n;
    Eigen::Matrix<T, 3, 1> ss, nn;

    // Radiation direction, x
    r = pos.norm();
    // Normal direction, z
    nn = pos.cross(vel);
    n = nn.norm();
    // Co-velocity direction, y
    ss = nn.cross(pos);
    s = ss.norm();

    Eigen::Matrix<T, 3, 3> mtxTmp;
    mtxTmp.row(0) = pos / r;
    mtxTmp.row(1) = ss / s;
    mtxTmp.row(2) = nn / n;

    mtx = mtxTmp.inverse();

    return true;
}

/// \brief Return transformation matrix from ICS to LVLH
/// 计算从直角坐标系到LVLH
/// 
/// The definition of LVLH coordinate system can be found in "LVLH2ICSMtx"
/// 
/// \Param[in]  pos     the position of vehicle 地心惯性系或地固系中的飞行器位置[m]
/// \Param[in]  vel     velocity of vehicle 地心惯性系或地固系中的飞行器速度[m/s]
/// \Param[out] mtx     ICS到LVLH的转移矩阵
/// \Return	            true = normal; false = error
template <typename T>
bool ICS2LVLHMtx(const Eigen::Matrix<T, 3, 1> & pos, 
    const Eigen::Matrix<T, 3, 1> & vel, Eigen::Matrix<T, 3, 3> & mtx) {

    T r, s, n;
    Eigen::Matrix<T, 3, 1> ss, nn;

    // Radiation direction, x
    r = pos.norm();
    // Normal direction, z
    nn = pos.cross(vel);
    n = nn.norm();
    // Co-velocity direction, y
    ss = nn.cross(pos);
    s = ss.norm();

    mtx.row(0) = pos / r;
    mtx.row(1) = ss / s;
    mtx.row(2) = nn / n;

    return true;
}

/// \brief Transformate the coordinates from ICS to LVLH
/// 从直角坐标系转换到LVLH坐标系
/// 
/// The definition of LVLH coordinate system can be found in "LVLHToICSMtx"
/// 
/// \Param[in]  Target      the coordinates of target in ICS (m, m/s)
/// \Param[in]  Chaser      the coordinates of chaser in ICS (m, m/s)
/// \Param[out] RelState    the coordinates of chaser in target's LVLH (m, m/s)
/// \Return	                true = normal; false = error
template <typename T>
bool ICS2LVLH(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> & Chaser, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & RelState) {

    Eigen::Matrix<T, 3, 3> Mo, Mw;                                // 坐标变换矩阵
    Eigen::Matrix<T, 3, 1> RTar, VTar, RCha, VCha, RRel, VRel, w;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(6, 6);
    Eigen::Matrix<T, Eigen::Dynamic, 1> state(6);

    RTar = Target.head(3);
    VTar = Target.tail(3);
    RCha = Chaser.head(3);
    VCha = Chaser.tail(3);

    ICS2LVLHMtx(RTar, VTar, Mo);
    w = Mo.row(2) * VTar.norm() / RTar.norm();
    Mw << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;

    M.topLeftCorner(3, 3) = Mo;
    M.topRightCorner(3, 3).fill(0);
    M.bottomLeftCorner(3, 3) = -Mo*Mw;
    M.bottomRightCorner(3, 3) = Mo;
    state = Chaser - Target;
    RelState = M*state;

    return true;
}

/// \brief Transformate the coordinates from LVLH to ICS
/// 从LVLH坐标系转换到直角坐标系
/// 
/// The definition of LVLH coordinate system can be found in "LVLHToICSMtx"
/// 
/// \Param[in]  Target      the coordinates of target in ICS (m, m/s)
/// \Param[in]  RelState    the coordinates of chaser in target's LVLH (m, m/s)
/// \Param[out] Chaser      the coordinates of chaser in ICS (m, m/s)
/// \Return	                true = normal; false = error
template <typename T>
bool LVLH2ICS(const Eigen::Matrix<T, Eigen::Dynamic, 1> &Target, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &RelState, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> &Chaser) {

    Eigen::Matrix<T, 3, 3> Mo, Mw;                                // 坐标变换矩阵
    Eigen::Matrix<T, 3, 1> RTar, VTar, RCha, VCha, RRel, VRel, w;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(6, 6), M_inv(6, 6);
    Eigen::Matrix<T, Eigen::Dynamic, 1> state(6);

    RTar = Target.head(3);
    VTar = Target.tail(3);
    RCha = Chaser.head(3);
    VCha = Chaser.tail(3);

    ICS2LVLHMtx(RTar, VTar, Mo);
    w = Mo.row(2) * VTar.norm() / RTar.norm();
    Mw << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;

    M.topLeftCorner(3, 3) = Mo;
    M.topRightCorner(3, 3).fill(0);
    M.bottomLeftCorner(3, 3) = -Mo*Mw;
    M.bottomRightCorner(3, 3) = Mo;
    M_inv = M.inverse();

    Chaser = M_inv*RelState + Target;

    return true;
}

/// Transformate Cartesian coordinates to classical elements
/// 计算从直角坐标到轨道根数
/// 
/// \Param[in]  cart    Cartesian coordinates (m, m/s)
/// \Param[out] elem    classical elements (m, rad)
/// \Return	            true = normal; false = error
template <typename T>
bool Cart2Elem(const Eigen::Matrix<T, Eigen::Dynamic, 1> & cart, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & elem) {

    double GM_km = GM * 1e-9;

    Eigen::Matrix<T, 3, 1> RVec, VVec; // 距离和速度矢量
    RVec = cart.head(3) / 1000; // km
    VVec = cart.tail(3) / 1000; // km/s
    T rr = RVec.norm();         // 地心距, km
    T vv = VVec.norm();

    // unit vector
    Eigen::Matrix<T, 3, 1> XVec, YVec, ZVec;
    XVec << 1, 0, 0;
    YVec << 0, 1, 0;
    ZVec << 0, 0, 1;

    // calculate angular-momentum and eccentricity
    Eigen::Matrix<T, 3, 1> HVec, NVec, EVec;
    HVec = RVec.cross(VVec);
    T hh = HVec.norm();
    NVec = ZVec.cross(HVec);
    T nn = NVec.norm();
    EVec = ((vv*vv - GM_km/rr)*RVec - (RVec.dot(VVec))*VVec)/GM_km;
    T Ecc = EVec.norm();

    // calculate semi-axis
    T SemiA, p, xi;
    xi = vv*vv/2.0 - GM_km/rr;
    if(abs(cons(Ecc)-1.0) < 1e-6) {
        SemiA = 0;
        p = hh*hh/GM_km;
        std::cout << "Parabolic orbit!!" << std::endl;
    } else {
        SemiA = -GM_km/(2.0*xi);
        p = SemiA * (1 - Ecc*Ecc);
    }

    // calculate inclination
    T Inc;
    Inc = acos(HVec.dot(ZVec)/hh);

    // calculate RAAN
    T RAAN;
    RAAN = acos(NVec.dot(XVec)/nn);
    if(cons(NVec.dot(YVec)) < 0.0) {
        RAAN = 2*ODS::Pi - RAAN;
    }

    // calculate omega
    T w;
    w = acos(NVec.dot(EVec)/(nn*Ecc));
    if(cons(EVec.dot(ZVec)) < 0.0) {
        w = 2*ODS::Pi - w;
    }

    // calculate omega
    T TrueA;
    TrueA = acos(EVec.dot(RVec)/(Ecc*rr));
    if(cons(RVec.dot(VVec)) < 0.0) {
        TrueA = 2*ODS::Pi - TrueA;
    }

    elem(0) = SemiA * 1000.0;
    elem(1) = Ecc;
    elem(2) = Inc;
    elem(3) = RAAN;
    elem(4) = w;
    elem(5) = TrueA;

    return true;
}

/// Transformate Cartesian coordinates to classical elements
/// 计算从轨道根数到直角坐标
/// 
/// \Param[in]  elem    classical elements (m, rad)
/// \Param[out] cart    Cartesian coordinates (m, m/s)
/// \Return	            true = normal; false = error
template <typename T>
bool Elem2Cart(const Eigen::Matrix<T, Eigen::Dynamic, 1> & elem, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & cart) {

    T SemiA = elem(0);
    T Ecc = elem(1);
    T Inc = elem(2);
    T RAAN = elem(3);
    T w = elem(4);
    T TrueA = elem(5);

    T p = SemiA * (1 - Ecc * Ecc); // 半通径
    T r = p / (1 + Ecc * cos(TrueA));

    // position and velocity vector in orbital frame
    Eigen::Matrix<T, 3, 1> r_orb, v_orb;
    r_orb << r * cos(TrueA), r * sin(TrueA), 0.0;
    v_orb << -sin(TrueA), Ecc + cos(TrueA), 0.0;
    v_orb = sqrt(GM/p) * v_orb;

    // orbital frame to inertial frame
    Eigen::Matrix<T, 3, 3> alpha;
    alpha(0, 0) =  cos(RAAN)*cos(w) - sin(RAAN)*sin(w)*cos(Inc);
    alpha(0, 1) = -cos(RAAN)*sin(w) - sin(RAAN)*cos(w)*cos(Inc);
    alpha(0, 2) =  sin(RAAN)*sin(Inc);
    alpha(1, 0) =  sin(RAAN)*cos(w) + cos(RAAN)*sin(w)*cos(Inc);
    alpha(1, 1) = -sin(RAAN)*sin(w) + cos(RAAN)*cos(w)*cos(Inc);
    alpha(1, 2) = -cos(RAAN)*sin(Inc);
    alpha(2, 0) = sin(w)*sin(Inc);
    alpha(2, 1) = cos(w)*sin(Inc);
    alpha(2, 2) = cos(Inc);

    cart.head(3) = alpha * r_orb;
    cart.tail(3) = alpha * v_orb;

    return true;
}

}
