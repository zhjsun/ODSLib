/// \file ODS_t.h
/// This file contains the implementations for the methods in "ODS.h"
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

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
    Eigen::Array<T, Eigen::Dynamic, 1> & elem) {

    double GM_km = GM * 1e-9;

    Eigen::Matrix<T, 3, 1> RVec, VVec; // 距离和速度矢量
    RVec = cart.head(3) / 1000; // km
    VVec = cart.tail(3) / 1000; // km/s
    T rr = RVec.norm();         // 地心距, km
    // 计算动量矩
    Eigen::Matrix<T, 3, 1> HVec;
    HVec = RVec.cross(VVec);
    T hh = HVec.norm();
    // 计算轨道倾角
    T Inc = acos(HVec(2) / hh);
    // 计算升交点赤经
    T RAAN = atan(-HVec(0) / HVec(1));
    if (cons(HVec(1)) > 0)
    {
        RAAN += Pi;
    }
    if (cons(RAAN) < 0)
    {
        RAAN += 2 * Pi;
    }
    // 计算偏心率矢量和大小
    Eigen::Matrix<T, 3, 1> EVec;
    EVec = VVec.cross(HVec) / GM_km - RVec / rr;
    T Ecc = EVec.norm();
    // 计算近地点纬度幅角
    T w = atan(EVec(2) / ((EVec(1) * sin(RAAN) + EVec(0) * cos(RAAN)) * sin(Inc)));
    if (cons(EVec(2)) > 0 && cons(w) < 0)
    {
        w += Pi;
    }
    else if (cons(EVec(2)) < 0 && cons(w) > 0)
    {
        w -= Pi;
    }
    // 计算半长轴
    T SemiA = hh * hh / (GM_km * (1 - Ecc * Ecc)) * 1000; // m
    // 计算真近点角
    T u = atan(RVec(2) / ((RVec(1) * sin(RAAN) + RVec(0) * cos(RAAN)) * sin(Inc)));
    if (cons(RVec(2)) > 0 && cons(u) < 0)
    {
        u += Pi;
    }
    else if (cons(RVec(2)) < 0 && cons(u) > 0)
    {
        u -= Pi;
    }
    T TrueA = u - w;
    if (cons(TrueA) < -Pi)
    {
        TrueA += 2 * Pi;
    }
    else if (cons(TrueA) > Pi)
    {
        TrueA -= 2 * Pi;
    }

    elem(0) = SemiA;
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
bool Elem2Cart(const Eigen::Array<T, Eigen::Dynamic, 1> & elem, 
    Eigen::Matrix<T, Eigen::Dynamic, 1> & cart) {

    T SemiA = elem(0);
    T Ecc = elem(1);
    T Inc = elem(2);
    T RAAN = elem(3);
    T w = elem(4);
    T TrueA = elem(5);

    T p = SemiA * (1 - Ecc * Ecc); // 半通径
    T h = sqrt(p * GM);
    T r = p / (1 + Ecc * cos(TrueA));

    Eigen::Matrix<T, 3, 3> M3NR, M1NI, M3Nw, M3NwPi2;
    RotationAxis(3, -RAAN, M3NR);
    RotationAxis(1, -Inc, M1NI);
    RotationAxis(3, -w, M3Nw);
    RotationAxis(3, -w - Pi / 2, M3NwPi2);
    Eigen::Matrix<T, 3, 1> RInPlane, UnitVec, ii0, jj0;
    RInPlane << r * cos(TrueA), r * sin(TrueA), 0;
    UnitVec << 1, 0, 0;
    ii0 = M3NR * M1NI * M3Nw * UnitVec;
    jj0 = M3NR * M1NI * M3NwPi2 * UnitVec;

    cart.head(3) = M3NR * M1NI * M3Nw * RInPlane;
    cart.tail(3) = -h / p * sin(TrueA) * ii0 + h / p * (Ecc + cos(TrueA)) * jj0;

    return true;
}

}