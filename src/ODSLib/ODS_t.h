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

/// Get the root of function for distance_point_ellipse by dichotomy
/// 点到椭圆距离涉及非线性方程的二分法求解
/// 
/// \Param[in]  r0
/// \Param[in]  z0
/// \Param[in]  z1
/// \Param[in]  g
/// \Return     s
template <typename T>
T GetRoot(T r0, T z0, T z1, T g) {
    T n0 = r0*z0;
    T s0 = z1 - 1;
    T s1;
    if (cons(g) < 0) {
        s1 = 0;
    } else {
        s1 = sqrt(Sqr(n0) + Sqr(z1)) - 1;
    }
    T s, ratio0, ratio1;
    do {
        s = (s0 + s1) / 2;
        if (s==s0||s==s1) {
            break;
        }
        ratio0 = n0 / (s + r0);
        ratio1 = z1 / (s + 1);
        g = Sqr(ratio0) + Sqr(ratio1) - 1;
        if (cons(g) > 0) {
            s0 = s;
        } else if (cons(g) < 0) {
            s1 = s;
        } else {
            break;
        }
    } while (fabs(cons(s1) - cons(s0) ) > 1e-6);
    return s;
}

/// Get the root of function for distance_point_ellipsoid by dichotomy
/// 点到椭球距离涉及非线性方程的二分法求解
/// 
/// \Param[in]  r0
/// \Param[in]  r1
/// \Param[in]  z0
/// \Param[in]  z1
/// \Param[in]  z2
/// \Param[in]  g
/// \Return     s
template <typename T>
T GetRoot(T r0, T r1, T z0, T z1, T z2, T g) {
    T n0 = r0*z0;
    T n1 = r1*z1;
    T s0 = z1 - 1;
    T s1;
    if (cons(g) < 0) {
        s1 = 0;
    } else {
        s1 = sqrt(Sqr(n0) + Sqr(n1) + Sqr(z2)) - 1;
    }
    T s, ratio0, ratio1, ratio2;
    do {
        s = (s0 + s1) / 2;
        if (s == s0 || s == s1) {
            break;
        }
        ratio0 = n0 / (s + r0);
        ratio1 = n1 / (s + r1);
        ratio2 = z2 / (s + 1);
        g = Sqr(ratio0) + Sqr(ratio1) + Sqr(ratio2) - 1;
        if (cons(g) > 0) {
            s0 = s;
        } else if (cons(g) < 0) {
            s1 = s;
        } else {
            break;
        }
    } while (fabs(cons(s1) - cons(s0) ) > 1e-6);
    return s;
}

/// Compute distance from a point to an ellipse
/// 点到椭圆距离求解
/// 
/// Distance from a Point to an Ellipse, an Ellipsoid, or a Hyperellipsoid
/// Geometric Tools, LLC, 
/// https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
/// 
/// \Param[in]  e0          椭圆第一个主轴，e0 >= e1
/// \Param[in]  e1          椭圆第二个主轴
/// \Param[in]  y0          点的第一个坐标，y0 >= 0
/// \Param[in]  y1          点的第二个坐标，y1 >= 0
/// \Param[out] x0          椭圆上最近点第一坐标
/// \Param[out] x1          椭圆上最近点第二坐标
/// \Return     distance    最近距离
template <typename T>
T DistancePointEllipse(T e0, T e1, T y0, T y1, T &x0, T &x1) {
    T distance;
    if (cons(y1) > 0) {
        if (cons(y0) > 0) {
            T z0 = y0 / e0;
            T z1 = y1 / e1;
            T g = Sqr(z0) + Sqr(z1) - 1;
            if (cons(g) != 0) {
                T r0 = Sqr(e0 / e1);
                T sbar = GetRoot(r0, z0, z1, g);
                x0 = r0*y0 / (sbar + r0);
                x1 = y1 / (sbar + 1);
                distance = sqrt(Sqr(x0 - y0) + Sqr(x1 - y1));
            } else {
                x0 = y0;
                x1 = y1;
                distance = 0;
            }
        } else {   // y0=0
            x0 = 0;
            x1 = e1;
            distance = fabs(y1 - e1);
        }
    }
    else {    // y1=0
        T numer0 = e0*y0;
        T denom0 = Sqr(e0) - Sqr(e1);
        if (cons(numer0) < cons(denom0) ) {
            T xde0 = numer0 / denom0;
            x0 = e0*xde0;
            x1 = e1*sqrt(1 - Sqr(xde0));
            distance = sqrt(Sqr(x0 - y0) + Sqr(x1));
        } else {
            x0 = e0;
            x1 = 0;
            distance = fabs(y0 - e0);
        }
    }
    return distance;
}

/// Compute distance from a point to an ellipsoid
/// 点到椭球距离求解
/// 
/// Distance from a Point to an Ellipse, an Ellipsoid, or a Hyperellipsoid
/// Geometric Tools, LLC, 
/// https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
/// 
/// \Param[in]  e0          椭球第一个主轴，e0 >= e1 >= e2 >0
/// \Param[in]  e1          椭球第二个主轴
/// \Param[in]  e2          椭球第三个主轴
/// \Param[in]  y0          点的第一个坐标，y0 >= 0
/// \Param[in]  y1          点的第二个坐标，y1 >= 0
/// \Param[in]  y2          点的第三个坐标，y2 >= 0
/// \Param[out] x0          椭圆上最近点第一坐标
/// \Param[out] x1          椭圆上最近点第二坐标
/// \Param[out] x2          椭圆上最近点第三坐标
/// \Return	    distance    最近距离
template <typename T>
T DistancePointEllipsoid(T e0, T e1, T e2, T y0, T y1, T y2, 
                         T &x0, T &x1, T &x2) {
    T distance;
    if (cons(y2) > 0) {
        if (cons(y1) > 0) {
            if (cons(y0) > 0) {
                T z0 = y0 / e0;
                T z1 = y1 / e1;
                T z2 = y2 / e2;
                T g = Sqr(z0) + Sqr(z1) + Sqr(z2) - 1;
                if (cons(g) != 0) {
                    T r0 = Sqr(e0 / e2);
                    T r1 = Sqr(e1 / e2);
                    T sbar = GetRoot(r0, r1, z0, z1, z2, g);
                    x0 = r0*y0 / (sbar + r0);
                    x1 = r1*y1 / (sbar + r1);
                    x2 = y2 / (sbar + 1);
                    distance = sqrt(Sqr(x0 - y0) + Sqr(x1 - y1) + Sqr(x2 - y2));
                } else {
                    x0 = y0;
                    x1 = y1;
                    x2 = y2;
                    distance = 0;
                }

            } else { // y0 =0
                x0 = 0;
                distance = DistancePointEllipse(e1, e2, y1, y2, x1, x2);
            }
        } else {    // y1 = 0
            if (cons(y0) > 0) {
                x1 = 0;
                distance = DistancePointEllipse(e0, e2, y0, y2, x0, x2);
            } else {    // y0 = 0
                x0 = 0;
                x1 = 0;
                x2 = e2;
                distance = fabs(y2 - e2);
            }
        }
    } else {   // y2 = 0
        T numer0 = e0*y0;
        T numer1 = e1*y1;
        T denom0 = Sqr(e0) - Sqr(e2);
        T denom1 = Sqr(e1) - Sqr(e2);
        bool computed = false;
        if ((numer0 < denom0) && (numer1 < denom1) ) {
            T xde0 = numer0 / denom0;
            T xde1 = numer1 / denom1;
            T discr = 1 - Sqr(xde0) - Sqr(xde1);
            if (discr > 0) {
                x0 = e0*xde0;
                x1 = e1*xde1;
                x2 = e2*sqrt(discr);
                distance = sqrt(Sqr(x0 - y0) + Sqr(x1 - y1) + Sqr(x2));
                computed = true;
            }
        }
        if (!computed) {
            x2 = 0;
            distance = DistancePointEllipse(e0, e1, y0, y1, x0, x1);
        }
    }
    return distance;
}

/// Compute closet distance from a point to an ellipsoid
/// 计算点和椭球之间的最近距离
/// 
/// \param[in]  pointPos        点的坐标，在椭球主轴坐标系中表示
/// \param[in]  ellipAx1        椭球的主轴长度x
/// \param[in]  ellipAx2        椭球的主轴长度y
/// \param[in]  ellipAx3        椭球的主轴长度z
/// \param[out] distanceVec     最小距离矢量，点到椭球
/// \param[out] posAtEllip      椭球面上最小距离对应的点，在椭球主轴坐标系中表示
/// \return true                计算成功
///         false               点不在椭球外部
template <typename T>
bool DistancePointEllipsoid(const Eigen::Matrix<T, 3, 1> &pointPos, 
                            T ellipAx1, T ellipAx2, T ellipAx3, 
                            Eigen::Matrix<T, 3, 1> &distanceVec, 
                            Eigen::Matrix<T, 3, 1> &posAtEllip) {
    Eigen::Matrix<T, 3, 1> s;          // 点的坐标符号变换标志
    Eigen::Matrix<T, 3, 1> e;          // 椭球主轴长度
    Eigen::Matrix<T, 3, 1> y;          // 点的坐标，变换到第一象限
    Eigen::Matrix<T, 3, 1> x;          // 椭球上最近点坐标

    // 将点变换到第一象限
    s.fill(1);
    for (int i = 0; i < 3; i++) {
        if (pointPos(i) < 0) {
            s(i) = -1;
        }
    }
    y = pointPos.cwiseProduct(s);
    e(0) = ellipAx1;
    e(1) = ellipAx2;
    e(2) = ellipAx3;

    // 判断点与椭球关系
    if ((Sqr(y(0)/e(0))+ Sqr(y(1) / e(1))+ Sqr(y(2) / e(2)))<=1) {
        return false;
    } else {
        // 椭球主轴降序排列
        T etmp, stmp, ytmp;
        for (int i = 0; i < 2; i++) {
            for (int j = i+1; j < 3; j++) {
                if (e(i)<e(j)) {
                    etmp = e(i);
                    e(i) = e(j);
                    e(j) = etmp;
                    stmp = s(i);
                    s(i) = s(j);
                    s(j) = stmp;
                    ytmp = y(i);
                    y(i) = y(j);
                    y(j) = ytmp;
                }
            }
        }

        DistancePointEllipsoid(e(0), e(1), e(2), y(0), y(1), y(2), x(0), x(1), x(2));

        posAtEllip = x.cwiseProduct(s);
        distanceVec = posAtEllip - pointPos;

        return true;
    }
}

/// Compute closet distance between two ellipsoids
/// 计算椭球和椭球之间的最近距离
/// 
/// \param[in]  ellipAx11           第一个椭球的主轴长度x
/// \param[in]  ellipAx12           第一个椭球的主轴长度y
/// \param[in]  ellipAx13           第一个椭球的主轴长度z
/// \param[in]  ellipAx21           第二个椭球的主轴长度x
/// \param[in]  ellipAx22           第二个椭球的主轴长度y
/// \param[in]  ellipAx23           第二个椭球的主轴长度z
/// \param[in]  mtx1To2             从第一个椭球的主轴坐标系转到第二个椭球的主轴坐标系的矩阵
/// \param[in]  ellip2AtEllip1Pos   第二个椭球中心在第一个椭球主轴坐标系中的位置，在第一个椭球主轴坐标系中表示
/// \param[out] distanceVec12       最小距离矢量，从第一个椭球到第二个椭球，在第一个椭球主轴坐标系中表示
/// \param[out] posAtEllip1         第一个椭球面上最小距离对应的点，在第一个椭球主轴坐标系中表示
/// \param[out] posAtEllip2         第二个椭球面上最小距离对应的点，在第二个椭球主轴坐标系中表示
/// \return     true                计算成功
///             false               椭球接触或相交，或者没有在规定次数内收敛
template <typename T>
bool DistanceEllipsoid(T ellipAx11, T ellipAx12, T ellipAx13, 
                       T ellipAx21, T ellipAx22, T ellipAx23,
                       const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mtx1To2, 
                       const Eigen::Matrix<T, 3, 1> &ellip2AtEllip1Pos, 
                       Eigen::Matrix<T, 3, 1> &distanceVec12, 
                       Eigen::Matrix<T, 3, 1> &posAtEllip1, 
                       Eigen::Matrix<T, 3, 1> &posAtEllip2) {
    //最大迭代次数50
    const int ITER = 50;
    const T EPSD = 1.0e-8;  T maxAx, minAx;
    Eigen::Matrix<T, 3, 1> point, distanceVec21;  //取第一个椭球上一点
    maxAx = std::max<T>(ellipAx11, ellipAx12);
    maxAx = std::max<T>(maxAx, ellipAx13);
    minAx = std::min<T>(ellipAx11, ellipAx12);
    minAx = std::min<T>(minAx, ellipAx13);
    point[0] = point[1] = sqrt(0.5*minAx*minAx);
    point[2] = ellipAx13*sqrt(fabs(1.0 - point[0] * point[0] / (ellipAx11*ellipAx11)
               - point[1] * point[1] / (ellipAx12*ellipAx12)));
    posAtEllip1 = point;

    for (int iter = 0; iter<ITER; iter++) {
        point = mtx1To2*(posAtEllip1 - ellip2AtEllip1Pos);
        if (!ODS::DistancePointEllipsoid(point, ellipAx21, ellipAx22, ellipAx23,
            distanceVec21, posAtEllip2)) {
            return false;
        }

        point = mtx1To2.transpose()*posAtEllip2 + ellip2AtEllip1Pos;
        if (!ODS::DistancePointEllipsoid(point, ellipAx11, ellipAx12, ellipAx13,
            distanceVec12, posAtEllip1)) {
            return false;
        }
        if (fabs(distanceVec12.norm() - distanceVec21.norm() ) < EPSD) {
            return true;
        }
    }
    return false;
}

/// Sort eigen values and vectors in descending order
/// 对特征值和特征向量按照降序排列
/// 
/// @Param[in,out]  eigenvalue
/// @Param[in,out]  eigenvector
template <typename T>
void EigenDesSort(Eigen::Matrix<T, Eigen::Dynamic, 1> &eigenvalue, 
                  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &eigenvector) {
    int N = eigenvalue.size();

    T tmp;
    Eigen::Matrix<T, Eigen::Dynamic, 1> Vtmp(N);

    for (int i = 0; i < N-1; i++) {
        for (int j = i+1; j < N; j++) {
            if (eigenvalue(i)<eigenvalue(j)) {
                tmp = eigenvalue(i);
                eigenvalue(i) = eigenvalue(j);
                eigenvalue(j) = tmp;
                Vtmp = eigenvector.col(i);
                eigenvector.col(i) = eigenvector.col(j);
                eigenvector.col(j) = Vtmp;
            }
        }
    }
}

}