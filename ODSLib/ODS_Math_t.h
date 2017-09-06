/// \file ODS_t.h
/// This file contains the implementations for the methods in "ODS.h"
/// 
/// \author Sun, Zhenjiang
/// \date 2016.Dec.20

#pragma once

namespace ODS {

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