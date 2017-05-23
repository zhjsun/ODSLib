// Orbit Dynamics and Safety Library
// 孙振江，2016.12.20
/////////////////////////////////////////////////////////////

#pragma once
#include <cstdlib>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "ODSConstant.h"
#include "ODSRightFun.h"

using namespace Eigen;

namespace ODS
{

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 坐标变换相关函数

    //********************************************************************
    /// 坐标旋转矩阵
    /// @Author	孙振江
    /// @Date	2016.12.22
    /// @Input	
    /// @Param	axis	制定旋转轴[1, 2, 3]
    /// @Param	alpha	旋转角度[rad]
    /// @Output	
    /// @Param	Mtx 	旋转矩阵
    /// @Return			true=计算正确; false=输入数据异常
    //********************************************************************
    template<typename T> bool RotationAxis(const int axis, const T alpha, Matrix<T, 3, 3> & Mtx)
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

    //********************************************************************
    /// 计算从VVLH坐标系到输入直角坐标系（地心惯性系或地固系）的转换矩阵
    ///  VVLH定义为：z轴指向地心，x轴与z轴垂直指向速度方向，y轴与其它两轴成右手坐标系.
    /// Return transformation matrix from VVLH(Vehicle Velocity Local Horizontal) to ICS.
    /// @Author	孙振江
    /// @Date	2016.12.20
    /// @Input	
    /// @Param	pos		the position of vehicle 地心惯性系或地固系中的飞行器位置[m]
    /// @Param	vel		velocity of vehicle 地心惯性系或地固系中的飞行器速度[m/s]
    /// @Output	
    /// @Param	mtx		VVLH到ICS的转移矩阵
    /// @Return			true=计算正确; false=输入数据异常
    //********************************************************************
    bool VVLHToICSMtx(const Vector3d& pos, const Vector3d& vel, Matrix3d& mtx);

    //********************************************************************
    /// 计算从直角坐标系到VVLH坐标系的转换矩阵
    /// VVLH定义为：z轴指向地心，x轴与z轴垂直指向速度方向，y轴与其它两轴成右手坐标系.
    /// Return transformation matrix from VVLH(Vehicle Velocity Local Horizontal) to ICS.
    /// @Author	孙振江
    /// @Date	2016.12.20
    /// @Input	
    /// @Param	pos		the position of vehicle 地心惯性系或地固系中的飞行器位置[m]
    /// @Param	vel		velocity of vehicle 地心惯性系或地固系中的飞行器速度[m/s]
    /// @Output	
    /// @Param	mtx		ICS到VVLH的转移矩阵
    /// @Return			true=计算正确; false=输入数据异常
    //********************************************************************
    bool ICSToVVLHMtx(const Vector3d& pos, const Vector3d& vel, Matrix3d& mtx);

    //********************************************************************
    /// 计算从直角坐标系到VVLH坐标系的转换矩阵
    /// VVLH定义为：z轴指向地心，x轴与z轴垂直指向速度方向，y轴与其它两轴成右手坐标系.
    /// Return transformation matrix from VVLH(Vehicle Velocity Local Horizontal) to ICS.
    /// @Author	孙振江
    /// @Date	2016.12.20
    /// @Input	
    /// @Param	pos		the position of vehicle 地心惯性系或地固系中的飞行器位置[m]
    /// @Param	vel		velocity of vehicle 地心惯性系或地固系中的飞行器速度[m/s]
    /// @Output	
    /// @Param	mtx		ICS到VVLH的转移矩阵
    /// @Return			true=计算正确; false=输入数据异常
    //********************************************************************
    bool ICSToVVLH(const VectorXd& Target, const VectorXd& Chaser, VectorXd& RelState);

    //********************************************************************
    /// 计算从直角坐标到轨道根数
    /// @Author	孙振江
    /// @Date	2016.12.21
    /// @Input	
    /// @Param	cart	绝对轨道直角坐标[m, m/s]
    /// @Output	
    /// @Param	elem	轨道根叔[m, rad]
    /// @Return			true=计算正确; false=输入数据异常
    //********************************************************************
    //bool Cart2Elem(const VectorXd & cart, Array<double, Dynamic, 1> & elem);
    //bool Cart2Elem(const Matrix<DA, Dynamic, 1> & cart, Array<DA, Dynamic, 1> & elem);
    template<typename T> bool Cart2Elem(const Matrix<T, Dynamic, 1> & cart, Array<T, Dynamic, 1> & elem)
    {
        double GM_km = GM*1e-9;

        Matrix<T, 3, 1> RVec, VVec;                                // 距离和速度矢量
        RVec = cart.head(3)/1000;                      // km
        VVec = cart.tail(3)/1000;                      // km/s
        T rr = RVec.norm();                            // 地心距, km
        // 计算动量矩
        Matrix<T, 3, 1> HVec;
        HVec = RVec.cross(VVec);
        T hh = HVec.norm();
        // 计算轨道倾角
        T Inc = acos(HVec(2) / hh);
        // 计算升交点赤经
        T RAAN = atan(-HVec(0) / HVec(1));
        if (cons(HVec(1))>0)
        {
            RAAN += Pi;
        }
        if (cons(RAAN) < 0)
        {
            RAAN += 2 * Pi;
        }
        // 计算偏心率矢量和大小
        Matrix<T, 3, 1> EVec;
        EVec = VVec.cross(HVec)/GM_km - RVec/rr;
        T Ecc = EVec.norm();
        // 计算近地点纬度幅角
        T w = atan(EVec(2) / ((EVec(1)*sin(RAAN) + EVec(0)*cos(RAAN)) * sin(Inc)));
        if (cons(EVec(2)) > 0 && cons(w) < 0)
        {
            w += Pi;
        }
        else if(cons(EVec(2)) < 0 && cons(w) > 0)
        {
            w -= Pi;
        }
        // 计算半长轴
        T SemiA = hh*hh / (GM_km*(1 - Ecc*Ecc)) * 1000;         // m
        // 计算真近点角
        T u = atan(RVec(2) / ((RVec(1)*sin(RAAN) + RVec(0)*cos(RAAN)) * sin(Inc)));
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

    //********************************************************************
    /// 计算从轨道根数到直角坐标
    /// @Author	孙振江
    /// @Date	2016.12.22
    /// @Input	
    /// @Param	elem	轨道根叔[m, rad]
    /// @Output	
    /// @Param	cart	绝对轨道直角坐标[m, m/s]
    /// @Return			true=计算正确; false=输入数据异常
    //********************************************************************
    template<typename T> bool Elem2Cart(const Array<T, Dynamic, 1> & elem, Matrix<T, Dynamic, 1> & cart)
    {
        T SemiA = elem(0);
        T Ecc = elem(1);
        T Inc = elem(2);
        T RAAN = elem(3);
        T w = elem(4);
        T TrueA = elem(5);
    
        T p = SemiA *(1 - Ecc*Ecc);            // 半通径
        T h = sqrt(p*GM);
        T r = p / (1 + Ecc*cos(TrueA));
    
        Matrix<T, 3, 3> M3NR, M1NI, M3Nw, M3NwPi2;
        RotationAxis(3, -RAAN, M3NR);
        RotationAxis(1, -Inc, M1NI);
        RotationAxis(3, -w, M3Nw);
        RotationAxis(3, -w-Pi/2, M3NwPi2);
        Matrix<T, 3, 1> RInPlane, UnitVec, ii0, jj0;
        RInPlane << r*cos(TrueA), r*sin(TrueA), 0;
        UnitVec << 1, 0, 0;
        ii0 = M3NR*M1NI*M3Nw*UnitVec;
        jj0 = M3NR*M1NI*M3NwPi2*UnitVec;
    
        cart.head(3) = M3NR*M1NI*M3Nw*RInPlane;
        cart.tail(3) = -h / p*sin(TrueA)*ii0 + h / p*(Ecc + cos(TrueA))*jj0;
    
        return true;
    }

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
    template<typename T> T Sqr(T x){ return x*x;}

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
    double DistancePointEllipse(double e0, double e1, double y0, double y1, double& x0, double& x1);

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
    double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2, double& x0, double& x1, double& x2);

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
    bool DistancePointEllipsoid(const Vector3d& pointPos, double ellipAx1, double ellipAx2, double ellipAx3, Vector3d& distanceVec, Vector3d& posAtEllip);

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
        const MatrixXd& mtx1To2, const Vector3d& ellip2AtEllip1Pos, Vector3d& distanceVec12, Vector3d& posAtEllip1, Vector3d& posAtEllip2);

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
    void EigenDesSort(VectorXd& eigenvalue, MatrixXd& eigenvector);

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
    template<typename T>
    void RK4Step(const CRightFun<T>& f, double t, double h, double eps, Matrix<T, Dynamic, 1>& y)
    {
        int n = y.size();
        int m,i,j,k;
        double hh,p,dt,x,tt,q,a[4];
        Matrix<T, Dynamic, 1> g, b, c, d, e;
        g.resize(n);
        b.resize(n);
        c.resize(n);
        d.resize(n);
        e.resize(n);
        hh=h; m=1; p=1.0+eps; x=t;
        c=y;
        while (p>=eps)
          { a[0]=hh/2.0; a[1]=a[0]; a[2]=hh; a[3]=hh;
            g=y;
            y=c;
            dt=h/m; t=x;
            for (j=0; j<=m-1; j++)
              { f(t,y,d);
                b = y;
                e = y;
                for (k=0; k<=2; k++)
                  {
                     y=e+a[k]*d;
                     b=b+a[k+1]*d/3.0;
                     tt=t+a[k];
                     f(tt,y,d);
                  }
                y=b+hh*d/6.0;
                t=t+dt;
              }
            p=0.0;
            for (i=0; i<=n-1; i++)
              { q=abs(cons(y[i])-cons(g[i]));
                if (q>p) p=q;
              }
            hh=hh/2.0; m=m+m;
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
    void RK78(const CRightFun<double>& f, double T0, double T1, VectorXd& X);

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
    double IntLRGS(const CFun11& f, double a, double b, double eps=1.0e-10);

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
    double IntLRGS3D(const CFun31& func, const CFun11& funcy1, const CFun11& funcy2, const CFun21& funcz1, const CFun21& funcz2,
        double x1, double x2, int js[3]);

}

