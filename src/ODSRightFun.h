// Orbit Dynamics and Safety Library
// 右函数类头文件
// 孙振江，2017.01.19
/////////////////////////////////////////////////////////////

#pragma once
#include <Eigen/Dense>

using namespace Eigen;

namespace ODS
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//右函数类，用于其它右函数计算类的基类
//
template <typename T>
class CRightFun
{
  public:
    CRightFun(){};
    virtual ~CRightFun(){};

    /// @Input
    /// @Param	t		自变量的值
    /// @Param	x		初始函数值
    /// @Output
    /// @Param	result	计算得到的函数值
    virtual void operator()(double t, const Matrix<T, Dynamic, 1> &x, Matrix<T, Dynamic, 1> &result) const = 0;
};

//
//右函数类，用于其它右函数计算类的基类
//double * version
//
class AsTRightFun
{
  public:
    AsTRightFun(){};
    virtual ~AsTRightFun(){};

    /// @Input
    /// @Param	t		自变量的值
    /// @Param	n		初始函数值的个数
    /// @Param	x		初始函数值
    /// @Output
    /// @Param	result	计算得到的函数值
    virtual void operator()(double t, int xN, const double *x, double *result) const = 0;
};

//
//输入1个变量,输出1个变量的函数类
//
class CFun11
{
  public:
    CFun11(){};
    virtual ~CFun11(){};

    /// @Input
    /// @Param	x		自变量的值
    /// @Output
    /// @Param	result	计算得到的函数值
    /// @Return	计算是否正常，上层程序检测到这个异常后应该退出
    virtual bool operator()(double x, double &result) const = 0;
};

//
//输入2个变量,输出1个变量的函数类
//
class CFun21
{
  public:
    CFun21(){};
    virtual ~CFun21(){};

    /// @Input
    /// @Param	x1		自变量的值
    /// @Param	x2		自变量的值
    /// @Output
    /// @Param	result	计算得到的函数值
    virtual void operator()(double x1, double x2, double &result) const = 0;
};

//
//输入3个变量,输出1个变量的函数类
//
class CFun31
{
  public:
    CFun31(){};
    virtual ~CFun31(){};

    /// @Input
    /// @Param	x1		自变量的值
    /// @Param	x2		自变量的值
    /// @Param	x3		自变量的值
    /// @Output
    /// @Param	result	计算得到的函数值
    virtual void operator()(double x1, double x2, double x3, double &result) const = 0;
};

//
//输入2个变量,输出2个变量的函数类
//
class CFun22
{
  public:
    CFun22(){};
    virtual ~CFun22(){};

    /// @Input
    /// @Param	x1		自变量的值
    /// @Param	x2		自变量的值
    /// @Output
    /// @Param	result1	计算得到的函数值
    /// @Param	result2	计算得到的函数值
    virtual void operator()(double x1, double x2, double &result1, double &result2) const = 0;
};

//
//输入多个变量,输出1个变量的函数类
//
class CFunN1
{
  public:
    CFunN1(){};
    virtual ~CFunN1(){};

    /// @Input
    /// @Param	xN		自变量的个数
    /// @Param	x		自变量的值
    /// @Output
    /// @Param	result	计算得到的函数值
    virtual void operator()(int xN, const double *x, double &result) const = 0;
};

//
//输入多个变量,输出多个变量的函数类
//
class CFunNN
{
  public:
    CFunNN(){};
    virtual ~CFunNN(){};

    /// @Input
    /// @Param	xN		自变量的个数
    /// @Param	x		自变量的值
    /// @Param	resultN	结果的个数
    /// @Output
    /// @Param	result	计算得到的函数值
    virtual bool operator()(int xN, const double *x, int resultN, double *result) const = 0;
};
}
