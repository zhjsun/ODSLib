// Orbit Dynamics and Safety Library
// 常数定义
// 孙振江，2016.01.19
/////////////////////////////////////////////////////////////

#pragma once

namespace ODS
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 常数定义
    const double Pi = 4 * atan(1.0);
    const double GM = 3.986004405e14;       // 地球引力常量，m^3/s^2
    const double EarthR = 6378137.0;        // 地球半径，m
    const double Earth_J2 = 1.0826261e-3;			///< 地球摄动J2项
    const double Earth_J3 = -2.54e-6;				///< 地球摄动J3项
    const double Earth_J4 = -1.61e-6;				///< 地球摄动J4项
}
