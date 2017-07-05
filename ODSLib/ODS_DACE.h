/// \file ODS_DACE.h
/// This file contains some functions in DACE,
/// to make ODS compatible with and without DACE.
/// 
/// \author Sun, Zhenjiang
/// \date 2017.Jun.16

#pragma once

namespace ODS {

template <typename T>
T cons(const T x) {return x;}

}