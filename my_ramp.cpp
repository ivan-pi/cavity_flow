#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" {
#include "incbeta.h"
};
#endif

#include "my_ramp.hpp"

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath> 
#include <bits/stdc++.h> 

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double my_ramp(double x, double p) {
    double f;

    double xsqr = abs(x)*abs(x);
    f = sgn(x)*incbeta(0.5,p+1,xsqr);
    return 0.5*(1 + f);
}
