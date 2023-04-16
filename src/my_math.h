#ifndef SRC_MY_MATH_H_
#define SRC_MY_MATH_H_
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define my_M_PI 3.141592653589793116
#define my_M_LN2 0.693147180559945309417232121458176568
#define my_EXP 2.71828182845904523536028747135266250
#define my_INFINITY 1.0 / 0.0
#define my_NAN 0.0 / 0.0
#define EPS_10 1e-10
#define EPS_6 1e-6

int my_abs(int x);
long double my_acos(double x);
long double my_asin(double x);
long double my_atan(double x);
long double my_ceil(double x);
long double my_cos(double x);
long double my_exp(double x);
long double my_fabs(double x);
long double my_floor(double x);
long double my_fmod(double x, double y);
long double my_log(double x);
long double my_pow(double base, double exp);
long double my_sin(double x);
long double my_sqrt(double x);
long double my_tan(double x);
long double my_sqrtl(long double x);
long double my_logl(long double x);
long double AGMl(long double x);
int my_isnan(double x);

#endif  // SRC_MY_MATH_H_
