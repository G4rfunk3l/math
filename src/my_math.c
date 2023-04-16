#include "my_math.h"

int my_abs(int x) {
    return x > 0 ? x : -x;
}


long double my_acos(double x) {
    long double res = my_NAN;
    if (-1 <= x && x <= 1) {
        long double a = 0, b = my_M_PI;
        res = 0.5 * (a + b);
        while (my_fabs(a - b) > 1e-7) {
            long double y = my_cos(res);
            if (y < x) b = res;
            if (y > x) a = res;
            res = 0.5 * (a + b);
        }
    }
    return res;
}


long double my_asin(double x) {
    long double res = my_NAN;
    if (-1e-7 < x && x < 1e-7) {
        res = 0.0;
    } else if (-1 <= x && x <= 1) {
        long double a = -my_M_PI / 2, b = my_M_PI / 2;
        res = 0.5 * (a + b);
        while (my_fabs(a - b) > 1e-7) {
            long double y = my_sin(res);
            if (y < x) a = res;
            if (y > x) b = res;
            res = 0.5 * (a + b);
        }
    }
    return res;
}


long double my_atan(double x) {
    long double res = 0.0;
    if (my_isnan(x)) {
        res = my_NAN;
    } else if (x <= -1e-7 || 1e-7 <= x) {
        long double a = -my_M_PI / 2, b = my_M_PI / 2;
        res = 0.5 * (a + b);
        while (my_fabs(a - b) > 1e-7) {
            long double y = my_tan(res);
            if (y < x) a = res;
            if (y > x) b = res;
            res = 0.5 * (a + b);
        }
    }
    return res;
}


long double my_ceil(double x) {
    long double result;
    double prm = x;
    int total = prm;
    if (x == -my_INFINITY || x == my_INFINITY) {
        result = x;
    } else if (prm > (long double)total) {
        result = total + 1;
    } else if (prm < (long double)total) {
        result = total;
    } else if (prm == total) {
        result = total;
    }
    return x != x ? my_NAN : result;
}


long double my_cos(double x) {
    x = my_fmod(x, 2.0 * my_M_PI);
    long double res = 0, last = 1;
    for (int k = 1; my_fabs(last) > EPS_10; ++k) {
        res += last;
        last *= -x * x /((2.0 * k - 1.0)*(2.0 * k));
    }
    return res;
}


long double my_exp(double x) {
    long double y = (long double)x;
    long double x1 = 1.0, sum = 0.0;
    int precision = 100000;

    for (int i = 0; i < precision; ) {
        sum += my_fabs(x1);
        x1  *= (y/++i);
    }

    sum = x < 0 ? 1/sum : sum;
    return sum;
}


long double my_fabs(double x) {
    return x > 0 ? x : -x;
}


long double my_floor(double x) {
    long double result;
    double prm = x;
    int total = prm;
    if (x == -my_INFINITY || x == my_INFINITY) {
        result = x;
    } else if (prm < (long double)total) {
        result = total - 1;
    } else if (prm > (long double)total) {
        result = total;
    } else if (prm == total) {
        result = total;
    }
    return x != x ? my_NAN : result;
}


long double my_fmod(double x, double y) {
    if (y <= -my_INFINITY || y >= my_INFINITY) return x;
    return y != 0 || (x && y) ? (x / y - (long int)(x / y)) * y : my_NAN;
}


long double my_log(double x) {
    return my_logl(x);
}


long double my_pow(double base, double exp) {
    long double powed = 0.0;
    int int_exp = (int) exp;
    if (base == 0.0 && exp != 0.0) {
        powed = base;
    } else if (exp == 0.0) {
        powed = 1.0;
    } else if (exp == int_exp) {
        if (int_exp > 0) {
            powed = base;
            int_exp--;
            while (int_exp) {
                powed *= base;
                int_exp--;
            }
        } else if (int_exp < 0) {
            powed = 1 / base;
            int_exp++;
            while (int_exp) {
                powed /= base;
                int_exp++;
            }
        }
    } else {
    powed = my_exp(exp * my_logl(base));
    }

    return powed;
}


long double my_sin(double x) {
    x = my_fmod(x, 2.0 * my_M_PI);
    long double res = x, last = -x * x * x / 6;
    for (int k = 2; my_fabs(last) > EPS_10; ++k) {
        res += last;
        last *= -x * x /((2.0 * k + 1.0)*(2.0 * k));
    }
    return res;
}


long double my_sqrt(double x) {
    return my_sqrtl(x);
}


long double my_tan(double x) {
    return (my_sin(x)/my_cos(x));
}



int my_isnan(double x) {
    return x != x;
}


long double my_sqrtl(long double x) {
    long double y = x;
    if (y < 0.0) {
        y = my_NAN;
    } else {
        long double y_temp = 0.0;
        while (y != y_temp) {
            y_temp = y;
            y = 0.5 *(y_temp + (x / y_temp));
        }
    }
    return y;
}


long double my_logl(long double x) {
    long double y;
    if (x < 0.0) {
        y = my_NAN;
    } else if (x == 0.0) {
        y = -my_INFINITY;
    } else {
        y = (my_M_PI / (2.0 * AGMl(x))) - (16.0 * my_M_LN2);
    }

    return y;
}


long double AGMl(long double x) {
    long double a = 1.0;
    long double b = 0.000061035156250 / x;
    long double e = 1e-17;
    while (a > (b + e) || a < (b - e)) {
        long double a_temp = a;
        long double b_temp = b;
        a = (a_temp + b_temp) / 2.0;
        b = my_sqrtl(a_temp * b_temp);
    }
    return a;
}
