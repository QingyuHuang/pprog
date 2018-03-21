#include <stdio.h>
#include <float.h>
#include <math.h>
#include"komplex.h"

// Define global precision parameters
#ifndef TAU
#define TAU 1e-6
#endif
#ifndef EPS
#define EPS 1e-6
#endif

const komplex komplex_I = { 0, 1 };


void komplex_print(char* s, komplex z) {
    printf("%s %g + %g i\n", s, z.re, z.im);
}


void komplex_set(komplex* z, double x, double y) {
    z->re = x;
    z->im = y;
    //    (*z).re = x;
    //    (*z).im = y;
}


komplex komplex_new(double x, double y) {
    komplex z;
    komplex_set(&z, x, y);
    return z;
}


komplex komplex_add(komplex a, komplex b) {
    double x = a.re + b.re;
    double y = a.im + b.im;
    komplex result = {.re = x, .im = y};
    return result;
}


komplex komplex_sub(komplex a, komplex b) {
    double x = a.re - b.re;
    double y = a.im - b.im;
    komplex result = {.re = x, .im = y};
    return result;
}


int double_equal(double a, double b) {
    if (fabs(a - b) < TAU)
        return 1;
    if (fabs(a - b) / (fabs(a) + fabs(b)) < (EPS / 2))
        return 1;
    else
        return 0;
}


int komplex_equal(komplex a, komplex b ) {
    int result = double_equal(a.re, b.re) && double_equal(a.im, b.im);
    return result;
}


komplex komplex_mul(komplex a, komplex b) {
    double x = (a.re * b.re) - (a.im * b.im);
    double y = (a.re * b.im) + (a.im * b.re);
    komplex result = {.re = x, .im = y};
    return result;
}


komplex komplex_div(komplex a, komplex b) {
    if (fabs(b.im) < fabs(b.re)) {
        double e = b.im / b.re;
        double f = b.re + b.im * e;
        komplex result =
        {.re = (a.re + a.im * e) / f,.im = (a.im - a.re * e) / f };
        return result;
    } else {
        double e = b.re / b.im;
        double f = b.im + b.re * e;
        komplex result =
        {.re = (a.im + a.re * e) / f,.im = (-a.re + a.im * e) / f };
        return result;
    }
}


komplex komplex_conjugate(komplex z) {
    komplex result = {z.re = z.re, z.im = -z.im};
    return result;
}


komplex komplex_abs(komplex z) {
    double abs;
    abs = sqrt(pow(z.re, 2) + pow(z.im, 2));
    komplex result = komplex_new(abs, 0);
    return result;
}


komplex komplex_exp(komplex a) {
    komplex z;
    z.re = cos(a.im) * exp(a.re);
    z.im = sin(a.im) * exp(a.re);
    return z;
}


komplex komplex_sin(komplex z) {
    komplex a;
    a.re = sin(z.re) * cosh(z.im);
    a.im = cos(z.re) * sinh(z.im);
    return a;
}


komplex komplex_cos(komplex z) {
    komplex a;
    a.re = cos(z.re) * cosh(z.im);
    a.im = - sin(z.re) * sinh(z.im);
    return a;
}


komplex komplex_sqrt(komplex z) {
    double r, theta;
    komplex w, sqrt_z;
    theta = atan(z.im / z.re) / 2;
    w = komplex_abs(z);
    r = sqrt(w.re);
    sqrt_z.re = r * cos(theta);
    sqrt_z.im = r * sin(theta);
    return sqrt_z;
}
