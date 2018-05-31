#ifndef ANN_H
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>

typedef struct
{
    int n;
    double (*f)(double);
    gsl_vector* data;
} ann;

ann* ann_alloc(int n, double(*act_fun)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist,double eps,double step);
int newton_broyden(double f(gsl_vector* x), gsl_vector* x, double dx, double eps);

#define ANN_h
#endif
