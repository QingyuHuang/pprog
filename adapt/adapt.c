#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include<assert.h>

double adapt24(double f(double x),double a,double b,double acc,double eps, double f2, double f3,int nrec,double *err);
double my_integrator(double f(double x),double a, double b, double acc, double eps, double *err);
double my_integrator_infinity(double f(double x),double a, double b, double acc, double eps,double *err);
double clenshaw_curtis(double f(double x),double a, double b, double acc, double eps, double *err);

double adapt24(double f(double x),double a,double b,double acc,double eps, double f2, double f3,int nrec,double *err)
{

    assert(nrec<1e6);

    double f1 = f(a+(b-a)/6);
    double f4 = f(a+5*(b-a)/6);
    double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
    double q = (f1+f2+f3+f4)/4*(b-a);
    double tol = acc+eps*fabs(Q);
    *err = fabs(Q-q);

    if(*err<tol)
    {
        return Q;
    }
    else
    {
        double Q1 = adapt24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1,err);
        double Q2 = adapt24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1,err);
        return Q1+Q2;
    }
}

double my_integrator(double f(double x),double a, double b, double acc, double eps, double *err)
{
    double f2 = f(a+2*(b-a)/6);
    double f3 = f(a+4*(b-a)/6);
    int nrec = 0;
    return adapt24(f,a,b,acc,eps,f2,f3,nrec,err);
}

double my_integrator_infinity(double f(double x),double a, double b, double acc, double eps,double *err)
{

    int acheck=isinf(-a);
    int bcheck=isinf(b);

    if(acheck && bcheck)
    {
        double g(double t)
        {
            return f(t/(1-t*t))*(1+t*t)/((1-t*t)*(1-t*t));
        };
        double a0 = -1, b0 = 1;
        double f2 = f(a0+2*(b0-a0)/6);
        double f3 = f(a0+4*(b0-a0)/6);
        int nrec = 0;
        return adapt24(g,a0,b0,acc,eps,f2,f3,nrec,err);

    }
    else if(acheck)
    {
        double g(double t)
        {
            return f(b-(1-t)/t)*1/(t*t);
        };
        double a0 = 0, b0 = 1;
        double f2 = f(a0+2*(b0-a0)/6);
        double f3 = f(a0+4*(b0-a0)/6);
        int nrec = 0;
        return adapt24(g,a0,b0,acc,eps,f2,f3,nrec,err);

    }

    else if(bcheck)
    {
        double g(double t)
        {
            return f(a+(1-t)/t)*1/(t*t);
        };
        double a0 = 0, b0 = 1;
        double f2 = f(a0+2*(b0-a0)/6);
        double f3 = f(a0+4*(b0-a0)/6);
        int nrec = 0;
        return adapt24(g,a0,b0,acc,eps,f2,f3,nrec,err);

    }

    else
    {

        double f2 = f(a+2*(b-a)/6);
        double f3 = f(a+4*(b-a)/6);
        int nrec = 0;
        return adapt24(f,a,b,acc,eps,f2,f3,nrec,err);
    }
}


double clenshaw_curtis(double f(double x),double a, double b, double acc, double eps, double *err)
{
    double g(double t)
    {
        return f((a+b)/2+(a-b)/2*cos(t) )*sin(t)*(b-a)/2;
    }
    return my_integrator_infinity(g,0,M_PI,acc,eps,err);
}


