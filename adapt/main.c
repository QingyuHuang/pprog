#include<math.h>
#include<assert.h>
#include<stdio.h>
double adapt24(double f(double x),double a,double b,double acc,double eps, double f2, double f3,int nrec,double *err);
double my_integrator(double f(double x),double a, double b, double acc, double eps, double *err);
double my_integrator_infinity(double f(double x),double a, double b, double acc, double eps,double *err);
double clenshaw_curtis(double f(double x),double a, double b, double acc, double eps, double *err);

int main() //uses gcc nested functions
{
    double err;
    int calls=0;
    // Part A.
    printf("===================================================================\n");
    printf("Part A:\n\n");
    double a=0,b=1,acc=1e-4,eps=1e-4;

    double f1(double x)
    {
        calls++;
        return sqrt(x);
    }; //nested function
    calls=0;
    double Q=my_integrator(f1,a,b,acc,eps,&err);
    double exact=2./3;
    printf("Recursive adaptive integrator: integrating sqrt(x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));


    double f2(double x)
    {
        calls++;
        return 1/sqrt(x);
    }; //nested function
    calls=0;
    Q=my_integrator(f2,a,b,acc,eps,&err);
    exact=2;
    printf("\nRecursive adaptive integrator: integrating 1/sqrt(x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    double f3(double x)
    {
        calls++;
        return log(x)/sqrt(x);
    }; //nested function
    calls=0;
    Q=my_integrator(f3,a,b,acc,eps,&err);
    exact=-4;
    printf("\nRecursive adaptive integrator: integrating log(x)/sqrt(x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    double f4(double x)
    {
        calls++;
        return 4.0*sqrt(1-(1-x)*(1-x));
    }; //nested function
    calls=0;
    Q=my_integrator(f4,a,b,acc,eps,&err);
    exact=M_PI;
    printf("\nRecursive adaptive integrator: integrating 4.0*sqrt(1-(1-x)*(1-x)) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    // Part B.
    printf("\n\n===================================================================\n");
    printf("Part B:\n");
    a=-INFINITY,b=INFINITY,acc=1e-4,eps=1e-4;
    calls=0;
    double f5(double x)
    {
        calls++;
        return 1/(1+x*x);
    }; //nested function
    Q=my_integrator_infinity(f5,a,b,acc,eps,&err);
    exact=M_PI;
    printf("\nInfinite limits: integrating 1/(1+x*x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));


    a=2.0/M_PI,b=INFINITY,acc=1e-4,eps=1e-4;
    calls=0;
    double f6(double x)
    {
        calls++;
        return 1.0/(x*x)*sin(1.0/x);
    }; //nested function
    Q=my_integrator_infinity(f6,a,b,acc,eps,&err);
    exact=1;
    printf("\nInfinite limits: integrating 1.0/(x*x)*sin(1.0/x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    a=0,b=INFINITY,acc=1e-4,eps=1e-4;
    calls=0;
    double f7(double x)
    {
        calls++;
        return x*log(x)/((1+x*x)*(1+x*x));
    }; //nested function
    Q=my_integrator_infinity(f7,a,b,acc,eps,&err);
    exact=0;
    printf("\nInfinite limits: integrating x*log(x)/((1+x*x)*(1+x*x)) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    a=01,b=2,acc=1e-4,eps=1e-4;
    calls=0;
    double f8(double x)
    {
        calls++;
        return x/sqrt(x-1);
    }; //nested function
    Q=my_integrator_infinity(f8,a,b,acc,eps,&err);
    exact=8.0/3;
    printf("\nInfinite limits: integrating x/sqrt(x-1) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    // Part C.
    printf("\n\n===================================================================\n");
    printf("Part C:\n");
    a=0,b=1,acc=1e-4,eps=1e-4;
    calls=0;
    Q=clenshaw_curtis(f1,a,b,acc,eps,&err);
    exact=2./3;
    printf("\nclenshaw_curtis: integrating sqrt(x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    calls=0;
    Q=clenshaw_curtis(f2,a,b,acc,eps,&err);
    exact=2;
    printf("\nclenshaw_curtis: integrating 1/sqrt(x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    calls=0;
    Q=clenshaw_curtis(f3,a,b,acc,eps,&err);
    exact=-4;
    printf("\nclenshaw_curtis: integrating log(x)/sqrt(x) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    calls=0;
    Q=clenshaw_curtis(f4,a,b,acc,eps,&err);
    exact=M_PI;
    printf("\nclenshaw_curtis: integrating 4.0*sqrt(1-(1-x)*(1-x)) from %g to %g\n",a,b);
    printf("acc=%g eps=%g\n",acc,eps);
    printf("              Q = %g\n",Q);
    printf("          exact = %g\n",exact);
    printf("          calls = %d\n",calls);
    printf("estimated error = %g\n",acc+fabs(Q)*eps);
    printf("   actual error = %g\n",fabs(Q-exact));

    return 0 ;
}
