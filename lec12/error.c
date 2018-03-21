#include <stdio.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int error_function(double t, double y[], double dydt[]) {
    dydt[0] = (2 / sqrt(M_PI)) * exp(- pow(t, 2));
    return GSL_SUCCESS;
}

double my_function(double x) {
    gsl_odeiv2_system sys;
    sys.function = error_function;
    sys.jacobian = NULL;
    sys.dimension = 1;
    sys.params = NULL;
    
    double epsabs = 1e-6, epsrel = 1e-6;
    double hstart = copysign(0.1, x);
    double t=0, y[1] = {0};
    
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                                              hstart, epsabs, epsrel);
    
    gsl_odeiv2_driver_apply(driver, &t, x, y);
    
    gsl_odeiv2_driver_free(driver);
    return y[0];
}

int main(int argc, char **argv) {
    double a = strtod(argv[1], NULL), b = strtod(argv[2], NULL), dt = strtod(argv[3], NULL);
    printf("t\t y(t)\n");
    for (double t = a; t <= b + dt; t += dt) {
        printf("%g\t %g\n", t, my_function(t));
    }
    return 0;
}

