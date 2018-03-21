#include <stdio.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int logistic_function(double t, double y[], double dydt[]) {
    dydt[0] =  y[0] * (1 - y[0]);
    // An ODE system should return GSL_SUCCES if the calculation was completed succesfully.
    return GSL_SUCCESS;
}

double my_function(double x) {
    gsl_odeiv2_system sys;
    sys.function = logistic_function;
    sys.jacobian = NULL;
    sys.dimension = 1;
    sys.params = NULL;
    
    double epsabs = 1e-6, epsrel = 1e-6;
    double hstart = copysign(0.1, x); // hstart determines the initial step size and copysign(a,b) is magnitude of a and sign of b.
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                                              hstart, epsabs, epsrel);
    
    double t=0, y[1] = {0.5};  // initial conditions
    gsl_odeiv2_driver_apply(driver, &t, x, y);
    
    gsl_odeiv2_driver_free(driver);
    return y[0];
}

double math_function(double t) {
    return 1.0 / (1 + exp(-t));
}

int main() {
    double start=0.0, end=3.0, dt = 0.1;
    printf("t\t y(t)\t exact\n");
    for (double t=start; t <= end + dt; t += dt) {
        printf("%g %g %g\n", t, my_function(t), math_function(t));
    }
    return 0;
}
