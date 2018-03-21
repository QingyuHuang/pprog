#include <stdio.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int motion(double phi, double u[], double dudt[], void *params) {
    double epsilon = *(double *) params;
    dudt[0] =  u[1];
    dudt[1] = 1 + epsilon * pow(u[0], 2) - u[0];
    return GSL_SUCCESS;
}

int main() {
    double epsabs = 1e-6, epsrel = 1e-6;
    double hstart = 1e-3;
    double epsilon = 0.01, dudt = -0.5;
    double phi_max = 90 * M_PI, delta_phi = 0.1;
    double t=0, u[2] = {1, dudt};  // initial conditions
    
    gsl_odeiv2_system sys;
    sys.function = motion;
    sys.jacobian = NULL;
    sys.dimension = 2;
    sys.params = (void *) &epsilon;
    
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
                                                              hstart, epsabs, epsrel);
    
    for (double phi = 0; phi < phi_max; phi += delta_phi) {
        gsl_odeiv2_driver_apply(driver, &t, phi, u);
        printf("%g %g\n", phi, u[0]);
    }
    gsl_odeiv2_driver_free(driver);
    return 0;
}
