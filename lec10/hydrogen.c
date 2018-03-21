#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>


int print_state (size_t iter, gsl_multiroot_fsolver *s) {
    fprintf (stderr, "%3lu\t %.8f\t%.8e\n",
             iter,
             gsl_vector_get (s->x, 0),
             gsl_vector_get (s->f, 0));
    
    return 0;
}

int diff_equation (double r, const double y[], double yprime[], void *params) {
    double e = *(double*) params;
    yprime[0] = y[1];
    yprime[1] = 2*(-1/r-e)*y[0];
    
    return GSL_SUCCESS;
}

double hydrogen_function (double e, double r) {
    assert (r>=0);
    
    const double rmin = 1e-3;
    if (r<rmin) {
        return r-pow(r, 2);
    }
    
    gsl_odeiv2_system solver;
    solver.function = diff_equation;
    solver.jacobian = NULL;
    solver.dimension = 2;
    solver.params = (void*) &e;
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new (&solver, gsl_odeiv2_step_rkf45, 1e-3, 1e-6, 1e-6);
    
    double x = rmin;
    double y[2] = {x-pow(x, 2), 1-2*x};
    int status = gsl_odeiv2_driver_apply (driver, &x, r, y);
    if(status!=GSL_SUCCESS) {
        fprintf (stderr, "odeiv2 error: %d\n", status);
    }
    
    gsl_odeiv2_driver_free (driver);
    
    return y[0];
}

int shooting_method (const gsl_vector *x, void *params, gsl_vector *f) {
    const double e = gsl_vector_get (x, 0);
    double rmax = *(double*) params;
    
    gsl_vector_set (f, 0, hydrogen_function (e, rmax));
    
    return GSL_SUCCESS;
}

int main () {
    const size_t dim =1;
    const double rmax = 12;
    
    gsl_multiroot_function f = {&shooting_method, dim, (void*) &rmax};
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, dim);
    gsl_vector *e = gsl_vector_alloc (dim);
    gsl_vector_set (e, 0, -1);
    gsl_multiroot_fsolver_set (s, &f, e);
    
    int iter = 0, status;
    const double acc = 1e-6;
    
    fprintf (stderr, "iter\te\tf(rmax)\n");
    print_state (iter, s);
    
    do {
        iter++;
        int flag = gsl_multiroot_fsolver_iterate (s);
        if (flag!=0) break;
        status = gsl_multiroot_test_residual (s->f, acc);
        print_state (iter, s);
    } while (status==GSL_CONTINUE && iter<999);
    
    double Emin = gsl_vector_get (s->x, 0);
    double radius, delta = 0.1;
    
    printf ("radius\tcalculated\tanalytical\n");
    for (radius=0.0; radius<=12.0; radius+=delta) {
        double analytical = radius*exp(-radius);
        printf("%g\t%g\t%g\n", radius, hydrogen_function (Emin, radius), analytical);
    }
    
    
    gsl_vector_free (e);
    gsl_multiroot_fsolver_free (s);
    
    return 0;
}
