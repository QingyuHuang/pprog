#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_vector.h>


int print_state (size_t iter, gsl_multiroot_fsolver *s) {
    printf ("iter=%3lu\tx=%.3f\ty =%.3f\t dxf=%.3e \t dyf=%.3f\n",
            iter,
            gsl_vector_get (s->x, 0),
            gsl_vector_get (s->x, 1),
            gsl_vector_get (s->f, 0),
            gsl_vector_get(s->f,1));
    
    return 0;
}

int rosenbrock (const gsl_vector *x, void *params, gsl_vector *f) {
    const double x1 = gsl_vector_get(x,0);
    const double x2 = gsl_vector_get(x,1);
    
    double rosenbrock1 = 2*(200*x1*(x1*x1-x2)+x1-1);
    double rosenbrock2 = 200*(x2-x1*x1);
    
    gsl_vector_set (f, 0, rosenbrock1);
    gsl_vector_set (f, 1, rosenbrock2);
    
    return GSL_SUCCESS;
}

int main () {
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    int test;
    int status;
    int iter = 0;
    
    const size_t dim = 2;
    
    gsl_multiroot_function f = {&rosenbrock, dim, NULL};
    
    double x_init[2] = {0.5, 0.6};
    
    gsl_vector *x = gsl_vector_alloc(dim);
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[0]);
    
    T = gsl_multiroot_fsolver_hybrids;
    
    s = gsl_multiroot_fsolver_alloc (T, dim);
    
    gsl_multiroot_fsolver_set(s, &f, x);
    
    print_state (iter, s);
    
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        print_state (iter, s);
        if (status != 0) break;
        
        test = gsl_multiroot_test_residual (s->f, 1e-6);
    } while(test==GSL_CONTINUE && iter<999);
    
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    
    return 0;
}
