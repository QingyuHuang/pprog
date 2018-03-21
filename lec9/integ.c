#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double integral (double x, void *params) {
    return log(x)/sqrt(x);
}

double integFun () {
    gsl_function f;
    f.function = integral;
    f.params = NULL;
    
    int limit = 100, key = 4;
    double a = 0, b = 1, acc = 1e-6, eps = 1e-6, result, err;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(limit);
    
    int status = gsl_integration_qag(&f, a, b, acc, eps, limit, key, workspace, &result, &err);
    
    gsl_integration_workspace_free(workspace);
    
    if(status != GSL_SUCCESS) {
        return NAN;
    }
    return result;
}

int main (void) {
    double result = integFun ();
    printf("Thus the result of the integral is: %g\n", result);
}

