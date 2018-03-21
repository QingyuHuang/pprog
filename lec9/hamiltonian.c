#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double normInteg (double x, void *params) {
    double a = *(double*) params;
    return exp(-a*pow(x,2));
}

double normIntegFun (double a) {
    gsl_function f;
    f.function = &normInteg;
    f.params = (void*) &a;
    
    int limit = 100;
    double acc = 1e-6, eps = 1e-6, result, err;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(limit);
    
    int status = gsl_integration_qagi(&f, acc, eps, limit, workspace, &result, &err);
    
    gsl_integration_workspace_free(workspace);
    
    if(status!=GSL_SUCCESS) {
        return NAN;
    }
    return result;
}

double hamInteg (double x, void *params) {
    double a = *(double*) params;
    return (-pow(a,2)*pow(x,2)/2 + a/2 + pow(x,2)/2) * exp(-a*pow(x,2));
}

double hamIntegFun (double a) {
    gsl_function f;
    f.function = &hamInteg;
    f.params = (void*) &a;
    
    int limit = 100;
    double acc = 1e-6, eps = 1e-6, result, err;
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(limit);
    
    int status = gsl_integration_qagi(&f, acc, eps, limit, workspace, &result, &err);
    
    gsl_integration_workspace_free(workspace);
    
    if(status!=GSL_SUCCESS) {
        return NAN;
    }
    return result;
}

int main (void) {
    for (double a=0.1; a<=2; a+=0.05) {
        double E = hamIntegFun(a)/normIntegFun(a);
        printf("%g \t %g\n", a, E);
    }
    
}
