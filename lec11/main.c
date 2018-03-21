#include <gsl/gsl_multimin.h>
#include <stdio.h>
#include <math.h>

double rosen_func(const gsl_vector *x, void *params){
    double x1 = gsl_vector_get(x,0);
    double x2 = gsl_vector_get(x,1);
    double f = pow(1-x1,2) + 100*pow(x2-x1*x1,2);
    return f;
}

struct experimental_data {int n; double *t, *y, *e;};

double f_to_min(const gsl_vector *x, void *params) {
    double A = gsl_vector_get(x,0);
    double T = gsl_vector_get(x,1);
    double B = gsl_vector_get(x,2);
    struct experimental_data *p = (struct experimental_data*) params;
    int n = p->n;
    double *t = p->t;
    double *y = p->y;
    double *e = p->e;
#define f(t) A*exp(-t/T)+B
    double sum = 0;
    for(int i = 0; i<n; i++) sum += pow( ( f(t[i])-y[i] )/e[i],2);
    return sum;
}

int main(){
    FILE* data = fopen("data.txt","w");
    //*********** Problem 1************//
    {
        size_t iter = 0;
        int status;
        double size;
        
        gsl_vector *x = gsl_vector_alloc(2);
        gsl_vector_set(x,0,2.0);
        gsl_vector_set(x,1,5.0);
        
        int dim = 2;
        gsl_multimin_function f;
        f.f = &rosen_func;
        f.n = dim;
        f.params = NULL;
        
        gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,dim);
        
        gsl_vector *step = gsl_vector_alloc(2);
        gsl_vector_set(step,0,0.01);
        gsl_vector_set(step,1,0.01);
        gsl_multimin_fminimizer_set(s,&f,x,step);
        
        do {
            iter++;
            status = gsl_multimin_fminimizer_iterate (s);
            
            if(status) break;
            
            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-4);
            
            fprintf (data,"%5ld\t%.5g\t%.5g\t%10.5g\n", iter,
                     gsl_vector_get (s->x, 0),
                     gsl_vector_get (s->x, 1),
                     s->fval);
            
            
        } while(status == GSL_CONTINUE && iter < 100);
        
        
    }
    //*********** Problem 2************//
    {
        fprintf(data,"\n\n");
        size_t iter = 0;
        int status;
        double size;
        
        gsl_vector *x = gsl_vector_alloc(3);
        gsl_vector_set(x,0,2.0);
        gsl_vector_set(x,1,5.0);
        gsl_vector_set(x,1,7.0);
        
        double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
        double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
        double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
        int n = sizeof(t)/sizeof(t[0]);
        
        struct experimental_data params = {n,t,y,e};
        
        int dim = 3;
        gsl_multimin_function f;
        f.f = &f_to_min;
        f.n = dim;
        f.params = &params;
        
        //gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,dim);
        
        gsl_vector *step = gsl_vector_alloc(3);
        gsl_vector_set(step,0,0.01);
        gsl_vector_set(step,1,0.01);
        gsl_vector_set(step,2,0.01);
        gsl_multimin_fminimizer_set(s,&f,x,step);
        
        do {
            iter++;
            status = gsl_multimin_fminimizer_iterate (s);
            
            if(status) break;
            
            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-3);
            
            fprintf (data,"%5ld\t%.5g\t%.5g\t%.5g\t%10.5g\n", iter,
                     gsl_vector_get (s->x, 0),
                     gsl_vector_get (s->x, 1),
                     gsl_vector_get (s->x, 2),
                     s->fval);
            
            
        } while(status == GSL_CONTINUE && iter < 1000);
        
        double A_min = gsl_vector_get (s->x, 0);
        double T_min = gsl_vector_get (s->x, 1);
        double B_min = gsl_vector_get (s->x, 2);
        printf("The minimum is at (%g,%g,%g)", A_min
               , T_min
               , B_min);
        
        fprintf(data,"\n\n");
        fprintf(data, "t\ty_exp\terror\n");
        for(int i = 0; i<n; i++) {
            fprintf(data, "%g\t%g\t%g\n", t[i], y[i], e[i]);
        }
        fprintf(data, "\n\n");
        fprintf(data, "t\ty_fit\n");
        for(double t = 0; t<10; t+=0.1){
            fprintf(data, "%g\t%g\n", t, A_min*exp(-t/T_min)+B_min);
        }
        gsl_multimin_fminimizer_free(s);
        gsl_vector_free(x);
        
    }
    return 0;
}
