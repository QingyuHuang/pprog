#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector* x), int N, double* result, double* error)
{
    int dim = a->size;

    double V = 1;
    for(int i=0; i<dim; i++)
    {
        V*=gsl_vector_get(b,i)-gsl_vector_get(a,i);
    }
    double sum =0, sum2 = 0, fx;
    gsl_vector* x = gsl_vector_alloc(dim);


    for(int j=0; j<N; j++)
    {
        int n = x->size;


        for(int i =0; i<n; i++)
        {
            double a_i = gsl_vector_get(a,i),
                   b_i = gsl_vector_get(b,i);

            gsl_vector_set(x,i,a_i+((double)rand()/RAND_MAX)*(b_i-a_i));
        }
        fx = f(x);
        sum+=fx;
        sum2+=fx*fx;
    }

    double avr = sum/N, var = sum2/N-avr*avr;
    *result = avr*V;
    *error = sqrt(var/N)*V;

    gsl_vector_free(x);
}

