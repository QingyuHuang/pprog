#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>


void plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector* x), int N,double* result, double* error);
double f2(gsl_vector* x);
double f3(gsl_vector* x);

int main()
{

    int dim = 3, N = 1e6;
    gsl_vector* a = gsl_vector_alloc(dim);
    gsl_vector* b = gsl_vector_alloc(dim);

    double result, error;

    printf("\n\nPart A:\n\n");

    gsl_vector_set_zero(a);
    gsl_vector_set(b,0,1);
    gsl_vector_set(b,1,2*M_PI);
    gsl_vector_set(b,2,M_PI);

    plainmc(a,b,f2,N,&result,&error);

    printf("Integrating r^2 with r from 0 to 1, phi from 0 to 2*pi and theta from 0 to pi with %i points, resulting in %lg with an error of %lg.\n\n",N,result,error);

    // Question 3
    gsl_vector_set_zero(a);
    gsl_vector_set(b,0,M_PI);
    gsl_vector_set(b,1,M_PI);
    gsl_vector_set(b,2,M_PI);

    plainmc(a,b,f3,N,&result,&error);

    printf("Integrating 1/pi^3*1/(1-cos(x)cos(y)cos(z) with x from 0 to pi, with y from 0 to pi and with z from 0 to pi with %i points, resulting in %lg with an error of %lg.\n\n",N,result,error);

    // Part B
    printf("\nPart B\nSee Graph\n\n");

    FILE* data = fopen("data.txt","w+");
    for(int i = 5; i<500; i+=2)
    {
        gsl_vector_set_zero(a);
        gsl_vector_set(b,0,1);
        gsl_vector_set(b,1,2*M_PI);
        gsl_vector_set(b,2,M_PI);

        plainmc(a,b,f2,i,&result,&error);

        fprintf(data,"%i\t%lg\n",i,error);
    }
    fclose(data);


    gsl_vector_free(a);
    gsl_vector_free(b);
    return 0;
}


double f2(gsl_vector* x)
{
    double x_0 = gsl_vector_get(x,0);
    return x_0*x_0;
}

double f3(gsl_vector* x)
{
    double x_0 = gsl_vector_get(x,0),
           y_0 = gsl_vector_get(x,1),
           z_0 = gsl_vector_get(x,2);
    return 1.0/(M_PI*M_PI*M_PI)*1.0/(1.0-cos(x_0)*cos(y_0)*cos(z_0));
}
