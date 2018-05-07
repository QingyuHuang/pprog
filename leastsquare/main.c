#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R,const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);
double fitfunctions(int i, double x);

void lsfit(int m, double f(int i,double x),gsl_vector* x, gsl_vector* y, gsl_vector* dy,gsl_vector* c, gsl_matrix* S);

double fit(double x, gsl_vector* c,int m);
double fit_plus(int i, double x,gsl_vector* c,int m, gsl_vector* dc);
double fit_minus(int i, double x,gsl_vector* c,int m, gsl_vector* dc);

void print_matrix(gsl_matrix* mat);
void print_vector(gsl_vector* vec);

int main()
{
    double x[] = {0.1,1.33,2.55,3.78,5,6.22,7.45,8.68,9.9};
    double y[] = {-15.3,0.32,2.45,2.75,2.27,1.35,0.157,-1.23,-2.75};
    double dy[] = {1.04,0.594,0.983,0.998,1.11,0.398,0.535,0.968,0.478};

    int n=sizeof(x)/sizeof(x[0]);

    for(int i=0; i<n; i++)
        printf("%g %g %g\n",x[i],y[i],dy[i]);
    printf("\n\n");

    gsl_vector* vx = gsl_vector_alloc(n);
    gsl_vector* vy = gsl_vector_alloc(n);
    gsl_vector* vdy = gsl_vector_alloc(n);
    for(int i=0; i<n; i++)
    {
        gsl_vector_set(vx,i,x[i]);
        gsl_vector_set(vy,i,y[i]);
        gsl_vector_set(vdy,i,dy[i]);
    }

    int m=3;

    gsl_vector* c = gsl_vector_alloc(m);
    gsl_matrix* S = gsl_matrix_alloc(m,m);
    lsfit(m,fitfunctions,vx,vy,vdy,c,S);

    gsl_vector* dc = gsl_vector_alloc(m);
    for(int k=0; k<m; k++)
    {
        double skk=gsl_matrix_get(S,k,k);
        gsl_vector_set(dc,k,sqrt(skk));
    }



    double z,dz=(x[n-1]-x[0])/90;
    for(int i=0; i<m; i++)
    {
        z=x[0]-dz/2;
        do
        {
            printf("%g %g %g %g\n",z,fit(z,c,m),fit_plus(i,z,c,m,dc),fit_minus(i,z,c,m,dc));
            z+=dz;
        }
        while(z<x[n-1]+dz);
        printf("\n\n");
    }

    return 0;
};



void print_matrix(gsl_matrix* mat)
{
    int n,m,i,j;
    n=mat->size1;
    m=mat->size2;
    for(i=0; i<n; i++)
    {
        for(j=0; j<m; j++)
        {
            printf("%10.4f ",gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(gsl_vector* vec)
{
    int n,i;
    n=vec->size;
    for(i=0; i<n; i++)
    {
        printf("%10.4f\n",gsl_vector_get(vec, i));
    }
    printf("\n");
}
