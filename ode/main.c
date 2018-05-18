#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void rkstep12 (double t, double h, gsl_vector *yt, void f(double t, gsl_vector *y, gsl_vector *dydt), gsl_vector *yth, gsl_vector *err);
void driver (double* tt, double b, double* hh, gsl_vector *yt, double acc, double eps, void stepper(double t, double h, gsl_vector *y, void f(double t, gsl_vector *y, gsl_vector *dydt), gsl_vector *yth, gsl_vector *err), void f(double t, gsl_vector *y, gsl_vector *dydt));
int driver_path (gsl_vector *tpath, gsl_matrix *ypath, double b, double h, gsl_vector *y, double acc, double eps, int max, void stepper(double t, double h, gsl_vector *y, void f(double t, gsl_vector *y, gsl_vector *dydt), gsl_vector *yth, gsl_vector *err), void f(double t, gsl_vector *y, gsl_vector *dydt));
void sine (double t, gsl_vector *y, gsl_vector *dydt);
void print_vector(gsl_vector* vec);

int main () {
	int dim = 2, max=1e7;
	double t =0, b = M_PI, h = 01, acc = 1e-6, eps = 1e-3;

	printf("PART A\n\n");
	gsl_vector *y = gsl_vector_alloc(dim);
	gsl_vector_set(y,0,0);
	gsl_vector_set(y,1,1);

	printf("Startpoint\ta = %g\nEndpoint\tb = %g\n", t, b);
	driver(&t,b,&h,y,acc,eps,rkstep12,sine);
	printf("Function\ty(b) = %g\nDifferential\tdy(b) = %g\n", gsl_vector_get(y,0), gsl_vector_get(y,1));

	printf("\nPART B\n\n");
	gsl_vector *tpath = gsl_vector_alloc(max);
	gsl_matrix *ypath = gsl_matrix_alloc(max,dim);
	gsl_vector_set(y,0,0);
	gsl_vector_set(y,1,1);

	printf("Startpoint\ta = %g\nEndpoint\tb = %g\n", t, b);
	int k = driver_path(tpath,ypath,b,h,y,acc,eps,max,rkstep12,sine);
	if (k<0) {
		printf("Maximum step count reached\n");
	}
	printf("Function\ty(b) = %g\nDifferential\tdy(b) = %g\n", gsl_vector_get(y,0), gsl_vector_get(y,1));

	printf("\n\n");
	for (int i=0; i<k; i++) {
		printf("%g %g %g %g %g\n", gsl_vector_get(tpath,i), gsl_matrix_get(ypath,i,0), gsl_matrix_get(ypath,i,1), sin(gsl_vector_get(tpath,i)), cos(gsl_vector_get(tpath,i)));
	}

	gsl_vector_free(y);
	gsl_vector_free(tpath);
	gsl_matrix_free(ypath);
	return 0;
}

void sine (double t, gsl_vector *y, gsl_vector *dydt) {
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-1*gsl_vector_get(y,0));
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
