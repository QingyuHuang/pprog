#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void rkstep12 (double t, double h, gsl_vector *yt, void f(double t, gsl_vector *y, gsl_vector *dydt), gsl_vector *yth, gsl_vector *err) {
	int n = yt->size;
	
	gsl_vector *k0 = gsl_vector_alloc(n);
	gsl_vector *k12 = gsl_vector_alloc(n);
	gsl_vector *ytt = gsl_vector_alloc(n);

	f(t,yt,k0);
	for (int i=0; i<n; i++) {
		gsl_vector_set(ytt,i,gsl_vector_get(yt,i) + h/2*gsl_vector_get(k0,i));
	}
	f(t+h/2,ytt,k12);
	for (int i=0; i<n; i++) {
		gsl_vector_set(yth,i,gsl_vector_get(yt,i) + h*gsl_vector_get(k12,i));
		gsl_vector_set(err,i,(gsl_vector_get(k0,i) - gsl_vector_get(k12,i))*h/2);
	}

	gsl_vector_free(k0);
	gsl_vector_free(k12);
	gsl_vector_free(ytt);
}

void driver (double* tt, double b, double* hh, gsl_vector *yt, double acc, double eps, void stepper(double t, double h, gsl_vector *y, void f(double t, gsl_vector *y, gsl_vector *dydt), gsl_vector *yth, gsl_vector *err), void f(double t, gsl_vector *y, gsl_vector *dydt)) {
	int n = yt->size;
	double t=*tt;
	double h=*hh;
	double a = t;

	gsl_vector *yth = gsl_vector_alloc(n);
	gsl_vector *dy = gsl_vector_alloc(n);

	while(t<b) {
		if (t+h>b) {
			h = b-t;
		}
		stepper(t,h,yt,f,yth,dy);
		double err = gsl_blas_dnrm2(dy);
		double normy = gsl_blas_dnrm2(yth);
		double tol = (normy*eps+acc)*sqrt(h/(b-a));
		if (err<tol) {
			t += h;
			for (int i=0; i<n; i++) {
				gsl_vector_set(yt,i,gsl_vector_get(yth,i));
			}
		}
		if (err>0) {
			h *= pow(tol/err,0.25)*0.95;
		}else {
			h *= 2;
		}
	}

	gsl_vector_free(yth);
	gsl_vector_free(dy);
}

int driver_path (gsl_vector *tpath, gsl_matrix *ypath, double b, double h, gsl_vector *y, double acc, double eps, int max, void stepper(double t, double h, gsl_vector *y, void f(double t, gsl_vector *y, gsl_vector *dydt), gsl_vector *yth, gsl_vector *err), void f(double t, gsl_vector *y, gsl_vector *dydt)) {
	int n = y->size;
	int k = 0;
	double a = gsl_vector_get(tpath,0);
	
	gsl_vector *yth = gsl_vector_alloc(n);
	gsl_vector *dy = gsl_vector_alloc(n);

	while(gsl_vector_get(tpath,k)<b) {
		double t = gsl_vector_get(tpath,k);
		if (t+h>b) {
			h = b-t;
		}
		stepper(t,h,y,f,yth,dy);
		double err = gsl_blas_dnrm2(dy);
		double normy = gsl_blas_dnrm2(yth);
		double tol = (normy*eps + acc)*sqrt(h/(b-a));
		if (err<tol) {
			k++;
			if (k>max-1) {
				return -k;
			}
			gsl_vector_set(tpath,k,t+h);
			for (int i=0; i<n; i++) {
				gsl_vector_set(y,i,gsl_vector_get(yth,i));
				gsl_matrix_set(ypath,k,i,gsl_vector_get(y,i));
			}
		}
		if (err>0) {
			h *= pow(tol/err,0.25)*0.95;
		}else {
			h *= 2;
		}
	}
	
	gsl_vector_free(yth);
	gsl_vector_free(dy);
	return k;
}