#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


//makes linear spline interpolation from a table {x[i], y[i]} at a given point z

typedef struct {int n; double *x, *y, *b, *c, *d;} cubicSpline;

cubicSpline * cubicSpline_alloc(int n, double *x, double *y) /* allocates and builds the quadratic spline */
{
	cubicSpline *cs = (cubicSpline*)malloc(sizeof(cubicSpline));
	cs->n = n;
	cs->x = (double*)malloc(n*sizeof(double));
	cs->y = (double*)malloc(n*sizeof(double));
	cs->b = (double*)malloc((n-1)*sizeof(double));
	cs->c = (double*)malloc((n-1)*sizeof(double));
	cs->d = (double*)malloc((n-1)*sizeof(double));
	for( int i=0; i<n; i++)
	{
		cs->x[i] = x[i];
		cs->y[i] = y[i];
	}
	double p[n-1], h[n-1];
	double u[n-2], v[n-2], z[n];
	for( int i=0; i<n-1; i++)
	{
		h[i] = x[i+1] - x[i];
		assert( h[i] > 0 );
		p[i] = (y[i+1] - y[i])/h[i];
	}
	u[0] = 2*(h[0] + h[1]);
	v[0] = 6*(p[1] - p[0]);
	for( int i=2; i<n-1; i++)
	{
		u[i-1] = 2*(h[i] + h[i-1]) - h[i-1]*h[i-1]/u[i-2];
		v[i-1] = 6*(p[i] - p[i-1]) - h[i-1]*v[i-2]/u[i-2];
	}
	z[n-1] = 0.0;
	for( int i=n-2; i>0 ; i--)
		z[i] = (v[i-1] - z[i+1]*h[i])/u[i-1];
	z[0] = 0.0;
	for(int i=0; i<n-1; i++)
	{
		cs->b[i] = -h[i]/6*z[i+1] - h[i]/3*z[i] + p[i]/h[i];
		cs->c[i] = z[i]/2;
		cs->d[i] = (z[i+1] - z[i])/6.0/h[i];
	}
	return cs;
}

double cubicSpline_evaluate( cubicSpline *s, double z)        /* evaluates the prebuilt spline at point z */
{
	int L,R,m; // position for z location
	L = 0;
	R = s->n -1;
	double xl;
	while ( L <= R )
	{
		m = (L+R)/2;
		if ( z < s->x[m] )
			R = m;
		else if ( s->x[m+1] > z )
		{
			xl = z - s->x[m];
			return s->y[m]+xl*( s->b[m] + xl*(s->c[m] + xl*s->d[m]));
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}

double cubicSpline_derivative( cubicSpline *s, double z) /* evaluates the derivative of the prebuilt spline at point z */
{
	int L,R,m; // position for z location
	L = 0;
	R = s->n -1;
	double xl;
	while ( L <= R )
	{
		m = (L+R)/2;
		if ( z < s->x[m] )
			R = m;
		else if ( s->x[m+1] > z )
		{
			xl = z - s->x[m];
			return s->b[m] + xl*(2*s->c[m] + xl*3*s->d[m]) ;
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}

double cubicSpline_integral( cubicSpline *s, double z)  /* evaluates the integral of the prebuilt spline from x[0] to z */
{
	int L,R,m; // position for z location
	L = 0;
	R = s->n-1;
	int i;
	double sum = 0.0;
	double xl,h;
	while ( L <= R )
	{
		m = (L+R)/2;
		if ( z < s->x[m] )
			R = m;
		else if ( s->x[m+1] > z )
		{
			xl = (z - s->x[m]);
			for(i=0; i<m; i++)
			{
				h = s->x[i+1] - s->x[i];
				sum += h*(s->y[i] + h*(s->b[i]/2 + h*(s->c[i]/3 + h*s->d[i])));
			}
			sum += xl*(s->y[m] + xl*(s->b[m]/2 + xl*(s->c[m]/3 + xl*s->d[m])));
			return sum;
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}
void cubicSpline_free( cubicSpline *s) /* free memory allocated in cubicSpline_alloc */
{
	free( s->x);
	free( s->y);
	free( s->b);
	free( s->c);
	free( s->d);
	free( s);
}

int main()
{
	int n;
	double *x, *y;
	int i;
	n = 5;
	x = malloc(n*sizeof(double));
	y = malloc(n*sizeof(double));
	for( i=0; i<n; i++)
	{
		x[i] = i+1;
		//y[i] = 1;
		//y[i] = i+1;
		y[i] = (i+1)*(i+1);
	}
	cubicSpline *s = cubicSpline_alloc(n, x, y);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *gsl_s = gsl_spline_alloc(gsl_interp_cspline, n);
	gsl_spline_init(gsl_s, x,y, n);

	double z;
	printf(" test of cubicSpline_evaluate : \n");
	for( z=1; z<=5; z+=0.1)
		printf("%f %f %f \n",z, cubicSpline_evaluate(s, z), gsl_spline_eval(gsl_s,z, acc) );
	printf(" test of cubicSpline_derivative : \n");
	for( z=1; z<=5; z+=0.1)
		printf("%f %f \n",z, cubicSpline_derivative(s, z) );
	printf(" test of cubicSpline_integral : \n");
	for( z=1; z<=5; z+=0.1)
		printf("%f %f \n",z, cubicSpline_integral(s, z) );
	cubicSpline_free(s);
	gsl_spline_free( gsl_s );
	gsl_interp_accel_free( acc);
	return 0;
}
