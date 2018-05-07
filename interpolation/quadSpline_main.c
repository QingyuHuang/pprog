#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double *x, double *y) /* allocates and builds the quadratic spline */
{
	qspline *qs = (qspline*)malloc(sizeof(qspline));
	qs->n = n;
	qs->x = (double*)malloc(n*sizeof(double));
	qs->y = (double*)malloc(n*sizeof(double));
	qs->b = (double*)malloc((n-1)*sizeof(double));
	qs->c = (double*)malloc((n-1)*sizeof(double));
	for( int i=0; i<n; i++)
	{
		qs->x[i] = x[i];
		qs->y[i] = y[i];
	}
	double p[n-1], h[n-1];
	for( int i=0; i<n-1; i++)
	{
		h[i] = x[i+1] - x[i];
		assert( h[i] > 0 );
		p[i] = (y[i+1] - y[i])/h[i];
	}
	qs->c[0] = 0.0;
	for( int i=0; i<n-2; i++)
		qs->c[i+1] = (p[i+1] - p[i] - qs->c[i]*h[i])/h[i+1];
	qs->c[n-2] /= 2.0;
	for ( int i=n-3; i>=0; i--)
		qs->c[i] = (p[i+1] - p[i] - qs->c[i+1]*h[i+1])/h[i];
	for( int i=0; i<n-1; i++)
		qs->b[i] = p[i] - qs->c[i]*h[i];
	return qs;
}

double qspline_evaluate( qspline *s, double z)        /* evaluates the prebuilt spline at point z */
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
			return s->y[m]+xl*( s->b[m] + xl*s->c[m]);
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}

double qspline_derivative( qspline *s, double z) /* evaluates the derivative of the prebuilt spline at point z */
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
			return s->b[m] + 2*xl*s->c[m] ;
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}

double qspline_integral( qspline *s, double z)  /* evaluates the integral of the prebuilt spline from x[0] to z */
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
				sum += h*s->y[i] + h*h*(s->b[i]/2 + h*s->c[i]/3);
			}
			sum += xl*s->y[i] + xl*xl*(s->b[i]/2 + xl*s->c[i]/3);
			return sum;
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}
void qspline_free( qspline *s) /* free memory allocated in qspline_alloc */
{
	free( s->x);
	free( s->y);
	free( s->b);
	free( s->c);
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
	qspline *s = qspline_alloc(n, x, y);
	double z;
	printf(" test of qspline_evaluate : \n");
	for( z=1; z<=5; z+=0.1)
		printf("%f %f \n",z, qspline_evaluate(s, z) );
	printf(" test of qspline_derivative : \n");
	for( z=1; z<=5; z+=0.1)
		printf("%f %f \n",z, qspline_derivative(s, z) );
	printf(" test of qspline_integral : \n");
	for( z=1; z<=5; z+=0.1)
		printf("%f %f \n",z, qspline_integral(s, z) );
	qspline_free(s);
	return 0;
}
