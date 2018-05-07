#include <stdio.h>
#include <stdlib.h>

//makes linear spline interpolation from a table {x[i], y[i]} at a given point z
double linterp(int n, double *x, double *y, double z)
{
	int L,R,m; // position for z location
	L = 0;
	R = n-1;
	double lambda;
	while ( L <= R )
	{
		m = (L+R)/2;
		if ( z < x[m] )
			R = m;
		else if ( x[m+1] > z )
		{
			lambda = (z - x[m])/(x[m+1] - x[m]);
			return lambda*y[m+1] + (1-lambda)*y[m];
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}

//calculates the intergral of the linear spline from the point x[0] to the given point z
double linterp_integ(int n, double *x, double *y, double z)
{
	int L,R,m; // position for z location
	L = 0;
	R = n-1;
	int i;
	double sum = 0.0;
	double lambda;
	while ( L <= R )
	{
		m = (L+R)/2;
		if ( z < x[m] )
			R = m;
		else if ( x[m+1] > z )
		{
			lambda = (z - x[m])/(x[m+1] - x[m]);
			for(i=0; i<m; i++)
				sum += (x[i+1] - x[i])*(y[i] + y[i+1])/2;
			sum += (z-x[m])*(y[m] + lambda*y[m+1] + (1-lambda)*y[m] )/2;
			return sum;
		}
		else
			L = m+1;
	}
	printf(" z in not in range of x \n");
	return 0;
}

int main()
{
	int n;
	double *x, *y;
	FILE *fp;
	int i;
	fp = fopen("data.inp","r");
	fscanf(fp, "%d", &n);
	x = malloc(n*sizeof(double));
	y = malloc(n*sizeof(double));
	for( i=0; i<n; i++)
	{
		fscanf(fp,"%lg %lg", &x[i], &y[i]);
	}
	fclose(fp);

	fp = fopen("data.out","w");
	double z;
	fprintf(fp,"# x , y_interp , y_ans, int_y_interp, int_y_ans \n");
	for( z=x[0]; z<=x[n-1]; z+=0.1)
		fprintf(fp,"%f %f %f  %f %f\n",z, 
							linterp(n, x, y, z), 
							2*z,
							linterp_integ(n, x, y, z), 
							z*z);
	fclose(fp);
	free(x);
	free(y);
	return 0;
}
