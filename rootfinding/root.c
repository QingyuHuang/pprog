#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void qr_givens(gsl_matrix* A);
void qr_solve_givens(gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart, double dx, double epsilon);
void newton_with_jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart, double epsilon);
void newton_with_Jacobian_refined(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart, double dx, double epsilon);

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart, double dx, double epsilon)
{
    int n,i,j;
    double scale;
    gsl_matrix* J;
    gsl_vector* Dx, *y, *fy, *fx, *df;

    n=xstart->size;
    J=gsl_matrix_alloc(n,n);
    Dx=gsl_vector_alloc(n);
    y=gsl_vector_alloc(n);
    fy=gsl_vector_alloc(n);
    fx=gsl_vector_alloc(n);
    df=gsl_vector_alloc(n);

    do
    {
        f(xstart,fx);
        for(j=0; j<n; j++)
        {
            gsl_vector_set(xstart,j,gsl_vector_get(xstart,j)+dx);
            f(xstart,df);
            gsl_vector_sub(df,fx);
            for(i=0; i<n; i++)
            {
                gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx);
            }
            gsl_vector_set(xstart,j,gsl_vector_get(xstart,j)-dx);
        }

        qr_givens(J);
        gsl_vector_scale(fx,-1.0);
        qr_solve_givens(J,fx,Dx);
        gsl_vector_scale(fx,-1.0);

        scale=2;
        do
        {
            scale/=2;
            gsl_vector_memcpy(y,Dx);
            gsl_vector_scale(y,scale);
            gsl_vector_add(y,xstart);
            f(y,fy);
        }
        while(gsl_blas_dnrm2(fy)>(1-scale/2)*gsl_blas_dnrm2(fx) && scale>0.02);

        gsl_vector_memcpy(xstart,y);
        gsl_vector_memcpy(fx,fy);

    }
    while(gsl_blas_dnrm2(Dx)>dx && gsl_blas_dnrm2(fx)>epsilon);

    gsl_matrix_free(J);
    gsl_vector_free(Dx);
    gsl_vector_free(y);
    gsl_vector_free(fy);
    gsl_vector_free(fx);
    gsl_vector_free(df);
}


void qr_givens(gsl_matrix* A)
{
    int m,n,i,j,k;
    double angle,x,y,dx,dy;

    m=A->size2;
    n=A->size1;

    for(i = 0; i < m; i++)
    {
        for(j = i+1; j< n; j++)
        {
            angle=atan2(gsl_matrix_get(A,j,i),gsl_matrix_get(A,i,i));

            for (k=i; k < m; k++)
            {
                x = gsl_matrix_get(A,i,k);
                y = gsl_matrix_get(A,j,k);
                dx=x*cos(angle)+y*sin(angle);
                dy=-x*sin(angle)+y*cos(angle);
                gsl_matrix_set(A,i,k,dx);
                gsl_matrix_set(A,j,k,dy);
            }

            gsl_matrix_set(A,j,i,angle);
        }
    }
}


void qr_solve_givens(gsl_matrix* QR, gsl_vector* b, gsl_vector* x)
{
    int m,n,i,j;
    double angle,x0,y0,dx,dy,sum;
    m=QR->size2;
    n=QR->size1;
    for(int i=0; i<m; i++)
    {
        for(int j=i+1; j<n; j++)
        {
            angle=gsl_matrix_get(QR,j,i);
            x0=gsl_vector_get(b,i);
            y0=gsl_vector_get(b,j);
            dx=x0*cos(angle)+y0*sin(angle);
            dy=-x0*sin(angle)+y0*cos(angle);
            gsl_vector_set(b,i,dx);
            gsl_vector_set(b,j,dy);
        }
    }

    for (i=m-1; i>=0; i--)
    {
        sum=0;
        for(j=i+1; j<m; j++)
        {
            sum+=gsl_matrix_get(QR,i,j)*gsl_vector_get(x,j);
        }
        gsl_vector_set(x,i,(gsl_vector_get(b,i)-sum)/gsl_matrix_get(QR,i,i));
    }
}




void newton_with_jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart, double epsilon)
{
    int n;
    double scale;
    gsl_matrix* J;
    gsl_vector* Dx, *y, *fy, *fx, *df;

    n=xstart->size;
    J=gsl_matrix_alloc(n,n);
    Dx=gsl_vector_alloc(n);
    y=gsl_vector_alloc(n);
    fy=gsl_vector_alloc(n);
    fx=gsl_vector_alloc(n);
    df=gsl_vector_alloc(n);

    do
    {
        f(xstart,fx,J);
        qr_givens(J);
        gsl_vector_scale(fx,-1.0);
        qr_solve_givens(J,fx,Dx);
        gsl_vector_scale(fx,-1.0);

        scale=2;
        do
        {
            scale/=2;
            gsl_vector_memcpy(y,Dx);
            gsl_vector_scale(y,scale);
            gsl_vector_add(y,xstart);
            f(y,fy,J);
        }
        while(gsl_blas_dnrm2(fy)>(1-scale/2)*gsl_blas_dnrm2(fx) && scale>0.02);
        gsl_vector_memcpy(xstart,y);
        gsl_vector_memcpy(fx,fy);
    }
    while(gsl_blas_dnrm2(fx)>epsilon);

    gsl_matrix_free(J);
    gsl_vector_free(Dx);
    gsl_vector_free(y);
    gsl_vector_free(fy);
    gsl_vector_free(fx);
    gsl_vector_free(df);
}


void newton_with_Jacobian_refined(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart, double dx, double epsilon)
{
    int n;
    double scale;
    gsl_matrix* J;
    gsl_vector* Dx, *y, *fy, *fx, *df;

    n=xstart->size;
    J=gsl_matrix_alloc(n,n);
    Dx=gsl_vector_alloc(n);
    y=gsl_vector_alloc(n);
    fy=gsl_vector_alloc(n);
    fx=gsl_vector_alloc(n);
    df=gsl_vector_alloc(n);

    do
    {
        f(xstart,fx,J);
        qr_givens(J);
        gsl_vector_scale(fx,-1.0);
        qr_solve_givens(J,fx,Dx);
        gsl_vector_scale(fx,-1.0);
        scale=0.99;
        do
        {
            gsl_vector_memcpy(y,Dx);
            gsl_vector_scale(y,scale);
            gsl_vector_add(y,xstart);
            f(y,fy,J);
            scale=gsl_blas_dnrm2(fx)*gsl_blas_dnrm2(fx)/(2.0*(0.5*gsl_blas_dnrm2(fy)*gsl_blas_dnrm2(fy)-0.5*gsl_blas_dnrm2(fx)*gsl_blas_dnrm2(fx)+gsl_blas_dnrm2(fx)*gsl_blas_dnrm2(fx)*scale)/scale/scale);
        }
        while(gsl_blas_dnrm2(fy)>(1-scale/2)*gsl_blas_dnrm2(fx) && scale>0.02);
        gsl_vector_memcpy(xstart,y);
        gsl_vector_memcpy(fx,fy);
    }
    while(gsl_blas_dnrm2(Dx)>dx && gsl_blas_dnrm2(fx)>epsilon);

    gsl_matrix_free(J);
    gsl_vector_free(Dx);
    gsl_vector_free(y);
    gsl_vector_free(fy);
    gsl_vector_free(fx);
    gsl_vector_free(df);
}
