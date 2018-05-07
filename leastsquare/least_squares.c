#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R)
{
    int n=A->size1;
    int m=A->size2;
    double val,data;

    for(int i=0; i < m; i++)
    {
        val = 0;
        for(int j = 0; j < n; j++)
        {
            val += gsl_matrix_get(A,j,i)*gsl_matrix_get(A,j,i);
        }
        val=sqrt(val);

        gsl_matrix_set(R, i, i, val);

        for(int x = 0; x < n; x++)
        {
            data=gsl_matrix_get(A,x,i)/val;
            gsl_matrix_set(A, x, i, data);
        }

        for(int j = i + 1; j<m; j++)
        {
            val = 0;

            for(int x = 0; x < n; x++)
            {
                val += gsl_matrix_get(A,x,i)*gsl_matrix_get(A,x,j);
            }

            gsl_matrix_set(R, i, j, val);

            for(int x = 0; x < n; x++)
            {
                data=gsl_matrix_get(A,x,j)-gsl_matrix_get(A,x,i)*val;
                gsl_matrix_set(A,x,j,data);
            }
        }

    }
}

void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R,const gsl_vector* b, gsl_vector* x)
{
    int n=Q->size1;
    int m=R->size2;
    double val,data;

    for(int i = 0; i < m; i++)
    {
        val = 0;
        for(int j = 0; j < n; j++)
        {
            val += gsl_matrix_get(Q,j,i)*gsl_vector_get(b,j);
        }
        gsl_vector_set(x,i,val);
    }

    for(int i = m-1; i >= 0; i--)
    {
        for(int j = m-1; j >= 0; j--)
        {
            if(i == j)
            {
                data=gsl_vector_get(x,j)/gsl_matrix_get(R,i,j);
            }
            else
            {
                data=gsl_vector_get(x,j)-gsl_vector_get(x,i)*gsl_matrix_get(R,j,i);
            }
            gsl_vector_set(x,j,data);
        }
    }
}


void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B)
{
    int m = Q->size2;
    double data;
    for(int k = 0; k<m; k++)
    {
        for(int i = 0; i < m; i++)
        {
            data=gsl_matrix_get(Q,k,i);
            gsl_matrix_set(B,i,k,data);
        }

        for(int i = m-1; i >= 0; i--)
        {
            for(int j = m-1; j >= 0; j--)
            {
                if(i == j)
                {
                    data=gsl_matrix_get(B,j,k)/gsl_matrix_get(R,i,j);
                }
                else
                {
                    data=gsl_matrix_get(B,j,k)-gsl_matrix_get(B,i,k)*gsl_matrix_get(R,j,i);
                }
                gsl_matrix_set(B,j,k,data);
            }
        }

    }

}

double fitfunctions(int i, double x)
{
    switch(i)
    {
    case 0:
        return log(x);
        break;
    case 1:
        return 1.0;
        break;
    case 2:
        return x;
        break;
    default:
    {
        fprintf(stderr,"funs: wrong i:%d",i);
        return NAN;
    }
    }
}

double fit(double x, gsl_vector* c,int m)
{
    double s=0;
    for(int k=0; k<m; k++)
        s+=gsl_vector_get(c,k)*fitfunctions(k,x);
    return s;
}

double fit_plus(int i, double x, gsl_vector* c, int m, gsl_vector* dc)
{
    return fit(x,c,m)+gsl_vector_get(dc,i)*fitfunctions(i,x);
}

double fit_minus(int i, double x, gsl_vector* c, int m, gsl_vector* dc)
{
    return fit(x,c,m)-gsl_vector_get(dc,i)*fitfunctions(i,x);
}

void lsfit(int m, double f(int i,double x), gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c, gsl_matrix* S)
{
    int n = x->size;

    gsl_matrix *A    = gsl_matrix_alloc(n,m);
    gsl_vector *b    = gsl_vector_alloc(n);
    gsl_matrix *R    = gsl_matrix_alloc(m,m);
    gsl_matrix* invR = gsl_matrix_alloc(m,m);
    gsl_matrix *I    = gsl_matrix_alloc(m,m);

    for(int i=0; i<n; i++)
    {
        double xi  = gsl_vector_get(x,i);
        double yi  = gsl_vector_get(y,i);
        double dyi = gsl_vector_get(dy,i);
        assert(dyi>0);
        gsl_vector_set(b,i,yi/dyi);
        for(int k=0; k<m; k++)
            gsl_matrix_set(A,i,k,f(k,xi)/dyi);
    }
    qr_gs_decomp(A,R);
    qr_gs_solve(A,R,b,c);

    gsl_matrix_set_identity(I);
    qr_gs_inverse(I,R,invR);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,invR,invR,0,S);

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_matrix_free(R);
    gsl_matrix_free(invR);
    gsl_matrix_free(I);
}
