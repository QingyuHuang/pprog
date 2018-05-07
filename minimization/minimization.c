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


int newton(double (*f)(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), void hessian(gsl_vector* x, gsl_matrix* H), gsl_vector* xstart, double eps)
{

    int n = xstart->size;
    int nsteps = 0;

    gsl_matrix* H = gsl_matrix_alloc(n, n);
    gsl_matrix* R = gsl_matrix_alloc(n, n);
    gsl_vector* grad_f = gsl_vector_alloc(n);
    gsl_vector* delta = gsl_vector_alloc(n);



    gradient(xstart, grad_f);
    double norm_f;
    gsl_blas_ddot (grad_f, grad_f, &norm_f);
    norm_f = sqrt(norm_f);

    do
    {

        hessian(xstart, H);


        qr_gs_decomp(H, R);
        qr_gs_solve(H, R, grad_f, delta);

        double lambda = 1.;

        double fx = f(xstart);

        gsl_blas_daxpy(-lambda, delta, xstart);

        double tmp;
        gsl_blas_ddot (delta, grad_f, &tmp);
        while (f(xstart) > (fx + 1e-2*lambda*tmp))
        {
            lambda /= 2.;
            gsl_blas_daxpy(lambda, delta, xstart);

        }

        gradient(xstart, grad_f);

        gsl_blas_ddot (grad_f, grad_f, &norm_f);
        norm_f = sqrt(norm_f);

        nsteps++;

    }
    while(norm_f > eps && nsteps < 1e6);

    if(nsteps >= 1e6)
    {
        return -1;
    }


    gsl_matrix_free(H);
    gsl_matrix_free(R);
    gsl_vector_free(delta);
    gsl_vector_free(grad_f);


    return nsteps;
}


int gradient(double (*f)(gsl_vector* x), gsl_vector* grad_f, gsl_vector* x, double dx)
{

    if(grad_f->size != x->size)
    {
        return -1;
    }

    double fx = f(x);
    for (int i = 0; i < grad_f->size; ++i)
    {
        double x_i = gsl_vector_get(x, i);
        gsl_vector_set(x, i, x_i + dx);
        double fx_dx = f(x);
        gsl_vector_set(x, i, x_i);
        double grad_f_i = (fx_dx - fx) / dx;
        gsl_vector_set(grad_f, i, grad_f_i);
    }

    return 0;
}


int hessian(gsl_matrix* H, gsl_matrix* Hnew, gsl_vector* y, gsl_vector* s, double lambda)
{

    if(H->size1 != H->size2 || Hnew->size1 != Hnew->size2 || Hnew->size1 != H->size2 || H->size1 != y->size || y->size != s->size)
    {
        return -1;
    }


    gsl_blas_dscal(-lambda,s);

    double y_trans_H_inv_s = 0.;
    for (int i = 0; i < y->size; ++i)
    {
        double y_i = gsl_vector_get(y, i);
        for (int j = 0; j < y->size; ++j)
        {
            double s_j = gsl_vector_get(y, j);
            double H_ij = gsl_matrix_get(H, i, j);
            y_trans_H_inv_s += y_i*s_j*H_ij;
        }
    }
    if (fabs(y_trans_H_inv_s) < 1e-3)
    {
        gsl_matrix_set_identity(H);
        return 0;
    }
    for (int i = 0; i < H->size1; ++i)
    {
        for (int j = 0; j < H->size1; ++j)
        {
            double s_i = gsl_vector_get(s, i);
            double result = 0.;
            for (int k = 0; k < H->size1; ++k)
            {
                double H_y_i = 0.;
                for (int q = 0; q < H->size1; ++q)
                {
                    double H_iq = gsl_matrix_get(H, i, q);
                    double y_q = gsl_vector_get(y, q);
                    H_y_i += H_iq*y_q;
                }
                double s_k = gsl_vector_get(s, k);
                double H_kj = gsl_matrix_get(H, k, j);
                result += s_i * s_k * H_kj - H_y_i * s_k * H_kj;
            }
            gsl_matrix_set(Hnew, i, j, result/y_trans_H_inv_s);
        }
    }

    gsl_matrix_memcpy(H,Hnew);

    return 0;
}


int quasi_newton(double (*f)(gsl_vector* x), gsl_vector* xstart, double dx, double eps)
{


    int n = xstart->size;
    int nsteps = 0;

    gsl_matrix* H = gsl_matrix_alloc(n, n);
    gsl_matrix* extra_mem = gsl_matrix_alloc(n, n);
    gsl_vector* grad_f = gsl_vector_alloc(n);
    gsl_vector* grad_f_dx = gsl_vector_alloc(n);
    gsl_vector* delta = gsl_vector_alloc(n);


    gsl_matrix_set_identity(H);

    int error = gradient(f, grad_f, xstart, dx);
    if (error == -1)
        return -1;


    do
    {
        error = gsl_blas_dgemv(CblasNoTrans, 1.0, H, grad_f, 0.0, delta);
        if (error == -1)
            return -1;

        double lambda = 1.;

        double fx = f(xstart);
        error =gsl_blas_daxpy(-lambda, delta, xstart);
        if (error == -1)
            return -1;

        double tmp;
        gsl_blas_ddot (delta, grad_f, &tmp);
        while (f(xstart) > (fx + 1e-2*lambda*tmp) && lambda > dx/10)
        {
            lambda /= 2.;

            error =  gsl_blas_daxpy(lambda, delta, xstart);
            if (error == -1)
                return -1;

        }

        error = gradient(f, grad_f_dx, xstart, dx);
        if (error == -1)
            return -1;

        double norm_f_dx;

        gsl_blas_ddot (grad_f_dx, grad_f_dx, &norm_f_dx);
        norm_f_dx = sqrt(norm_f_dx);
        nsteps++;


        if (norm_f_dx < eps)
        {
            break;
        }
        else
        {
            gsl_blas_daxpy(-1, grad_f_dx, grad_f);
            gsl_blas_dscal(-2,grad_f);
            hessian(H, extra_mem, grad_f, delta, lambda);
            gsl_vector_memcpy(grad_f, grad_f_dx);
        }


    }
    while(nsteps < 1e6);

    if(nsteps >= 1e6)
    {
        return -1;
    }


    gsl_matrix_free(H);
    gsl_matrix_free(extra_mem);
    gsl_vector_free(delta);
    gsl_vector_free(grad_f);
    gsl_vector_free(grad_f_dx);

    return nsteps;
}


