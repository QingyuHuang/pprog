#include "ann.h"

ann* ann_alloc(int n, double(*act_fun)(double))
{
    ann* a = malloc(sizeof(ann));
    a->n = n;
    a->f=act_fun;
    a->data = gsl_vector_alloc(3*n);
    return a;
}

void ann_free(ann* network)
{
    gsl_vector_free(network->data);
    free(network);
}

double ann_feed_forward(ann* network, double x)
{
    int n = network->n;
    double sum=0;
    for(int i=0; i<n; i++)
    {
        double a_i = gsl_vector_get(network->data,0*network->n+i);
        double b_i = gsl_vector_get(network->data,1*network->n+i);
        double w_i = gsl_vector_get(network->data,2*network->n+i);

        double arg = (x+a_i)/b_i;
        double result = network->f(arg)*w_i;
        sum+=result;
    }
    return sum;
}

void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist,double eps,double step)
{
    if (xlist->size != ylist->size)
        return;

    int N=xlist->size;

    double func(gsl_vector* p)
    {
        gsl_vector_memcpy(network->data,p);
        double sum=0;
        for(int k=0; k<N; k++)
        {
            double x_k = gsl_vector_get(xlist,k);
            double y_k = gsl_vector_get(ylist,k);
            double F_k = ann_feed_forward(network,x_k);
            sum+=(F_k-y_k)*(F_k-y_k);
        }
        return sum/N;
    }

    gsl_vector* p = gsl_vector_alloc(network->data->size);
    gsl_vector_memcpy(p,network->data);
    newton_broyden(func,p,step,eps);
    gsl_vector_memcpy(network->data,p);
    gsl_vector_free(p);
}

void grad(double f(gsl_vector* x), gsl_vector* x, gsl_vector* df, double dx)
{
    double dfi, fi_1, fi_2, xi;

    for(int i=0; i<x->size; i++)
    {
        fi_1 = f(x);
        xi = gsl_vector_get(x,i);
        gsl_vector_set(x,i,xi+dx);
        fi_2=f(x);
        gsl_vector_set(x,i,xi-dx);
        dfi=(fi_2-fi_1)/dx;
        gsl_vector_set(df,i,dfi);
    }
}


int newton_broyden(double f(gsl_vector* x), gsl_vector* x, double dx, double eps)
{
    int n=x->size;
    gsl_vector* x_p = gsl_vector_alloc(n);
    gsl_vector* delta_x = gsl_vector_alloc(n);
    gsl_vector* df = gsl_vector_alloc(n);
    gsl_vector* df_p = gsl_vector_alloc(n);
    gsl_vector* y = gsl_vector_alloc(n);
    gsl_vector* H_inv_y = gsl_vector_alloc(n);
    gsl_vector* delta_H_inv_y = gsl_vector_alloc(n);
    gsl_matrix* d_H_inv = gsl_matrix_alloc(n,n);
    gsl_matrix* H_inv = gsl_matrix_alloc(n,n);

    gsl_matrix_set_identity(H_inv);
    double delta_x_norm, df_norm, fx, lambda, DxTdf, fx_p, H_inv_update, alpha = 0.5,update_num, update_denum;
    int n_step =0;
    do
    {
        n_step++;
        fx = f(x);
        grad(f,x,df,dx);
        gsl_vector_memcpy(x_p,x);
        gsl_vector_add(x_p,delta_x);
        fx_p = f(x_p);
        grad(f,x_p,df_p,dx);

        if(n_step>1)
        {
            gsl_vector_memcpy(y,df_p);
            gsl_vector_sub(y,df);
            gsl_blas_dgemv(CblasNoTrans,1.0,H_inv,y,0.0,H_inv_y);
            gsl_vector_memcpy(delta_H_inv_y,delta_x);
            gsl_vector_sub(delta_H_inv_y,H_inv_y);
            gsl_blas_ddot(delta_H_inv_y,delta_x,&update_num);
            gsl_blas_dgemv(CblasNoTrans,1.0,H_inv,delta_x,0.0,H_inv_y);
            gsl_blas_ddot(y,H_inv_y,&update_denum);
            H_inv_update = update_num/update_denum;
            gsl_matrix_memcpy(d_H_inv,H_inv);
            gsl_matrix_scale(d_H_inv,H_inv_update);
            gsl_matrix_add(H_inv,d_H_inv);
        }

        gsl_blas_dgemv(CblasNoTrans,1.0,H_inv,df,0.0,delta_x);
        gsl_vector_scale(delta_x,-1.0);
        lambda=1.0;

        gsl_vector_memcpy(x_p,delta_x);
        gsl_vector_scale(x_p,lambda);
        gsl_vector_add(x_p,x);

        fx_p = f(x_p);
        grad(f,x_p,df_p,dx);

        gsl_blas_ddot(delta_x,df,&DxTdf);
        delta_x_norm = gsl_blas_dnrm2(delta_x);

        do
        {
            lambda /=2.0;
            gsl_vector_memcpy(x_p,delta_x);
            gsl_vector_scale(x_p,lambda);
            gsl_vector_add(x_p,x);
            fx_p=f(x_p);
            grad(f,x_p,df_p,dx);
        }
        while(fx_p>fx+alpha*lambda*DxTdf && lambda > 1.0/64.0);

        if(gsl_blas_dnrm2(df)<gsl_blas_dnrm2(df_p))
        {
            gsl_matrix_set_identity(H_inv);
        }

        gsl_vector_scale(delta_x,lambda);
        gsl_vector_add(x,delta_x);
        fx=f(x);
        grad(f,x,df,dx);
        df_norm = gsl_blas_dnrm2(df);
    }
    while(df_norm>eps && delta_x_norm>dx);

    gsl_vector_free(x_p);
    gsl_vector_free(delta_x);
    gsl_vector_free(df);
    gsl_vector_free(df_p);
    gsl_vector_free(y);
    gsl_vector_free(H_inv_y);
    gsl_vector_free(delta_H_inv_y);
    gsl_matrix_free(d_H_inv);
    gsl_matrix_free(H_inv);
    return n_step;
}
