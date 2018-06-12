#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void gkl_bidiag(gsl_matrix* A,gsl_matrix* U,gsl_matrix* B,gsl_matrix* V);
void print_matrix(const gsl_matrix* mat);
int check_matrix_equal(const gsl_matrix * a, const gsl_matrix * b);


void gkl_bidiag(gsl_matrix* A,gsl_matrix* U,gsl_matrix* B,gsl_matrix* V)
{
    int m=A->size1;
    int n=A->size2;
    int p=(m<n)? m:n;

    gsl_vector* start=gsl_vector_alloc(m);
    gsl_vector* alpha=gsl_vector_alloc(p);
    gsl_vector* beta=gsl_vector_alloc(p);
    gsl_vector* tmp=gsl_vector_alloc(n);
    gsl_vector* Vj=gsl_vector_alloc(n);
    gsl_vector* Uj=gsl_vector_alloc(m);

    gsl_vector_set(start,0,1);
    gsl_vector_set(beta,0,gsl_blas_dnrm2(start));

    gsl_vector_scale(start,1.0/gsl_blas_dnrm2(start));
    gsl_matrix_set_col(U,0,start);
    gsl_blas_dgemv (CblasTrans, 1.0, A, start, 0, tmp);

    gsl_vector_set(alpha,0,gsl_blas_dnrm2(tmp));
    gsl_vector_scale(tmp,1.0/gsl_blas_dnrm2(tmp));
    gsl_matrix_set_col(V,0,tmp);



    for(int j=0; j<p-1; j++)
    {
        gsl_matrix_get_col(Vj, V, j);
        gsl_matrix_get_col(Uj, U, j);

        gsl_blas_dgemv (CblasNoTrans, 1.0, A, Vj, (-1)*gsl_vector_get(alpha,j), Uj);

        gsl_vector_set(beta,j+1,gsl_blas_dnrm2(Uj));
        gsl_vector_scale(Uj,1.0/gsl_blas_dnrm2(Uj));

        gsl_matrix_set_col(U,j+1,Uj);



        gsl_blas_dgemv (CblasTrans, 1.0, A, Uj, (-1)*gsl_vector_get(beta,j+1), Vj);
        gsl_vector_set(alpha,j+1,gsl_blas_dnrm2(Vj));
        gsl_vector_scale(Vj,1.0/gsl_blas_dnrm2(Vj));

        gsl_matrix_set_col(V,j+1,Vj);

    }

    for(int i=0; i<p; i++)
    {
        gsl_matrix_set(B,i,i,gsl_vector_get(alpha,i));
    }
    for(int i=1; i<p; i++)
    {
        gsl_matrix_set(B,i,i-1,gsl_vector_get(beta,i));
    }

    gsl_vector_free(start);
    gsl_vector_free(alpha);
    gsl_vector_free(beta);
    gsl_vector_free(tmp);
    gsl_vector_free(Vj);
    gsl_vector_free(Uj);

}


void print_matrix(const gsl_matrix* mat)
{
    int n=mat->size1;
    int m=mat->size2;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            printf("%10.4f ",gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
    printf("\n");
}



int check_matrix_equal(const gsl_matrix * a, const gsl_matrix * b)
{
    int n1=a->size1;
    int m1=a->size2;
    int n2=b->size1;
    int m2=b->size2;
    double eps=1e-10;

    if(n1!=n2||m1!=m2)
        return 0;

    for(int i=0; i<n1; i++)
    {
        for(int j=0; j<m1; j++)
        {
            if(fabs(gsl_matrix_get(a,i,j)-gsl_matrix_get(b,i,j))>eps)
                return 0;
        }
    }
    return 1;
}


