#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void gkl_bidiag(gsl_matrix* A,gsl_matrix* U,gsl_matrix* B,gsl_matrix* V);
void print_matrix(const gsl_matrix* mat);
int check_matrix_equal(const gsl_matrix * a, const gsl_matrix * b);

int main()
{
    printf("\n1. Number of rows m equal to number of columns n, m = n:\n");
    int m1=4;
    int n1=4;
    int p1=(m1<n1)? m1:n1;

    int mat1[4][4]=
    {
        {1,3,4,5},
        {1,2,3,4},
        {0,0,1,2},
        {3,3,2,4}
    };


    gsl_matrix* A1=gsl_matrix_alloc(m1,n1);
    gsl_matrix* U1=gsl_matrix_alloc(m1,p1);
    gsl_matrix* B1=gsl_matrix_alloc(p1,p1);
    gsl_matrix* V1=gsl_matrix_alloc(n1,p1);
    gsl_matrix* UA1=gsl_matrix_alloc(p1,n1);
    gsl_matrix* UAV1=gsl_matrix_alloc(p1,p1);

    for(int i=0; i<m1; i++)
    {
        for(int j=0; j<n1; j++)
        {
            gsl_matrix_set(A1,i,j,mat1[i][j]);
        }
    }

    gkl_bidiag(A1,U1,B1,V1);

    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, U1, A1, 0.0, UA1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, UA1, V1, 0.0, UAV1);


    printf("\nOrigin matrix A:\n");
    print_matrix(A1);
    printf("Matrix U:\n");
    print_matrix(U1);
    printf("Bidiagnal Matrix B:\n");
    print_matrix(B1);
    printf("Matrix V:\n");
    print_matrix(V1);

    printf("Matrix U'*A*V:\n");
    print_matrix(UAV1);

    if (check_matrix_equal(B1,UAV1))
        printf("B = U'*A*V is True.\n\n");
    else
        printf("B = U'*A*V is not True.\n\n");





    printf("\n2. Number of rows m equal to number of columns n, m < n:\n");
    int m2=4;
    int n2=5;
    int p2=(m2<n2)? m2:n2;

    int mat2[4][5]=
    {
        {1,3,4,5,6},
        {1,2,3,4,1},
        {0,0,1,2,6},
        {3,3,2,4,3}
    };


    gsl_matrix* A2=gsl_matrix_alloc(m2,n2);
    gsl_matrix* U2=gsl_matrix_alloc(m2,p2);
    gsl_matrix* B2=gsl_matrix_alloc(p2,p2);
    gsl_matrix* V2=gsl_matrix_alloc(n2,p2);
    gsl_matrix* UA2=gsl_matrix_alloc(p2,n2);
    gsl_matrix* UAV2=gsl_matrix_alloc(p2,p2);

    for(int i=0; i<m2; i++)
    {
        for(int j=0; j<n2; j++)
        {
            gsl_matrix_set(A2,i,j,mat2[i][j]);
        }
    }

    gkl_bidiag(A2,U2,B2,V2);

    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, U2, A2, 0.0, UA2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, UA2, V2, 0.0, UAV2);


    printf("\nOrigin matrix A:\n");
    print_matrix(A2);
    printf("Matrix U:\n");
    print_matrix(U2);
    printf("Bidiagnal Matrix B:\n");
    print_matrix(B2);
    printf("Matrix V:\n");
    print_matrix(V2);

    printf("Matrix U'*A*V:\n");
    print_matrix(UAV2);

    if (check_matrix_equal(B2,UAV2))
        printf("B = U'*A*V is True.\n\n");
    else
        printf("B = U'*A*V is not True.\n\n");


    printf("\n3. Number of rows m equal to number of columns n, m > n:\n");
    int m3=5;
    int n3=4;
    int p3=(m3<n3)? m3:n3;

    int mat3[5][4]=
    {
        {1,1,0,3},
        {3,2,0,3},
        {4,3,1,2},
        {5,4,2,4},
        {6,1,6,3}
    };


    gsl_matrix* A3=gsl_matrix_alloc(m3,n3);
    gsl_matrix* U3=gsl_matrix_alloc(m3,p3);
    gsl_matrix* B3=gsl_matrix_alloc(p3,p3);
    gsl_matrix* V3=gsl_matrix_alloc(n3,p3);
    gsl_matrix* UA3=gsl_matrix_alloc(p3,n3);
    gsl_matrix* UAV3=gsl_matrix_alloc(p3,p3);

    for(int i=0; i<m3; i++)
    {
        for(int j=0; j<n3; j++)
        {
            gsl_matrix_set(A3,i,j,mat3[i][j]);
        }
    }

    gkl_bidiag(A3,U3,B3,V3);

    gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, U3, A3, 0.0, UA3);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, UA3, V3, 0.0, UAV3);


    printf("\nOrigin matrix A:\n");
    print_matrix(A3);
    printf("Matrix U:\n");
    print_matrix(U3);
    printf("Bidiagnal Matrix B:\n");
    print_matrix(B3);
    printf("Matrix V:\n");
    print_matrix(V3);

    printf("Matrix U'*A*V:\n");
    print_matrix(UAV3);

    if (check_matrix_equal(B3,UAV3))
        printf("B = U'*A*V is True.\n\n");
    else
        printf("B = U'*A*V is not True.\n\n");



    gsl_matrix_free(A1);
    gsl_matrix_free(U1);
    gsl_matrix_free(B1);
    gsl_matrix_free(V1);
    gsl_matrix_free(UA1);
    gsl_matrix_free(UAV1);

    gsl_matrix_free(A2);
    gsl_matrix_free(U2);
    gsl_matrix_free(B2);
    gsl_matrix_free(V2);
    gsl_matrix_free(UA2);
    gsl_matrix_free(UAV2);

    gsl_matrix_free(A3);
    gsl_matrix_free(U3);
    gsl_matrix_free(B3);
    gsl_matrix_free(V3);
    gsl_matrix_free(UA3);
    gsl_matrix_free(UAV3);
    return 0;
}

