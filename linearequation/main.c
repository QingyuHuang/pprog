#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R,const gsl_vector* b,gsl_vector* x);

void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);

void qr_gr(gsl_matrix* A);

void print_matrix(gsl_matrix* mat);
void print_vector(gsl_vector* vec);

int check_upper_triangular(gsl_matrix* mat);
int check_matrix_equal(const gsl_matrix * a, const gsl_matrix * b);
int check_vector_equal(const gsl_vector * a, const gsl_vector * b);

int main()
{
    // Part A.
    printf("\nSolving linear equations using QR-decomposition by modified Gram-Schmidt orthogonalization:\n");
    printf("\nA. 1. Implement \"void qr_gs_decomp(matrix* A, matrix* R)\":\n\n");
    //generate a random tall (n>m) matrix A (of a modest size);
    int n=5,m=4;
    gsl_matrix *A = gsl_matrix_alloc(n,m);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            gsl_matrix_set (A, i, j, i+2*j);
        }
    }
    printf("generate a (n>m) matrix A:\n");
    print_matrix(A);

    //factorize it into QR;
    gsl_matrix *Q = gsl_matrix_alloc(n,m);
    gsl_matrix *R = gsl_matrix_alloc(m,m);
    gsl_matrix_memcpy (Q, A);
    qr_gs_decomp(Q, R);

    printf("matrix Q:\n");
    print_matrix(Q);

    printf("matrix R:\n");
    print_matrix(R);

    // check that R is upper triangular
    if (check_upper_triangular(R))
        printf("R is upper triangular.\n\n");
    else
        printf("R is not upper triangular.\n\n");

    // check that QR=A;
    gsl_matrix *QR = gsl_matrix_alloc(n,m);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, Q, R,
                    0.0, QR);
    if (check_matrix_equal(QR,A))
        printf("QR=A is True.\n\n");
    else
        printf("QR=A is not True.\n\n");

    gsl_matrix_free(A);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_matrix_free(QR);


    printf("\nA. 2. Implement \"void qr_gs_solve(const matrix* Q, const matrix* R,const vector* b,vector* x)\":\n\n");
    // generate a random square matrix A (of a modest size);
    n=2,m=2;
    A = gsl_matrix_alloc(n,m);
    int arr[2][2]= {{4,7},{2,6}};
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            gsl_matrix_set (A, i, j, arr[i][j]);
        }
    }
    printf("generate square matrix A:\n");
    print_matrix(A);
    //generate a random vector b (of the same size);
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector_set_all(b,1.0);
    printf("generate vector b:\n");
    print_vector(b);

    //factorize A into QR;
    Q = gsl_matrix_alloc(n,m);
    R = gsl_matrix_alloc(m,m);
    gsl_matrix_memcpy (Q, A);
    qr_gs_decomp(Q, R);

    printf("matrix Q:\n");
    print_matrix(Q);

    printf("matrix R:\n");
    print_matrix(R);

    //solve QRx=b;
    gsl_vector *x = gsl_vector_alloc(m);
    qr_gs_solve(Q,R,b,x);
    printf("vector x:\n");
    print_vector(x);

    //check that Ax=b;
    gsl_vector *Ax = gsl_vector_alloc(n);
    gsl_blas_dgemv (CblasNoTrans,
                    1.0, A, x,
                    0.0, Ax);
    printf("vector Ax:\n");
    print_vector(Ax);

    if (check_vector_equal(Ax,b))
        printf("Ax=b is True.\n\n");
    else
        printf("Ax=b is not True.\n\n");



    // Part B.
    printf("B. Matrix inverse by Gram-Schmidt QR factorization\n\n");
    printf("Implement \"void qr_gs_inverse(const matrix* Q, const matrix* R, matrix* B)\"\n\n");
    // generate a random square matrix A (of a modest size);
    printf("generate square matrix A:\n");
    print_matrix(A);
    // factorize A into QR;
    printf("matrix Q:\n");
    print_matrix(Q);

    printf("matrix R:\n");
    print_matrix(R);
    // calculate the inverse B;
    gsl_matrix* B = gsl_matrix_alloc(n,m);
    qr_gs_inverse(Q, R, B);

    printf("Inverse matrix B:\n");
    print_matrix(B);
    // check that AB=I, where I is the identity matrix;
    gsl_matrix* AB = gsl_matrix_alloc(n,m);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, A, B,
                    0.0, AB);

    gsl_matrix* I = gsl_matrix_alloc(n,m);
    gsl_matrix_set_identity(I);

    if (check_matrix_equal(AB,I))
        printf("AB=I is True.\n\n");
    else
        printf("AB=I is not True.\n\n");

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(AB);
    gsl_matrix_free(I);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_vector_free(b);
    gsl_vector_free(Ax);
    gsl_vector_free(x);

    // Part C.
    printf("C. QR-decomposition by Givens rotations\n\n");

    //generate a random tall (n>m) matrix A (of a modest size);
    n=5,m=4;
    A = gsl_matrix_alloc(n,m);
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            gsl_matrix_set (A, i, j, i+2*j);
        }
    }
    printf("generate a (n>m) matrix A:\n");
    print_matrix(A);

    //factorize it into QR;
    qr_gr(A);

    printf("matrix QR:\n");
    print_matrix(A);
    gsl_matrix_free(A);

    return 0;
}

