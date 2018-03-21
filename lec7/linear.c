#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

int main(){
    gsl_vector *x,*x2,*v,*v2,*a,*t;
    gsl_matrix *M,*M2;
    gsl_permutation *p;
    int *sig;
    size_t n=3;
    size_t i;
    //-----------------------------------------------------------------------------------------------------------------------------
    printf("Allocating memory.\n");
    sig=malloc(n*sizeof(int));
    a=gsl_vector_calloc(n+1);
    x=gsl_vector_calloc(n);
    v=gsl_vector_calloc(n);
    t=gsl_vector_calloc(n);
    M=gsl_matrix_calloc(n,n);
    x2=gsl_vector_calloc(n);
    v2=gsl_vector_calloc(n);
    M2=gsl_matrix_calloc(n,n);
    p=gsl_permutation_calloc(n);
    //-----------------------------------------------------------------------------------------------------------------------------
    printf("Initializing vector v and matrix M.\n");
    for(i=0;i<3;i++){
        gsl_vector_fscanf(stdin,a);
        gsl_vector_set(v,i,gsl_vector_get(a,3));
        gsl_matrix_set(M,i,0,gsl_vector_get(a,0));
        gsl_matrix_set(M,i,1,gsl_vector_get(a,1));
        gsl_matrix_set(M,i,2,gsl_vector_get(a,2));
        
    }
    
    printf("v=\n");
    gsl_vector_fprintf(stdout,v,"%g");
    
    printf("Copying vector v into v2 and matrix M to M2 for calculations.\n");
    gsl_vector_memcpy(v2,v);
    gsl_matrix_memcpy(M2,M);

    gsl_linalg_HH_solve(M2,v2,x);
    printf("The value found for x was:\n");
    gsl_vector_fprintf(stdout,x,"%g");
    printf("Verifying the result. Be the judge of this method yourself:\n");
    gsl_blas_dgemv(CblasNoTrans,1.0,M,x,0.0,x2);
    printf("Printing Ax.\n");
    gsl_vector_fprintf(stdout,x2,"%.10g");
    gsl_vector_memcpy(v2,v);
    gsl_matrix_memcpy(M2,M);
    //-----------------------------------------------------------------------------------------------------------------------------
    printf("Solving by the LU Decomposition:\n");
    gsl_linalg_LU_decomp(M2,p,sig);
    gsl_linalg_LU_solve(M2,p,v2,x);
    printf("The value found for x was:\n");
    gsl_vector_fprintf(stdout,x,"%g");
    printf("Verifying the result. Be the judge of this method yourself:\n");
    gsl_blas_dgemv(CblasNoTrans,1.0,M,x,0.0,x2);
    printf("Printing Ax.\n");
    gsl_vector_fprintf(stdout,x2,"%.10g");
    gsl_vector_memcpy(v2,v);
    gsl_matrix_memcpy(M2,M);
    //-----------------------------------------------------------------------------------------------------------------------------
    printf("Solving by the QR Decomposition:\n");
    gsl_linalg_QR_decomp(M2,t);
    gsl_linalg_QR_solve(M2,t,v2,x);
    printf("The value found for x was:\n");
    gsl_vector_fprintf(stdout,x,"%g");
    printf("Verifying the result. Be the judge of this method yourself:\n");
    gsl_blas_dgemv(CblasNoTrans,1.0,M,x,0.0,x2);
    printf("Printing Ax.\n");
    gsl_vector_fprintf(stdout,x2,"%.10g");
    gsl_vector_memcpy(v2,v);
    gsl_matrix_memcpy(M2,M);
    gsl_vector_free(a);
    gsl_vector_free(x);
    gsl_vector_free(v);
    gsl_vector_free(t);
    gsl_matrix_free(M);
    gsl_vector_free(x2);
    gsl_vector_free(v2);
    gsl_matrix_free(M2);
    gsl_permutation_free(p);
    free(sig);
    return 0;
}
