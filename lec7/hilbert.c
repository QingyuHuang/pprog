#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

int main(){
    gsl_vector *eval;
    gsl_vector *test;
    gsl_vector_view eveci;
    gsl_matrix *evec;
    gsl_matrix *H,*H2;
    gsl_matrix *Id,*Id2;
    gsl_permutation *p;
    gsl_eigen_symmv_workspace *w;
    double det;
    size_t i,j,n=4;
    int sig;
    printf("Allocating memory.\n");
    //sig=malloc(sizeof(int));
    Id=gsl_matrix_calloc(n,n);
    Id2=gsl_matrix_calloc(n,n);
    H=gsl_matrix_calloc(n,n);
    H2=gsl_matrix_calloc(n,n);
    test=gsl_vector_calloc(n);
    eval=gsl_vector_calloc(n);
    evec=gsl_matrix_calloc(n,n);
    p=gsl_permutation_calloc(n);
    w=gsl_eigen_symmv_alloc(n);
    printf("Creating Hilbert matrix. And identity matrix.\n");
    for(i=0;i<n;i++){
        gsl_matrix_set(Id,i,i,1.0);
        for(j=0;j<n;j++){
            gsl_matrix_set(H,i,j,1.0/(1+i+j));
        }
    }
    gsl_matrix_memcpy(Id2,Id);
    gsl_matrix_memcpy(H2,H);
    printf("Solving The eigenvalues/vectors problem.\n");
    gsl_eigen_symmv(H2,eval,evec,w);
    printf("\t 1. Displaying eigenvalues and eigenvectors:\n");
    for(i=0;i<n;i++){
        eveci=gsl_matrix_column(evec,i);
        printf("\tEigenvalue: %g\n",gsl_vector_get(eval,i));
        printf("\tCorresponding eigenvector:\n");
        gsl_vector_fprintf(stdout,&eveci.vector,"%g");
    }
    printf("Verifying the result.\n");
    printf("To evaluate the determinant we will use the functions defined in the linear algebra header\n");
    printf("that depend on the LU decomposition of the matrix whose determinant we want.\n");
    printf("\t1. Verifying determinants.\n");
    gsl_matrix_memcpy(H2,H);
    for(i=0;i<n;i++){
        printf("\tVerifying eiganvalue %g.\n",gsl_vector_get(eval,i));
        printf("\t\t1. Creating H-evalI.\n");
        gsl_matrix_scale(Id2,gsl_vector_get(eval,i));
        gsl_matrix_sub(H2,Id2);
        printf("\t\t2. LU decomposing.\n");
        gsl_linalg_LU_decomp(H2,p,&sig);
        printf("\t\t3. Calculating determinant.\n");
        det=gsl_linalg_LU_det(H2,sig);
        printf("\t\t|A-evalI|=%g.\n",det);
        gsl_matrix_memcpy(H2,H);
        gsl_matrix_memcpy(Id2,Id);
    }
    printf("\t2. Verifying eigenvectors:\n");
    printf("In order to do that, we'll use the fact that Aevec=eval*evec.\n");
    for(i=0;i<n;i++){
        eveci=gsl_matrix_column(evec,i);
        gsl_blas_dgemv(CblasNoTrans,1.0,H2,&eveci.vector,0.0,test);
        gsl_vector_scale(test,1.0/gsl_vector_get(eval,i));
        printf("\tObtained vector for eigenvalue %g.\n",gsl_vector_get(eval,i));
        gsl_vector_fprintf(stdout,test,"%g");
    }
    gsl_matrix_free(Id);
    gsl_matrix_free(Id2);
    gsl_matrix_free(H);
    gsl_matrix_free(H);
    gsl_vector_free(eval);
    gsl_vector_free(test);
    gsl_matrix_free(evec);
    gsl_permutation_free(p);
    gsl_eigen_symmv_free(w);
    return 0;
}
