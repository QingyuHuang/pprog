#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void qr_givens(gsl_matrix* A);
void qr_solve_givens(gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart, double dx, double epsilon);
void newton_with_jacobian(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart, double epsilon);
void newton_with_Jacobian_refined(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J), gsl_vector* xstart, double dx, double epsilon);

void print_matrix(gsl_matrix* mat);
void print_vector(gsl_vector* vec);

void equation (gsl_vector *p, gsl_vector *fx);
void rosenbrock (gsl_vector *p, gsl_vector *fx);
void himmelblau (gsl_vector *p, gsl_vector *fx);

void equation_jacobian (gsl_vector *p, gsl_vector *fx, gsl_matrix* J);
void rosenbrock_jacobian (gsl_vector *p, gsl_vector *fx, gsl_matrix* J);
void himmelblau_jacobian (gsl_vector *p, gsl_vector *fx, gsl_matrix* J);

int main()
{
    double dx,epsilon;
    dx=1e-6;
    epsilon=1e-3;
    gsl_vector *x = gsl_vector_alloc(2);

    printf("Part A:\n");
    printf("\nSolve the system of equations:\n");
    gsl_vector_set(x,0,1);
    gsl_vector_set(x,1,10);

    newton(equation,x,dx,epsilon);
    print_vector(x);

    printf("\nFind the minimum of the Rosenbrock's valley function:\n");
    gsl_vector_set(x,0,5);
    gsl_vector_set(x,1,5);

    newton(rosenbrock,x,dx,epsilon);
    print_vector(x);

    printf("Find the minimum of the Himmelblau's function:\n");
    gsl_vector_set(x,0,5);
    gsl_vector_set(x,1,5);

    newton(himmelblau,x,dx,epsilon);
    print_vector(x);


    printf("Part B:\n");
    printf("\nSolve the system of equations with analytic Jacobian:\n");
    gsl_vector_set(x,0,1);
    gsl_vector_set(x,1,10);

    newton_with_jacobian(equation_jacobian,x,epsilon);
    print_vector(x);

    printf("\nFind the minimum of the Rosenbrock's valley function with analytic Jacobian:\n");
    gsl_vector_set(x,0,5);
    gsl_vector_set(x,1,5);

    newton_with_jacobian(rosenbrock_jacobian,x,epsilon);
    print_vector(x);

    printf("Find the minimum of the Himmelblau's function with analytic Jacobian:\n");
    gsl_vector_set(x,0,5);
    gsl_vector_set(x,1,5);

    newton_with_jacobian(himmelblau_jacobian,x,epsilon);
    print_vector(x);

    printf("Part C:\n");
    printf("\nSolve the system of equations with refined linesearch:\n");
    gsl_vector_set(x,0,1);
    gsl_vector_set(x,1,10);

    newton_with_Jacobian_refined(equation_jacobian,x,dx,epsilon);
    print_vector(x);

    printf("\nFind the minimum of the Rosenbrock's valley function with refined linesearch:\n");
    gsl_vector_set(x,0,5);
    gsl_vector_set(x,1,5);

    newton_with_Jacobian_refined(rosenbrock_jacobian,x,dx,epsilon);
    print_vector(x);

    printf("Find the minimum of the Himmelblau's function with refined linesearch:\n");
    gsl_vector_set(x,0,5);
    gsl_vector_set(x,1,5);

    newton_with_Jacobian_refined(himmelblau_jacobian,x,dx,epsilon);
    print_vector(x);


    gsl_vector_free(x);

    return 0;
}

void equation (gsl_vector *p, gsl_vector *fx)
{
    double x, y, f_x, f_y;
    x = gsl_vector_get(p,0);
    y = gsl_vector_get(p,1);
    f_x = 10000*x*y - 1;
    f_y = exp(-x) + exp(-y) - 1 - 1.0/10000;
    gsl_vector_set(fx,0,f_x);
    gsl_vector_set(fx,1,f_y);
}

void rosenbrock (gsl_vector *p, gsl_vector *fx)
{
    double x, y, f_x, f_y;
    x = gsl_vector_get(p,0);
    y = gsl_vector_get(p,1);
    f_x = 2*(1-x)*(-1)+100*2*(y-x*x)*(-2)*x;
    f_y = 100*2*(y-x*x);
    gsl_vector_set(fx,0,f_x);
    gsl_vector_set(fx,1,f_y);
}

void himmelblau (gsl_vector *p, gsl_vector *fx)
{
    double x, y, f_x, f_y;
    x = gsl_vector_get(p,0);
    y = gsl_vector_get(p,1);
    f_x = 4*x*(x*x+y-11) + 2*(x+y*y-7);
    f_y = 2*(x*x+y-11) + 4*y*(x+y*y-7);
    gsl_vector_set(fx,0,f_x);
    gsl_vector_set(fx,1,f_y);
}


void equation_jacobian (gsl_vector *p, gsl_vector *fx, gsl_matrix* J)
{
    double J00, J11, J01, J10;
    double A = 10000;
    double x = gsl_vector_get(p,0);
    double y = gsl_vector_get(p,1);
    double f_x = A*x*y - 1;
    double f_y = exp(-1*x) + exp(-1*y) - 1 - 1./A;;
    gsl_vector_set(fx,0,f_x);
    gsl_vector_set(fx,1,f_y);

    J00 = A*y;
    J01 = A*x;
    J10 = -1*exp(-1*x);
    J11 = -1*exp(-1*y);
    gsl_matrix_set(J,0,0,J00);
    gsl_matrix_set(J,0,1,J01);
    gsl_matrix_set(J,1,0,J10);
    gsl_matrix_set(J,1,1,J11);
}

void rosenbrock_jacobian (gsl_vector *p, gsl_vector *fx, gsl_matrix* J)
{
    double x, y, f_x, f_y;
    double J00, J11, J01, J10;
    x = gsl_vector_get(p,0);
    y = gsl_vector_get(p,1);
    f_x = 2*(1-x)*(-1)+100*2*(y-x*x)*(-1)*2*x;
    f_y = 100*2*(y-x*x);
    gsl_vector_set(fx,0,f_x);
    gsl_vector_set(fx,1,f_y);

    J00 = 2-400*(y-3*x*x);
    J01 = -400*x;
    J10 = -400*x;
    J11 = 200;
    gsl_matrix_set(J,0,0,J00);
    gsl_matrix_set(J,0,1,J01);
    gsl_matrix_set(J,1,0,J10);
    gsl_matrix_set(J,1,1,J11);
}

void himmelblau_jacobian (gsl_vector *p, gsl_vector *fx, gsl_matrix* J)
{
    double x, y, f_x, f_y;
    double J00, J11, J01, J10;
    x = gsl_vector_get(p,0);
    y = gsl_vector_get(p,1);
    f_x = 4*x*(x*x+y-11) + 2*(x+y*y-7);
    f_y = 2*(x*x+y-11) + 4*y*(x+y*y-7);
    gsl_vector_set(fx,0,f_x);
    gsl_vector_set(fx,1,f_y);
    J00 = 12*x*x + 4*y - 42;
    J01 = 4*(x+y);
    J10 = 4*(x+y);
    J11 = 12*y*y + 4*x - 26;
    gsl_matrix_set(J,0,0,J00);
    gsl_matrix_set(J,0,1,J01);
    gsl_matrix_set(J,1,0,J10);
    gsl_matrix_set(J,1,1,J11);
}

void print_matrix(gsl_matrix* mat)
{
    int n,m,i,j;
    n=mat->size1;
    m=mat->size2;
    for(i=0; i<n; i++)
    {
        for(j=0; j<m; j++)
        {
            printf("%10.4f ",gsl_matrix_get(mat, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(gsl_vector* vec)
{
    int n,i;
    n=vec->size;
    for(i=0; i<n; i++)
    {
        printf("%10.4f\n",gsl_vector_get(vec, i));
    }
    printf("\n");
}
