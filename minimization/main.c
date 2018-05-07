#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R,const gsl_vector* b, gsl_vector* x);
void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);

int newton(double (*f)(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), void hessian(gsl_vector* x, gsl_matrix* H), gsl_vector* xstart, double eps);
int gradient(double (*f)(gsl_vector* x), gsl_vector* grad_f, gsl_vector* x, double dx);
int hessian(gsl_matrix* H, gsl_matrix* H_new, gsl_vector* y, gsl_vector* s, double lambda);
int quasi_newton(double (*f)(gsl_vector* x), gsl_vector* xstart, double dx, double eps);
double rosenbrock(gsl_vector* x);
void rosenbrock_grad(gsl_vector* x, gsl_vector* fx);
void rosenbrock_hessian(gsl_vector* x, gsl_matrix* H);
double himmelblau(gsl_vector* x);
void himmelblau_grad(gsl_vector* x, gsl_vector* fx);
void himmelblau_hessian(gsl_vector* x, gsl_matrix* H);
double fit_fun(gsl_vector* x);

void print_matrix(gsl_matrix* mat);
void print_vector(gsl_vector* vec);


int main()
{
    gsl_vector* xstart = gsl_vector_alloc(2);

    gsl_vector_set(xstart, 0, 0);
    gsl_vector_set(xstart, 1, 0);

    int nsteps = newton(&rosenbrock, &rosenbrock_grad, &rosenbrock_hessian, xstart, 1e-4);
    if(nsteps == -1)
        return -1;


    printf("Newton's minimization method:\n");
    printf("\nFind the minima of the Rosenbrock's valley function f(x,y) = (1-x)^2+100*(y-x^2)^2 \n");
    printf("The minima [x;y] = \n");
    print_vector(xstart);
    printf(" %i steps it takes for the algorithm to reach the minimum.\n", nsteps);




    gsl_vector_set(xstart, 0, -3);
    gsl_vector_set(xstart, 1, 3);

    nsteps = newton(&himmelblau, &himmelblau_grad, &himmelblau_hessian, xstart, 1e-4);
    if(nsteps == -1)
        return -1;


    printf("\nFind the minima of the Himmelblau's function f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2 \n");
    printf("The minima [x;y] = \n");
    print_vector(xstart);
    printf(" %i steps it takes for the algorithm to reach the minimum.\n", nsteps);



    gsl_vector_set(xstart, 0, 0);
    gsl_vector_set(xstart, 1, 0);

    nsteps = quasi_newton(&rosenbrock, xstart, 1e-6, 1e-4);
    if(nsteps == -1)
        return -1;

    printf("\n\nQuasi-Newton method with Broyden's update:\n");
    printf("\nFind the minima of the Himmelblau's function f(x,y) = (1-x)^2+100*(y-x^2)^2 \n");
    printf("The minima [x;y] = \n");
    print_vector(xstart);
    printf("%i steps it takes for the algorithm to reach the minimum.\n", nsteps);




    gsl_vector_set(xstart, 0, -3);
    gsl_vector_set(xstart, 1, 3);

    nsteps = quasi_newton(&himmelblau, xstart, 1e-6, 1e-4);
    if(nsteps == -1)
        return -1;


    printf("\nFind the minima of the Himmelblau's function f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2 \n");
    printf("The minima [x;y] = \n");
    print_vector(xstart);
    printf(" %i steps it takes for the algorithm to reach the minimum.\n", nsteps);


    gsl_vector_free(xstart);


    xstart = gsl_vector_alloc(3);
    gsl_vector_set(xstart, 0, 1);
    gsl_vector_set(xstart, 1, 1);
    gsl_vector_set(xstart, 2, 1);

    nsteps = quasi_newton(&fit_fun, xstart, 1e-6, 1e-4);
    if(nsteps == -1)
        return -1;


    printf("\nNon-linear lieast-squares fitting:\n\n");
    printf("The minima [x;y] = \n");
    print_vector(xstart);
    printf(" %i steps it takes for the algorithm to reach the minimum.\n", nsteps);


    FILE* file = fopen("data.txt", "w");

    double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
    double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
    double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
    int N = sizeof(t)/sizeof(t[0]);

    for (int i = 0; i < N; ++i)
    {
        fprintf(file, "%g %g %g\n", t[i], y[i], e[i]);
    }
    fprintf(file, "\n\n");
    double a = gsl_vector_get(xstart, 0);
    double b = gsl_vector_get(xstart, 1);
    double c = gsl_vector_get(xstart, 2);
    for (double x = t[0]; x < t[N-1]+1e-5; x+=0.1)
    {
        fprintf(file, "%g %g\n", x, a*exp(-x/b) + c);
    }


    gsl_vector_free(xstart);

    return 0;
};



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

double rosenbrock(gsl_vector* x)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    return pow(1 - a, 2) + 100*pow(b - a*a, 2);
}

void rosenbrock_grad(gsl_vector* x, gsl_vector* fx)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double f1 = -2*(1-a) - 400*(b-a*a)*a;
    double f2 = 200*(b-a*a);
    gsl_vector_set(fx, 0, f1);
    gsl_vector_set(fx, 1, f2);
}


void rosenbrock_hessian(gsl_vector* x, gsl_matrix* H)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double H11 = 2 - 400*(b - a*a) + 800*a*a;
    double H12 = -400*a;
    double H21 = -400*a;
    double H22 = 200;
    gsl_matrix_set(H, 0, 0, H11);
    gsl_matrix_set(H, 0, 1, H12);
    gsl_matrix_set(H, 1, 0, H21);
    gsl_matrix_set(H, 1, 1, H22);
}

double himmelblau(gsl_vector* x)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    return pow(a*a + b -11, 2) + pow(a + b*b -7, 2);
}

void himmelblau_grad(gsl_vector* x, gsl_vector* fx)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double f1 = 4*(a*a + b - 11)*a + 2*(a + b*b - 7);
    double f2 = 2*(a*a + b - 11) + 4*(a + b*b - 7)*b;
    gsl_vector_set(fx, 0, f1);
    gsl_vector_set(fx, 1, f2);
}


void himmelblau_hessian(gsl_vector* x, gsl_matrix* H)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double H11 = 4 * (a * a + b - 11) + 8 * a * a + 2;
    double H12 = 4 * a + 4 * b;
    double H21 = 4 * a + 4 * b;
    double H22 = 4 * (a + b * b - 7) + 8 * b * b + 2;
    gsl_matrix_set(H, 0, 0, H11);
    gsl_matrix_set(H, 0, 1, H12);
    gsl_matrix_set(H, 1, 0, H21);
    gsl_matrix_set(H, 1, 1, H22);
}


double fit_fun(gsl_vector* x)
{
    double a = gsl_vector_get(x, 0);
    double b = gsl_vector_get(x, 1);
    double c = gsl_vector_get(x, 2);

    double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
    double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
    double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
    int N = sizeof(t)/sizeof(t[0]);

    double F = 0.;

    for (int i = 0; i < N; ++i)
    {
        F += pow((a*exp(-t[i]/b) + c - y[i])/e[i], 2);
    }

    return F;
}
