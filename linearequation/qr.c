#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


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


void qr_gr(gsl_matrix* A)
{
    int m=A->size2;
    double angle,x,y,dx,dy;

    for(int i = 0; i < m; i++)
    {
        for(int j = i+1; j< m; j++)
        {
            angle=atan2(gsl_matrix_get(A,j,i),gsl_matrix_get(A,i,i));

            for (int k=i; k < m; k++)
            {
                x = gsl_matrix_get(A,i,k);
                y = gsl_matrix_get(A,j,k);
                dx=x*cos(angle)+y*sin(angle);
                dy=-x*sin(angle)+y*cos(angle);
                gsl_matrix_set(A,i,k,dx);
                gsl_matrix_set(A,j,k,dy);
            }

            gsl_matrix_set(A,j,i,angle);
        }
    }
}



void print_matrix(gsl_matrix* mat)
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

void print_vector(gsl_vector* vec)
{
    int n=vec->size;
    for(int i=0; i<n; i++)
    {
        printf("%10.4f\n",gsl_vector_get(vec, i));
    }
    printf("\n");
}

int check_upper_triangular(gsl_matrix* mat)
{
    int n=mat->size1;
    int m=mat->size2;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            if(i>j && gsl_matrix_get(mat, i, j)!=0)
                return 0;
        }
    }
    return !gsl_matrix_isnull(mat);
}

int check_matrix_equal(const gsl_matrix * a, const gsl_matrix * b)
{
    int n1=a->size1;
    int m1=a->size2;
    int n2=b->size1;
    int m2=b->size2;
    if(n1!=n2||m1!=m2)
        return 0;

    for(int i=0; i<n1; i++)
    {
        for(int j=0; j<m1; j++)
        {
            if(fabs(gsl_matrix_get(a,i,j)-gsl_matrix_get(b,i,j))>1e-6)
                return 0;
        }
    }
    return 1;
}

int check_vector_equal(const gsl_vector * a, const gsl_vector * b)
{
    int n1=a->size;
    int n2=b->size;

    if(n1!=n2)
        return 0;

    for(int i=0; i<n1; i++)
    {
        if(fabs(gsl_vector_get(a,i)-gsl_vector_get(b,i))>1e-6)
            return 0;
    }
    return 1;
}
