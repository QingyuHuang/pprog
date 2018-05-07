#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>

// difine the max matrix dimension for 100 
#define MATRIXMAXDIMENSION 100
// define the min matrix dimension for 10
#define MATRIXMINDEMENSION 100

// the scope of the element in the matrix
#define MATRIXELEMENTSCOPE 10.0f

// the eps for double data
#define EPS 1.0e-10f

void generateMatrix(double **matrix, int *dimension);
void jacobi(double *matrix, int dimension, double *eigenvalue, double *eigenvector);
void findMaxOffDiagonalElement(double *matrix, int dimension, int *row, int *col);
void matrixMul(double *matrix1, double *matrix2, double *matrixRes, int dimension);
bool checkMatrixJacobi(double *matrix1, int dimension, double *eigenvalue, double *eigenvector);
void printMatrix(double *matrix, int dimension);
void printEigenValueAndEigenVector(double *eigenvalue, double *eigenvector, int dimension);
void printSeperate(int dimension);

int main()
{
    double *matrix = NULL;
    int dimension = 0;
    generateMatrix(&matrix, &dimension);

    printSeperate(dimension);
    printf("original matrix is: \n");
    printMatrix(matrix, dimension);

    double *eigenvalue = (double *)malloc(sizeof(double) * dimension);
    memset(eigenvalue, 0, sizeof(double)*dimension);
    double *eigenvector = (double *)malloc(sizeof(double) * dimension * dimension);
    memset(eigenvector, 0, sizeof(double)*dimension*dimension);
    jacobi(matrix, dimension, eigenvalue, eigenvector);
    bool checkStatus = false;
    printSeperate(dimension);
    printf("chech V'AV==D.\n");
    checkStatus = checkMatrixJacobi(matrix, dimension, eigenvalue, eigenvector);
    if(checkStatus)
    {
        printSeperate(dimension);
        printf("check OK!\n");
    }
    else
    {
        printSeperate(dimension);
        printf("check failed!\n");
    }
    printEigenValueAndEigenVector(eigenvalue, eigenvector, dimension);

    if(matrix != NULL)
    {
        free(matrix);
        matrix = NULL;
    }
    free(eigenvalue);
    free(eigenvector);
    return 0;
}

void generateMatrix(double **matrix, int *dimension)
{
    int i=0, j=0;
    srand(time(NULL));
    *dimension = rand() % (MATRIXMAXDIMENSION - MATRIXMINDEMENSION + 1) + MATRIXMINDEMENSION;

    printf("matrix dimension is %d.\n", *dimension);

    *matrix = (double *)malloc((*dimension) * (*dimension) * sizeof(double));
    memset(*matrix, 0, sizeof(double)*(*dimension)*(*dimension));
    for(i=0; i<(*dimension); i++)
    {
        for(j=i; j<(*dimension); j++)
        {
            (*matrix)[i*(*dimension) + j] = (rand() / (double)(RAND_MAX) - 0.5f) * MATRIXELEMENTSCOPE * 2;
            (*matrix)[j*(*dimension) + i] = (*matrix)[i*(*dimension) + j];
        }
    }

    return;
}

void jacobi(double *matrix, int dimension, double *eigenvalue, double *eigenvector)
{
    int i=0;
    int cycle = 0;
    bool finishing = false;
    int row=-1, col=-1;
    double angle = 0.0f;
    double cosAngle = 0.0f;
    double sinAngle = 0.0f;
    double elementii = 0.0f;
    double elementjj = 0.0f;
    double elementij = 0.0f;
    double elementik = 0.0f;
    double elementjk = 0.0f;
    double *matrixTmp = (double *)malloc(sizeof(double)*dimension*dimension);
    memcpy(matrixTmp, matrix, sizeof(double)*dimension*dimension);
    for(i=0; i<dimension; i++)
    {
        eigenvector[i*dimension+i] = 1;
    }
    while(!finishing)
    {
        findMaxOffDiagonalElement(matrixTmp, dimension, &row, &col);
        angle = 0.5*atan2(2*matrixTmp[row*dimension+col], matrixTmp[col*dimension+col]-matrixTmp[row*dimension+row]);
        cosAngle = cos(angle);
        sinAngle = sin(angle);
        elementii = matrixTmp[row*dimension+row];
        elementjj = matrixTmp[col*dimension+col];
        elementij = matrixTmp[row*dimension+col];

        // update
        matrixTmp[row*dimension+row] = cosAngle*cosAngle*elementii - 2*cosAngle*sinAngle*elementij + sinAngle*sinAngle*elementjj;
        matrixTmp[col*dimension+col] = sinAngle*sinAngle*elementii + 2*cosAngle*sinAngle*elementij + cosAngle*cosAngle*elementjj;
        matrixTmp[row*dimension+col] = (cosAngle*cosAngle - sinAngle*sinAngle)*elementij + cosAngle*sinAngle*(elementii-elementjj);
        matrixTmp[col*dimension+row] = matrixTmp[row*dimension+col];

        for(i=0; i<dimension; i++)
        {
            if(i!=row && i!=col)
            {
                elementik = matrixTmp[i*dimension+row];
                elementjk = matrixTmp[i*dimension+col];
                matrixTmp[i*dimension+row] = cosAngle*elementik - sinAngle*elementjk;
                matrixTmp[row*dimension+i] = matrixTmp[i*dimension+row];
                matrixTmp[i*dimension+col] = sinAngle*elementik + cosAngle*elementjk;
                matrixTmp[col*dimension+i] = matrixTmp[i*dimension+col];
            }
        }

        // eigenvector
        for(i=0; i<dimension; i++)
        {
            elementik = eigenvector[i*dimension+row];
            elementjk = eigenvector[i*dimension+col];
            eigenvector[i*dimension+row] = elementik*cosAngle - elementjk*sinAngle;
            eigenvector[i*dimension+col] = elementik*sinAngle + elementjk*cosAngle;
        }

        // check for finishing
        if(cycle > 0)
        {
            finishing = true;
            for(i=0; i<dimension; i++)
            {
                if(fabs(eigenvalue[i]-matrixTmp[i*dimension+i]) > EPS)
                {
                    finishing = false;
                }
            }
        }

        for(i=0; i<dimension; i++)
        {
            eigenvalue[i] = matrixTmp[i*dimension+i];
        }

        cycle++;
    }

    printSeperate(dimension);
    printf("run cycle counts is %d.\n", cycle);
    free(matrixTmp);
    return;
}

void findMaxOffDiagonalElement(double *matrix, int dimension, int *row, int *col)
{
    int i=0, j=0;
    double minElement = DBL_MIN;
    for(i=0; i<dimension; i++)
    {
        for(j=0; j<dimension; j++)
        {
            if(i!=j)
            {
                if(fabs(matrix[i*dimension+j]) > minElement)
                {
                    *row = i;
                    *col = j;
                    minElement = fabs(matrix[i*dimension+j]);
                }
            }
        }
    }
    return;
}

void matrixMul(double *matrix1, double *matrix2, double *matrixRes, int dimension)
{
    memset(matrixRes, 0, sizeof(double)*dimension*dimension);
    int i=0, j=0, k=0;
    for(i=0; i<dimension; i++)
    {
        for(j=0; j<dimension; j++)
        {
            for(k=0; k<dimension; k++)
            {
                matrixRes[i*dimension +j] += matrix1[i*dimension + k] * matrix2[k*dimension +j];
            }
        }
    }
    return;
}

// check V'AV == D
bool checkMatrixJacobi(double *matrix1, int dimension, double *eigenvalue, double *eigenvector)
{
    int i=0, j=0;
    double *matrixTranspose = (double *)malloc(sizeof(double)*dimension*dimension);
    memset(matrixTranspose, 0, sizeof(double)*dimension*dimension);
    for(i=0; i<dimension; i++)
    {
        for(j=0; j<dimension; j++)
        {
            matrixTranspose[j*dimension + i] = eigenvector[i*dimension +j];
        }
    }
    
    double *matrixMulTmp = (double *)malloc(sizeof(double)*dimension*dimension);
    memset(matrixMulTmp, 0, sizeof(double)*dimension*dimension);
    matrixMul(matrixTranspose, matrix1, matrixMulTmp, dimension);
    matrixMul(matrixMulTmp, eigenvector, matrixTranspose, dimension);

    bool isEqual = true;
    for(i=0; i<dimension; i++)
    {
        if(fabs(matrixTranspose[i*dimension+i]-eigenvalue[i])>EPS)
        {
            isEqual = false;
            break;
        }
    }

    free(matrixTranspose);
    free(matrixMulTmp);
    return isEqual;
}

void printMatrix(double *matrix, int dimension)
{
    printSeperate(dimension);

    int i=0, j=0;
    for(i=0; i<dimension; i++)
    {
        for(j=0; j<dimension; j++)
        {
            printf("%7.4lf ",matrix[i*dimension + j]);
        }
        printf("\n");
    }
    return;
}

void printEigenValueAndEigenVector(double *eigenvalue, double *eigenvector, int dimension)
{
    printSeperate(dimension);
    printf("eigenvalue is :\n");
    printSeperate(dimension);

    int i=0;
    for(i=0; i<dimension; i++)
    {
        printf("%7.4lf ", eigenvalue[i]);
    }
    printf("\n");

    printSeperate(dimension);
    printf("eigenvector is :\n");
    printMatrix(eigenvector, dimension);
    return;
}

void printSeperate(int dimension)
{
    int i=0;
    for(i=0; i<dimension; i++)
    {
        printf("-----");
    }
    printf("\n");
    return;
}
