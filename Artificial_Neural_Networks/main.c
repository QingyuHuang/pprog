#include"ann.h"

double fit_fun(double x)
{
    return cos(5*x)*sin(2*x);
}
double act_fun(double x)
{
    return x*exp(-x*x);
}

int main()
{
    printf("\nPart A:\nSee Graph.\n");

    int n = 20;

    gsl_vector* x_data = gsl_vector_alloc(n);
    gsl_vector* y_data = gsl_vector_alloc(n);

    FILE* data = fopen("data.txt","w+");

    double a=-2, b=2;
    for(int i=0; i<n; i++)
    {
        double x = a+(b-a)*i/(n-1);
        double y = fit_fun(x);

        gsl_vector_set(x_data,i,x);
        gsl_vector_set(y_data,i,y);

        fprintf(data,"%lg\t%lg\n",x,y);
    }
    fprintf(data,"\n\n");

    int neurons = 10;
    double eps = 1e-6, step = 1e-6;

    ann* network = ann_alloc(neurons,&act_fun);
    for(int i=0; i<neurons; i++)
    {
        gsl_vector_set(network->data,0*network->n+i,a+(b-a)*i/(network->n-1));
        gsl_vector_set(network->data,1*network->n+i,1);
        gsl_vector_set(network->data,2*network->n+i,1);
    }

    ann_train(network,x_data,y_data,eps,step);

    for(double i=a; i<b; i+=fabs(b-a)/1000)
    {
        double result = ann_feed_forward(network,i);
        fprintf(data,"%lg\t%lg\n",i,result);
    }
    fprintf(data,"\n\n");


    gsl_vector_free(x_data);
    gsl_vector_free(y_data);
    ann_free(network);
    return 0;
}
