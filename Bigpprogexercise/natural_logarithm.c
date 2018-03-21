#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double natural_integrand(double x, void* params){
	return 1.0/x;
}

double natural_logarithm(double b){

	gsl_function f;
	f.function = natural_integrand;
	f.params = 0;

    int limit = 1e3, key=6, status1,status2;
    double a=1,acc=1e-8,eps=1e-8,result,err,result2,err2;
    int n;
    gsl_integration_workspace * workspace =gsl_integration_workspace_alloc(limit);

    if(b>0 && b<=0.5)
    {
        b=1.0/b; // change 0<b<1/2 to 2<b<inf
        n=0;
        while(b>2.0)
        {
            b=b/2.0;
            n++;
        }
        status1=gsl_integration_qag(&f,a,b,acc,eps,limit,key,workspace,&result,&err);
        status2=gsl_integration_qag(&f,a,2.0,acc,eps,limit,key,workspace,&result2,&err2);
        result+=n*result2;
        result*=-1;
    }
    else if(b>0.5 && b<=1.0)
    {
        b=1.0/b; // change 1/2<b<1 to 1<b<2
        status1=gsl_integration_qag(&f,a,b,acc,eps,limit,key,workspace,&result,&err);
        status2=GSL_SUCCESS;
        result*=-1.0;
    }
    else if(b>=1.0 && b<2.0)
    {
        status1=gsl_integration_qag(&f,a,b,acc,eps,limit,key,workspace,&result,&err);
        status2=GSL_SUCCESS;
    }
    else
    {
        n=0;
        while(b>2.0)
        {
            b=b/2.0;
            n++;
        }
        status1=gsl_integration_qag(&f,a,b,acc,eps,limit,key,workspace,&result,&err);
        status2=gsl_integration_qag(&f,a,2.0,acc,eps,limit,key,workspace,&result2,&err2);
        result+=n*result2;
    }

	
	gsl_integration_workspace_free(workspace);
	if(status1!=GSL_SUCCESS||status2!=GSL_SUCCESS) return NAN;
	else return result;
}
