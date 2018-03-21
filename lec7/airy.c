#include<stdio.h>
#include<gsl/gsl_sf_airy.h>

int main(){
    int i,N=300;
    double x[2]={-3,3},prec=0.00001;
    double dx,airyA[N],r[300],airyB[N];
    dx=(x[1]-x[0])/N;
    r[0]=x[0];
    airyA[0]=gsl_sf_airy_Ai(r[0],prec);
    airyB[0]=gsl_sf_airy_Bi(r[0],prec);
    printf("%lg\t%lg\t%lg\n",r[0],airyA[0],airyB[0]);
    for(i=1;i<N;i++){
        r[i]=r[i-1]+dx;
        airyA[i]=gsl_sf_airy_Ai(r[i],prec);
        airyB[i]=gsl_sf_airy_Bi(r[i],prec);
        printf("%lg\t%lg\t%lg\n",r[i],airyA[i],airyB[i]);
    }
    return 0;
}

