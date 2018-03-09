#include<stdio.h>
#include<math.h>
#include<complex.h>
#include<float.h>

int main(void){
    double e = exp(1);
    double pi = 3.14159265;
    
    double g = tgamma(5);
    printf("Gamma(5) = %g\n", g);
    
    double j = j0(5);
    printf("J_0(5) = %g\n", j);
    
    complex double sq = cpow(-2, 0.5);
    double isq = cimag(sq);
    printf("sqrt(-2) = %gi\n", isq);
    
    complex double ei = cexp(I);
    double rei = creal(ei);
    double iei = cimag(ei);
    printf("e^i = %g + %gi\n", rei, iei);
    
    complex double eip = cexp(pi*I);
    double reip = creal(eip);
    printf("e^i\u03c0 = %g\n", reip);
    
    complex double ie = cpow(I, e);
    double rie = creal(ie);
    double iie = cimag(ie);
    printf("e^i = %g + %gi\n\n", rie, iie);
    
    float f = 0.11111111111111111111111111111111;
    double d = 0.11111111111111111111111111111111;
    long double l;
    l = 0.11111111111111111111111111111111L;
    
    printf("Float:        %.25g\n", f);
    printf("Double:       %.25lg\n", d);
    printf("Long Double:  %.25Lg\n", l);
    
    return 0;
}
