#include<stdio.h>
#include<complex.h>
#include<tgmath.h>
#include"komplex.h"
#include <stdlib.h>
#define RND (double)rand()/RAND_MAX

#define KOMPLEX(z) komplex_new(creal(z),cimag(z));

int main()
{
    komplex a = { 0, 2 };
    komplex b = { 0, 2 };
    komplex z, w;
    int eq1, eq2, eq1_c, eq2_c;
    complex A = a.re + a.im * I;
    complex B = b.re + b.im * I;
    
    printf("\nTesting komplex_set...\n");
    komplex_print("a=", a);
    printf("z = a\n");
    komplex_set(&z, a.re, a.im);
    komplex_print("z=", z);
    if (komplex_equal(z, a))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_new...\n");
    komplex_print("a =", a);
    printf("x = %g, y = %g\n", a.re, a.im);
    printf("z = x+iy\n");
    z = komplex_new(a.re, a.im);
    komplex_print("z =", z);
    if (komplex_equal(z, a))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_add...\n");
    komplex_print("a =", a);
    komplex_print("b =", b);
    z = komplex_add(a, b);
    w = KOMPLEX(A + B);
    komplex_print("a+b should   =", w);
    komplex_print("a+b actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_sub...\n");
    komplex_print("a =", a);
    komplex_print("b =", b);
    z = komplex_sub(a, b);
    w = KOMPLEX(A - B);
    komplex_print("a-b should   =", w);
    komplex_print("a-b actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_equal...\n");
    komplex_print("z = a =", a);
    z = a;
    eq1 = komplex_equal(a, z);
    printf("Is a = z? Answer = %s\n", (eq1 ? "true" : "false"));
    eq1_c = (z.re == a.re && z.im == a.im);
    printf("Actually %s\n", (eq1_c ? "true" : "false"));
    komplex_print("z = b =", b);
    z = b;
    eq2 = komplex_equal(a, z);
    printf("Is a = z? Answer = %s\n", (eq2 ? "true" : "false"));
    eq2_c = (z.re == a.re && z.im == a.im);
    printf("Actually %s\n", (eq2_c ? "true" : "false"));
    if (eq1 == eq1_c && eq2 == eq2_c)
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_mul...\n");
    komplex_print("a=", a);
    komplex_print("b=", b);
    z = komplex_mul(a, b);
    w = KOMPLEX(A * B);
    komplex_print("a*b should   =", w);
    komplex_print("a*b actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_div...\n");
    komplex_print("a=", a);
    komplex_print("b=", b);
    z = komplex_div(a, b);
    w = KOMPLEX(A / B);
    komplex_print("a/b should   =", w);
    komplex_print("a/b actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_conjugate...\n");
    komplex_print("a=", a);
    z = komplex_conjugate(a);
    w = komplex_new(a.re, -a.im);
    komplex_print("a* should   =", w);
    komplex_print("a* actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_abs...\n");
    komplex_print("a=", a);
    z = komplex_abs(a);
    w = KOMPLEX(cabs(A));
    komplex_print("abs(a) should   =", w);
    komplex_print("abs(a) actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("testing komplex_exp...\n");
    komplex_print("a=", a);
    z = komplex_exp(a);
    w = KOMPLEX(exp(A));
    komplex_print("exp(a) should   =", w);
    komplex_print("exp(a) actually =", z);
    if (komplex_equal(w, z))
        printf("test passed\n\n");
    else
        printf("test falied\n\n");
    
    printf("Testing komplex_sin...\n");
    komplex_print("a=", a);
    z = komplex_sin(a);
    w = KOMPLEX(csin(A));
    komplex_print("sin(a) should   =", w);
    komplex_print("sin(a) actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_cos...\n");
    komplex_print("a=", a);
    z = komplex_cos(a);
    w = KOMPLEX(ccos(A));
    komplex_print("cos(a) should   =", w);
    komplex_print("cos(a) actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    printf("Testing komplex_sqrt...\n");
    komplex_print("a=", a);
    z = komplex_sqrt(a);
    w = KOMPLEX(csqrt(A));
    komplex_print("sqrt(a) should   =", w);
    komplex_print("sqrt(a) actually =", z);
    if (komplex_equal(w, z))
        printf("Test passed\n\n");
    else
        printf("Test failed\n\n");
    
    return 0;
}
