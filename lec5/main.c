#include <stdio.h>
#include <stdlib.h>
#include "nvector.h"
#define RND (double)rand()/RAND_MAX

int main() {
    int n = 5;
    
    printf("Main: Testing nvector_alloc ...\n");
    nvector* v = nvector_alloc(n);
    if (v == NULL)
        printf("Test failed\n");
    else
        printf("Test passed\n");
    
    printf("\nMain: now test nvector_set and nvector_get ...\n");
    double value = RND;
    int i = n / 2;
    nvector_set(v, i, value);
    double vi = nvector_get(v, i);
    if (double_equal(vi, value))
        printf("Test passed\n");
    else
        printf("Test failed\n");
    
    printf("\nMain: now testing nvector_dot_product ...\n");
    nvector* a = nvector_alloc(n);
    nvector* dp = nvector_alloc(n);
    for (int i = 0; i < n; i++) {
        double x = RND, y = RND;
        nvector_set(a, i, x);
        nvector_set(v, i, y);
        nvector_set(dp, i, x*y);
    }
    double sum;
    double dot = nvector_dot_product(v, a);
    for (int j=0; j < n; j++) {
        sum += nvector_get(dp, j);
    }
    printf("a * b should = %g\n", sum);
    printf("a * b actually = %g\n", dot);
    if (double_equal(sum, dot))
        printf("Test passed\n");
    else
        printf("Test failed\n");
    
    printf("\nMain: now testing nvector_add ...\n");
    nvector* b = nvector_alloc(n);
    nvector* c = nvector_alloc(n);
    for (int i = 0; i < n; i++) {
        double x = RND, y = RND;
        nvector_set(a, i, x);
        nvector_set(b, i, y);
        nvector_set(c, i, x + y);
    }
    nvector_add(a, b);
    nvector_print("a + b should   = ", c);
    nvector_print("a + b actually = ", a);
    
    if (nvector_equal(c, a))
        printf("Test passed\n");
    else
        printf("Test failed\n");
    
    printf("\nMain: now testing nvector_sub ...\n");
    for (int i = 0; i < n; i++) {
        double x = RND, y = RND;
        nvector_set(a, i, x);
        nvector_set(b, i, y);
        nvector_set(c, i, x - y);
    }
    nvector_sub(a, b);
    nvector_print("a - b should   = ", c);
    nvector_print("a - b actually = ", a);
    
    if (nvector_equal(c, a))
        printf("Test passed\n");
    else
        printf("Test failed\n");
    
    printf("\nMain: now testing nvector_scale ...\n");
    double x = RND;
    printf("x is %g\n", x);
    for (int i = 0; i < n; i++) {
        double y = RND;
        nvector_set(a, i, y);
        nvector_set(b, i, x*y);
    }
    nvector_scale(a, x);
    nvector_print("a * x should   = ", b);
    nvector_print("a * x actually = ", a);
    
    if (nvector_equal(b, a))
        printf("Test passed\n");
    else
        printf("Test failed\n");
    
    nvector_free(v);
    nvector_free(a);
    nvector_free(dp);
    nvector_free(b);
    nvector_free(c);
    
    return 0;
}
