#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "nvector.h"

nvector* nvector_alloc(int n){
    nvector* v = malloc(sizeof(nvector));
    (*v).size = n;
    (*v).data = malloc(n * sizeof(double));
    if( v==NULL )
        fprintf(stderr,"Error in nvector_alloc\n");
    return v;
}

void nvector_free(nvector* v){
    // Free the memory v points to. You need a free statement per malloc. Unused pointers are set to NULL in order to protect against dangling pointer bugs. If v is accessed after it is freed, you may read or overwrite random memory, but if it is a null pointer, most systems would crash, telling you what the error is.
    free(v->data);
    free(v);
    v=NULL;
}

void nvector_set(nvector* v, int i, double value){
    // make sure that i is within the length of (*v).data
    assert(0 <= i && i < (*v).size);
    (*v).data[i]=value;
}

double nvector_get(nvector* v, int i){
    assert(0 <= i && i < (*v).size);
    return (*v).data[i];
}

double nvector_dot_product(nvector* u, nvector* v){
    assert(v->size == u->size);
    double sum;
    for (int i = 0; i < v->size; i++){
        double c = nvector_get(u, i) * nvector_get(v, i);
        sum += c;
    }
    return sum;
}

void nvector_print(char* s, nvector* v){
    printf("%s", s);
    for (int i = 0; i < v->size; i++)
        printf("%g", v->data[i]);
    printf("\n");
}

void nvector_set_zero(nvector* v){
    for (int i = 0; i < v->size; i++){
        nvector_set(v, i, 0.0);
    }
}

int double_equal(double a, double b){
    double TAU = 1e-6, EPS = 1e-6;
    if (fabs(a - b) < TAU)
        return 1;
    if (fabs(a - b) / (fabs(a) + fabs(b)) < EPS / 2)
        return 1;
    return 0;
}

int nvector_equal(nvector* a, nvector* b){
    assert(a->size == b->size);
    for (int i = 0; i < a->size; i++){
        if (!double_equal(a->data[i], b->data[i]))
            return 0;
    }
    return 1;
}

void nvector_add(nvector* a, nvector* b){
    assert(a->size == b->size);
    for (int i = 0; i < a->size; i++){
        double c = nvector_get(a, i) + nvector_get(b, i);
        nvector_set(a, i, c);
    }
}

void nvector_sub(nvector* a, nvector* b){
    assert(a->size == b->size);
    for (int i = 0; i < a->size; i++){
        double c = nvector_get(a, i) - nvector_get(b, i);
        nvector_set(a, i, c);
    }
}

void nvector_scale(nvector* a, double x){
    for (int i = 0; i < a->size; i++){
        double c = nvector_get(a, i) * x;
        nvector_set(a, i, c);
    }
}
