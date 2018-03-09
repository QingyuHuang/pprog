#include"stdio.h"
#include"float.h"
#include"limits.h"
#include"math.h"

int equal(double a, double b, double tau, double epsilon) {
    if(abs(a-b)<tau || abs(a-b)/(abs(a)+abs(b))<epsilon/2){return 1;}
    return 0;
}

