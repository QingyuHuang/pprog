#include<stdio.h>
#include<math.h>
#include<limits.h>

int main () {
    int max = INT_MAX/3;
    
    float sum_up_float = 0.0;
    for (int i=1; i<=max; i++) {
        sum_up_float += 1.0f/i;
    }
    
    float sum_down_float = 0.0;
    for (int i=max; i>0; i--) {
        sum_down_float += 1.0f/i;
    }
    
    printf("Summing up and down with floats:\n");
    printf("sum_up_float: %g\n", sum_up_float);
    printf("sum_down_float: %g\n", sum_down_float);
    printf("The reason for this difference is that when summing down we are initiating the float to have a high precision whereas another method we begin at a low precision and the high precision terms at the end will be ignored. The series converges as a function of max.");
    printf("The sum has already converged by MAX_INT/3 since we are already beyond the precision of the float.");
    double sum_up_double = 0.0;
    for (int i=1; i<=max; i++) {
        sum_up_double += 1.0/i;
    }
    
    double sum_down_double = 0.0;
    for (int i=max; i>0; i--) {
        sum_down_double += 1.0/i;
    }
    
    printf("Summing up and down with doubles:\n");
    printf("sum_up_double: %g\n", sum_up_double);
    printf("sum_down_double: %g\n\n", sum_down_double);
    printf("Since double can handle more decimels we do not reach that point and get matching results for summing up and down with double.\n\n");
    
    
    return 0;
}

