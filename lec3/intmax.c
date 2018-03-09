#include<stdio.h>
#include<limits.h>

int main () {
    
    int realMax = INT_MAX;
    printf("max integer: %d\n", realMax);
    
    /* while loop */
    int max = 0;
    int k = 0;
    while(k+1 > k) {
        k++;
        max = k;
    }
    printf("max integer with while loop: %d\n", max);
    
    /* for loop */
    int max2 = 0;
    for(int l=0; l+1>l; l++) {
        max2 = l+1;
    }
    printf("max integer with for loop: %d\n", max2);
    
    /* do while loop */
    int max3 = 0;
    int m = 0;
    do {
        m++;
        max3 = m;
    } while(m+1 > m);
    printf("max integer with do while loop: %d\n \n", max3);
    
   
    return 0;
}
