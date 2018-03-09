#include<stdio.h>
#include<limits.h>

int main () {
    
    int realMin = INT_MIN;
    printf("min integer: %d\n", realMin);
    
    /* while loop */
    int min = 0;
    int k = 0;
    while(k-1 < k) {
        k--;
        min = k;
    }
    printf("min integer with while loop: %d\n", min);
    
    /* for loop */
    int min2 = 0;
    for(int l = 0;  l-1 < l; l--) {
        min2 = l-1;
    }
    printf("min integer with for loop: %d\n", min2);
    
    /* do while loop */
    int min3 = 0;
    int m = 0;
    do {
        m--;
        min3 = m;
    } while(m-1 < m);
    printf("min integer with do while loop: %d\n \n", min3);
    
    return 0;
}

