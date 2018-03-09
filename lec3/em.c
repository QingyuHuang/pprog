#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<float.h>
int main(){
/* float */
float g = 1;
while (1 + g != 1){
    g/=2;
};
g*=2;
printf("\nThe machine epsilon for loats are:         %g\n", g);

/* doubles */
double h = 1.0;
while (1.0 + h != 1.0){
    h/=2;
};
h*=2;
printf("The machine epsilon for double is:           %g\n", h);

/* long double */
long double i = 1.0;
while (1.0 + i != 1.0){
    i/=2;
};
i*=2;
printf("The machine epsilon for long double is:      %Lg\n", i);


    return 0;
}
