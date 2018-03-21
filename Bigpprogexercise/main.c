#include<stdio.h>
#include<math.h>
double natural_logarithm(double);

int main(){
	double a=0.1, b=10.0, dx=0.2;
	for(double x=a;x<=b;x+=dx)
		printf("%.10g %.10g %.10g\n",x,natural_logarithm(x),log(x));

return 0;
}
