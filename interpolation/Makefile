CFLAG = -std=gnu99
all : linSpline quadSpline cubicSpline a.png

linSpline : linSpline_main.c
	gcc -$(CFLAG) linSpline_main.c -o linSpline 

quadSpline : quadSpline_main.c
	gcc -$(CFLAG) quadSpline_main.c -o quadSpline 

cubicSpline : cubicSpline_main.c
	gcc -$(CFLAG) cubicSpline_main.c -o cubicSpline -lgsl -lgslcblas

a.png: plt data.out
	gnuplot plt
