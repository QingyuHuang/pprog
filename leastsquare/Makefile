CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

all:plot.svg

plot.svg: out.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "least-squares fit of F_c(x)=c_0*log(x)+c_1+c_2*x to some data" ;\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 using 1:2 with lines title "F_c(x)"\
	,"$<" index 1 using 1:3 with lines title "F_C(x)_{c0+dc0}"\
	,"$<" index 1 using 1:4 with lines title "F_C(x)_{c0-dc0}"\
	,"$<" index 2 using 1:3 with lines title "F_C(x)_{c1+dc1}"\
	,"$<" index 2 using 1:4 with lines title "F_C(x)_{c1-dc1}"\
	,"$<" index 3 using 1:3 with lines title "F_C(x)_{c2+dc2}"\
	,"$<" index 3 using 1:4 with lines title "F_C(x)_{c2-dc2}"\
	'|gnuplot

out.txt: main
	./$< > $@

main: main.o least_squares.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o least_squares.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

least_squares.o : least_squares.c            
	$(CC) $(CFLAGS) -c least_squares.c -o least_squares.o -Wall -lm

clean :
	$(RM) *.o *.txt *.svg $(objects) main 