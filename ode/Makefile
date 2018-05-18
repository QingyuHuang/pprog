CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

all:plot.svg

plot.svg: out.txt makefile
	echo '\
	set term svg background rgb "white";\
	set output "plot.svg";\
	plot \
	"out.txt" every ::15 using 1:4 with points title "sin(x)" \
	,"out.txt" every ::15 using 1:5 with points title "cos(x)" \
	,"out.txt" every ::15 using 1:2 with lines title "y(x)" \
	,"out.txt" every ::15 using 1:3 with lines title "dy(x)" \
	'|gnuplot

out.txt: main
	./$< > $@

main: main.o ode.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o ode.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

ode.o : ode.c            
	$(CC) $(CFLAGS) -c ode.c -o ode.o -Wall -lm

clean :
	$(RM) *.o *.txt *.svg $(objects) main 
