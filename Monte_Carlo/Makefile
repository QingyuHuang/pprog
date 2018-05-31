CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

all:out.txt plot.svg
	cat $<


plot.svg: data.txt
	echo ' \
	set term $(subst .,,$(suffix $@)) background rgb "white"; \
	set out "$@"; \
	set grid; \
	set xlabel "n"; \
	set ylabel "error"; \
	set title "error vs points"; \
	f(x) = a1/sqrt(x); \
	FIT_LIMIT = 1e-8; \
	fit f(x) "$<" using 1:2 via a1; \
	plot "$<" index 0 with points title "error" \
	,a1/sqrt(x) with lines title "1/sqrt(x)" ; \
	' | gnuplot

out.txt: main
	./$< > $@

main: main.o monte_carlo.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o monte_carlo.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

monte_carlo.o : monte_carlo.c            
	$(CC) $(CFLAGS) -c monte_carlo.c -o monte_carlo.o -Wall -lm

clean :
	$(RM) *.o *.txt *.svg *.log $(objects) main 
