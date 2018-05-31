CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

all:out.txt plot.svg
	cat $<


plot.svg: data.txt
	@echo ' \
	set term $(subst .,,$(suffix $@)) background rgb "white"; \
    set out "$@"; \
    set grid; \
    set xlabel "x"; \
    set ylabel "y"; \
    set title "ANN interpolation"; \
    plot "$<" index 0 with points title "points" \
    ,"$<" index 1 with lines title "interpolation" ; \
    ' | gnuplot

out.txt: main
	./$< > $@

main: main.o ann.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o ann.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

ann.o : ann.c            
	$(CC) $(CFLAGS) -c ann.c -o ann.o -Wall -lm

clean :
	$(RM) *.o *.txt *.svg *.log $(objects) main 
