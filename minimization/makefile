CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

all: out.txt plot.svg

plot.svg: plot.gpi data.txt; gnuplot plot.gpi

out.txt data.txt: main
	./$< > $@

main: main.o minimization.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o minimization.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

minimization.o : minimization.c            
	$(CC) $(CFLAGS) -c minimization.c -o minimization.o -Wall -lm

clean :
	$(RM) *.o *.txt *.svg $(objects) main 
