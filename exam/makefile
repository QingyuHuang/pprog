CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

all: out.txt 

out.txt : main
	./$< > $@

main: main.o gkl_bidiag.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o gkl_bidiag.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

gkl_bidiag.o : gkl_bidiag.c            
	$(CC) $(CFLAGS) -c gkl_bidiag.c -o gkl_bidiag.o -Wall -lm

clean :
	$(RM) *.o *.txt  $(objects) main 
