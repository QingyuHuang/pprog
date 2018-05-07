CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags) -ggdb
LDLIBS = -lm $$(gsl-config --libs)
CC = gcc

out.txt: main
	./$< > $@

main: main.o root.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o root.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

root.o : root.c            
	$(CC) $(CFLAGS) -c root.c -o root.o -Wall -lm

clean :
	$(RM) *.o *.txt $(objects) main 