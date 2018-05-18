CFLAGS = -Wall -std=gnu11
LDLIBS = -lm
CC = gcc

all:out.txt
	cat $<

out.txt:main
	./$< > $@

main: main.o adapt.o
	$(CC) $(CFLAGS) $(LDFLAGS) main.o adapt.o -o main $(LDLIBS) 

main.o : main.c
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall

adapt.o : adapt.c            
	$(CC) $(CFLAGS) -c adapt.c -o adapt.o -Wall -lm

clean :
	$(RM) *.o *.txt $(objects) main 