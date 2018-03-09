CFLAGS = -std=gnu99 -Wall -O 
LDLIBS = -lm              
CC = gcc

.PHONEY: default clean
default: out.txt       	
	cat out.txt     

out.txt: main  
	./main > out.txt   

main: intmax.o intmin.o em.o sum.o equal.o
	$(CC) $(CFLAGS) $(LDFLAGS) intmax.o intmin.o em.o sum.o equal.o -o main $(LDLIBS)    

intmax.o : intmax.c             
	$(CC) $(CFLAGS) -c intmax.c -o intmax.o -Wall             

intmin.o : intmin.c             
	$(CC) $(CFLAGS) -c intmin.c -o intmin.o -Wall   

em.o : em.c             
	$(CC) $(CFLAGS) -c em.c -o em.o -Wall   

sum.o : sum.c
	$(CC) $(CFLAGS) -c sum.c -o sum.o -Wall 

equal.o : equal.c            
	$(CC) $(CFLAGS) -c equal.c -o equal.o -Wall -lm  

clean:
	rm *.o main