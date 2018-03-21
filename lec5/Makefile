CFLAGS = -std=gnu99 -Wall 
LDLIBS = -lm              
CC = gcc

.PHONEY: default
default: out.txt        
	cat out.txt     

out.txt: main  
	./main > out.txt   

main: main.o nvector.o        	
	$(CC) $(CFLAGS) $(LDFLAGS) main.o nvector.o -o main $(LDLIBS)    

main.o : main.c             
	$(CC) $(CFLAGS) -c main.c -o main.o -Wall             

nvector.o : nvector.c            
	$(CC) $(CFLAGS) -c nvector.c -o nvector.o -Wall -lm            
