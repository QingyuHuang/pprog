CC=gcc
CFLAGS=-Wall -g
LIBS=-lm

all :
	make a
	make b
	make c
	echo "spend time" >> timeA.txt
	\time -f "%U" -ao timeA.txt ./a.out >> timeA.txt
	echo "spend time" >> timeB.txt
	\time -f "%U" -ao timeB.txt ./b.out > timeB.txt
	echo "spend time" >> timeC.txt
	\time -f "%U" -ao timeC.txt ./c.out > timeC.txt

a : questionA.o
	${CC} ${CFLAGS} -o a.out questionA.o ${LIBS}

b : questionB.o
	${CC} ${CFLAGS} -o b.out questionB.o ${LIBS}

c : questionC.o
	${CC} ${CFLAGS} -o c.out questionC.o ${LIBS}

questionA.o : questionA.c
	${CC} ${CFLAGS} -c questionA.c

questionB.o : questionB.c
	${CC} ${CFLAGS} -c questionB.c

questionC.o : questionC.c
	${CC} ${CFLAGS} -c questionC.c

.PHONY :
clean:
	rm *.o *.out
