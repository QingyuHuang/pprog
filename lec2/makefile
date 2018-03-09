out.txt : math
	./math > out.txt

math : math.o
	cc math.o -o math -lm

math.o : math.c
	cc -c math.c -o math.o

clean:
	rm *.o