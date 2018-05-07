CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = -lm `gsl-config --libs`

out.txt: main
	./$< > $@

main: main.o qr.o

clean:
	$(RM) main *.o 