CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`

plot.pdf: plot.gpi out.dat
	gnuplot $<

out.dat: main
	./$< > $@

main: main.o natural_logarithm.o

clean:
	$(RM) main *.o *.pdf *.dat