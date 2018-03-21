CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


.PHONEY: all
all: plot.pdf rosenbrock.txt

plot.pdf: plot.gpi hydrogen.txt
	gnuplot $<

hydrogen.txt: hydrogen
	./$< >$@

rosenbrock.txt: rosenbrock
	./$< >$@

.PHONEY: clean
clean:
	$(RM) *.txt *.o rosenbrock hydrogen *.pdf