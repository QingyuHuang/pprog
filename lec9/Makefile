CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)


.PHONEY: all
all: integ.out.txt hamiltonian.out.txt plot.pdf

plot.pdf: plot.gpi hamiltonian.out.txt
	gnuplot $<

hamiltonian.out.txt: hamiltonian
	./$< >$@

integ.out.txt: integ
	./$< >$@


.PHONEY: clean
clean:
	$(RM) *.txt *.o integ
