CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`
CC = gcc
SCRIPT = log.gpi
FORMAT = pdf
PLOTTER = gnuplot

plot.pdf : plot.gpi data.txt
	$(PLOTTER) $<


out.txt data.txt: main
	./$< > out.txt

PHONY : clean
clean :
	$(RM) main *.o *.txt
