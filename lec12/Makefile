CFLAGS = -Wall -std=gnu99 
LDLIBS = -lm 
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
FORMAT = pdf
SCRIPT = log
PLOTTER = gnuplot
SOURCES = error.c
OBJECTS = $(SOURCES: .c=.o)
EXECUTEABLE = error

.PHONEY: all plot.$(FORMAT) clean 

all : plot.$(FORMAT)

plot.$(FORMAT) : out_error.txt Makefile
	echo 'set term $(FORMAT)' >> $(SCRIPT)
	echo 'set title "Error function"' >> $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set size ratio -1' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "u(x)"' >> $(SCRIPT)
	echo 'plot \' >> $(SCRIPT)
	echo '"$(word 1, $^)" using "t":"y(t)" with lines notitle \' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)

out_error.txt : $(EXECUTEABLE)
	./$^ -3.0 3.0 0.2  > $@

$(EXECUTEABLE) : $(SOURCES)
	cc $^ -o $@ $(CFLAGS) $(LDLIBS)
clean : 
	find ./ -type f -executable -delete
	$(RM) *.txt $(SCRIPT) log* plot*