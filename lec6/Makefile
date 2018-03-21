CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
SCRIPT = plot.gpi
FORMAT = png
PLOTTER = gnuplot

.PHONY: all out.txt clean

all: main-stdio main-cmdline plot.$(FORMAT)

plot.$(FORMAT) : out.io.txt out.cmd.txt Makefile
	echo 'set term $(FORMAT)' >> $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \' >> $(SCRIPT)
	echo '"$(word 1,$^)" with points title "stdio" \' >> $(SCRIPT)
	echo ',"$(word 2,$^)" with points title "cmdline" \' >> $(SCRIPT)
	echo ',sin(x) with lines \' >> $(SCRIPT)
	echo ',cos(x) with lines' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)

main-stdio : main-stdio.c
	cc $^ -o $@ $(CFLAGS) $(LDLIBS)

main-cmdline : main-cmdline.c
	cc $^ -o $@ $(CFLAGS) $(LDLIBS)

out.io.txt: in.dat main-stdio
	./$(word 2, $^) < $< > $@

out.cmd.txt: in.dat main-cmdline
	./$(word 2, $^) `cat $<` > $@

in.dat : Makefile
	cat /dev/null > $@
	for x in `seq 0 0.2 5`; do echo $$x >> $@; done

clean :                                                            
	find ./ -type f -executable -delete
	$(RM) *.dat *.txt $(SCRIPT) *.pdf log*
