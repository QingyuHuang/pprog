CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
.PHONY: all clean

all : out.txt
	cat $<

out.txt : main
	./$< > $@

main : main.o komplex.o

main.o komplex.o : komplex.h

clean :                                                            
	find ./ -type f -executable -delete
	$(RM) *.txt $(SCRIPT) *.pdf log*