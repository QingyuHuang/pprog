CFLAGS += -Wall -std=gnu1x `gsl-config --cflags`
LDLIBS += $$(gsl-config --libs)
CC = gcc

#test:
#	echo $(CC)

.PHONEY: all
all: answers.txt results.txt plot.svg clean
answers.txt: linear hilbert inputlinear.txt
	echo The answer to question 1 is plot.svg > $@
	echo Solving question 2 >> $@
	echo The first 3 columns of inputlinear.txt are the matrix and the fourth colum is the result of Ax >> $@
	cat inputlinear.txt | ./$< >> $@
	echo Solving optional question 3 >>$@
	./hilbert >> $@
results.txt: airy
	./$< > $@
linear: linear.c
hilbert: hilbert.c
airy: airy.c
plot.svg: plot.gpi
	gnuplot $<
clean:
	rm airy linear hilbert results.txt