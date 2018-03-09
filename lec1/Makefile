CFLAGS = -std=gnu99 -Wall -O

.PHONEY:default
default: out.txt
	cat out.txt

out.txt: hello
	./hello > out.txt

hello: hello.c
	cc $(CFLAGS) hello.c -o hello

.PHONEY:clean
clean:
	$(RM) hello out.txt