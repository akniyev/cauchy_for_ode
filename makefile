CC = gcc

default: libexamples.a

libexamples.a: examples.o
	ar rcs $@ $^

examples.o: example.c example.h
	$(CC) -c $<

clean:
	rm *.o *.a