
CC=cc




#CFLAGS=-O2 -std=gnu11 -Wall -Wextra
#CFLAGS=-O3 -Wall -Wextra -std=c99 -Wno-deprecated-declarations
CFLAGS=-O2 -Wall -Wextra -Wno-deprecated-declarations


all:	
	$(CC) $(CFLAGS) -c countime.c
	$(CC) $(CFLAGS) -c function.c 
	$(CC) $(CFLAGS) -c fft.c
	$(CC) $(CFLAGS) -c main.c
	$(CC) $(CFLAGS) -o test countime.o function.o fft.o main.o 

test: all
	./test

clean:
	rm -f *.o
	rm -f test

prettyprint:
	astyle --style=java --indent=tab --pad-header --pad-oper --align-pointer=name --align-reference=name --suffix=none *.c *.h
