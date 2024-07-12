CC=gcc
CFLAGS= -O3 -Wall -Wextra -march=native -fomit-frame-pointer

SOURCES= KeccakHash.c KeccakSponge.c KeccakP-1600-opt64.c
OBJECTS= KeccakHash.o KeccakSponge.o KeccakP-1600-opt64.o
SHAKELIB= libshake.a

all: $(SOURCES) $(SHAKELIB)

$(SHAKELIB): $(OBJECTS)
	ar rcs $@ $^

clean:
	rm -f *.o
	rm -f libshake.a
