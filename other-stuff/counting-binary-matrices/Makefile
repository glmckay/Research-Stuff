# Compiler
CC=gcc
# Compiler flags
CFLAGS=-O3 -Wall
# GMP library root
GMPPATH=gmp-6.1.0
# GMP library file/thing
GMPLIB=$(GMPPATH)/libgmp.la


all: binMatCount

binMatCount: binMatCount.c
	$(GMPPATH)/libtool --mode=link $(CC) $(CFLAGS) binMatCount.c -o bmcount $(GMPLIB)

clean:
	rm bmcount .libs/stable-parts .libs/.DS_Store
	rmdir .libs