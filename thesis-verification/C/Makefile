# Comipler
CC=gcc
# Compiler flags
CFLAGS=-O3 -Wall
# GMP library root
GMPPATH=gmp-6.1.0
# GMP library file/thing
GMPLIB=$(GMPPATH)/libgmp.la


all: log_conc_check

log_conc_check: log_conc_check.c
	$(GMPPATH)/libtool --mode=link $(CC) $(CFLAGS) log_conc_check.c -o log_conc_check $(GMPLIB)

clean:
	rm log_conc_check .libs/log_conc_check .libs/.DS_Store
	rmdir .libs