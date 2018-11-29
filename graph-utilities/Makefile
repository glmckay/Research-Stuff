# Comipler
CC=gcc
# Compiler flags
CFLAGS=-O3 -Wall


all: g6compl g6conn

g6compl: g6compl.c
	$(CC) $(CFLAGS) g6compl.c -o g6compl

g6conn: g6connected.c
	$(CC) $(CFLAGS) g6connected.c -o g6connected

clean:
	rm g6compl
	rm g6connected