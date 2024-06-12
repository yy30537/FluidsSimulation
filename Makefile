# Compiler and linker
CC = gcc

# Compiler flags
CFLAGS = -Wall -std=c99

# Libraries
LIBS = -lGL -lGLU -lglut -lm

# Targets and dependencies
all: demo

demo: demo.o solver.o vorticity.o
	$(CC) -o demo demo.o solver.o vorticity.o $(LIBS)

demo.o: demo.c
	$(CC) $(CFLAGS) -c demo.c

solver.o: solver.c
	$(CC) $(CFLAGS) -c solver.c

vorticity.o: vorticity.c
	$(CC) $(CFLAGS) -c vorticity.c

clean:
	rm -f *.o demo
