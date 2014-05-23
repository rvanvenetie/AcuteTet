CC=gcc -fopenmp 
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic-errors
CFLAGS=  -g -O3 -std=c11 $(WARNING_FLAGS) -DINLINE_MACROS
LDFLAGS= 
.c.o:
	$(CC) -c $(CFLAGS) $<

all: main
main: main.o vector.o triangle.o tetraeder.o combinations.o mem_list.o triangulate.o
	$(CC) -o main main.o vector.o triangle.o tetraeder.o combinations.o mem_list.o triangulate.o $(LDFLAGS)

test: test.o vector.o triangle.o tetraeder.o combinations.o mem_list.o triangulate.o
	$(CC) -o test test.o vector.o triangle.o tetraeder.o combinations.o mem_list.o triangulate.o $(LDFLAGS)
	
run: main
	main
clean:
	rm -f *.o
	rm -f main

combinations.o: combinations.h combinations.c
vector.o   : vector.h vector.c
triangle.o : vector.h triangle.h triangle.c
tetraeder.o: vector.h triangle.h tetraeder.h tetraeder.c 
triangulate.o: vector.h triangle.h tetraeder.h triangulate.h triangulate.c
mem_list.o : vector.h triangle.h mem_list.h mem_list.c
#datastructures.o: vector.h triangle.h tetraeder.h datastructures.h datastructures.c
main.o     : vector.h triangle.h tetraeder.h datastructures.h mem_list.h triangulate.h main.c
test.o     : vector.h triangle.h tetraeder.h datastructures.h mem_list.h triangulate.h test.c
