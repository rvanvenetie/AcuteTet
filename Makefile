CC=gcc -fopenmp 
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS=  -g -O3 -std=gnu99 $(WARNING_FLAGS) -DINLINE_MACROS
LDFLAGS= -lm 
.c.o:
	$(CC) -c $(CFLAGS) $<

all: main
main: main.o vector.o triangle.o tetraeder.o  mem_list.o tri_list.o
	$(CC) -o main main.o vector.o triangle.o tetraeder.o  mem_list.o tri_list.o  $(LDFLAGS)

test: test.o vector.o triangle.o tetraeder.o  mem_list.o tri_list.o triangulate.o
	$(CC) -o test test.o vector.o triangle.o tetraeder.o  mem_list.o tri_list.o triangulate.o  $(LDFLAGS)
	
run: main
	main
clean:
	rm -f *.o
	rm -f main

vector.o   : vector.h vector.c
triangle.o : vector.h triangle.h triangle.c
tetraeder.o: vector.h triangle.h tetraeder.h tetraeder.c 
triangulate.o: vector.h triangle.h tetraeder.h tri_list.h triangulate.h triangulate.c
mem_list.o : vector.h triangle.h mem_list.h mem_list.c
tri_list.o : mem_list.h tri_list.h tri_list.c 
#datastructures.o: vector.h triangle.h tetraeder.h datastructures.h datastructures.c
main.o     : vector.h triangle.h tetraeder.h mem_list.h tri_list.h main.c
test.o     : vector.h triangle.h tetraeder.h mem_list.h  test.c
