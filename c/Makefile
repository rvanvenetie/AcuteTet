CC=gcc -fopenmp 
WARNING_FLAGS=-Wall -Wextra -Werror-implicit-function-declaration -Wshadow -Wstrict-prototypes -pedantic
CFLAGS=  -g -O3 -std=c11 $(WARNING_FLAGS) -DINLINE_MACROS
LDFLAGS= -lm 
OBJ= vector.o triangle.o tetraeder.o mem_list.o tri_list.o triangulation.o
DEPS=$(OBJ:.o=.h)
MAIN_OBJ= main.o $(OBJ)
TEST_OBJ= test.o $(OBJ)
TRIANG_OBJ= triangulate.o edge_list.o $(OBJ)

EDGE_OBJ= vector.o triangle.o edge_list.o edge.o

%.o : %.c
	$(CC) -c $(CFLAGS) $<

all: main triang test edge
edge: $(EDGE_OBJ)
	$(CC) -o edge $(EDGE_OBJ) $(LDFLAGS)
main: $(MAIN_OBJ) 
	$(CC) -o main $(MAIN_OBJ) $(LDFLAGS)

test: $(TEST_OBJ) 
	$(CC) -o test $(TEST_OBJ) $(LDFLAGS)

triang: $(TRIANG_OBJ)
	$(CC) -o triang $(TRIANG_OBJ) $(LDFLAGS)
	
run: main
	main
clean:
	rm -f *.o
	rm -f main

vector.o     : vector.h 
edge_list.o  : edge_list.h
edge.o       : edge.h
triangle.o   : vector.h triangle.h 
mem_list.o   : vector.h triangle.h mem_list.h
tri_list.o   : vector.h triangle.h mem_list.h tri_list.h 
tetraeder.o  : vector.h triangle.h mem_list.h tri_list.h tetraeder.h
triangulate.o: $(DEPS)
main.o       : $(DEPS)
test.o       : $(DEPS)
