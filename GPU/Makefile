NVCC		= nvcc
CC		= gcc
#CU_FLAGS	= -O3 -g -deviceemu
CU_FLAGS	= -O3 -g -DINLINE_MACROS 
CC_FLAGS	= -O3 -Wall -std=c11

CU_SOURCES	= main.cu
CC_SOURCES	= ../vector.c ../tetraeder.c ../mem_list.c ../triangle.c ../combinations.c

CU_OBJECTS	= $(CU_SOURCES:%.cu=%.o)
CU_PTX		= $(CU_SOURCES:%.cu=%.ptx)
CC_OBJECTS	= $(CC_SOURCES:%.c=%.o)

%.o:		%.cu
		$(NVCC) $(CU_FLAGS) -c $< -o $@

%.o:		%.c
		$(CC) $(CC_FLAGS) -c $< -o $@

%.ptx:		%.cu
		$(NVCC) $(CU_FLAGS) --ptx $< -o $@

main:	$(CU_OBJECTS) $(CC_OBJECTS)
		$(NVCC) $^ -o $@ -lgomp

ptx:		$(CU_PTX) 

clean:
		rm -f *.o main *~
