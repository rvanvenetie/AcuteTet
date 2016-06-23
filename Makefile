UNAME = $(shell echo $$USER)
CXX      = g++

OMP_ENABLE = 0
ifeq ($(OMP_ENABLE), 1)
OMP_FLAGS = -fopenmp -DUSE_OMP
else
OMP_FLAGS =
endif

ifeq ($(UNAME), raymond)
CXXFLAGS = -g  -march=native $(OMP_FLAGS)  -std=c++14 -Wall -Ofast    -Winline
else
CXXFLAGS = -g -march=haswell  $(OMP_FLAGS) -std=c++14 -Wall -Ofast    -Winline
endif
LDFLAGS  = -lboost_program_options

# srun -v -p fatq -t 500:00:00 -N 1 main '/local/rvveneti/fund_33.19_jun.fund'

TARGET = main
IDIR   = #-I/usr/lib/gcc/x86_64-linux-gnu/5/include
SRCS   = $(wildcard *.cpp)
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.depends)

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(IDIR) -c $< -o $@

%.depends: %.cpp
	$(CXX) -M $(CXXFLAGS) $(IDIR)  $< > $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)
	echo $(UNAME)

-include $(DEPS)
