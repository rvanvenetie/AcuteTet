CXX      = g++
CXXFLAGS = -g  -march=native -fopenmp -std=c++14 -Wall -Ofast -Winline
LDFLAGS  =

TARGET = main
IDIR   = -I/usr/lib/gcc/x86_64-linux-gnu/5/include
SRCS   = $(wildcard *.cpp)
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.depends)

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(IDIR) -c $< -o $@

%.depends: %.cpp
	$(CXX) -M $(CXXFLAGS) $(IDIR)  $< > $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
