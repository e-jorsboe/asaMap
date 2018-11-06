FLAGS=-O3 

LIBS= -lz

## whether to use EIGEN or not
WITH_EIGEN=1
EIGEN_PATH=/home/emil/software/eigen-eigen-5a0156e40feb

ifdef WITH_EIGEN
FLAGS += -I$(realpath $(EIGEN_PATH))
## for declaring variable for compiling
FLAGS += -DEIGEN
endif

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

PRG=asaMap

all: $(PRG)

.PHONY: clean

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)


%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d
%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d


asaMap: $(OBJ)
	$(CXX) $(FLAGS) -o asaMap *.o $(LIBS)

## with valgrind
##asaMap: $(OBJ)
##	$(CXX) $(FLAGS)  -o asaMap *.o -lz -lpthread -ggdb -static -g


clean:
	rm  -f *.o *.d asaMap *~
