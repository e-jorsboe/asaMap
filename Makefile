FLAGS=-O3 


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
	$(CXX) $(FLAGS)  -o asaMap *.o -lz -lpthread

## with valgrind
##asaMap: $(OBJ)
##	$(CXX) $(FLAGS)  -o asaMap *.o -lz -lpthread -ggdb -static -g


clean:
	rm  -f *.o *.d asaMap *~
