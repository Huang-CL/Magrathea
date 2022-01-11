#DEBUG=-g
CC=g++ #-fopenmp

CFLAGS= -std=c++0x -O3 -Wall -I ~/Downloads/gsl/include
LDFLAGS= -L/usr/local/lib -lgsl -lgslcblas -lm -L ~/Downloads/gsl/lib
# HEADDIR=src
SRCDIR = src
LIBDIR = src

BIN = planet
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJ := $(patsubst $(SRCDIR)/%,%,$(SOURCES))
OBJ := $(patsubst %.cpp,%.o,$(OBJ))
OBJ := $(addprefix ./$(LIBDIR)/,$(OBJ))

#all:
#	echo $(OBJ)
#	$(BIN)

planet: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(DEBUG)

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp  
	$(CC) -o $@ -c $< $(CFLAGS) $(DEBUG)

# planet: main.o EOS.o EOSlist.o phase.o hydro.o EOSmodify.o

# 	$(CC) -o $@ $^ $(LDFLAGS) $(DEBUG) 

# %.o: ./src/%.cpp
# 	$(CC) $(CFLAGS) $(DEBUG) -c $^

clean:
	rm -f *.o ./src/*.o planet Magrathea
