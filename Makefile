#DEBUG=-g
CC=g++ #-fopenmp

CFLAGS= -std=c++0x -O3 -Wall #-I ~/Downloads/gsl/include
LDFLAGS= -L/usr/local/lib -lgsl -lgslcblas -lm -L #~/Downloads/gsl/lib

SRCDIR = src
LIBDIR = src

BIN = planet
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJ := $(patsubst $(SRCDIR)/%,%,$(SOURCES))
OBJ := $(patsubst %.cpp,%.o,$(OBJ))
OBJ := $(addprefix ./$(LIBDIR)/,$(OBJ))

planet: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(DEBUG)

$(LIBDIR)/%.o: $(SRCDIR)/%.cpp  
	$(CC) -o $@ -c $< $(CFLAGS) $(DEBUG)

clean:
	rm -f *.o ./src/*.o planet
