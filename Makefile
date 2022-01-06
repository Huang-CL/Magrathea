#DEBUG=-g
CC=g++ -fopenmp
CFLAGS= -std=c++0x -O3 -Wall #-I ~/include

LDFLAGS= -L/usr/local/lib -lgsl -lgslcblas -lm #-L ~/lib

planet: main.o EOS.o EOSlist.o phase.o hydro.o EOSmodify.o

	$(CC) -o $@ $^ $(LDFLAGS) $(DEBUG) 

%.o: %.cpp
	$(CC) $(CFLAGS) $(DEBUG) -c $^

clean:
	rm -f *.o planet
