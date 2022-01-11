#DEBUG=-g
CC=g++ -fopenmp
CFLAGS= -std=c++0x -O3 -Wall #-I ~/include

LDFLAGS= -L/usr/local/lib -lgsl -lgslcblas -lm #-L ~/lib

planet: ./src/main.o ./src/EOS.o ./src/EOSlist.o ./src/phase.o ./src/hydro.o ./src/EOSmodify.o

	$(CC) -o $@ $^ $(LDFLAGS) $(DEBUG) 

./src/%.o: %.cpp
	$(CC) $(CFLAGS) $(DEBUG) -c $^

clean:
	rm -f *.o ./src/*.o planet
