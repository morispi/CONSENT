CC = g++ -std=c++11
CFLAGS  = -Wall -O3 -std=c++11
LDFLAGS = -lpthread

all: utils.o LRSelfCorrection.o main.o LRSelfCorrection

LRSelfCorrection: src/main.o src/LRSelfCorrection.o
	$(CC) -o bin/LRSelfCorrection src/main.o src/utils.o src/LRSelfCorrection.o $(LDFLAGS)

LRSelfCorrection.o: src/LRSelfCorrection.cpp src/LRSelfCorrection.h src/Alignment.h src/utils.h
	$(CC) -o src/LRSelfCorrection.o -c src/LRSelfCorrection.cpp $(CFLAGS)

utils.o: src/utils.cpp
	$(CC) -o src/utils.o -c src/utils.cpp $(CFLAGS)

main.o: src/main.cpp src/LRSelfCorrection.h
	$(CC) -o src/main.o -c src/main.cpp

clean:
	rm src/*.o bin/LRSelfCorrection
