CC = g++ -std=c++11
CFLAGS  = -Wall -O3 -std=c++11
LDFLAGS = -lpthread

all: alignmentPiles.o reverseComplement.o localAlignment.o kMersProcessing.o LRSelfCorrection.o main.o LRSelfCorrection

LRSelfCorrection: src/main.o src/LRSelfCorrection.o
	$(CC) -o bin/LRSelfCorrection src/main.o src/localAlignment.o src/kMersProcessing.o src/reverseComplement.o src/alignmentPiles.o src/LRSelfCorrection.o $(LDFLAGS)

LRSelfCorrection.o: src/LRSelfCorrection.cpp src/LRSelfCorrection.h src/alignmentPiles.h src/localAlignment.h src/kMersProcessing.h
	$(CC) -o src/LRSelfCorrection.o -c src/LRSelfCorrection.cpp $(CFLAGS)

localAlignment.o: src/localAlignment.cpp
	$(CC) -o src/localAlignment.o -c src/localAlignment.cpp $(CFLAGS)

reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS)

alignmentPiles.o: src/alignmentPiles.cpp src/Alignment.h src/reverseComplement.h
	$(CC) -o src/alignmentPiles.o -c src/alignmentPiles.cpp $(CFLAGS)

kMersProcessing.o: src/kMersProcessing.cpp
	$(CC) -o src/kMersProcessing.o -c src/kMersProcessing.cpp $(CFLAGS)

main.o: src/main.cpp src/LRSelfCorrection.h
	$(CC) -o src/main.o -c src/main.cpp

clean:
	rm src/*.o bin/LRSelfCorrection
