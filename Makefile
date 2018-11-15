CC = g++ -std=c++11
CFLAGS  = -Wall -O3 -std=c++11
LDFLAGS = -lpthread

all: alignmentPiles.o reverseComplement.o localAlignment.o kMersProcessing.o LRSelfCorrection.o DBG.o main.o LRSelfCorrection

LRSelfCorrection: main.o LRSelfCorrection.o utils.o
	$(CC) -o bin/LRSelfCorrection src/main.o src/localAlignment.o src/kMersProcessing.o src/reverseComplement.o src/alignmentPiles.o src/LRSelfCorrection.o src/DBG.o BMEAN/bmean.o BMEAN/utils.o BMEAN/BOA/align_lpo2.o  BMEAN/BOA/align_lpo_po2.o  BMEAN/BOA/align_score.o  BMEAN/BOA/black_flag.o  BMEAN/BOA/buildup_lpo.o  BMEAN/BOA/create_seq.o  BMEAN/BOA/fasta_format.o  BMEAN/BOA/heaviest_bundle.o  BMEAN/BOA/lpo_format.o  BMEAN/BOA/lpo.o   BMEAN/BOA/msa_format.o  BMEAN/BOA/numeric_data.o  BMEAN/BOA/remove_bundle.o  BMEAN/BOA/seq_util.o  BMEAN/BOA/stringptr.o BMEAN/Complete-Striped-Smith-Waterman-Library/src/*.o $(LDFLAGS)

LRSelfCorrection.o: src/LRSelfCorrection.cpp src/LRSelfCorrection.h src/alignmentPiles.h src/localAlignment.h src/kMersProcessing.h src/DBG.h
	$(CC) -o src/LRSelfCorrection.o -c src/LRSelfCorrection.cpp $(CFLAGS) -IBMEAN/BOA/

localAlignment.o: src/localAlignment.cpp
	$(CC) -o src/localAlignment.o -c src/localAlignment.cpp $(CFLAGS)

reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS)

alignmentPiles.o: src/alignmentPiles.cpp src/Alignment.h src/reverseComplement.h
	$(CC) -o src/alignmentPiles.o -c src/alignmentPiles.cpp $(CFLAGS)

kMersProcessing.o: src/kMersProcessing.cpp
	$(CC) -o src/kMersProcessing.o -c src/kMersProcessing.cpp $(CFLAGS)

DBG.o: src/DBG.cpp
	$(CC) -o src/DBG.o -c src/DBG.cpp $(CFLAGS)

#BMEAN.o: BMEAN/bmean.cpp
#	$(CC) -o BMEAN/bmean.o -c BMEAN/bmean.cpp $(CFLAGS) -IBMEAN/BOA/

utils.o: BMEAN/utils.cpp
	$(CC) -o BMEAN/utils.o -c BMEAN/utils.cpp $(CFLAGS)

main.o: src/main.cpp src/LRSelfCorrection.h
	$(CC) -o src/main.o -c src/main.cpp

clean:
	rm src/*.o bin/LRSelfCorrection
