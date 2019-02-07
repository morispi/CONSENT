CC = g++ -std=c++11
CFLAGS  = -Wall -O3 -std=c++11
LDFLAGS = -lpthread

all: alignmentPiles.o reverseComplement.o kMersProcessing.o CONSENT.o DBG.o main.o CONSENT

CONSENT: main.o CONSENT.o utils.o
	$(CC) -o bin/CONSENT src/main.o src/kMersProcessing.o src/reverseComplement.o src/alignmentPiles.o src/CONSENT.o src/DBG.o BMEAN/bmean.o BMEAN/utils.o BMEAN/BOA/align_lpo2.o  BMEAN/BOA/align_lpo_po2.o  BMEAN/BOA/align_score.o  BMEAN/BOA/black_flag.o  BMEAN/BOA/buildup_lpo.o  BMEAN/BOA/create_seq.o  BMEAN/BOA/fasta_format.o  BMEAN/BOA/heaviest_bundle.o  BMEAN/BOA/lpo_format.o  BMEAN/BOA/lpo.o   BMEAN/BOA/msa_format.o  BMEAN/BOA/numeric_data.o  BMEAN/BOA/remove_bundle.o  BMEAN/BOA/seq_util.o  BMEAN/BOA/stringptr.o BMEAN/Complete-Striped-Smith-Waterman-Library/src/*.o $(LDFLAGS)

CONSENT.o: src/CONSENT.cpp src/CONSENT.h src/alignmentPiles.h src/kMersProcessing.h src/DBG.h
	$(CC) -o src/CONSENT.o -c src/CONSENT.cpp $(CFLAGS) -IBMEAN/BOA/

reverseComplement.o: src/reverseComplement.cpp
	$(CC) -o src/reverseComplement.o -c src/reverseComplement.cpp $(CFLAGS)

alignmentPiles.o: src/alignmentPiles.cpp src/Alignment.h src/reverseComplement.h
	$(CC) -o src/alignmentPiles.o -c src/alignmentPiles.cpp $(CFLAGS)

kMersProcessing.o: src/kMersProcessing.cpp
	$(CC) -o src/kMersProcessing.o -c src/kMersProcessing.cpp $(CFLAGS)


DBG.o: src/DBG.cpp src/reverseComplement.h
	$(CC) -o src/DBG.o -c src/DBG.cpp $(CFLAGS)

#BMEAN.o: BMEAN/bmean.cpp
#	$(CC) -o BMEAN/bmean.o -c BMEAN/bmean.cpp $(CFLAGS) -IBMEAN/BOA/

utils.o: BMEAN/utils.cpp
	$(CC) -o BMEAN/utils.o -c BMEAN/utils.cpp $(CFLAGS)

main.o: src/main.cpp src/CONSENT.h
	$(CC) -o src/main.o -c src/main.cpp

clean:
	rm src/*.o bin/CONSENT
