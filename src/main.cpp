#include "LRSelfCorrection.h"

int main(int argc, char* argv[]) {
	if (argc < 2) {
		fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-d RawLongReadsDir] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-j threadsNb] \n\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	std::string alignmentFile, readsDir, readsFile, proofFile;
	unsigned minSupport, maxSupport, windowSize, nbThreads, opt, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, nbReads;

	readsDir =  "RawLongReads/";
	minSupport = 4;
	maxSupport = 100;
	windowSize = 500;
	merSize = 9;
	commonKMers = 8;
	minAnchors = 10;
	solidThresh = 4;
	windowOverlap = 50;
	nbThreads = 1;
	nbReads = 0;


	while ((opt = getopt(argc, argv, "a:A:d:k:s:S:l:f:e:p:c:m:j:w:m:r:R:n:")) != -1) {
        switch (opt) {
			case 'a':
				alignmentFile = optarg;
				break;
			case 'd':
				readsDir = optarg;
				break;
			case 's':
				minSupport = atoi(optarg);
				break;
			case 'S':
				maxSupport = atoi(optarg);
				break;
			case 'l':
				windowSize = atoi(optarg);
				break;
			case 'k':
				merSize = atoi(optarg);
				break;
			case 'c':
				commonKMers = atoi(optarg);
				break;
			case 'A':
				minAnchors = atoi(optarg);
				break;
			case 'f':
				solidThresh = atoi(optarg);
				break;
			case 'm':
				windowOverlap = atoi(optarg);
				break;
			case 'r':
				readsFile = optarg;
				break;
			case 'R':
				proofFile = optarg;
				break;
			case 'n':
				nbReads = atoi(optarg);
				break;
			case 'j':
				nbThreads = atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-d RawLongReadsDir] [-k merSize] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-f freqThresholdForKMers] [-e maxError] [-p freqThresholdForKPersFreqs] [-c freqThresholdForKPersCons] [-m mode (0 for regions, 1 for cluster)] [-j threadsNb] \n\n", argv[0]);
				exit(EXIT_FAILURE);
        }
    }
    
	runCorrection(alignmentFile, readsDir, minSupport, maxSupport, windowSize, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, nbThreads, readsFile, proofFile);

	return EXIT_SUCCESS;
}