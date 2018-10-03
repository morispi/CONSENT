#include "LRSelfCorrection.h"

int main(int argc, char* argv[]) {
	if (argc < 2) {
		fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-d RawLongReadsDir] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-j threadsNb] \n\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	std::string alignmentFile, readsDir;
	unsigned minSupport, windowSize, nbThreads, opt, merSize, commonKMers, solidThresh, windowOverlap;

	readsDir =  "RawLongReads/";
	minSupport = 10;
	windowSize = 50;
	merSize = 8;
	commonKMers = 2;
	solidThresh = 10;
	windowOverlap = 10;
	nbThreads = 1;


	while ((opt = getopt(argc, argv, "a:d:k:s:l:f:e:p:c:m:j:w:m:")) != -1) {
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
			case 'l':
				windowSize = atoi(optarg);
				break;
			case 'k':
				merSize = atoi(optarg);
				break;
			case 'c':
				commonKMers = atoi(optarg);
				break;
			case 'f':
				solidThresh = atoi(optarg);
				break;
			case 'm':
				windowOverlap = atoi(optarg);
				break;
			case 'j':
				nbThreads = atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-d RawLongReadsDir] [-k merSize] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-f freqThresholdForKMers] [-e maxError] [-p freqThresholdForKPersFreqs] [-c freqThresholdForKPersCons] [-m mode (0 for regions, 1 for cluster)] [-j threadsNb] \n\n", argv[0]);
				exit(EXIT_FAILURE);
        }
    }
    
	runCorrection(alignmentFile, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap, nbThreads);

	return EXIT_SUCCESS;
}