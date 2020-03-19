#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "Alignment.h"
#include "reverseComplement.h"
#include "robin_hood.h"

unsigned* getCoverages(std::vector<Alignment>& alignments);

std::vector<std::pair<unsigned, unsigned>> getAlignmentPilesPositions(unsigned tplLen, std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, int overlappingWindows);

robin_hood::unordered_map<std::string, std::string> getSequencesunordered_maps(std::vector<Alignment>& alignments, std::string readsDir);

std::vector<std::string> getAlignmentPileSeq(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, robin_hood::unordered_map<std::string, std::string>& sequences, unsigned beg, unsigned end, unsigned merSize, unsigned maxSupport, unsigned commonKMers);

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentPiles(std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned windowOverlap, robin_hood::unordered_map<std::string, std::string> sequences, unsigned merSize, unsigned commonKMers);
