#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include "Alignment.h"
#include "reverseComplement.h"

unsigned* getCoverages(std::vector<Alignment>& alignments);

std::vector<std::pair<unsigned, unsigned>> getAlignmentPilesPositions(unsigned tplLen, std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, int overlappingWindows);

std::unordered_map<std::string, std::string> getSequencesunordered_maps(std::vector<Alignment>& alignments, std::string readsDir);

std::vector<std::string> getAlignmentPileSeq(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, std::unordered_map<std::string, std::string>& sequences, unsigned beg, unsigned end, unsigned merSize, unsigned maxSupport);

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentPiles(std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned windowOverlap, std::unordered_map<std::string, std::string> sequences, unsigned merSize);