#include <vector>
#include <string>
#include <map>
#include <fstream>
#include "Alignment.h"
#include "reverseComplement.h"

std::vector<std::pair<unsigned, unsigned>> getAlignmentPilesPositions(unsigned tplLen, std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, int overlappingWindows);

std::map<std::string, std::string> getSequencesMaps(std::vector<Alignment>& alignments, std::string readsDir);

std::vector<std::vector<std::string>> getAlignmentPiles(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, std::string readsDir);