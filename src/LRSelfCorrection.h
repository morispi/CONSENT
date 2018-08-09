#include <mutex>
#include <future>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm> 
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <utility>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <unistd.h>
#include "Alignment.h"

std::pair<std::string, std::vector<std::string>> computeConsensus(std::vector<std::string>& sequences, std::string readsDir, int minSupport, int merSize, int commonKMers, int solidThresh, int windowOverlap);

std::vector<std::pair<std::string, std::vector<std::string>>> processRead(std::vector<Alignment>& alignments, std::string readsDir, int minSupport, int windowSize, int merSize, int commonKMers, int solidThresh, int windowOverlap);

void processReads(std::vector<std::vector<std::string>>& reads, std::string readsDir, int minSupport, int windowSize, int merSize, int commonKMers, int solidThresh, int windowOverlap);

Alignment getAlignmentFromString(std::string al);

void runCorrection(std::string alignmentFile, std::string readsDir, int minSupport, int windowSize, int merSize, int commonKMers, int solidThresh, int windowOverlap, int nbThreads);

