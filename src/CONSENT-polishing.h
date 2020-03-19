#include <mutex>
#include <future>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <set>
#include <algorithm> 
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <unistd.h>
#include <map>

// Polishing

std::pair<std::string, robin_hood::unordered_map<kmer, unsigned>> computeConsensuses(int id, std::string& readId, std::vector<std::string>& piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path, unsigned nbThreads);

std::pair<std::string, std::string> processRead(std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned maxMSA, std::string path, unsigned nbThreads);

