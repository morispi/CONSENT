#include <utility>
#include <string>
#include <vector>
#include <unordered_map>
#include "../BMEAN/utils.h"

std::pair<std::string, std::unordered_map<kmer, unsigned>> computeConsensusReadCorrection(std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path);

std::pair<std::string, std::unordered_map<kmer, unsigned>> computeConsensusAssemblyPolishing(int id, std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path, unsigned nbThreads);