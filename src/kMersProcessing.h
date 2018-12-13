#include <string>
#include <unordered_map>
#include <vector>

// Get positions of the k-mers of the sequence
std::unordered_map<std::string, std::vector<unsigned>> getKMersPos(std::string sequence, unsigned merSize);

std::unordered_map<std::string, unsigned> getKMersCounts(std::vector<std::string>& sequences, unsigned merSize, unsigned solidThresh);

std::unordered_map<std::string, std::vector<unsigned>> getKMersOccs(std::vector<std::string>& sequences, unsigned merSize);