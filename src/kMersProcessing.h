#include <string>
#include <map>
#include <vector>

// Get positions of non-repeated k-mers in sequence (1-based positions)
// Positions of repeated k-mers are set to -1
std::map<std::string, int> getKMersPos(std::string sequence, unsigned merSize);

std::map<std::string, int> getKMersCounts(std::vector<std::string> sequences, unsigned merSize);