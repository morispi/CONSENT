#include <string>
#include <map>
#include <vector>

// Get positions of the k-mers of the sequence
std::map<std::string, std::vector<int>> getKMersPos(std::string sequence, unsigned merSize);

std::map<std::string, int> getKMersCounts(std::vector<std::string> sequences, unsigned merSize);