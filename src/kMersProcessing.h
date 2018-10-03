#include <string>
#include <map>
#include <vector>

// Get positions of the k-mers of the sequence
std::map<std::string, std::vector<unsigned>> getKMersPos(std::string sequence, unsigned merSize);

std::map<std::string, unsigned> getKMersCounts(std::vector<std::string> sequences, unsigned merSize);

std::map<std::string, std::vector<unsigned>> getKMersOccs(std::vector<std::string> sequences, unsigned merSize);