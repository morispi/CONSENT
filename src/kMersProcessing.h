#include <string>
#include <vector>
#include "robin_hood.h"

// Get positions of the k-mers of the sequence
robin_hood::unordered_map<std::string, std::vector<unsigned>> getKMersPos(std::string sequence, unsigned merSize);

robin_hood::unordered_map<std::string, unsigned> getKMersCounts(std::vector<std::string>& sequences, unsigned merSize, unsigned solidThresh);

robin_hood::unordered_map<std::string, std::vector<unsigned>> getKMersOccs(std::vector<std::string>& sequences, unsigned merSize);
