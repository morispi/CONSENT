#include <vector>
#include <string>
#include <unordered_map>
#include "Overlap.h"
#include "utils.h"

std::vector<Overlap> getNextReadPile(std::ifstream& f, unsigned maxSupport);

std::unordered_map<std::string, std::string> getSequencesMap(std::vector<Overlap>& alignments, std::unordered_map<std::string, std::vector<bool>>& readIndex);