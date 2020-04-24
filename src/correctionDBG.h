#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>
#include "../BMEAN/utils.h"
#include "utils.h"
#include "DBG.h"

std::string polishCorrection(std::string correctedRead, std::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh);
