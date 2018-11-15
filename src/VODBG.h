#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

int link(std::unordered_map<int, std::unordered_map<std::string, unsigned>> unordered_mapMerCounts, std::string srcSeed, std::string tgtSeed, unsigned curK, std::set<std::string> &visited, unsigned* curBranches, unsigned dist, std::string curExt, std::string &missingPart, unsigned merSize, unsigned LRLen, unsigned maxBranches, unsigned solidThresh, unsigned minOrder);