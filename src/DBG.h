#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include "../BMEAN/utils.h"
#include "reverseComplement.h"
using namespace std;

unsigned extendLeft(std::unordered_map<kmer, unsigned> merCounts, unsigned curK, unsigned extLen, string &LR, unsigned solidThresh);

unsigned extendRight(std::unordered_map<kmer, unsigned> merCounts, unsigned curK, unsigned extLen, string &LR, unsigned solidThresh);

int link(std::unordered_map<kmer, unsigned> mapMerCounts, std::string srcSeed, std::string tgtSeed, unsigned curK, std::set<std::string> &visited, unsigned* curBranches, unsigned dist, std::string curExt, std::string &missingPart, unsigned merSize, unsigned LRLen, unsigned maxBranches, unsigned solidThresh, unsigned minOrder);