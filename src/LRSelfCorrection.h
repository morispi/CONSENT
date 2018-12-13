#include <mutex>
#include <future>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <set>
#include <algorithm> 
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <unistd.h>
#include <map>
#include "alignmentPiles.h"
#include "../BMEAN/utils.h"

struct POASeq {
	std::string seq;
	int beg;
	int end;

	bool operator<(const POASeq& s2) const {
		if (beg < s2.beg) {
			return true;
		} else if (beg == s2.beg and end < s2.end) {
			return true;
		} else {
			return false;
		}
	}

	POASeq() {
		
	}

	POASeq(std::string s, int b, int e) {
		seq = s;
		beg = b;
		end = e;
	}
};

std::string polishCorrection(std::string correctedRead, std::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh);

// std::vector<std::pair<std::string, std::string>> polishCorrection(std::string correctedRead, std::vector<std::pair<std::pair<int, int>, int>>& corPosPiles, std::vector<std::vector<std::string>>& piles, std::unordered_map<std::string, unsigned>& pilesMers, unsigned merSize, int solidThresh, int minGap, int maxGap);

// std::vector<std::pair<std::string, std::string>> polishCorrection(std::string correctedRead, std::vector<std::pair<std::pair<int, int>, int>>& corPosPiles, std::vector<std::vector<std::string>>& piles, unsigned merSize, int solidThresh, int minGap, int maxGap);

void removeBadSequences(std::vector<std::string>& sequences, std::string tplSeq, std::unordered_map<std::string, unsigned>& merCounts, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize);

std::string computeConsensuses(std::string& readId, std::vector<std::string>& piles, std::pair<unsigned, unsigned>& pilesPos, std::string& readsDir, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& solidThresh, unsigned& windowSize);

std::string alignConsensuses(std::string rawRead, std::string sequence, std::vector<std::string>& consensuses, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::vector<std::string>>& piles, int startPos, unsigned windowSize, unsigned windowOverlap);

void processRead(std::vector<Alignment>& alignments, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap);

void processReads(std::vector<std::vector<std::string>>& reads, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap);

void runCorrection(std::string alignmentFile, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads);

