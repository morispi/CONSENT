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
#include "alignmentPiles.h"

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

void removeBadSequences(std::vector<std::string>& sequences, std::string tplSeq, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize);

std::vector<std::pair<std::string, std::string>> computeConsensuses(std::string tplId, std::vector<std::vector<std::string>>& piles, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::string readsDir, unsigned minSupport, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize);

std::vector<std::pair<std::string, std::string>> alignConsensuses(std::vector<Alignment>& alignments, std::string readsDir, std::vector<std::pair<std::string, std::string>>& consensuses);

void processRead(std::vector<Alignment>& alignments, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap);

void processReads(std::vector<std::vector<std::string>>& reads, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap);

void runCorrection(std::string alignmentFile, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads);

