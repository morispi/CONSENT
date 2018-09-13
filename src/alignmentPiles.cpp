#include "alignmentPiles.h"
#include <unistd.h>
#include <iostream>

std::vector<std::pair<unsigned, unsigned>> getAlignmentPilesPositions(unsigned tplLen, std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, int overlappingWindows) {
	unsigned* coverages = new unsigned[tplLen];
	unsigned i;
	for (i = 0; i < tplLen; i++) {
		coverages[i] = 1;
	}
	unsigned beg, end;

	for (Alignment al : alignments) {
		beg = al.qStart;
		end = al.qEnd;

		// std::cerr << beg << " ; " << end << std::endl;

		for (i = beg; i <= end; i++) {
			coverages[i]++;
		}
	}

	// std::cerr << "alignments.size() = " << alignments.size() << std::endl;
	// std::cerr << "coverages : " << std::endl;
	// for (i = 0; i < tplLen; i++) {
	// 	std::cerr << coverages[i] << " ";
	// }
	// std::cerr << std::endl;

	std::vector<std::pair<unsigned, unsigned>> pilesPos;

	unsigned curLen = 0;
	beg = 0;
	end = 0;
	i = 0;
	while (i < tplLen) {
		if (curLen >= windowSize) {
			pilesPos.push_back(std::make_pair(beg, beg + windowSize - 1));
			if (overlappingWindows) {
				i = i - overlappingWindows;
			}
			beg = i;
			curLen = 0;
		}
		if (coverages[i] < minSupport) {
			curLen = 0;
			i++;
			beg = i;
		} else {
			curLen++;
			i++;
		}
	}

	delete [] coverages;

	return pilesPos;
}

std::map<std::string, std::string> getSequencesMaps(std::vector<Alignment>& alignments, std::string readsDir) {
	std::map<std::string, std::string> sequences;
	std::string header, seq;

	// Insert template sequence
	std::ifstream f(readsDir + alignments.begin()->qName);
	getline(f, header);
	header.erase(0, 1);
	getline(f, seq);
	sequences[header] = seq;
	f.close();

	// Insert aligned sequences
	for (Alignment al : alignments) {
		std::ifstream f(readsDir + al.tName);
		getline(f, header);
		header.erase(0, 1);
		getline(f, seq);
		sequences[header] = seq;
		f.close();
	}

	return sequences;
}

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentPiles(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, std::string readsDir) {
	int tplLen = alignments.begin()->qLength;

	std::vector<std::pair<unsigned, unsigned>> pilesPos = getAlignmentPilesPositions(tplLen, alignments, minSupport, windowSize, windowOverlap);
	std::map<std::string, std::string> sequences = getSequencesMaps(alignments, readsDir);

	int beg, end, length, shift;
	std::vector<std::string> curPile;

	unsigned curPos = 0;
	unsigned prevPos = 0;
	int passed = 0;
	std::vector<std::vector<std::string>> piles;

	std::vector<std::pair<unsigned, unsigned>> resPilesPos;

	for (std::pair<int, int> p : pilesPos) {
		curPile.clear();
		beg = p.first;
		end = p.second;
		length = end - beg + 1;
		
		// Insert template sequence
		curPile.push_back(sequences[alignments.begin()->qName].substr(beg, length));

		// Insert aligned sequences
		curPos = prevPos;
		while (curPos < alignments.size()) {
			Alignment al = alignments[curPos];
			// For alignments spanning the query window
			// if (al.qStart <= beg and end <= al.qEnd and al.tStart + beg - al.qStart + length - 1 <= al.tEnd) {
			// For all alignments than span, or begin/end in the query window
			if ( ((al.qStart <= beg and al.qEnd > beg) or (end <= al.qEnd and al.qStart < end)) and al.tStart + beg - al.qStart + length - 1 <= al.tEnd) {
				if (beg > al.qStart) {
					shift = beg - al.qStart;
				} else {
					shift = 0;
				}
				std::string tmpSeq = sequences[al.tName].substr(al.tStart, al.tEnd - al.tStart + 1);
				if (al.strand) {
					tmpSeq = rev_comp::run(tmpSeq);
				}
				tmpSeq = tmpSeq.substr(shift, length);
				// Insert aligned sequence only if it is not already present (duplicates can cause POA to crash)
				if (std::find(curPile.begin(), curPile.end(), tmpSeq) == curPile.end()) {
					curPile.push_back(tmpSeq);
				}
				passed++;
			}
			curPos++;
		}

		resPilesPos.push_back(std::make_pair(beg, end));
		piles.push_back(curPile);
	}

	return std::make_pair(resPilesPos, piles);
}