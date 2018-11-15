#include "alignmentPiles.h"
#include <unistd.h>
#include <iostream>

unsigned* getCoverages(std::vector<Alignment>& alignments) {
	unsigned tplLen = alignments.begin()->qLength;
	unsigned* coverages = new unsigned[tplLen];
	unsigned i;
	for (i = 0; i < tplLen; i++) {
		coverages[i] = 1;
	}
	unsigned beg, end;

	for (Alignment al : alignments) {
		beg = al.qStart;
		end = al.qEnd;

		for (i = beg; i <= end; i++) {
			coverages[i]++;
		}
	}

	return coverages;
}

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

		for (i = beg; i <= end; i++) {
			coverages[i]++;
		}
	}

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

std::unordered_map<std::string, std::string> getSequencesunordered_maps(std::vector<Alignment>& alignments, std::string readsDir) {
	std::unordered_map<std::string, std::string> sequences;
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

std::vector<std::string> getAlignmentPileSeq(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, std::unordered_map<std::string, std::string>& sequences, unsigned qBeg, unsigned end) {
	std::vector<std::string> curPile;
	unsigned length, shift;
	length = end - qBeg + 1;
	unsigned curPos = 0;
	unsigned tBeg, tEnd;

	// Insert template sequence
	if (qBeg + length - 1 >= sequences[alignments.begin()->qName].length()) {
		return curPile;
	}

	curPile.push_back(sequences[alignments.begin()->qName].substr(qBeg, length));

	Alignment al;
	std::string tmpSeq;
	// Insert aligned sequences
	while (curPos < alignments.size()) {
		al = alignments[curPos];
		tBeg = al.tStart;
		tEnd = al.tEnd;
		length = end - qBeg + 1;
		if (qBeg > al.qStart) {
			shift = qBeg - al.qStart;
		} else {
			shift = 0;
		}
		// For alignments spanning the query window
		// if (al.qStart <= qBeg and end <= al.qEnd and al.tStart + qBeg - al.qStart + length - 1 <= al.tEnd) {
		// if (al.qStart <= qBeg and end <= al.qEnd and al.tStart + shift <= al.tEnd) {
		// if (al.qStart <= qBeg and end <= al.qEnd and al.tStart + qBeg - al.qStart + length - 1 < sequences[al.tName].length()) {
		
		// For all alignments than span, or begin/end in the query window
		// if ( ((al.qStart <= qBeg and al.qEnd > qBeg) or (end <= al.qEnd and al.qStart < end)) and al.tStart + shift + length - 1 <= al.tEnd) {
		if ( ((al.qStart <= qBeg and al.qEnd > qBeg) or (end <= al.qEnd and al.qStart < end)) and al.tStart + shift <= al.tEnd) {

			if (qBeg < al.qStart and al.qEnd < end) {
				shift = 0;
				tBeg = std::max(0, (int) al.tStart - ((int) al.qStart - (int) qBeg));
				tEnd = std::min((int) al.tLength - 1, (int) al.tEnd + ((int) end - (int) al.qEnd));
				length = tEnd - tBeg + 1;
			} else if (qBeg < al.qStart) {
				shift = 0;
				tBeg = std::max(0, (int) al.tStart - ((int) al.qStart - (int) qBeg));
				// tEnd = std::min((int) al.tLength - 1, (int) tBeg + (int) length - 1);
				length = std::min((int) length, std::min((int) al.tLength - 1, (int) tBeg + (int) length - 1) - (int) tBeg + 1);
			} else if (al.qEnd < end) {
				tEnd = std::min((int) al.tLength - 1, (int) al.tEnd + ((int) end - (int) al.qEnd));
				// // tBeg = std::max(0, (int) tEnd - (int) length + 1);
				length = std::min((int) length, (int) tEnd - std::max(0, (int) tEnd - (int) length + 1) + 1);
				// length = tEnd - std::max(0, (int) tEnd - (int) length + 1) + 1;
			}

			// tmpSeq = sequences[al.tName].substr(al.tStart, al.tEnd - al.tStart + 1);
			tmpSeq = sequences[al.tName].substr(tBeg, tEnd - tBeg + 1);
			if (al.strand) {
				tmpSeq = rev_comp::run(tmpSeq);
			}

			// length = std::min(length + (int) windowSize, al.qEnd - qBeg + 1 + (int) 1 * (int) windowSize);
			// length = std::min(length, al.qEnd - qBeg + 1);
			tmpSeq = tmpSeq.substr(shift, length);
			
			if (tmpSeq.length() >= 0.25 * windowSize) {
				curPile.push_back(tmpSeq);
			}
		}
		curPos++;
	}

	return curPile;
}

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentPiles(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, std::unordered_map<std::string, std::string> sequences) {
	unsigned tplLen = alignments.begin()->qLength;

	std::vector<std::pair<unsigned, unsigned>> pilesPos = getAlignmentPilesPositions(tplLen, alignments, minSupport, windowSize, windowOverlap);

	std::vector<std::vector<std::string>> piles;

	for (std::pair<int, int> p : pilesPos) {
		piles.push_back(getAlignmentPileSeq(alignments, minSupport, windowSize, windowOverlap, sequences, p.first, p.second));
	}

	return std::make_pair(pilesPos, piles);
}