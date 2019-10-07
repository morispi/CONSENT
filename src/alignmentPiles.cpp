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

	// std::cerr << "coverages" << std::endl;
	// for (i = 0; i < tplLen; i++) {
	// 	std::cerr << coverages[i] << " ";
	// }
	// std::cerr << std::endl;

	return coverages;
}

std::vector<std::pair<unsigned, unsigned>> getAlignmentPilesPositions(unsigned tplLen, std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, int overlappingWindows) {
	unsigned* coverages = getCoverages(alignments);
	unsigned i;
	unsigned beg = 0;
	unsigned end = tplLen - 1;
	unsigned meanCov = 0;

	// while (beg < tplLen and coverages[beg] < minSupport) {
	// 	beg++;
	// }

	// while (end > 0 and coverages[end] < minSupport) {
	// 	end--;
	// }

	// = new unsigned[tplLen];
	// unsigned i;
	// for (i = 0; i < tplLen; i++) {
	// 	coverages[i] = 1;
	// }
	// unsigned beg, end;

	// for (Alignment al : alignments) {
	// 	beg = al.qStart;
	// 	end = al.qEnd;

	// 	for (i = beg; i <= end; i++) {
	// 		coverages[i]++;
	// 	}
	// }

	std::vector<std::pair<unsigned, unsigned>> pilesPos;

	unsigned curLen = 0;
	beg = 0;
	end = 0;
	i = 0;
	while (i < tplLen) {
		if (curLen >= windowSize) {
			// std::cerr << "pushed : " << beg << " ; " << beg + curLen - 1 << std::endl;
			pilesPos.push_back(std::make_pair(beg, beg + curLen - 1));
			if (overlappingWindows) {
				i = i - overlappingWindows;
			}
			beg = i;
			curLen = 0;
		}
		// if (coverages[i] < minSupport or coverages[i] > maxSupport) {
		if (coverages[i] < minSupport) {
			curLen = 0;
			i++;
			beg = i;
		} else {
			curLen++;
			i++;
		}
	}

	// Special case for the last window
	int pushed = 0;
	beg = 0;
	end = tplLen - 1;
	curLen = 0;
	i = tplLen - 1;
	while (i > 0 and !pushed) {
		if (curLen >= windowSize) {
			// std::cerr << "pushed : " << end - curLen + 1 << " ; " << end << std::endl;
			pilesPos.push_back(std::make_pair(end - curLen + 1, end));
			pushed = 1;
			end = i;
			curLen = 0;
		}
		// if (coverages[i] < minSupport or coverages[i] > maxSupport) {
		if (coverages[i] < minSupport) {
			curLen = 0;
			i--;
			end = i;
		} else {
			curLen++;
			i--;
		}
	}

	// int totalCov = 0;
	// for (int k = 0; k < tplLen; k++) {
	// 	totalCov += coverages[i];
	// }
	// std::cerr << alignments[0].qName << " : " << totalCov / tplLen << std::endl;

	if (pilesPos.size() > 0) {
		for (i = pilesPos[0].first; i <= pilesPos[pilesPos.size() - 1].second; i++) {
			if (coverages[i] >= minSupport) {
				meanCov++;
			}
		}

		if ((float) meanCov / (pilesPos[pilesPos.size() - 1].second - pilesPos[0].first + 1) < 0.5) {
			// std::cerr << "drop read : " << alignments[0].qName << " : " << std::endl;
			// std::cerr << "support was : " << (float) meanCov / (pilesPos[pilesPos.size() - 1].second - pilesPos[0].first + 1) << std::endl;
			// std::cerr << "meanCov was : " << meanCov << std::endl;
			// std::cerr << "length was : " << pilesPos[pilesPos.size() - 1].second - pilesPos[0].first + 1 << std::endl;
			// std::cerr << "beg : " << pilesPos[0].first << std::endl;
			// std::cerr << "end : " << pilesPos[pilesPos.size() - 1].second << std::endl;
			pilesPos.clear();
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

std::vector<std::string> getAlignmentPileSeq(std::vector<Alignment>& alignments, unsigned minSupport, unsigned windowSize, unsigned windowOverlap, std::unordered_map<std::string, std::string>& sequences, unsigned qBeg, unsigned end, unsigned merSize, unsigned maxSupport, unsigned commonKMers) {
	int bmeanSup;
	bmeanSup = std::min((int) commonKMers, (int) alignments.size() / 2);
	int MAX = 100;
	// std::vector<std::string> curPile(MAX);
	// std::vector<unsigned> curScore(MAX);
	std::vector<std::string> curPile;
	std::vector<unsigned> curScore;
	unsigned length, shift;
	length = end - qBeg + 1;
	unsigned curPos = 0;
	unsigned tBeg, tEnd;
	unsigned nbElems = 0;
	unsigned minScore, posMin;

	if (qBeg + length - 1 >= sequences[alignments.begin()->qName].length()) {
		return curPile;
	}

	// Insert template sequence
	curPile.push_back(sequences[alignments.begin()->qName].substr(qBeg, length));
	
	// If list of best overlaps
	// curPile[nbElems] = sequences[alignments.begin()->qName].substr(qBeg, length);
	// curScore[nbElems] = alignments.begin()->qLength;
	// nbElems++;
	// minScore = curScore[nbElems];
	
	posMin = 0;

	Alignment al;
	std::string tmpSeq;
	// Insert aligned sequences
	while (curPos < alignments.size()) {// and curPile.size() < maxSupport) {
		al = alignments[curPos];
		tBeg = al.tStart;
		tEnd = al.tEnd;
		length = end - qBeg + 1;
		if (qBeg > al.qStart) {
			shift = qBeg - al.qStart;
		} else {
			shift = 0;
		}
		
		// For all alignments than span, or begin/end in the query window
		if ( ((al.qStart <= qBeg and al.qEnd > qBeg) or (end <= al.qEnd and al.qStart < end)) and al.tStart + shift <= al.tEnd) {

			if (qBeg < al.qStart and al.qEnd < end) {
				shift = 0;
				tBeg = std::max(0, (int) al.tStart - ((int) al.qStart - (int) qBeg));
				tEnd = std::min((int) al.tLength - 1, (int) al.tEnd + ((int) end - (int) al.qEnd));
				length = tEnd - tBeg + 1;
			} else if (qBeg < al.qStart) {
				shift = 0;
				tBeg = std::max(0, (int) al.tStart - ((int) al.qStart - (int) qBeg));
				length = std::min((int) length, std::min((int) al.tLength - 1, (int) tBeg + (int) length - 1) - (int) tBeg + 1);
			} else if (al.qEnd < end) {
				tEnd = std::min((int) al.tLength - 1, (int) al.tEnd + ((int) end - (int) al.qEnd));
				length = std::min((int) length, (int) tEnd - std::max(0, (int) tEnd - (int) length + 1) + 1);
			}

			tmpSeq = sequences[al.tName].substr(tBeg, tEnd - tBeg + 1);
			if (al.strand) {
				tmpSeq = rev_comp::run(tmpSeq);
			}

			tmpSeq = tmpSeq.substr(shift, length);

			// Default
			if (tmpSeq.length() >= merSize) {
				curPile.push_back(tmpSeq);
			}

			// Unsorted list of MAX best overlaps
			// if (nbElems >= MAX) {
			// 	if (al.resMatches > minScore) {
			// 		curPile[posMin] = tmpSeq;
			// 		curScore[posMin] = al.resMatches;
			// 	}

			// 	for (int i = 0; i < nbElems; i++) {
			// 		if (curScore[i] < minScore) {
			// 			minScore = curScore[i];
			// 			posMin = i;
			// 		}
			// 	}
			// } else {
			// 	curPile[nbElems] = tmpSeq;
			// 	curScore[nbElems] = al.resMatches;
			// 	if (al.resMatches < minScore) {
			// 		minScore = al.resMatches;
			// 		posMin = nbElems;
			// 	}
			// 	nbElems++;
			// }

			// Sorted list of MAX best overlaps
			// if (tmpSeq.length() >= merSize) {
			// 	int i = nbElems;
			// 	while (i > 1 and al.resMatches > curScore[i-1]) {
			// 		i--;
			// 	}
			// 	for (int j = std::min(nbElems, (unsigned) MAX - 1); j > i; j--) {
			// 		curPile[j] = curPile[j-1];
			// 		curScore[j] = curScore[j-1];
			// 	}
			// 	if (i >= 1 and i < MAX) {
			// 		curPile[i] = tmpSeq;
			// 		curScore[i] = al.resMatches;
			// 		if (nbElems < MAX) {
			// 			nbElems++;
			// 		}
			// 	}
			// }


		}
		curPos++;
	}

	// for (int i = 0; i < nbElems; i++) {
	// 	std::cerr << curScore[i] << " ; ";
	// }
	// std::cerr << std::endl << std::endl;
	// std::cerr << curPile.size() << std::endl;

	return curPile;
}

std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> getAlignmentPiles(std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned windowOverlap, std::unordered_map<std::string, std::string> sequences, unsigned merSize, unsigned commonKMers) {
	unsigned tplLen = alignments.begin()->qLength;

	std::vector<std::pair<unsigned, unsigned>> pilesPos = getAlignmentPilesPositions(tplLen, alignments, minSupport, maxSupport, windowSize, windowOverlap);

	std::vector<std::vector<std::string>> piles;

	for (std::pair<int, int> p : pilesPos) {
		piles.push_back(getAlignmentPileSeq(alignments, minSupport, windowSize, windowOverlap, sequences, p.first, p.second, merSize, maxSupport, commonKMers));
	}

	return std::make_pair(pilesPos, piles);
}
