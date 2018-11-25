#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include "LRSelfCorrection.h"
#include "localAlignment.h"
#include "kMersProcessing.h"
#include "DBG.h"
#include "../BMEAN/bmean.h"
#include "../BMEAN/utils.h"

std::mutex outMtx;

bool isUpperCase(char c) {
	return 'A' <= c and c <= 'Z';
}

// std::string trimRead(std::string correctedRead, unsigned merSize) {
// 	unsigned beg, end, n;
// 	int i;
// 	i = 0;
// 	n = 0;
// 	while (i < correctedRead.length() and n < merSize) {
// 		if (isUpperCase(correctedRead[i])) {
// 			n++;
// 		} else {
// 			n = 0;
// 		}
// 		i++;
// 	}
// 	beg = i - merSize;// + 1;

// 	i = correctedRead.length() - 1;
// 	n = 0;
// 	while (i >= 0 and n < merSize) {
// 		if (isUpperCase(correctedRead[i])) {
// 			n++;
// 		} else {
// 			n = 0;
// 		}
// 		i--;
// 	}
// 	end = i + merSize;// - 1;

// 	if (end > beg) {
// 		return correctedRead.substr(beg, end - beg + 1);
// 	} else {
// 		return "";
// 	}
// }

int nbCorBases(std::string correctedRead) {
	int n = 0;
	for (unsigned i = 0; i < correctedRead.length(); i++) {
		if ('A' <= correctedRead[i] && correctedRead[i] <= 'Z') {
			n++;
		}
	}

	return n;
}

bool dropRead(std::string correctedRead) {
	return (float) nbCorBases(correctedRead) / correctedRead.length() < 0.75;
}

std::vector<std::string> trimRead(std::string correctedRead, unsigned merSize) {
	std::vector<std::string> res;
	unsigned beg, end, n;
	beg = 0;
	end = 0;
	int i = 0;
	n = 0;

	while (i < correctedRead.length()) {
		while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
			i++;
		}
		beg = i;
		n = 0;

		while (i < correctedRead.length() and n < merSize) {
			if (!isUpperCase(correctedRead[i])) {
				n++;
			} else {
				n = 0;
			}
			i++;
		}

		end = i - n - 1;
		if (end >= beg) {
			std::string split = correctedRead.substr(beg, end - beg + 1);
			std::cerr << "split : " << split << std::endl;
			if (!dropRead(split)) {
				res.push_back(split);
			}
		}
	}

	// res.push_back(correctedRead);

	return res;
}

bool compareLen(const std::string& a, const std::string& b) {
    return (a.size() > b.size()); 
}

int weakMersNumbers(std::string str, std::unordered_map<std::string, unsigned>& merCounts, unsigned merSize, unsigned solidThresh) {
	int i = 0;
	int weakMers = 0;
	while (i < str.length() - merSize + 1) {
		if (merCounts[str.substr(i, merSize)] < solidThresh) {
			weakMers++;
		}
		i++;
	}

	return weakMers;
}

std::vector<std::string> removeBadSequencesPrev(std::vector<std::string>& sequences, std::string tplSeq, std::unordered_map<std::string, unsigned>& merCounts, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize) {
	std::unordered_map<std::string, std::vector<unsigned>> kMers = getKMersPos(tplSeq, merSize);
	// std::cerr << sequences.size() << std::endl;
	// std::cerr << "remove pre : " << merCounts.size() << std::endl;
	// merCounts = getKMersCounts(sequences, merSize);
	// std::cerr << "remove post: " << merCounts.size() << std::endl;
	std::string curSeq;
	unsigned i, j, c;
	int pos, tplBegPos, tplEndPos, curSeqBegPos, curSeqEndPos;
	std::set<std::string> anchored;
	std::string mer;
	bool found;
	unsigned p;
	int beg, end;
	std::unordered_map<int, std::vector<std::string>> unordered_mapCommonMers;
	std::vector<std::string> newSequences;
	newSequences.push_back(sequences[0]);

	// Iterate through the sequences of the alignment window, and only keep those sharing enough solid k-mers with the template
	i = 1;
	while (i < sequences.size()) {
		curSeq = sequences[i];
		j = 0;
		c = 0;
		anchored.clear();
		pos = -1;
		tplBegPos = -1;
		tplEndPos = -1;
		curSeqBegPos = -1;
		curSeqEndPos = -1;

		std::vector<std::string> smallerMSAs;

		// Allow overlapping k-mers
		// while (j < curSeq.length() - merSize + 1) {
		while (j < curSeq.length() - merSize + 1 and c < commonKMers) {
			mer = (curSeq.substr(j, merSize));
			found = false;
			p = 0;
			// Check if the current k-mer (of the current sequence) appears in the template sequence after the current position
			while (!found and p < kMers[mer].size()) {
				found = pos == -1 or kMers[mer][p] > pos;
				p++;
			}
			// Allow repeated k-mers
			if (merCounts[mer] >= solidThresh and found) {
			// Non-repeated k-mers only
			// if (merCounts[mer] >= solidThresh and found and anchored.find(mer) == anchored.end() and kMers[mer].size() == 1) {
				pos = kMers[mer][p-1];
				if (tplBegPos == -1) {
					tplBegPos = pos;
					curSeqBegPos = j;
				}

				c++;
				j += 1;
				anchored.insert(mer);
				tplEndPos = pos + merSize - 1;
				curSeqEndPos = j + merSize - 2;
			} else {
				j += 1;
			}
		}

		// Non-overlapping k-mers only
		// while (j < curSeq.length() - merSize + 1) {
		// while (j < curSeq.length() - merSize + 1 and c < commonKMers) {
		// 	mer = (curSeq.substr(j, merSize));
		// 	bool found = false;
		// 	unsigned p = 0;
		// 	while (!found and p < kMers[mer].size()) {
		// 		found = pos == -1 or kMers[mer][p] >= pos + merSize;
		// 		p++;
		// 	}
		// 	// Allow repeated k-mers
		// 	// if (merCounts[mer] >= solidThresh and found) {
		// 	// Non-repeated k-mers only
		// 	if (merCounts[mer] >= solidThresh and found and anchored.find(mer) == anchored.end() and kMers[mer].size() == 1) {
		// 		pos = kMers[mer][p-1];
		// 		if (tplBegPos == -1) {
		// 			tplBegPos = pos;
		// 			curSeqBegPos = j;
		// 		}

		// 		c++;
		// 		j += merSize;
		// 		anchored.insert(mer);
		// 		tplEndPos = pos + merSize - 1;
		// 		curSeqEndPos = j - 1;
		// 	} else {
		// 		j += 1;
		// 	}
		// }

		if (c >= commonKMers) {
			// Also get previous / following nts, as they can still correct the template
            beg = std::max(0, curSeqBegPos - tplBegPos);
            end = std::min(curSeqEndPos + tplSeq.length() - tplEndPos - 1, curSeq.length() - 1);
            // beg = curSeqBegPos;
            // end = curSeqEndPos - curSeqBegPos;
            // std::cerr << "beg : " << beg << " ; " << "end : " << end << std::endl;
            // std::string tmpSeq = sequences[i].substr(beg, end + 1);
			// newSequences.push_back(tmpSeq);
			// newSequences.push_back(sequences[i].substr(beg, end + 1));
			// if (end - beg + 1 >= 0.25 * windowSize) {
            	// std::string fctSeq = sequences[i].substr(beg, end - beg + 1);
            	// int weakMers = weakMersNumbers(fctSeq, merCounts, merSize, solidThresh);
            	// int totalMers = fctSeq.length() - merSize + 1;
				// std::cerr << "weak k-mers : " << (float) weakMers / (float) totalMers * 1.0 << std::endl;
				// if ((float) weakMers / (float) totalMers * 1.0 < 0.8) {
					newSequences.push_back(sequences[i].substr(beg, end - beg + 1));
				// }
			// }

			// Only add the substring from the leftmost anchor to the rightmost anchor to the MSA (anchor = shared k-mer wth the template)
			// unordered_mapCommonMers[c].push_back(tmpSeq);
			// Add the whole sequence to the MSA
			// unordered_mapCommonMers[c].push_back(sequences[i]);

		}

		i++;
	}

	return newSequences;
}

void toUpperCase(std::string& s, int beg, int end) {
	std::locale loc;
	for (int i = beg; i < end; i++) {
		s[i] = std::toupper(s[i], loc);
	}
}

void toLowerCase(std::string& s, int beg, int end) {
	std::locale loc;
	for (int i = beg; i < end; i++) {
		s[i] = std::tolower(s[i], loc);
	}
}



std::string weightConsensus(std::string& consensus, std::vector<std::string>& pile, std::unordered_map<std::string, unsigned>& merCounts, unsigned merSize, unsigned windowSize, unsigned solidThresh) {
	std::vector<std::string> splits;
	std::string curSplit;

	std::string header = "";
	std::string sequence = "";
	std::string curFct;

	unsigned i = 0;
	while (i < consensus.length() - merSize + 1) {
		curFct = consensus.substr(i, merSize);
		toUpperCase(curFct, 0, merSize);
		if (merCounts[curFct] >= solidThresh) {
			toUpperCase(consensus, i, i + merSize - 1);
		} else {
			toLowerCase(consensus, i, i + merSize - 1);
		}
		i++;
	}

	return consensus;
}

std::vector<std::pair<std::string, std::string>> computeConsensuses(std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, std::unordered_map<std::string, unsigned>& pilesMers, std::string& readsDir, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& solidThresh, unsigned& windowSize) {
	auto start_antoine = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<std::string, std::string>> res(1);

	for (std::string s : piles) {
		std::cerr << s << std::endl;
	}

	std::vector<std::vector<std::string>> result = MSABMAAC(piles, 7, 7);
	auto end_antoine = std::chrono::high_resolution_clock::now();
	std::cerr << "antoine took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_antoine - start_antoine).count() << " ms\n";

	auto c_start = std::chrono::high_resolution_clock::now();
	std::string consSeq = result[0][0];
	// Align the computed consensus to the template, to retrieve the (raw, corrected) pairs of subsequences
	std::string corTpl = consSeq;
	std::string rawTpl = piles[0];
	if (!consSeq.empty()) {
		std::cerr << ">consSeq" << std::endl << consSeq << std::endl;
		int i = 0;
		std::pair<std::pair<int, int>, std::pair<int, int>> posPair = NeedlemanWunschLocalAlignments(piles[0], consSeq);
		std::pair<int, int> rawPos = posPair.first;
		std::pair<int, int> corPos = posPair.second;
		corTpl = consSeq.substr(corPos.first, corPos.second - corPos.first + 1);
		rawTpl = piles[0].substr(rawPos.first, rawPos.second - rawPos.first + 1);
		// std::cerr << "on template : " << rawPos.first << " - " << rawPos.second << std::endl;
		// std::cerr << "on consensu : " << corPos.first << " - " << corPos.second << std::endl;
		// std::cerr << "size before : " << consSeq.size() << std::endl;
		// std::cerr << "size after : " << corTpl.size() << std::endl;

		if (corTpl.length() >= merSize) {
			corTpl = weightConsensus(corTpl, piles, pilesMers, merSize, windowSize, solidThresh);
		}

		std::cerr << ">corTpl" << std::endl << corTpl << std::endl;
		// std::cerr << ">newCorTpl" << std::endl << newCorTpl << std::endl;
		res[0] = std::make_pair(rawTpl, corTpl);
		i++;
	}
	// cons.close();
	auto c_end = std::chrono::high_resolution_clock::now();
	std::cerr << "voting took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

	// c_end = std::chrono::high_resolution_clock::now();

	return res;
}

std::pair<std::string, std::vector<std::pair<std::pair<int, int>, int>>> alignConsensuses(std::string rawRead, std::string sequence, std::vector<std::pair<std::string, std::string>>& consensuses, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::vector<std::string>>& piles, int startPos) {
	std::vector<std::pair<std::string, std::string>> correctedReads;
	int beg, end;
	std::string outSequence;
	outSequence = sequence;
	std::transform(outSequence.begin() + startPos, outSequence.end(), outSequence.begin() + startPos, ::tolower);
	std::vector<std::pair<std::pair<int, int>, int>> corPosPiles(consensuses.size());

	// Anchor the corrected templates on the read
	std::string corWindow;
	unsigned i = 0;
	std::pair<std::string, std::string> c;
	std::string tmpSequence;
	for (i = 0; i < consensuses.size(); i++) {
		c = consensuses[i];
		// Don't proceed if no consensus was produced by POA
		if (!c.second.empty()) {
			// Replace the template by its correction on the read
			std::transform(c.first.begin(), c.first.end(), c.first.begin(), ::tolower);
			tmpSequence = outSequence;
			std::transform(tmpSequence.begin(), tmpSequence.end(), tmpSequence.begin(), ::tolower);
			beg = (int) tmpSequence.find(c.first);
			end = beg + c.first.length() - 1;
			if ((int) beg != -1) {
				outSequence.replace(beg, c.first.length(), c.second);
			}
			end = beg + c.second.length() - 1;
		} else {
			// TODO: may cause bugs
			beg = -1;
			end = -1;
		}
		// Specify that positions beg..end of the corrected read correspond to the current pile
		corPosPiles[i] = std::make_pair(std::make_pair(beg, end), i);
	}

	return std::make_pair(outSequence, corPosPiles);
}

int getNextSrc(std::string correctedRead, unsigned beg, unsigned merSize) {
	unsigned nb = 0;
	unsigned i = beg;

	while (i < correctedRead.length() and (isUpperCase(correctedRead[i]) or nb < merSize)) {
		if (isUpperCase(correctedRead[i])) {
			nb++;
		} else {
			nb = 0;
		}
		i++;
	}

	return nb >= merSize ? i - 1 : -1;
}

int getNextDst(std::string correctedRead, unsigned beg, unsigned merSize) {
	unsigned nb = 0;
	unsigned i = beg;

	while (i < correctedRead.length() and nb < merSize) {
		if (isUpperCase(correctedRead[i])) {
			nb++;
		} else {
			nb = 0;
		}
		i++;
	}

	return nb >= merSize ? i - 1 : -1;
}

std::pair<std::vector<int>, int> getRegPiles(std::vector<std::pair<std::pair<int, int>, int>>& corPosPiles, std::vector<std::vector<std::string>>& piles, int beg, int end, int lastPile) {
	int pileBeg, pileEnd;
	std::vector<int> res;

	unsigned i = lastPile;
	while (i < corPosPiles.size() and corPosPiles[i].first.first <= end) {
		pileBeg = corPosPiles[i].first.first;
		pileEnd = corPosPiles[i].first.second;

		if ((pileBeg <= beg and beg <= pileEnd) or (pileBeg <= end and end <= pileEnd)) {
			// for (unsigned j = 0; j < piles[corPosPiles[i].second].size(); j++) {
			// 	// Push current pile, as it spans the required positions
			// 	res.push_back(piles[corPosPiles[i].second][j]);

			// 	// Push previous and following piles, as they can contain information
			// 	// if (j - 1 >= 0) {
			// 	// 	res.push_back(piles[corPosPiles[i].second][j-1]);
			// 	// }
			// 	// if (j + 1 < piles[corPosPiles[i].second].size()) {
			// 	// 	res.push_back(piles[corPosPiles[i].second][j+1]);
			// 	// }
			// }
			res.push_back(corPosPiles[i].second);
		}
		i++;
	}

	return std::make_pair(res, i - 1);
}

std::vector<std::string> removeUnanchoredPiles(std::vector<std::string>& regPiles, std::string src, std::string dst, unsigned merSize) {
	std::vector<std::string> res;
	std::unordered_map<std::string, std::vector<unsigned>> mersPos;
	std::vector<unsigned> posSrc;
	std::vector<unsigned> posDst;

	for (std::string s : regPiles) {
		mersPos = getKMersPos(s, merSize);

		posSrc = mersPos[src];
		posDst = mersPos[dst];

		if (!posSrc.empty() or !posDst.empty()) {
			res.push_back(s);
		}

		// bool found = false;
		// int p, q;			
		// if (!posSrc.empty() and !posDst.empty()) {
		// 	p = 0;
		// 	found = false;
		// 	while (!found and p < posSrc.size()) {
		// 		q = 0;
		// 		while (!found and q < posDst.size()) {
		// 			found = posDst[q] > posSrc[p];
		// 			q++;
		// 		}
		// 		p++;
		// 	}

		// 	if (found) {
		// 		// res.push_back(s.substr(posSrc[p-1], posDst[q-1] - posSrc[p-1] + merSize));
		// 		res.push_back(s);
		// 	}
		// }
	}

	return res;
}

std::pair<int, int> getNextBadRegion(std::string correctedRead, unsigned beg, unsigned merSize) {
	int posBeg, posEnd;
	unsigned i = beg;

	while (i < correctedRead.length() and isUpperCase(correctedRead[i])) {
		i++;
	}
	posBeg = i;

	int nb = 0;
	while (i < correctedRead.length() and nb < merSize) {
		if (isUpperCase(correctedRead[i])) {
			nb++;
		} else {
			nb = 0;
		}
		i++;
	}
	posEnd = i - nb - 1;

	return std::make_pair(posBeg, posEnd);
}


// Anchors without repeated k-mers
std::vector<std::string> getAnchors(std::unordered_map<std::string, unsigned>& merCounts, std::string region, unsigned merSize, int nb) {
	std::vector<std::string> res;

	// Get the counts of the k-mers of the pile
	// std::unordered_map<std::string, unsigned> merCounts = getKMersCounts(pile, merSize);

	std::unordered_map<std::string, std::vector<unsigned>> mersPos = getKMersPos(region, merSize);

	// Consider all k-mers of the region as potential anchors
	std::vector<std::string> candidates(region.size() - merSize + 1);
	for (unsigned i = 0; i < region.size() - merSize + 1; i++) {
		candidates[i] = region.substr(i, merSize);
	}

	// Sort the k-mers set in ascending order of number of occurrences
	std::sort(candidates.begin(), candidates.end(), 
		[&merCounts](std::string& n1, std::string& n2) {
			int occ1 = merCounts[n1];
			int occ2 = merCounts[n2];
			return  occ1 > occ2;
		}
	);

	// Add the nb most frequent k-mers to the result set
	unsigned i = 0;
	int n = 0;
	while (i < candidates.size() and n < nb) {
		if (mersPos[candidates[i]].size() == 1) {
			res.push_back(candidates[i]);
			n++;
		}
		i++;
	}

	return res;
}

std::vector<std::pair<std::string, std::string>> polishCorrection(std::string correctedRead, std::vector<std::pair<std::pair<int, int>, int>>& corPosPiles, std::vector<std::vector<std::string>>& piles, std::vector<std::unordered_map<std::string, unsigned>>& pilesMers, unsigned merSize, int solidThresh, int minGap, int maxGap) {
	auto main_start = std::chrono::high_resolution_clock::now();
	std::unordered_map<std::string, unsigned> merCounts;
	std::set<std::string> visited;
	unsigned curBranches;
	unsigned dist;
	std::string curExt;
	std::string correctedRegion;
	unsigned maxSize;
	unsigned maxBranches = 25;
	std::vector<std::pair<std::string, std::string>> corList;
	int minOrder = merSize;
	int zone = 1;
	int srcBeg, srcEnd, dstBeg, dstEnd;//, k;
	int tmpSrcBeg, tmpSrcEnd, tmpDstBeg, tmpDstEnd;
	std::string src, dst;
	int lastPile = 0;
	std::pair<int, int> pos;
	std::unordered_map<std::string, unsigned> solidMerCounts;
	int curCount;
	std::vector<std::pair<std::string, std::string>> anchors;
	std::pair<std::vector<int>, int> p;
	std::vector<int> regPilesId;
	unsigned anchorNb;
	std::string srcZone, dstZone;

	// Skip uncorrected head of the read
	auto skip_start = std::chrono::high_resolution_clock::now();
	// std::cerr << "declaration took " << std::chrono::duration_cast<std::chrono::milliseconds>(skip_start - main_start).count() << " ms\n";
	unsigned i = 0;
	while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
		i++;
	}
	auto skip_end = std::chrono::high_resolution_clock::now();
	// std::cerr << "skipping took " << std::chrono::duration_cast<std::chrono::milliseconds>(skip_end - skip_start).count() << " ms\n";
	// search for poorly supported regions bordered by solid corrected regions
	while (i < correctedRead.length()) {
		auto c_start = std::chrono::high_resolution_clock::now();
		srcEnd = getNextSrc(correctedRead, i, merSize + zone);
		dstEnd = getNextDst(correctedRead, srcEnd + 1, merSize + zone);
		// std::cerr << "stats : " << i << " " << srcEnd << " " << dstEnd << std::endl;
		srcBeg = srcEnd - merSize - zone + 1;
		dstBeg = dstEnd - merSize - zone + 1;
		// std::cerr << correctedRead.substr(srcBeg, dstEnd - srcBeg + 1) << std::endl;
		auto c_end = std::chrono::high_resolution_clock::now();
		// std::cerr << "srcdst took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

		// k = zone;
		tmpSrcBeg = srcBeg + zone;
		tmpSrcEnd = tmpSrcBeg + merSize - 1;
		tmpDstBeg = dstBeg;
		tmpDstEnd = tmpDstBeg + merSize - 1;

		auto main_loop_beg = std::chrono::high_resolution_clock::now();
		// Polish the poorly supported region region if 2 anchors were found
		if (srcEnd != -1 and dstEnd != -1) {
			if (dstBeg - srcEnd - 1 <= maxGap) {
				auto c_start = std::chrono::high_resolution_clock::now(); //A	
				// Retrieve sequences of the window overlapping to the region to polish
				// auto c_start = std::chrono::high_resolution_clock::now();
				// auto clean_end = std::chrono::high_resolution_clock::now();
				p = getRegPiles(corPosPiles, piles, srcBeg, dstEnd, lastPile);
				// auto c_end = std::chrono::high_resolution_clock::now();
				// std::cerr << "getRegPiles took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
				regPilesId = p.first;
				lastPile = p.second;
				// std::cerr << "regPilesSizes" << regPilesId.size() << std::endl;

				merCounts = pilesMers[regPilesId[0]];
				for (int i = 1; i < regPilesId.size(); i++) {
					for (std::pair<std::string, int> pp : pilesMers[regPilesId[1]]) {
						merCounts[pp.first] += pp.second;
					}
				}
				auto c_end = std::chrono::high_resolution_clock::now(); //B
				auto clean_beg = std::chrono::high_resolution_clock::now(); 

				// Get the k-mers of the window's sequences
				// Only keep the solid k-mers
				solidMerCounts.clear();
				for (std::pair<std::string, unsigned> p : merCounts) {
					if (p.second >= solidThresh) {
						solidMerCounts[p.first] = p.second;
					}
				}

				auto clean_end = std::chrono::high_resolution_clock::now(); //C
				
				// std::cerr << "adding took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
				// std::cerr << "cleaning took " << std::chrono::duration_cast<std::chrono::milliseconds>(clean_end - clean_beg).count() << " ms\n";
				// std::cerr << "total took " << std::chrono::duration_cast<std::chrono::milliseconds>(clean_end - c_start).count() << " ms\n";

				c_start = std::chrono::high_resolution_clock::now();
				correctedRegion = "";
				// if (correctedRegion.empty()) {
					// Compute frequent anchors to try to link again if first linking fails
				srcZone = correctedRead.substr(srcBeg, merSize + zone);
				dstZone = correctedRead.substr(dstBeg, merSize + zone);
					std::vector<std::string> srcList = getAnchors(merCounts, srcZone, merSize, zone + 1);
					std::vector<std::string> dstList = getAnchors(merCounts, dstZone, merSize, zone + 1);
					// std::cerr << "srcList.size() : " << srcList.size() << std::endl;
					anchors.clear();				
					for (std::string s : srcList) {
						for (std::string d : dstList) {
							// if (s != src or d != dst) {
								anchors.push_back(std::make_pair(s, d));
							// }
						}
					}
					std::unordered_map<std::string, std::vector<unsigned>> srcPos = getKMersPos(srcZone, merSize);
					std::unordered_map<std::string, std::vector<unsigned>> dstPos = getKMersPos(dstZone, merSize);
					c_end = std::chrono::high_resolution_clock::now();
					// std::cerr << "anchors took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

					// Attempt to link frequent anchors
					anchorNb = 0;
					auto loop_start = std::chrono::high_resolution_clock::now();
					// while (anchorNb < anchors.size() and minGap <= tmpDstBeg - tmpSrcEnd - 1 and tmpDstBeg - tmpSrcEnd - 1 <= maxGap and correctedRegion.empty()) {
					// std::cerr << "got : " << tmpDstBeg << "(" << dstBeg << ") ; " << tmpSrcEnd << "(" << srcEnd << ")" << std::endl;
					while (anchorNb < anchors.size() and correctedRegion.empty() and tmpDstBeg - tmpSrcEnd - 1 <= maxGap) {
						auto selec_start = std::chrono::high_resolution_clock::now();
						src = anchors[anchorNb].first;
						dst = anchors[anchorNb].second;
						tmpSrcBeg = srcBeg + srcPos[src][0];
						tmpSrcEnd = tmpSrcBeg + merSize - 1;
						tmpDstBeg = dstBeg + dstPos[dst][0];
						tmpDstEnd = tmpDstBeg + merSize - 1;
						auto selec_end = std::chrono::high_resolution_clock::now();
						// std::cerr << "select took " << std::chrono::duration_cast<std::chrono::milliseconds>(selec_end - selec_start).count() << " ms\n";
						
						auto c_start = std::chrono::high_resolution_clock::now();
						if (src != dst) {
							// visited.clear();
							curBranches = 0;
							dist = 0;
							curExt = src;
							correctedRegion = "";
							maxSize = 15.0 / 100.0 * 2.0 * (tmpDstBeg - tmpSrcEnd - 1) + (tmpDstBeg - tmpSrcEnd - 1) + merSize;
							link(solidMerCounts, src, dst, merSize, visited, &curBranches, dist, curExt, correctedRegion, merSize, maxSize, maxBranches, solidThresh, minOrder);
							// std::cerr << "linked" << std::endl;
						}
						anchorNb++;
						auto c_end = std::chrono::high_resolution_clock::now();
						std::cerr << "linking2 took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
					}
					auto loop_end = std::chrono::high_resolution_clock::now();
					std::cerr << "looping took " << std::chrono::duration_cast<std::chrono::milliseconds>(loop_end - loop_start).count() << " ms\n";
				// }

				c_start = std::chrono::high_resolution_clock::now();
				if (!correctedRegion.empty()) {
					corList.push_back(std::make_pair(correctedRead.substr(tmpSrcBeg, tmpDstEnd - tmpSrcBeg + 1), correctedRegion));
				}
				c_end = std::chrono::high_resolution_clock::now();
				std::cerr << "pushing took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
			}
			i = tmpDstBeg;
		} else {
			// break;
			auto c_start = std::chrono::high_resolution_clock::now();
			i = correctedRead.length();
			auto c_end = std::chrono::high_resolution_clock::now();
			std::cerr << "else took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
		}
		auto main_loop_end = std::chrono::high_resolution_clock::now();
		std::cerr << "mainloop took " << std::chrono::duration_cast<std::chrono::milliseconds>(main_loop_end - main_loop_beg).count() << " ms\n";
	}

	// std::cerr << "return" << std::endl;
	auto main_end = std::chrono::high_resolution_clock::now();
	std::cerr << "function took " << std::chrono::duration_cast<std::chrono::milliseconds>(main_end - main_start).count() << " ms\n";
	return corList;
}

void processRead(std::vector<Alignment>& alignments, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::string readId = alignments.begin()->qName;

	// Compute alignment piles
	std::unordered_map<std::string, std::string> sequences = getSequencesunordered_maps(alignments, readsDir);
	std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> pairPiles = getAlignmentPiles(alignments, minSupport, windowSize, windowOverlap, sequences);
	std::vector<std::vector<std::string>> piles = pairPiles.second;
	std::vector<std::pair<unsigned, unsigned>> pilesPos = pairPiles.first;

	auto c_start = std::chrono::high_resolution_clock::now();
	std::vector<std::unordered_map<std::string, unsigned>> pilesMers(piles.size());
	for (unsigned i = 0; i < piles.size(); i++) {
		pilesMers[i] = getKMersCounts(piles[i], merSize);
	}
	auto c_end = std::chrono::high_resolution_clock::now();
	std::cerr << "counting took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

	// Remove sequences sharing too few solid k-mers with their associated template
	auto c_startR = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<std::string>> oldPiles = piles;
	std::vector<std::unordered_map<std::string, unsigned>> oldPilesMers = pilesMers;
	unsigned i = 0;
	while (i < piles.size()) {
		piles[i] = removeBadSequencesPrev(piles[i], piles[i][0], pilesMers[i], merSize, commonKMers, solidThresh, windowSize);
		if (piles[i].size() < minSupport) {
			piles.erase(piles.begin() + i);
			pilesPos.erase(pilesPos.begin() + i);
			pilesMers.erase(pilesMers.begin() + i);
		} else {
			i++;
		}
	}
	auto c_endR = std::chrono::high_resolution_clock::now();
	std::cerr << "removing took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_endR - c_startR).count() << " ms\n";

	// Compute consensuses for all the piles
	auto c_start1 = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<std::string, std::string>> consensuses;
	std::vector<std::pair<std::string, std::string>> curCons;
	
	for (i = 0; i < piles.size(); i++) {
		curCons = computeConsensuses(readId, piles[i], pilesPos[i], pilesMers[i], readsDir, minSupport, merSize, commonKMers, solidThresh, windowSize);
		if (!curCons.empty()) {
			consensuses.push_back(curCons[0]);
		}
	}
	auto c_end1 = std::chrono::high_resolution_clock::now();
	std::cerr << "consensus took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end1 - c_start1).count() << " ms\n";

	// std::cerr << "computed all consensuses" << std::endl;
	
	// Align computed consensuses to the read and retrieve poorly supported regions windows positions
	c_start = std::chrono::high_resolution_clock::now();
	std::pair<std::string, std::vector<std::pair<std::pair<int, int>, int>>> p = alignConsensuses(readId, sequences[alignments[0].qName], consensuses, pilesPos, piles, 0);
	std::string correctedRead = p.first;
	std::vector<std::pair<std::pair<int, int>, int>> corPosPiles = p .second;
	c_end = std::chrono::high_resolution_clock::now();
	std::cerr << "anchoring1 took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

	// Drop read if it contains too many poorly supported bases
	// if (!dropRead(correctedRead)) {
	if (1) {
		// Polish poorly supported regions with local DBGs
		std::vector<std::pair<std::string, std::string>> corList, newList;

		// Polish small regions with a small k-mer size
		c_start = std::chrono::high_resolution_clock::now();
		corList = polishCorrection(correctedRead, corPosPiles, oldPiles, oldPilesMers, merSize, solidThresh, 1, 100);
		c_end = std::chrono::high_resolution_clock::now();
		std::cerr << "polishing took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

		// Anchor the polished regions to the read (replace the original region by its polished version). Non-optimal, should be done online.
		c_start = std::chrono::high_resolution_clock::now();
		for (std::pair<std::string, std::string> p : corList) {
			std::string r = p.first;
			std::string c = p.second;
			int b, l;
			b = correctedRead.find(r);
			l = r.length();
			if ((int) b != -1) {
				correctedRead.replace(b, l, c);
			}
		}
		c_end = std::chrono::high_resolution_clock::now();
		std::cerr << "anchoring took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

		// Trim the read (trims until a long sketch of corrected bases if found. Ex: aaaaCCggagtAttagGGACTTACGATCGATCGATCa => GGACTTACGATCGATCGATC)
		c_start = std::chrono::high_resolution_clock::now();
		std::vector<std::string> correctedSplits = trimRead(correctedRead, 100);
		int nbSplit = 0;
		while (nbSplit < correctedSplits.size()) {
			outMtx.lock();
			std::cout << ">" << readId << "_" << nbSplit + 1 << std::endl << correctedSplits[nbSplit] << std::endl;
			outMtx.unlock();
			nbSplit++;
		}
		// if (!correctedRead.empty()) {
		// 	outMtx.lock();
		// 	std::cout << ">" << readId << std::endl << correctedRead << std::endl;
		// 	outMtx.unlock();
		// }
		c_end = std::chrono::high_resolution_clock::now();
		std::cerr << "trimming took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
	}
}


void processReads(std::vector<std::vector<Alignment>>& reads, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::vector<std::pair<std::string, std::string>> consensuses;

	for (std::vector<Alignment> alignments : reads) {
		auto c_start = std::chrono::high_resolution_clock::now();
		processRead(alignments, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
		auto c_end = std::chrono::high_resolution_clock::now();
		std::cerr << "processing " << alignments[0].qName << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
	}

}

// Multithread ok, but running with GNU parallel for now.
// Multithreading seems to consume lots of resources when processing large number of reads.
void runCorrection(std::string alignmentFile, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads) {	
	std::ifstream f(alignmentFile);
	std::vector<Alignment> curReadAlignments;
	Alignment al;
	std::string curRead, line;
	curRead = "";
	int readNumber = 0;

	// Init threads
	std::vector<std::vector<std::vector<Alignment>>> reads(nbThreads);
	for (unsigned i = 0; i < nbThreads; i++) {
		reads[i] = std::vector<std::vector<Alignment>>();
	}

	getline(f, line);
	while(line.length() > 0 or !curReadAlignments.empty()) {
		if (line.length() > 0) {
			al = Alignment(line);
		}
		if (line.length() > 0 and (curRead == "" or al.qName == curRead)) {
			curRead = al.qName;
			curReadAlignments.push_back(al);
			getline(f, line);
		} else {
			std::sort(curReadAlignments.begin(), curReadAlignments.end());
			reads[readNumber % nbThreads].push_back(curReadAlignments);	
			readNumber++;
			curReadAlignments.clear();
			curRead = "";
		}
	}

	// Launch threads
	std::vector<std::future<void>> threads(nbThreads);
	for (unsigned i = 0 ; i < nbThreads ; i++) {
		std::vector<std::vector<Alignment>> als = reads[i];
		threads[i] = async(std::launch::async, [als, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap]() mutable {
			processReads(als, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
		});
	}
	
	// Get threads results
	for (std::future<void> &t: threads) {
		t.get();
	}
}
