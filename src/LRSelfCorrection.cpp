#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include "LRSelfCorrection.h"
#include "localAlignment.h"
#include "kMersProcessing.h"

std::mutex outMtx;
int totalTime = 0;
int totalIterations = 0;
int totalComputeConsensus = 0;

bool compareLen(const std::string& a, const std::string& b) {
    return (a.size() > b.size()); 
}

void removeBadSequences(std::vector<std::string>& sequences, std::string tplSeq, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize) {
	unsigned tplLen = tplSeq.length();
	std::map<std::string, std::vector<int>> kMers = getKMersPos(tplSeq, merSize);
	std::map<std::string, int> merCounts = getKMersCounts(sequences, merSize);
	std::string curSeq;
	unsigned i, j, c;
	int pos, tplBegPos, tplEndPos, curSeqBegPos, curSeqEndPos;
	std::set<std::string> anchored;
	std::string mer;
	std::map<int, std::vector<std::string>> mapCommonMers;

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
		while (j < curSeq.length() - merSize + 1) {
			mer = (curSeq.substr(j, merSize));
			bool found = false;
			unsigned p = 0;
			// Check if the current k-mer (of the current sequence) appears in the template sequence after the current position
			while (!found and p < kMers[mer].size()) {
				found = pos == -1 or kMers[mer][p] > pos;
				p++;
			}
			// Allow repeated k-mers
			// if (merCounts[mer] >= solidThresh and found) {
			// Non-repeated k-mers only
			if (merCounts[mer] >= solidThresh and found and anchored.find(mer) == anchored.end() and kMers[mer].size() == 1) {
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
			int beg, end;
            beg = std::max(0, curSeqBegPos - tplBegPos);
            end = std::min(curSeqEndPos + tplSeq.length() - tplEndPos - 1, curSeq.length() - 1);
            std::string tmpSeq = sequences[i].substr(beg, end + 1); 

			// Only add the substring from the leftmost anchor to the rightmost anchor to the MSA (anchor = shared k-mer wth the template)
			// if (tmpSeq.length() > windowSize / 2) {
				mapCommonMers[c].push_back(tmpSeq);
			// }
			// Add the whole sequence to the MSA
			// mapCommonMers[c].push_back(sequences[i]);
		}

		i++;
	}

	// Insert sequences from the alignment pile into the result vector according to their number of shared k-mers with the template sequence
	// (useful to compute consensus quicker with POA). If two sequences have the same numer of shared k-mer, sort them according to their length.
	// Both sort are done in descending order.
	std::vector<std::string> newSequences;
	newSequences.push_back(sequences[0]);
	std::vector<std::string> tmpNewSeqs;
	for (std::map<int, std::vector<std::string>>::reverse_iterator it = mapCommonMers.rbegin(); it != mapCommonMers.rend(); it++) {
		tmpNewSeqs.clear();
		for (std::string s : it->second) {
			tmpNewSeqs.push_back(s);
		}
		std::sort(tmpNewSeqs.begin(), tmpNewSeqs.end(), compareLen);
		for (std::string s : tmpNewSeqs) {
			newSequences.push_back(s);
		}
	}

	sequences.clear();
	sequences = newSequences;
}

std::vector<std::string> splitConsensus(std::string consensus, std::vector<std::string> pile, int merSize) {
	std::vector<std::string> splits;
	std::string curSplit;
	std::map<std::string, int> merCounts = getKMersCounts(pile, merSize);

	std::string header = "";
	std::string sequence = "";

	int i = 0;
	while (i < consensus.size() - merSize + 1) {
		header += std::to_string(merCounts[consensus.substr(i, merSize)]) + "_";
		if (sequence.length() == 0) {
			sequence = consensus.substr(i, merSize);
		} else {
			sequence += consensus[i + merSize - 1];
		}
		if (merCounts[consensus.substr(i, merSize)] < 2) {
			if (curSplit.length() >= 2 * merSize) {
				splits.push_back(curSplit);
			}
			curSplit = "";
			// i += merSize - 1;
		} else {
			if (curSplit.length() == 0) {
				curSplit = consensus.substr(i, merSize);
			} else {
				curSplit += consensus[i + merSize - 1];
			}
		}
		i++;
	}

	if (curSplit.length() >= 2 * merSize) {
		splits.push_back(curSplit);
	}

	// std::cerr << ">" << header << std::endl << sequence << std::endl;
	std::cerr << "splits size : " << splits.size() << std::endl;

	for (std::string s : splits) {
		std::cerr << s << std::endl;
	}

	return splits;
}

std::vector<std::pair<std::string, std::string>> computeConsensuses(std::string tplId, std::vector<std::vector<std::string>>& piles, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::string readsDir, unsigned minSupport, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowSize) {
	std::vector<std::pair<std::string, std::string>> res;

	// Remove sequences sharing too few solid k-mers with their associated template
	unsigned i = 0;
	while (i < piles.size()) {
		removeBadSequences(piles[i], piles[i][0], merSize, commonKMers, solidThresh, windowSize);
		if (piles[i].size() < minSupport) {
			piles.erase(piles.begin() + i);
			pilesPos.erase(pilesPos.begin() + i);
		} else {
			i++;
		}
	}

	// Prepare POA input file
	std::string POAinFile = readsDir + "../POAInput/" + tplId;
	const char* cPOAinFile = POAinFile.c_str();
	std::ofstream inPOA(POAinFile);
	int nbReg = 0;
	for (std::vector<std::string> p : piles) {
		// int mySize = std::min(3, (int) p.size());
		int mySize = (int) p.size();
		inPOA << mySize << std::endl;
		std::cerr << mySize << " poaFileSize " << std::endl;
		for (unsigned i = 0; i < mySize; i++) {
			std::string s = p[i];
			inPOA << ">" << nbReg << std::endl << s << std::endl;
			// std::cerr << ">" << nbReg << std::endl << s << std::endl;
			nbReg++;
		}
	}
	inPOA.close();

	// Run POA
	std::string POAoutFile = readsDir + "../POAOutput/" + tplId;
	const char* cPOAOutFile = POAoutFile.c_str();
	std::string cmd = "poaV2/poa -read_fasta " + POAinFile + " -lower -hb -best -pir " + POAoutFile + " poaV2/blosum80.mat > /dev/null 2> /dev/null";
	const char *ccmd = cmd.c_str();
	auto c_start = std::chrono::high_resolution_clock::now();
	system(ccmd);
	auto c_end = std::chrono::high_resolution_clock::now();
	remove(cPOAinFile);

	// Store templates corrections in the result vector
	res.clear();
	std::ifstream cons(POAoutFile);
	std::string consHeader = "";
	std::string consSeq = "";
	i = 0;
	while (getline(cons, consHeader)) {
		getline(cons, consSeq);
		// Align the computed consensus to the template, to retrieve the (raw, corrected) pair of subsequences
		std::cerr << "consSeq size : " << consSeq.length() << std::endl;
		std::pair<std::pair<int, int>, std::pair<int, int>> posPair = NeedlemanWunschLocalAlignments(piles[i][0], consSeq);
		std::pair<int, int> rawPos = posPair.first;
		std::pair<int, int> corPos = posPair.second;
		std::string corTpl = consSeq.substr(corPos.first, corPos.second - corPos.first + 1);
		// std::string rawTpl = piles[i][0].substr(rawPos.first, rawPos.second - rawPos.first + 1);
		std::string rawTpl = piles[i][0];

		// std::cerr << "on raw : " << rawPos.first << " ; " << rawPos.second << std::endl;
		// std::cerr << "on cor : " << corPos.first << " ; " << corPos.second << std::endl;
		std::cerr << ">";
		std::map<std::string, int> merCounts = getKMersCounts(piles[i], merSize);
		for (int j = 0; j < corTpl.size() - merSize + 1; j++) {
			std::cerr << merCounts[corTpl.substr(j, merSize)] << "_";
		}
		std::cerr << std::endl << corTpl << std::endl;
		// Only store the corrected template (storing the whole consensus would cause overcorrection of some regions)
		res.push_back(std::make_pair(rawTpl, corTpl));

		// Split the corrected template if it contains k-mers that do not appear in the pile (as they have a high chance to be erroneous)
		// std::vector<std::string> splits = splitConsensus(corTpl, piles[i], merSize);

		// // Align the splits of the corrected template to the original template to get their raw sequences
		// for (std::string s : splits) {
		// 	// std::cerr << corSplit.length() << std::endl;
		// 	std::pair<std::pair<int, int>, std::pair<int, int>> posPair = NeedlemanWunschLocalAlignments(s, rawTpl);
		// 	std::pair<int, int> rawPos = posPair.second;
		// 	std::pair<int, int> corPos = posPair.first;
		// 	std::string rawSplit = rawTpl.substr(rawPos.first, rawPos.second - rawPos.first + 1);
		// 	// std::string corSplit = s.substr(corPos.first, corPos.second - corPos.first + 1);
		// 	std::string corSplit = s;

		// 	if (rawPos.second < windowSize - 1) {
		// 		rawTpl = rawTpl.substr(rawPos.second + 1);
		// 	} else {
		// 		rawTpl = "";
		// 	}

		// 	std::cerr << "beg : " << rawPos.first << " ; " << "end : " << rawPos.second << std::endl;

		// 	// Store pairs of (rawSplit, corSplit) subsequences to the result vector
		// 	res.push_back(std::make_pair(rawSplit, corSplit));
		// 	std::cerr << "rawSplit : " << rawSplit << std::endl;
		// 	std::cerr << "corSplit : " << corSplit << std::endl;
		// 	std::cerr << ">corSplit" << std::endl << corSplit << std::endl;
		// }
		
		i++;
	}
	cons.close();
	remove(cPOAOutFile);

	c_end = std::chrono::high_resolution_clock::now();
	totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count();

	return res;
}

std::vector<std::pair<std::string, std::string>> alignConsensuses(std::vector<Alignment>& alignments, std::string readsDir, std::vector<std::pair<std::string, std::string>>& consensuses, std::vector<std::pair<unsigned, unsigned>> pilesPos) {
	std::vector<std::pair<std::string, std::string>> correctedReads;

	// for (unsigned k = 0; k < alignments.size(); k++) {
	for (unsigned k = 0; k < 1; k++) {
 		std::ifstream curRead(readsDir + alignments[k].qName);
		std::string header, sequence;
		getline(curRead, header);
		getline(curRead, sequence);
		curRead.close();
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::tolower);

		// Anchor the corrected templates on the read
		std::string corWindow;
		unsigned i = 0;
		for (i = 0; i < consensuses.size(); i++) {
			std::pair<std::string, std::string> c = consensuses[i];
			// Don't proceed if no consensus was produced by POA
			if (!c.second.empty()) {
				// Replace the template by its correction on the read
				int beg, end;
				std::transform(c.first.begin(), c.first.end(), c.first.begin(), ::tolower);
				std::string tmpSequence = sequence;
				std::transform(tmpSequence.begin(), tmpSequence.end(), tmpSequence.begin(), ::tolower);
				beg = (int) tmpSequence.find(c.first);
				end = beg + c.first.length() - 1;
				sequence.replace(beg, end - beg + 1, c.second);

				// Anchor the correction on the read by aligning
				// std::transform(c.second.begin(), c.second.end(), c.second.begin(), ::tolower);
				// std::string tmpSequence = sequence;
				// std::transform(tmpSequence.begin(), tmpSequence.end(), tmpSequence.begin(), ::tolower);
				// std::pair<std::pair<int, int>, std::pair<int, int>> posPair = NeedlemanWunschLocalAlignments(c.second, tmpSequence);
				// std::pair<int, int> rawPos = posPair.second;
				// std::pair<int, int> corPos = posPair.first;
				// std::cerr << "raw : " << rawPos.first << " ; " << rawPos.second << std::endl;
				// std::cerr << "cor : " << corPos.first << " ; " << corPos.second << std::endl;
				// sequence.replace(rawPos.first, rawPos.second - rawPos.first + 1, c.second);
				// std::cerr << std::endl;
			}
		}

		correctedReads.push_back(std::make_pair(header, sequence));
	}

	return correctedReads;
}

unsigned int edit_distance(const std::string& s1, const std::string& s2) {
	const std::size_t len1 = s1.size(), len2 = s2.size();
	std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

	d[0][0] = 0;
	for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
	for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

	for(unsigned int i = 1; i <= len1; ++i)
		for(unsigned int j = 1; j <= len2; ++j)
                      // note that std::min({arg1, arg2, arg3}) works only in C++11,
                      // for C++98 use std::min(std::min(arg1, arg2), arg3)
                      d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
	return d[len1][len2];
}

void patchConsensuses(std::vector<std::pair<std::string, std::string>>& consensuses, std::vector<std::vector<std::string>>& piles, int merSize) {
	int i = 0;
	while (i < consensuses.size()) {
	// for (int i = 0; i < consensuses.size(); i++) {
		std::map<std::string, int> merCounts = getKMersCounts(piles[i], merSize);
		std::string curCons = consensuses[i].second;
		
		if (!curCons.empty()) {
			int j = 0;
			bool allSolid = true;
			while (j < curCons.size() - merSize + 1 and allSolid) {
				std::string curMer = curCons.substr(j, merSize);
				if (merCounts[curMer] < 1) {
					allSolid = false;
				}
				j++;
			}

			if (allSolid) {
				outMtx.lock();
				std::cerr << ">";
				for (int j = 0; j < curCons.size() - merSize + 1; j++) {
					std::cerr << merCounts[curCons.substr(j, merSize)] << "_";
				}
				std::cerr << std::endl;
				std::cerr << curCons << std::endl;
				outMtx.unlock();
				i++;
			} else {
				// consensuses.erase(consensuses.begin() + i);
				consensuses[i].second = "";
			}
		} else {
			i++;
		}
	}
}

void saveConsensuses(std::vector<std::pair<std::string, std::string>> correctedReads, std::string readsDir) {
	// for (int i = 0; i < correctedReads.size(); i++) {
	for (int i = 0; i < 1; i++) {
		correctedReads[i].first.erase(0,1);
		std::ofstream f(readsDir + correctedReads[i].first);
		f << ">" << correctedReads[i].first << std::endl << correctedReads[i].second << std::endl;
		f.close();
	}
}

void processRead(std::vector<Alignment>& alignments, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> pairPiles = getAlignmentPiles(alignments, minSupport, windowSize, windowOverlap, readsDir);
	std::vector<std::vector<std::string>> piles = pairPiles.second;
	std::vector<std::pair<unsigned, unsigned>> pilesPos = pairPiles.first;

	std::vector<std::pair<std::string, std::string>> consensuses = computeConsensuses(alignments.begin()->qName, piles, pilesPos, readsDir, minSupport, merSize, commonKMers, solidThresh, windowSize);

	// patchConsensuses(consensuses, piles, merSize);

	// std::cerr << "piles.size() = " << piles.size() << std::endl;
	// std::cerr << "consensuses.size() = " << consensuses.size() << std::endl;
	// std::cerr << "pilesPos.size() = " << pilesPos.size() << std::endl;

	auto c_start = std::chrono::high_resolution_clock::now();
	std::vector<std::pair<std::string, std::string>>correctedReads = alignConsensuses(alignments, readsDir, consensuses, pilesPos);
	auto c_end = std::chrono::high_resolution_clock::now();
	// std::cerr << "aligning consensuses " << alignments[0].qName << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

	// saveConsensuses(correctedReads, readsDir);
	outMtx.lock();
	std::cout << correctedReads[0].first << std::endl << correctedReads[0].second << std::endl;
	outMtx.unlock();
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

void runCorrection(std::string alignmentFile, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads) {	
	std::ifstream f(alignmentFile);
	std::vector<Alignment> curReadAlignments;
	Alignment al;
	std::string curRead, line;
	curRead = "";
	int readNumber = 0;

	// Init threads
	std::vector<std::vector<std::vector<Alignment>>> reads;
	for (unsigned i = 0; i < nbThreads; i++) {
		reads.push_back(std::vector<std::vector<Alignment>>());
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
			reads[readNumber % nbThreads].push_back(curReadAlignments);	
			readNumber++;
			curReadAlignments.clear();
			curRead = "";
		}
	}

	// Launch threads
	std::vector<std::future<void>> threads;
	for (unsigned i = 0 ; i < nbThreads ; i++) {
		std::vector<std::vector<Alignment>> als = reads[i];
		threads.push_back(async(std::launch::async, [als, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap]() mutable {
			processReads(als, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
		}));
	}
	
	// Get threads results
	for (std::future<void> &t: threads) {
		t.get();
	}
}
