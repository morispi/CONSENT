#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include "LRSelfCorrection.h"
#include "localAlignment.h"
#include "kMersProcessing.h"
// #include "reverseComplement.h"
// #include "alignmentPiles.h"

std::mutex outMtx;
int totalTime = 0;
int totalIterations = 0;
int totalComputeConsensus = 0;

void removeBadSequences(std::vector<std::string>& sequences, std::string tplSeq, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::map<std::string, int> kMers = getKMersPos(tplSeq, merSize);
	std::map<std::string, int> merCounts = getKMersCounts(sequences, merSize);
	std::string curSeq;
	unsigned i, j, c;
	int pos, tplBegPos, tplEndPos, curSeqBegPos, curSeqEndPos;
	std::set<std::string> anchored;
	std::string mer;
	std::map<int, std::vector<std::string>> mapCommonMers;

	// Discriminate by number of share solid k-mers with template sequence
	i = 1;
	while (i < sequences.size()) {
		curSeq = sequences[i];
		j = 0;
		c = 0;
		anchored.clear();
		pos = 0;
		tplBegPos = -1;
		tplEndPos = -1;
		curSeqBegPos = -1;
		curSeqEndPos = -1;


		while (j < curSeq.length() - merSize + 1) {
			mer = (curSeq.substr(j, merSize));
			if (merCounts[mer] >= solidThresh && anchored.find(mer) == anchored.end() and kMers[mer] > 0 and (pos == 0 or kMers[mer] >= pos + merSize)) {
				pos = kMers[mer];
				if (tplBegPos == -1) {
					tplBegPos = pos - 1;
					curSeqBegPos = j;
				}
				c++;
				j += merSize;
				anchored.insert(mer);
				tplEndPos = pos - 1 + merSize - 1;
				curSeqEndPos = j - 1;
			} else {
				j += 1;
			}
		}

		if (c >= commonKMers) {
			int beg, end;
			beg = std::max(0, curSeqBegPos - tplBegPos);
			end = std::min(curSeqEndPos + tplSeq.length() - tplEndPos - 1, curSeq.length() - 1);
			std::string tmpSeq = sequences[i].substr(beg, end + 1); 

			mapCommonMers[c].push_back(tmpSeq);
		}

		i++;
	}

	// Insert sequences from the alignment pile into the result vector according to their number of shared k-mers with the template sequence
	// (useful to compute consensus quicker with POA)
	std::vector<std::string> newSequences;
	newSequences.push_back(sequences[0]);
	for (std::map<int, std::vector<std::string>>::reverse_iterator it = mapCommonMers.rbegin(); it != mapCommonMers.rend(); it++) {
		for (std::string s : it->second) {
			newSequences.push_back(s);
		}
	}

	sequences.clear();
	sequences = newSequences;

}

std::vector<std::pair<std::string, std::string>> computeConsensuses(std::string tplId, std::vector<std::vector<std::string>>& piles, std::string readsDir, unsigned minSupport, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::vector<std::pair<std::string, std::string>> res;

	// Remove sequences with too few shared solid k-mers with their associated window
	unsigned i = 0;
	while (i < piles.size()) {
		removeBadSequences(piles[i], piles[i][0], merSize, commonKMers, solidThresh, windowOverlap);
		if (piles[i].size() < minSupport) {
			piles.erase(piles.begin() + i);
		} else {
			i++;
		}
	}

	// Prepare POA input file
	std::string POAinFile = readsDir + "../TMPRegions/" + tplId;
	const char* cPOAinFile = POAinFile.c_str();
	std::ofstream inPOA(POAinFile);
	int nbReg = 0;
	for (std::vector<std::string> p : piles) {
		inPOA << p.size() << std::endl;
		for (std::string s : p) {
			inPOA << ">" << nbReg << std::endl << s << std::endl;
			nbReg++;
		}
	}
	inPOA.close();

	// Run POA
	std::string POAoutFile = readsDir + "../TMPConsensus/" + tplId;
	const char* cPOAOutFile = POAoutFile.c_str();

	std::string cmd = "poaV2/poa -read_fasta " + POAinFile + " -lower -preserve_seqorder -hb -best -pir " + POAoutFile + " poaV2/blosum80.mat > /dev/null 2> /dev/null";
	const char *ccmd = cmd.c_str();
	auto c_start = std::chrono::high_resolution_clock::now();
	system(ccmd);
	auto c_end = std::chrono::high_resolution_clock::now();

	remove(cPOAinFile);

	// Write POA results to a file
	res.clear();
	std::ifstream cons(POAoutFile);
	std::string consHeader, consSeq;
	i = 0;
	while (getline(cons, consHeader)) {
		getline(cons, consSeq);
		res.push_back(std::make_pair(piles[i][0], consSeq));
		i++;
	}
	cons.close();
	remove(cPOAOutFile);

	c_end = std::chrono::high_resolution_clock::now();
	totalTime += std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count();

	return res;
}

std::pair<std::string, std::string> alignConsensuses(std::vector<Alignment>& alignments, std::string readsDir, std::vector<std::pair<std::string, std::string>>& consensuses) {
	std::ifstream curRead(readsDir + alignments[0].qName);
	std::string header, sequence;
	getline(curRead, header);
	getline(curRead, sequence);
	curRead.close();
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::tolower);

	// Align consensus to window, correct window, align corrected window to read, replace mapped zone in the read by corrected window
	std::string corWindow;
	unsigned i = 0;
	for (i = 0; i < consensuses.size(); i++) {
		std::pair<std::string, std::string> c = consensuses[i];
		// Align consensus to window to get start and end position of the substring to extract from the consensus, as the corrected window
		std::pair<int, int> positions = NeedlemanWunschLocalAlignments(c.first, c.second);
		int beg, end;
		beg = positions.first;
		end = positions.second;
		corWindow = c.second.substr(beg, end - beg + 1);

		// Align corrected window to read, to get start and end position of the substring to replace in the read
		std::string tmpCorWindow, tmpSequence;
		tmpCorWindow = corWindow;
		tmpSequence = sequence;
		std::transform(tmpSequence.begin(), tmpSequence.end(), tmpSequence.begin(), ::tolower);
		std::transform(tmpCorWindow.begin(), tmpCorWindow.end(), tmpCorWindow.begin(), ::tolower);
		positions = NeedlemanWunschLocalAlignments(tmpCorWindow, tmpSequence);
		beg = positions.first;
		end = positions.second;
		sequence.replace(beg, end - beg + 1, corWindow);
	}

	return std::make_pair(header, sequence);
}

void processRead(std::vector<Alignment>& alignments, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::vector<std::vector<std::string>> piles = getAlignmentPiles(alignments, minSupport, windowSize, windowOverlap, readsDir);

	std::vector<std::pair<std::string, std::string>> consensuses = computeConsensuses(alignments.begin()->qName, piles, readsDir, minSupport, merSize, commonKMers, solidThresh, windowOverlap);

	std::pair<std::string, std::string> correctedRead = alignConsensuses(alignments, readsDir, consensuses);

	outMtx.lock();
	std::cout << correctedRead.first << std::endl << correctedRead.second << std::endl;
	outMtx.unlock();

	// Treat all alignments at once (ie exact read, aligned regions of other reads, and feed this to POA)
	// curPile.clear();
	// curPile.push_back(sequences[alignments.begin()->qName]);
	// for (Alignment al : alignments) {
	// 	std::string tmpSeq = sequences[al.tName].substr(al.tStart, al.tEnd - al.tStart + 1);
	// 	if (al.strand) {
	// 		tmpSeq = reverseComplement(tmpSeq);
	// 	}
	// 	curPile.push_back(tmpSeq);
	// }
	// std::pair<std::string, std::vector<std::string>> cons = computeConsensus(curPile, readsDir, minSupport, merSize, commonKMers, solidThresh, windowOverlap);
	// for (std::string s : cons.second) {
	// 	std::cerr << cons.first << std::endl << s << std::endl << std::endl;
	// }
	// consensuses.push_back(cons);
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