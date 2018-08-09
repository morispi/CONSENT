#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include "LRSelfCorrection.h"
#include "utils.h"

std::mutex outMtx;

// Get positions of non-repeated k-mers in sequence (1-based positions)
// Positions of repeated k-mers are set to -1
std::map<std::string, int> getKMersPos(std::string sequence, int merSize) {
	std::map<std::string, int> mers;

	for (int i = 0; i < sequence.length() - merSize + 1; i++) {
		// mers.insert((sequence.substr(i, merSize)));
		if (mers[sequence.substr(i, merSize)] == 0) {
			mers[sequence.substr(i, merSize)] = i + 1;
		} else {
			mers[sequence.substr(i, merSize)] = -1;
		}
	}

	return mers;
}

std::map<std::string, int> getKMersCounts(std::vector<std::string> sequences, int merSize) {
	std::map<std::string, int> merCounts;
	int i;

	for (std::string s : sequences) {
		i = 0;
		while (s.length() >= merSize && i < s.length() - merSize + 1) {
			merCounts[(s.substr(i, merSize))]++;
			i++;
		}
	}

	return merCounts;
}

void removeBadSequences(std::vector<std::string>& sequences, std::string tplSeq, int merSize, int commonKMers, int solidThresh, int windowOverlap) {
	// Redefine solid threshold according to the size of the alignment pile
	// solidThresh = sequences.size() / solidThresh + 1;
	std::map<std::string, int> kMers = getKMersPos(tplSeq, merSize);
	std::map<std::string, int> merCounts = getKMersCounts(sequences, merSize);
	std::string curSeq;
	int i, j, c, pos, tplBegPos, tplEndPos, curSeqBegPos, curSeqEndPos;
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

std::pair<std::string, std::vector<std::string>> computeConsensus(std::vector<std::string>& sequences, std::string readsDir, int minSupport, int merSize, int commonKMers, int solidThresh, int windowOverlap) {
	std::string tplSeq = sequences[0];
	std::vector<std::string> res;

	removeBadSequences(sequences, tplSeq, merSize, commonKMers, solidThresh, windowOverlap);

	if (sequences.size() < minSupport) {
		return std::make_pair(tplSeq, res);
	}


	std::string POAinFile = readsDir + "../TMPRegions/" + std::tmpnam(nullptr) + ".fasta";
	// Repeat while filename exist, to avoid writing over an already used file
	while (fileExists(POAinFile)) {
		POAinFile = readsDir + "../TMPRegions/" + std::tmpnam(nullptr) + ".fasta";
	}
	const char* cPOAinFile = POAinFile.c_str();
	std::ofstream inPOA(POAinFile);

	std::string POAoutFile = readsDir + "../TMPConsensus/" + std::tmpnam(nullptr) + ".fasta";
	// Repeat while filename exist, to avoid writing over an already used file
	while (fileExists(POAoutFile)) {
		POAoutFile = readsDir + "../TMPConsensus/" + std::tmpnam(nullptr) + ".fasta";
	}
	const char* cPOAOutFile = POAoutFile.c_str();

	int regLen;
	int nbReg = 0;
	for (std::string s : sequences) {
		inPOA << ">" << nbReg << std::endl << s << std::endl;
		regLen = s.length();
		nbReg++;
	}
	inPOA.close();

	std::string cmd = "poaV2/poa -read_fasta " + POAinFile + " -lower -preserve_seqorder -hb -best -pir " + POAoutFile + " poaV2/blosum80.mat > /dev/null 2> /dev/null";
	const char *ccmd = cmd.c_str();
	auto c_start = std::chrono::high_resolution_clock::now();
	system(ccmd);
	auto c_end = std::chrono::high_resolution_clock::now();

	remove(cPOAinFile);

	res.clear();
	std::ifstream cons(POAoutFile);
	std::string consHeader, consSeq;
	while (getline(cons, consHeader)) {
		getline(cons, consSeq);
		res.push_back(consSeq);
	}
	cons.close();
	remove(cPOAOutFile);

	return std::make_pair(tplSeq, res);
}

std::vector<std::pair<int, int>> getAlignmentPilesPositions(int tplLen, std::vector<Alignment>& alignments, int minSupport, int windowSize, int overlappingWindows) {
	int* coverages = new int[tplLen];
	int i;
	for (i = 0; i < tplLen; i++) {
		coverages[i] = 1;
	}
	int beg, end;

	for (Alignment al : alignments) {
		beg = al.qStart;
		end = al.qEnd;

		for (i = beg; i <= end; i++) {
			coverages[i]++;
		}
	}

	std::vector<std::pair<int, int>> pilesPos;
	int prev = -1;

	int curLen = 0;
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

std::vector<std::pair<std::string, std::vector<std::string>>> processRead(std::vector<Alignment>& alignments, std::string readsDir, int minSupport, int windowSize, int merSize, int commonKMers, int solidThresh, int windowOverlap) {
	std::vector<std::pair<std::string, std::vector<std::string>>> consensuses;
	int min, max, sup, myStartPos, tplLen;
	min = -1;
	max = -1;
	sup = 0;
	tplLen = alignments.begin()->qLength;

	std::vector<std::pair<int, int>> pilesPos = getAlignmentPilesPositions(tplLen, alignments, minSupport, windowSize, windowOverlap);
	std::map<std::string, std::string> sequences = getSequencesMaps(alignments, readsDir);

	int beg, end, length, start, shift;
	std::vector<std::string> pilesSeqs;
	std::vector<std::string> curPile;

	for (std::pair<int, int> p : pilesPos) {
		pilesSeqs.clear();
		beg = p.first;
		end = p.second;
		length = end - beg + 1;
		
		// Insert template sequence
		pilesSeqs.push_back(sequences[alignments.begin()->qName].substr(beg, length));
		
		// Insert aligned sequences
		for (Alignment al : alignments) {
			// Check if the alignment spans the query window, and if the target has enough bases
			if (al.qStart <= beg and end <= al.qEnd and al.tStart + beg - al.qStart + length - 1 <= al.tEnd) {
				shift = beg - al.qStart;
				std::string tmpSeq = sequences[al.tName].substr(al.tStart, al.tEnd - al.tStart + 1);
				int rcd = 0;
				if (al.strand) {
					tmpSeq = reverseComplement(tmpSeq);
					rcd = 1;
				}
				tmpSeq = tmpSeq.substr(shift, length);
				pilesSeqs.push_back(tmpSeq);
			}
		}

		// Compute consensus for current pile
		// TODO: useless to use pilesSeqs ?
		curPile.clear();
		for (std::string s : pilesSeqs) {
			curPile.push_back(s);
		}
		if (curPile.size() >= minSupport) {
			std::pair<std::string, std::vector<std::string>> cons = computeConsensus(curPile, readsDir, minSupport, merSize, commonKMers, solidThresh, windowOverlap);
			if (cons.second.size() > 0) {
				consensuses.push_back(cons);
			}
		}
	}

	// Treat all alignments as one (ie exact read, aligned regions of other read, and fee this to POA)
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

	return consensuses;
}


void processReads(std::vector<std::vector<Alignment>>& reads, std::string readsDir, int minSupport, int windowSize, int merSize, int commonKMers, int solidThresh, int windowOverlap) {
	std::vector<std::pair<std::string, std::vector<std::string>>> consensuses;
	int i;

	for (std::vector<Alignment> alignments : reads) {
		auto c_start = std::chrono::high_resolution_clock::now();
		consensuses = processRead(alignments, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);

		std::ifstream curRead(readsDir + alignments[0].qName);
		std::string header, sequence;
		getline(curRead, header);
		getline(curRead, sequence);
		curRead.close();
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::tolower);

		// Align consensus to window, correct window, align corrected window to read, replace mapped zone in the read by corrected window
		std::string corWindow;
		i = 1;
		for (std::pair<std::string, std::vector<std::string>> c : consensuses) {
			// TODO: see what we do if multiple consensus were created
			std::string r = c.second[0];
			// Align consensus to window to get start and end position of the substring to extract from the consensus, as the corrected window
			std::pair<int, int> positions = NeedlemanWunschLocalAlignments(c.first, r);
			int beg, end;
			beg = positions.first;
			end = positions.second;
			corWindow = r.substr(beg, end - beg + 1);

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

		outMtx.lock();
		std::cout << header << std::endl << sequence << std::endl;
		outMtx.unlock();

		auto c_end = std::chrono::high_resolution_clock::now();
		std::cerr << "processing " << alignments[0].qName << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

	}

}

void runCorrection(std::string alignmentFile, std::string readsDir, int minSupport, int windowSize, int merSize, int commonKMers, int solidThresh, int windowOverlap, int nbThreads) {	
	std::ifstream f(alignmentFile);
	std::vector<Alignment> curReadAlignments;
	Alignment al;
	std::string curRead, line;
	curRead = "";
	int readNumber = 0;

	// Init threads
	std::vector<std::vector<std::vector<Alignment>>> reads;
	for (int i = 0; i < nbThreads; i++) {
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
	for (int i = 0 ; i < nbThreads ; i++) {
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