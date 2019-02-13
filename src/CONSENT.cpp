#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <chrono>
#include "CONSENT.h"
#include "kMersProcessing.h"
#include "DBG.h"
#include "../BMEAN/bmean.h"
#include "../BMEAN/utils.h"
#include "../BMEAN/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h"

std::mutex outMtx;
std::map<std::string, std::vector<bool>> readIndex;
bool doTrimRead = true;

std::vector<bool> fullstr2num(const string& str){
  std::vector<bool> res;
  for(uint i(0);i<str.size();i++){
    switch (str[i]){
      case 'A':res.push_back(false);res.push_back(false);break;
      case 'C':res.push_back(false);res.push_back(true);break;
      case 'G':res.push_back(true);res.push_back(false);break;
      default:res.push_back(true);res.push_back(true);break;
    }
  }
  return res;
}

std::string fullnum2str(vector<bool> num){
  string str(num.size()/2, 'N');
  uint j = 0;
  for(uint i(0);i<num.size();i+=2){
    if(num[i]){
      if(num[i+1]){
      	str[j] = 'T';
      }else{
        str[j] = 'G';
      }
    }else{
      if(num[i+1]){
        str[j] = 'C';
      }else{
        str[j] = 'A';
      }
    }
    j++;
  }
  return str;
}

bool isUpperCase(char c) {
	return 'A' <= c and c <= 'Z';
}

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
	return (float) nbCorBases(correctedRead) / correctedRead.length() < 0.1;
}


std::string toUpperCase(std::string& s, int beg, int end) {
	std::string res = s;
	std::locale loc;
	for (int i = beg; i <= end; i++) {
		res[i] = std::toupper(res[i], loc);
	}

	return res;
}

std::string toLowerCase(std::string& s, int beg, int end) {
	std::string res = s;
	std::locale loc;
	for (int i = beg; i <= end; i++) {
		res[i] = std::tolower(res[i], loc);
	}

	return res;
}

std::string trimRead(std::string correctedRead, unsigned merSize) {
	unsigned beg, end, n;
	unsigned i;
	i = 0;
	n = 0;
	while (i < correctedRead.length() and n < merSize) {
		if (isUpperCase(correctedRead[i])) {
			n++;
		} else {
			n = 0;
		}
		i++;
	}
	beg = i - merSize;

	i = correctedRead.length() - 1;
	n = 0;
	while (i >= 0 and n < merSize) {
		if (isUpperCase(correctedRead[i])) {
			n++;
		} else {
			n = 0;
		}
		i--;
	}
	end = i + merSize;

	if (end > beg) {
		return correctedRead.substr(beg, end - beg + 1);
	} else {
		return "";
	}
}

std::string weightConsensus(std::string& consensus, std::vector<std::string>& pile, std::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, unsigned windowSize, unsigned solidThresh) {
	std::vector<std::string> splits;
	std::string curSplit;

	std::string header = "";
	std::string sequence = "";
	std::string curFct;

	unsigned i = 0;
	while (i < consensus.length() - merSize + 1) {
		curFct = consensus.substr(i, merSize);
		curFct = toUpperCase(curFct, 0, merSize);
		if (merCounts[str2num(curFct)] >= solidThresh) {
			consensus = toUpperCase(consensus, i, i + merSize - 1);
		} else {
			consensus = toLowerCase(consensus, i, i + merSize - 1);
		}
		i++;
	}

	return consensus;
}

std::pair<std::string, std::unordered_map<kmer, unsigned>> computeConsensuses(std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& minAnchors, unsigned& solidThresh, unsigned& windowSize, unsigned maxMSA, std::string path) {
	int bmeanSup;
	bmeanSup = std::min((int) commonKMers, (int) piles.size() / 2);
	std::pair<std::vector<std::vector<std::string>>, std::unordered_map<kmer, unsigned>> rOut = MSABMAAC(piles, merSize, bmeanSup, solidThresh, minAnchors, maxMSA, path);

	if (rOut.first.size() == 0) {
		return std::make_pair("", rOut.second);
	}
	auto result = rOut.first;
	auto merCounts = rOut.second;
	std::string corTpl = result[0][0];

	// Polish the consensus
	std::vector<std::pair<std::string, std::string>> corList;
	if (corTpl.length() >= merSize) {
		corTpl = weightConsensus(corTpl, piles, merCounts, merSize, windowSize, solidThresh);
		corTpl = polishCorrection(corTpl, merCounts, merSize, solidThresh);
	}

	return std::make_pair(corTpl, merCounts);
}

int nbSolidMers(std::string seq, std::unordered_map<kmer, unsigned> merCounts, unsigned merSize, unsigned solidThresh) {
	int nb = 0;
	for (unsigned i = 0; i < seq.length() - merSize + 1; i++) {
		if (merCounts[str2num(seq.substr(i, merSize))] >= solidThresh) {
			nb++;
		}
	}

	return nb;
}

int nbUpperCase(std::string s) {
	int nb = 0;
	for (unsigned i = 0; i < s.length(); i++) {
		if (isUpperCase(s[i])) {
			nb++;
		}
	}

	return nb;
}

std::pair<int, int> getIndels(std::string cigar){
	int ins = 0;
	int del = 0;
	int current = 0;
	for(unsigned i = 0; i < cigar.length(); i++){
		if('0' <= cigar[i] && cigar[i] <= '9'){
			current = (current * 10) + (cigar[i] - '0');
		} else {
			if (cigar[i] == 'I') {
				ins += current;
			} else if (cigar[i] == 'D') {
				del += current;
			}
			current = 0;
		}
	}
	return std::make_pair(ins, del);
}

std::string alignConsensuses(std::string rawRead, std::string sequence, std::vector<std::string>& consensuses, std::vector<std::unordered_map<kmer, unsigned>>& merCounts, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::string>& templates, int startPos, unsigned windowSize, unsigned windowOverlap, unsigned solidThresh, unsigned merSize) {
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	StripedSmithWaterman::Alignment subAlignment;	
	int32_t maskLen = 15;
	unsigned beg, end, oldEnd;
	oldEnd = 0;
	std::string outSequence;
	outSequence = sequence;
	std::transform(outSequence.begin(), outSequence.end(), outSequence.begin(), ::tolower);

	std::string corWindow;
	unsigned i = 0;
	std::string tmpSequence, consUp;
	int curPos = startPos;
	int sizeAl;
	std::string curCons, oldCons;
	std::unordered_map<kmer, unsigned> oldMers;
	std::unordered_map<kmer, unsigned> curMers;
	unsigned overlap;
	std::string seq1, seq2;
	int solidMersSeq1, solidMersSeq2;
	std::pair<int, int> indels;
	unsigned ins, del;

	for (i = 0; i < consensuses.size(); i++) {
		if (consensuses[i].length() < merSize) {
			curCons = templates[i];
			std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::tolower);
		} else {
			curCons = consensuses[i];
		}
		curMers = merCounts[i];
		
		curPos = std::max(0, (int) curPos - (int) windowOverlap);
		if (curPos + windowSize + 2 * windowOverlap >= outSequence.length()) {
			sizeAl = outSequence.length() - curPos;
		} else {
			sizeAl = windowSize + 2 * windowOverlap;
		}

		aligner.Align(curCons.c_str(), outSequence.c_str() + curPos, sizeAl, filter, &alignment, maskLen);
		beg = alignment.ref_begin + curPos;
		end = alignment.ref_end + curPos;
		curCons = curCons.substr(alignment.query_begin, alignment.query_end - alignment.query_begin + 1);

		// Check if alignment positions overlap the previous window. If they do, chose the best subsequence
		if (i != 0 and oldEnd >= beg) {
 			overlap = oldEnd - beg + 1;
			if (consensuses[i].length() >= merSize and oldCons.length() >= overlap and curCons.length() >= overlap) {
				seq1 = oldCons.substr(oldCons.length() - 1 - overlap + 1, overlap);
				seq2 = curCons.substr(0, overlap);
				if (toUpperCase(seq1, 0, seq1.length() - 1) != toUpperCase(seq2, 0, seq2.length() - 1)) {
					if (overlap >= merSize) {
						solidMersSeq1 = nbSolidMers(seq1, oldMers, merSize, solidThresh);
						solidMersSeq2 = nbSolidMers(seq2, curMers, merSize, solidThresh);
					} else {
						solidMersSeq1 = nbUpperCase(seq1);
						solidMersSeq2 = nbUpperCase(seq2);
					}
					if (solidMersSeq1 > solidMersSeq2) {
						aligner.Align(seq1.c_str(), seq2.c_str(), std::min(seq1.length(), seq2.length()), filter, &subAlignment, maskLen);
						indels = getIndels(subAlignment.cigar_string);
						ins = indels.first;
						del = indels.second;
						if (overlap - ins + del < curCons.length()) {
							curCons = seq1 + curCons.substr(overlap - ins + del);
						} else {
							curCons = "";
						}
					}
				}
			}
		}

		if (curCons != "") {
			if (consensuses[i].length() >= merSize){//} and curCons.length() >= merSize) {
				consUp = curCons;
				std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::toupper);
				outSequence.replace(beg, end - beg + 1, consUp);
			}
			curPos = beg + curCons.length() ;
			oldCons = curCons;
			oldMers = merCounts[i];
			oldEnd = beg + curCons.length() - 1;
		}
	}

	return outSequence;
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


// Anchors without repeated k-mers
std::vector<std::pair<std::string, std::string>> getAnchors(std::unordered_map<kmer, unsigned>& merCounts, std::string srcZone, std::string dstZone, unsigned merSize, unsigned nb) {
	std::vector<std::pair<std::string, std::string>> res;
	unsigned i;

	std::unordered_map<std::string, std::vector<unsigned>> mersPosSrc = getKMersPos(srcZone, merSize);
	std::unordered_map<std::string, std::vector<unsigned>> mersPosDst = getKMersPos(dstZone, merSize);

	// Consider all k-mers of the src zone as potential anchors
	std::vector<std::string> candidatesSrc(srcZone.size() - merSize + 1);
	for (i = 0; i < srcZone.size() - merSize + 1; i++) {
		candidatesSrc[i] = srcZone.substr(i, merSize);
	}
	// Same with the dst zone
	std::vector<std::string> candidatesDst(dstZone.size() - merSize + 1);
	for (i = 0; i < dstZone.size() - merSize + 1; i++) {
		candidatesDst[i] = dstZone.substr(i, merSize);
	}

	// Add the anchors pairs to the result vector, without allowing repeated k-mers
	for (std::string csrc : candidatesSrc) {
		if (mersPosSrc[csrc].size() == 1) {
			for (std::string cdst : candidatesDst) {
				if (mersPosDst[cdst].size() == 1) {
					res.push_back(std::make_pair(csrc, cdst));
				}
			}
		}
	}

	// Sort the anchors vector in ascending order of the number of occurrences of (src + dst)
	std::sort(res.begin(), res.end(),
		[&merCounts](std::pair<std::string, std::string>& r1, std::pair<std::string, std::string>& r2) {
			int occ1 = merCounts[str2num(r1.first)] + merCounts[str2num(r1.second)];
			int occ2 = merCounts[str2num(r2.first)] + merCounts[str2num(r2.second)];
			return occ1 > occ2;
		}
	);

	std::vector<std::pair<std::string, std::string>> finalRes;
	for (i = 0; i < nb and i < res.size(); i++) {
		finalRes.push_back(res[i]);
	}

	return finalRes;
}

std::string polishCorrection(std::string correctedRead, std::unordered_map<kmer, unsigned>& merCounts, unsigned merSize, int solidThresh) {
	std::set<std::string> visited;
	unsigned curBranches;
	unsigned dist;
	std::string curExt;
	std::string correctedRegion;
	unsigned maxSize;
	unsigned maxBranches = 50;
	std::vector<std::pair<std::string, std::string>> corList;
	int zone = 3;
	int srcBeg, srcEnd, dstBeg, dstEnd;
	unsigned tmpSrcBeg = 0, tmpSrcEnd = 0, tmpDstBeg = 0, tmpDstEnd = 0;
	std::string src, dst;
	std::pair<int, int> pos;
	std::vector<std::pair<std::string, std::string>> anchors;
	unsigned anchorNb;
	std::string srcZone, dstZone;
	std::unordered_map<std::string, std::vector<unsigned>> srcPos, dstPos;
	std::string oldCorrectedRead;
	int b, l;
	std::string r, c;

	// Skip uncorrected head of the read
	unsigned i = 0;
	while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
		i++;
	}

	if (i > 0 and i < correctedRead.length() and correctedRead.length() - i >= merSize) {
		int extLen = i;
		oldCorrectedRead = correctedRead;
		correctedRead = correctedRead.substr(i);
		int extSize = extendLeft(merCounts, merSize, extLen, correctedRead, solidThresh);
		if (extSize < extLen) {
			correctedRead = oldCorrectedRead.substr(0, extLen - extSize) + correctedRead;
			i = i - (extLen - extSize);
		}
	}

	// Search for poorly supported regions bordered by solid corrected regions
	while (i < correctedRead.length()) {
		srcEnd = getNextSrc(correctedRead, i, merSize + zone);
		dstEnd = getNextDst(correctedRead, srcEnd + 1, merSize + zone);
		srcBeg = srcEnd - merSize - zone + 1;
		dstBeg = dstEnd - merSize - zone + 1;

		// Polish the poorly supported region region if 2 anchors were found
		if (srcEnd != -1 and dstEnd != -1) {
			correctedRegion = "";
			srcZone = correctedRead.substr(srcBeg, merSize + zone);
			dstZone = correctedRead.substr(dstBeg, merSize + zone);
			anchors = getAnchors(merCounts, srcZone, dstZone, merSize, 5);
			srcPos = getKMersPos(srcZone, merSize);
			dstPos = getKMersPos(dstZone, merSize);

			// Attempt to link frequent anchors
			anchorNb = 0;
			while (anchorNb < anchors.size() and correctedRegion.empty()) {
				src = anchors[anchorNb].first;
				dst = anchors[anchorNb].second;
				tmpSrcBeg = srcBeg + srcPos[src][0];
				tmpSrcEnd = tmpSrcBeg + merSize - 1;
				tmpDstBeg = dstBeg + dstPos[dst][0];
				tmpDstEnd = tmpDstBeg + merSize - 1;
				
				if (src != dst) {
					curBranches = 0;
					dist = 0;
					curExt = src;
					correctedRegion = "";
					maxSize = 15.0 / 100.0 * 2.0 * (tmpDstBeg - tmpSrcEnd - 1) + (tmpDstBeg - tmpSrcEnd - 1) + merSize;
					link(merCounts, src, dst, merSize, visited, &curBranches, dist, curExt, correctedRegion, merSize, maxSize, maxBranches, solidThresh, merSize);
				}
				anchorNb++;
			}

			if (!correctedRegion.empty()) {
				// Anchor the correction to the read
				r = correctedRead.substr(tmpSrcBeg, tmpDstEnd - tmpSrcBeg + 1);
				c = correctedRegion;
				b = correctedRead.find(r);
				l = r.length();
				if ((int) b != -1) {
					correctedRead.replace(b, l, c);
					i = b;
				} else {
					i = tmpDstBeg > i ? tmpDstBeg : dstBeg;	
				}
			} else {
				i = tmpDstBeg > i ? tmpDstBeg : dstBeg;
			}
		} else {
			i = correctedRead.length();	
		}
	}

	i = correctedRead.length() - 1;
	while (i > 0 and !isUpperCase(correctedRead[i])) {
		i--;
	}

	if (i > 0 and i < correctedRead.length() - 1 and i + 1 >= merSize) {
		int extLen = correctedRead.length() - 1 - i;
		oldCorrectedRead = correctedRead;
		correctedRead = correctedRead.substr(0, i + 1);
		int extSize = extendRight(merCounts, merSize, extLen, correctedRead, solidThresh);
		if (extSize < extLen) {
			correctedRead = correctedRead + oldCorrectedRead.substr(oldCorrectedRead.length() - (extLen - extSize), extLen - extSize);
		}
	}

	return correctedRead;
}

std::unordered_map<std::string, std::string> getSequencesMap(std::vector<Alignment>& alignments) {
	std::unordered_map<std::string, std::string> sequences;
	std::string header, seq;

	// Insert template sequence
	sequences[alignments.begin()->qName] = fullnum2str(readIndex[alignments.begin()->qName]);

	// Insert aligned sequences
	for (Alignment al : alignments) {
		if (sequences[al.tName] == "") {
			sequences[al.tName] = fullnum2str(readIndex[al.tName]);
		}
	}

	return sequences;
}

std::pair<std::string, std::string> processRead(int id, std::vector<Alignment>& alignments, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors,unsigned solidThresh, unsigned windowOverlap, unsigned maxMSA, std::string path) {
	std::string readId = alignments.begin()->qName;
	std::unordered_map<std::string, std::string> sequences = getSequencesMap(alignments);
	std::vector<std::pair<unsigned, unsigned>> pilesPos = getAlignmentPilesPositions(alignments.begin()->qLength, alignments, minSupport, maxSupport, windowSize, windowOverlap);
	if (pilesPos.size() == 0) {
		return std::make_pair(readId, "");
	}
	unsigned i = 0;

	// Compute consensuses for all the piles
	std::pair<std::string, std::unordered_map<kmer, unsigned>> resCons;
	std::vector<std::string> consensuses(pilesPos.size());
	std::vector<std::unordered_map<kmer, unsigned>> merCounts(pilesPos.size()); 
	std::vector<std::string> curPile;
	std::vector<std::string> templates(pilesPos.size());
	for (i = 0; i < pilesPos.size(); i++) {
		curPile = getAlignmentPileSeq(alignments, minSupport, windowSize, windowOverlap, sequences, pilesPos[i].first, pilesPos[i].second, merSize, maxSupport);
		templates[i] = curPile[0];
		resCons = computeConsensuses(readId, curPile, pilesPos[i], minSupport, merSize, commonKMers, minAnchors, solidThresh, windowSize, maxMSA, path);
		if (resCons.first.length() < merSize) {
			consensuses[i] = resCons.first;
		} else {
			consensuses[i] = resCons.first;
		}
		merCounts[i] = resCons.second;
	}

	// Align computed consensuses to the read
	std::string correctedRead = alignConsensuses(readId, sequences[alignments[0].qName], consensuses, merCounts, pilesPos, templates, pilesPos[0].first, windowSize, windowOverlap, solidThresh, merSize);

	// Trim read if need (ie when performing correction), and drop it if it contains too many uncorrected bases
	if (doTrimRead) {
		if (!dropRead(correctedRead)) {
			correctedRead = trimRead(correctedRead, 1);
			return std::make_pair(readId, correctedRead);
		} else {
			return std::make_pair(readId, "");
		}
	} else {
		return std::make_pair(readId, correctedRead);
	}
}

void indexReads(std::map<std::string, std::vector<bool>>& index, std::string readsFile) {
	std::ifstream f(readsFile);
	std::string header, sequence;

	getline(f, header);
	while (header.length() > 0) {
		header.erase(0, 1);
		getline(f, sequence);
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		index[header] = fullstr2num(sequence);
		getline(f, header);
	}
}

std::vector<Alignment> getNextReadPile(std::ifstream& f) {
	std::vector<Alignment> curReadAlignments;
	Alignment al;
	std::string line, curRead;

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
			if (!f.eof()) {
				f.seekg(-line.length()-1, f.cur);
			}
			return curReadAlignments;
		}
	}

	return curReadAlignments;
}

void runCorrection(std::string alignmentFile, unsigned minSupport, unsigned maxSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned minAnchors, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile, std::string proofFile, unsigned maxMSA, std::string path) {
	std::ifstream f(alignmentFile);
	std::vector<Alignment> curReadAlignments;
	std::string curRead, line;
	curRead = "";

	indexReads(readIndex, readsFile);
	if (proofFile != "") {
		indexReads(readIndex, proofFile);
		doTrimRead = false;
	}


	int poolSize = 1000;
	ctpl::thread_pool myPool(nbThreads);
	int jobsToProcess = 100000000;
	int jobsLoaded = 0;
	int jobsCompleted = 0;

	// Load the first jobs
	vector<std::future<std::pair<std::string, std::string>>> results(poolSize);
    for(int i = 0; i < poolSize && !f.eof() && jobsLoaded < jobsToProcess; i++) {
        curReadAlignments = getNextReadPile(f);
        while (curReadAlignments.size() == 0 and !f.eof()) {
        	curReadAlignments = getNextReadPile(f);
        }
        results[i] = myPool.push(processRead, curReadAlignments, minSupport, maxSupport, windowSize, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, maxMSA, path);
        jobsLoaded++;
	}

	// Load the remaining jobs as other jobs terminate
	int curJob = 0;
    std::pair<std::string, std::string> curRes;
    while(!f.eof() && jobsLoaded < jobsToProcess) {
    	// Get the job results
        curRes = results[curJob].get();
        if (curRes.second.length() != 0) {
	        std::cout << ">" << curRes.first << std::endl << curRes.second << std::endl;
	    }
        jobsCompleted++;
        
        // Load the next job
        curReadAlignments = getNextReadPile(f);
        while (curReadAlignments.size() == 0 and !f.eof()) {
        	curReadAlignments = getNextReadPile(f);
        }
        results[curJob] = myPool.push(processRead, curReadAlignments, minSupport, maxSupport, windowSize, merSize, commonKMers, minAnchors, solidThresh, windowOverlap, maxMSA, path);
        jobsLoaded++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

	// Wait for the remaining jobs to terminate
	while(jobsCompleted < jobsLoaded) {
        // Get the job results
        curRes = results[curJob].get();
        if (curRes.second.length() != 0) {
	        std::cout << ">" << curRes.first << std::endl << curRes.second << std::endl;
	    }
        jobsCompleted++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

}
