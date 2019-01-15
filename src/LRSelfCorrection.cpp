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
#include "../BMEAN/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h"

std::mutex outMtx;
std::map<std::string, std::vector<bool>> readIndex;

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
        // str.push_back('T');
      	str[j] = 'T';
      }else{
        // str.push_back('G');
        str[j] = 'G';
      }
    }else{
      if(num[i+1]){
        // str.push_back('C');
        str[j] = 'C';
      }else{
        // str.push_back('A');
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
	return (float) nbCorBases(correctedRead) / correctedRead.length() < 0.75;
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

// std::vector<std::string> trimRead(std::string correctedRead, unsigned merSize) {
// 	return {correctedRead};
// 	std::vector<std::string> res;
// 	unsigned beg, end, n;
// 	beg = 0;
// 	end = 0;
// 	unsigned i = 0;
// 	n = 0;

// 	while (i < correctedRead.length()) {
// 		while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
// 			i++;
// 		}
// 		beg = i;
// 		n = 0;

// 		while (i < correctedRead.length() and n < merSize) {
// 			if (!isUpperCase(correctedRead[i])) {
// 				n++;
// 			} else {
// 				n = 0;
// 			}
// 			i++;
// 		}

// 		end = i - n - 1;
// 		if (end >= beg) {
// 			std::string split = correctedRead.substr(beg, end - beg + 1);
// 			// std::cerr << "split : " << split << std::endl;
// 			if (!dropRead(split)) {
// 				res.push_back(split);
// 			}
// 		}
// 	}

// 	return res;
// }

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
	beg = i - merSize + 1;

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
	end = i + merSize - 1;

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

std::pair<std::string, std::unordered_map<kmer, unsigned>> computeConsensuses(std::string& readId, std::vector<std::string> & piles, std::pair<unsigned, unsigned>& pilesPos, std::string& readsDir, unsigned& minSupport, unsigned& merSize, unsigned& commonKMers, unsigned& solidThresh, unsigned& windowSize) {
	// return piles[0];
	// auto start_antoine = std::chrono::high_resolution_clock::now();
	// std::cerr << "go MSABMAAC" << std::endl;
	// std::cerr << "in : " << piles[0].length() <<  std::endl;
	int bmeanSup;
	// if (piles.size() > commonKMers) {
	// 	bmeanSup = commonKMers;
	// } else {
	// 	bmeanSup = piles.size() / 2;
	// }
	bmeanSup = std::min((int) commonKMers, (int) piles.size() / 2);
	std::pair<std::vector<std::vector<std::string>>, std::unordered_map<kmer, unsigned>> rOut = MSABMAAC(piles, merSize, bmeanSup, solidThresh);
	// std::cerr << "out : " << rOut.first[0][0].length() << std::endl;
	// std::cerr << "end MSABMAAC" << std::endl;
	auto result = rOut.first;
	auto merCounts = rOut.second;
	// auto end_antoine = std::chrono::high_resolution_clock::now();
	// outMtx.lock();
	// std::cerr << "antoine took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_antoine - start_antoine).count() << " ms\n";
	// outMtx.unlock();
	// auto c_start = std::chrono::high_resolution_clock::now();
	std::string corTpl = result[0][0];

	// Polish the consensus
	std::vector<std::pair<std::string, std::string>> corList;
	// auto c_start = std::chrono::high_resolution_clock::now();
	if (corTpl.length() >= merSize) {
		corTpl = weightConsensus(corTpl, piles, merCounts, merSize, windowSize, solidThresh);
	}
	// std::cerr << "go polish" << std::endl;
	corTpl = polishCorrection(corTpl, merCounts, merSize, solidThresh);
	// corTpl = trimRead(corTpl, merSize);
	// std::cerr << ">corTpl" << std::endl << corTpl << std::endl;
	// std::cerr << "end polish" << std::endl;
	// auto c_end = std::chrono::high_resolution_clock::now();
	// outMtx.lock();
	// std::cerr << "polishing took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
	// outMtx.unlock();
	// c_end = std::chrono::high_resolution_clock::now();
	// std::cerr << "voting took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";

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

// std::pair<int, int> mapConsensus(std::string read, std::string consensus, int merSize, int solidThresh, int beg, int end) {
// 	unsigned i, j;
// 	int resBeg = 0;
// 	int resEnd = 0;
// 	std::vector<localisation> posList;

// 	std::vector<std::string> vectRead(1);
// 	vectRead[0] = read;
// 	std::vector<std::string> vectCons(1);
// 	vectCons[0] = consensus;

// 	kmer2localisation read_index;
// 	std::unordered_map<kmer, unsigned> merCountsRead;
// 	kmer2localisation consensus_index;
// 	std::unordered_map<kmer, unsigned> merCountsConsensus;

// 	fill_index_kmers(vectRead, read_index, merSize, merCountsRead, solidThresh);
// 	fill_index_kmers(vectCons, consensus_index, merSize, merCountsConsensus, solidThresh);

// 	std::vector<kmer> kMersConsensus(consensus.length() - merSize + 1);
// 	for (i = 0; i < consensus.length() - merSize + 1; i++) {
// 		kMersConsensus[i] = str2num(consensus.substr(i, merSize));
// 	}

// 	// std::cerr << "cons are : " << std::endl;
// 	// std::cerr << read.substr(0, merSize) << std::endl;
// 	// std::cerr << consensus.substr(0, merSize) << std::endl;
// 	// Find head anchor
// 	i = 0;
// 	while (i < kMersConsensus.size() and resBeg == 0) {
// 		posList = read_index[kMersConsensus[i]];
// 		if (posList.size() != 0) {
// 			j = 0;
// 			while (j < posList.size() and resBeg == 0) {
// 				if (posList[j].position >= beg) {//} and posList[j].position <= beg + merSize) {
// 					resBeg = posList[j].position;
// 					// std::cerr << "set resBeg : " << posList[j].position << std::endl;
// 				} else {
// 					// std::cerr << "update j beg" << std::endl;
// 					j++;
// 				}
// 			}
// 		}
// 		i++;
// 	}

// 	// Find tail anchor
// 	i = kMersConsensus.size() - 1;
// 	while (i >= 0 and resEnd == 0) {
// 		posList = read_index[kMersConsensus[i]];
// 		if (posList.size() != 0) {
// 			j = 0;
// 			while (j < posList.size() and resEnd == 0) {
// 				if (posList[j].position + merSize - 1 <= end){//} and posList[j].position + merSize - 1 >= end - merSize) {
// 					// std::cerr << "set resEnd : " << posList[j].position << std::endl;
// 					resEnd = posList[j].position + merSize - 1;
// 				} else {
// 					// std::cerr << "update j end" << std::endl;
// 					j++;
// 				}
// 			}
// 		}
// 		i--;
// 	}

// 	return std::make_pair(resBeg, resEnd);
// }

std::string alignConsensuses(std::string rawRead, std::string sequence, std::vector<std::string>& consensuses, std::vector<std::unordered_map<kmer, unsigned>>& merCounts, std::vector<std::pair<unsigned, unsigned>>& pilesPos, std::vector<std::vector<std::string>>& piles, int startPos, unsigned windowSize, unsigned windowOverlap, unsigned solidThresh, unsigned merSize) {
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	StripedSmithWaterman::Alignment subAlignment;	
	int32_t maskLen = 15;
	unsigned beg, end, oldEnd;
	oldEnd = 0;
	std::string outSequence;
	outSequence = sequence;
	std::transform(outSequence.begin() + startPos, outSequence.end(), outSequence.begin() + startPos, ::tolower);

	std::string corWindow;
	unsigned i = 0;
	std::string tmpSequence, consUp;
	int curPos = std::max(0, (int) startPos - (int) windowOverlap);
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
		curCons = consensuses[i];
		if (curCons.size() < merSize) {
			continue;
		}
		// std::cerr << ">cons" << std::endl << curCons << std::endl;
		// std::cerr << "antoine : " << curCons << std::endl;
		// std::cerr << "cons size : " << curCons.length() << std::endl;
		curMers = merCounts[i];
		if (curPos + windowSize + 2 * windowOverlap >= outSequence.length()) {
			sizeAl = outSequence.length() - curPos;
		} else {
			sizeAl = windowSize + 2 * windowOverlap;
		}
		// auto al_start = std::chrono::high_resolution_clock::now();
		aligner.Align(curCons.c_str(), outSequence.c_str() + curPos, sizeAl, filter, &alignment, maskLen);
		// auto al_end = std::chrono::high_resolution_clock::now();
		// outMtx.lock();
		// std::cerr << "aligning took " << std::chrono::duration_cast<std::chrono::milliseconds>(al_end - al_start).count() << " ms\n";
		// outMtx.unlock();

		// std::string mapSequence = outSequence;
		// std::transform(mapSequence.begin(), mapSequence.end(), mapSequence.begin(), ::toupper);
		
		// auto map_start = std::chrono::high_resolution_clock::now();
		// std::pair<int, int> newRes = mapConsensus(mapSequence, curCons, merSize, solidThresh, curPos, curPos + sizeAl - 1);
		// auto map_end = std::chrono::high_resolution_clock::now();
		// outMtx.lock();
		// std::cerr << "mapping took " << std::chrono::duration_cast<std::chrono::milliseconds>(map_end - map_start).count() << " ms\n";
		// outMtx.unlock();
		// int newBeg = newRes.first;
		// int newEnd = newRes.second;

		beg = alignment.ref_begin + curPos;
		end = alignment.ref_end + curPos;
		// beg = newBeg;
		// end = newEnd;

		// std::cerr << "begs : " << beg << " ; " << newBeg << std::endl;
		// std::cerr << "ends : " << end << " ; " << newEnd << std::endl;
		// std::cerr << std::endl;

		// std::cerr << "1" << std::endl;
		curCons = curCons.substr(alignment.query_begin, alignment.query_end - alignment.query_begin + 1);
		// std::cerr << "ok" << std::endl;
		// std::cerr << "curCons : " << curCons << std::endl;
		// Check if windows overlap, and if they do, chose the best subsequence
		if (i != 0) {
		// if (0) {
			if (oldEnd >= beg and oldCons.length() >= oldEnd - beg + 1) {
 				overlap = oldEnd - beg + 1;
 				// std::cerr << "2" << std::endl;
				seq1 = oldCons.substr(oldCons.length() - 1 - overlap + 1, overlap);
				seq2 = curCons.substr(0, overlap);
				// std::cerr << "ok" << std::endl;
				if (toUpperCase(seq1, 0, seq1.length() - 1) != toUpperCase(seq2, 0, seq2.length() - 1)) {
					// std::cerr << "overlap : " << overlap << std::endl;
					if (overlap >= merSize) {
						// std::cerr << "ici" << std::endl;
						solidMersSeq1 = nbSolidMers(seq1, oldMers, merSize, solidThresh);
						solidMersSeq2 = nbSolidMers(seq2, curMers, merSize, solidThresh);
						// std::cerr << "ok" << std::endl;
					} else {
						// std::cerr << "la" << std::endl;
						solidMersSeq1 = nbUpperCase(seq1);
						solidMersSeq2 = nbUpperCase(seq2);
						// std::cerr << "ok" << std::endl;
					}
					if (solidMersSeq1 > solidMersSeq2) {
						aligner.Align(seq1.c_str(), seq2.c_str(), std::min(seq1.length(), seq2.length()), filter, &subAlignment, maskLen);
						indels = getIndels(subAlignment.cigar_string);
						ins = indels.first;
						del = indels.second;
						if (overlap - ins + del < curCons.length()) {
							// std::cerr << "3" << std::endl;
							curCons = seq1 + curCons.substr(overlap - ins + del);
							// std::cerr << "ok" << std::endl;
						} else {
							curCons = "";
						}
					}
				}
			}
		}

		if (curCons != "") {
			consUp = curCons;
			std::transform(consUp.begin(), consUp.end(), consUp.begin(), ::toupper);
			outSequence.replace(beg, end - beg + 1, consUp);
			curPos = (beg + consUp.length() - 1) - windowOverlap;
			oldCons = curCons;
			oldMers = merCounts[i];
			// oldBeg = beg;
			oldEnd = beg + consUp.length() - 1;
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
		// std::cerr << "candidatesSrc : " << srcZone.substr(i, merSize) << std::endl;
	}
	// Same with the dst zone
	std::vector<std::string> candidatesDst(dstZone.size() - merSize + 1);
	for (i = 0; i < dstZone.size() - merSize + 1; i++) {
		candidatesDst[i] = dstZone.substr(i, merSize);
		// std::cerr << "candidatesDst : " << dstZone.substr(i, merSize) << std::endl;
	}

	// Add the anchors pairs to the result vector, without allowing repeated k-mers
	for (std::string csrc : candidatesSrc) {
		if (mersPosSrc[csrc].size() == 1) {
			for (std::string cdst : candidatesDst) {
				if (mersPosDst[cdst].size() == 1) {
					res.push_back(std::make_pair(csrc, cdst));
					// std::cerr << "push back : " << csrc << " ; " << cdst << std::endl;
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

	// std::cerr << "getAnchors, res.size() = " << res.size() << std::endl;

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
	// int zone = 4;
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
	// std::cerr << "go skip" << std::endl;
	while (i < correctedRead.length() and !isUpperCase(correctedRead[i])) {
		i++;
	}
	// std::cerr << "end skip" << std::endl;

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
		// std::cerr << "go src dst" << std::endl;
		srcEnd = getNextSrc(correctedRead, i, merSize + zone);
		dstEnd = getNextDst(correctedRead, srcEnd + 1, merSize + zone);
		srcBeg = srcEnd - merSize - zone + 1;
		dstBeg = dstEnd - merSize - zone + 1;
		// std::cerr << "end src dst" << std::endl;
		// std::cerr << "srcEnd : " << srcEnd << std::endl;
		// std::cerr << "dstEnd : " << dstEnd << std::endl;

		// Polish the poorly supported region region if 2 anchors were found
		if (srcEnd != -1 and dstEnd != -1) {
			correctedRegion = "";
			srcZone = correctedRead.substr(srcBeg, merSize + zone);
			dstZone = correctedRead.substr(dstBeg, merSize + zone);
			// std::cerr << "go anchor pos" << std::endl;
			// std::cerr << "srcZone : " << srcZone << std::endl;
			// std::cerr << "dstZone : " << dstZone << std::endl;
			anchors = getAnchors(merCounts, srcZone, dstZone, merSize, 5);
			srcPos = getKMersPos(srcZone, merSize);
			dstPos = getKMersPos(dstZone, merSize);
			// std::cerr << "end anchor pos" << std::endl;
			// std::cerr << "anchors size : " << anchors.size() << std::endl;
			// std::cerr << "correctedRegion : " << correctedRegion << std::endl;

			// Attempt to link frequent anchors
			anchorNb = 0;
			while (anchorNb < anchors.size() and correctedRegion.empty()) {
				// std::cerr << "go other part" << std::endl;
				src = anchors[anchorNb].first;
				dst = anchors[anchorNb].second;
				tmpSrcBeg = srcBeg + srcPos[src][0];
				tmpSrcEnd = tmpSrcBeg + merSize - 1;
				tmpDstBeg = dstBeg + dstPos[dst][0];
				tmpDstEnd = tmpDstBeg + merSize - 1;
				
				// auto c_start = std::chrono::high_resolution_clock::now();
				if (src != dst) {
					curBranches = 0;
					dist = 0;
					curExt = src;
					correctedRegion = "";
					maxSize = 15.0 / 100.0 * 2.0 * (tmpDstBeg - tmpSrcEnd - 1) + (tmpDstBeg - tmpSrcEnd - 1) + merSize;
					// std::cerr << "go link" << std::endl;
					link(merCounts, src, dst, merSize, visited, &curBranches, dist, curExt, correctedRegion, merSize, maxSize, maxBranches, solidThresh, merSize);
					// std::cerr << "end link" << std::endl;
				}
				anchorNb++;
				// auto c_end = std::chrono::high_resolution_clock::now();
				// std::cerr << "linking2 took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
			}

			if (!correctedRegion.empty()) {
				// corList.push_back(std::make_pair(correctedRead.substr(tmpSrcBeg, tmpDstEnd - tmpSrcBeg + 1), correctedRegion));

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

	// Anchor the polished region to the corrected read
	// int b, l;
	// std::string r, c;
	// for (std::pair<std::string, std::string> p : corList) {
	// 	r = p.first;
	// 	c = p.second;
	// 	b = correctedRead.find(r);
	// 	l = r.length();
	// 	if ((int) b != -1) {
	// 		correctedRead.replace(b, l, c);
	// 	}
	// }
	// std::transform(correctedRead.begin(), correctedRead.end(), correctedRead.begin(), ::toupper);
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

std::pair<std::string, std::string> processRead(int id, std::vector<Alignment>& alignments, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
	std::string readId = alignments.begin()->qName;

	// Compute alignment piles
	// auto c_start = std::chrono::high_resolution_clock::now();
	// std::unordered_map<std::string, std::str/ing> sequences = getSequencesunordered_maps(alignments, readsDir);
	std::unordered_map<std::string, std::string> sequences = getSequencesMap(alignments);
	// std::cerr << "got maps" << std::endl;
	std::pair<std::vector<std::pair<unsigned, unsigned>>, std::vector<std::vector<std::string>>> pairPiles = getAlignmentPiles(alignments, minSupport, windowSize, windowOverlap, sequences);
	// std::cerr << "got piles : " << pairPiles.first.size() << std::endl;
	std::vector<std::vector<std::string>> piles = pairPiles.second;
	std::vector<std::pair<unsigned, unsigned>> pilesPos = pairPiles.first;
	// auto c_end = std::chrono::high_resolution_clock::now();
	// std::cerr << "init took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
	unsigned i = 0;

	// Compute consensuses for all the piles
	// auto c_start1 = std::chrono::high_resolution_clock::now();
	std::pair<std::string, std::unordered_map<kmer, unsigned>> resCons;
	std::vector<std::string> consensuses(piles.size());
	std::vector<std::unordered_map<kmer, unsigned>> merCounts(piles.size()); 
	// std::vector<std::pair<std::string, std::string>> curCons;
	for (i = 0; i < piles.size(); i++) {
		// std::cerr << "compute cons : " << i << std::endl;
		resCons = computeConsensuses(readId, piles[i], pilesPos[i], readsDir, minSupport, merSize, commonKMers, solidThresh, windowSize);
		consensuses[i] = resCons.first;
		merCounts[i] = resCons.second;
		// std::cerr << "got it " << std::endl;
	}
	// auto c_end1 = std::chrono::high_resolution_clock::now();
	// std::cerr << "consensus took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end1 - c_start1).count() << " ms\n";

	// Align computed consensuses to the read
	// auto c_start = std::chrono::high_resolution_clock::now();
	// std::cerr << "go align" << std::endl;
	std::string correctedRead = alignConsensuses(readId, sequences[alignments[0].qName], consensuses, merCounts, pilesPos, piles, 0, windowSize, windowOverlap, solidThresh, merSize);
	// std::cerr << "end align" << std::endl;
	// auto c_end = std::chrono::high_resolution_clock::now();
	// outMtx.lock();
	// std::cerr << "anchoring took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
	// outMtx.unlock();

	// Drop read if it contains too many poorly supported bases
	// if (!dropRead(correctedRead)) {
	if (1) {
		// Split the read if it contains uncorrected windows
		// c_start = std::chrono::high_resolution_clock::now();
		// std::vector<std::string> correctedSplits = trimRead(correctedRead, windowSize);
		// unsigned nbSplit = 0;
		// while (nbSplit < correctedSplits.size()) {
		// 	outMtx.lock();
		// 	std::cout << ">" << readId << "_" << nbSplit + 1 << std::endl << correctedSplits[nbSplit] << std::endl;
		// 	outMtx.unlock();
		// 	nbSplit++;
		// }
		// std::cerr << "processed " << readId << std::endl;
		return std::make_pair(readId, correctedRead);

		// outMtx.lock();
		// std::cout << ">" << readId << std::endl << correctedRead << std::endl;
		// outMtx.unlock();
	}
}


// void processReads(std::vector<std::vector<Alignment>>& reads, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap) {
// 	std::vector<std::pair<std::string, std::string>> consensuses;

// 	for (std::vector<Alignment> alignments : reads) {
// 		auto c_start = std::chrono::high_resolution_clock::now();
// 		std::pair<std::string, std::string> correctedRead = processRead(alignments, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
// 		outMtx.lock();
// 		std::cout << ">" << correctedRead.first << std::endl << correctedRead.second << std::endl;
// 		outMtx.unlock();
// 		auto c_end = std::chrono::high_resolution_clock::now();
// 		std::cerr << "processing " << alignments[0].qName << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count() << " ms\n";
// 	}
// }

std::map<std::string, std::vector<bool>> indexReads(std::string readsFile) {
	std::ifstream f(readsFile);
	std::string header, sequence;
	std::map<std::string, std::vector<bool>> res;

	getline(f, header);
	while (header.length() > 0) {
		header.erase(0, 1);
		getline(f, sequence);
		std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		res[header] = fullstr2num(sequence);
		getline(f, header);
	}
	return res;
}

// Multithread ok, but running with GNU parallel for now as BOA still has memory leaks
// void runCorrection(std::string alignmentFile, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile) {	
// 	std::ifstream f(alignmentFile);
// 	std::vector<Alignment> curReadAlignments;
// 	Alignment al;
// 	std::string curRead, line;
// 	curRead = "";
// 	int readNumber = 0;

// 	readIndex = indexReads(readsFile);

// 	// Init threads
// 	std::vector<std::vector<std::vector<Alignment>>> reads(nbThreads);
// 	for (unsigned i = 0; i < nbThreads; i++) {
// 		reads[i] = std::vector<std::vector<Alignment>>();
// 	}

// 	getline(f, line);
// 	while(line.length() > 0 or !curReadAlignments.empty()) {
// 		if (line.length() > 0) {
// 			al = Alignment(line);
// 		}
// 		if (line.length() > 0 and (curRead == "" or al.qName == curRead)) {
// 			curRead = al.qName;
// 			curReadAlignments.push_back(al);
// 			getline(f, line);
// 		} else {
// 			std::sort(curReadAlignments.begin(), curReadAlignments.end());
// 			reads[readNumber % nbThreads].push_back(curReadAlignments);	
// 			readNumber++;
// 			curReadAlignments.clear();
// 			curRead = "";
// 		}
// 	}

// 	// Launch threads
// 	std::vector<std::future<void>> threads(nbThreads);
// 	for (unsigned i = 0 ; i < nbThreads ; i++) {
// 		std::vector<std::vector<Alignment>> als = reads[i];
// 		threads[i] = async(std::launch::async, [als, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap]() mutable {
// 			processReads(als, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
// 		});
// 	}
	
// 	// Get threads results
// 	for (std::future<void> &t: threads) {
// 		t.get();
// 	}
// }

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
			// std::cerr << (float) al.resMatches / al.alBlockLen << std::endl;
			curReadAlignments.push_back(al);
			getline(f, line);
		} else {
			// std::cerr << curRead << " : " << curReadAlignments.size() << std::endl;
			// std::sort(curReadAlignments.begin(), curReadAlignments.end(), cmpById);
			// std::vector<Alignment> tmp;
			// for (int i = 0; i < curReadAlignments.size() and i < 50; i++) {
			// 	tmp.push_back(curReadAlignments[i]);
			// }
			// curReadAlignments = tmp;
			std::sort(curReadAlignments.begin(), curReadAlignments.end());
			if (!f.eof()) {
				f.seekg(-line.length()-1, f.cur);
			}
			// std::cerr << "curReadAlignments : " << curReadAlignments.size() << std::endl;
			return curReadAlignments;
			// reads[readNumber % nbThreads].push_back(curReadAlignments);	
			// readNumber++;
			// curReadAlignments.clear();
			// curRead = "";
		}
	}

	// std::cerr << "curReadAlignments : " << curReadAlignments.size() << std::endl;
	return curReadAlignments;
}

void runCorrection(std::string alignmentFile, std::string readsDir, unsigned minSupport, unsigned windowSize, unsigned merSize, unsigned commonKMers, unsigned solidThresh, unsigned windowOverlap, unsigned nbThreads, std::string readsFile, unsigned nbReads) {
	std::ifstream f(alignmentFile);
	std::vector<Alignment> curReadAlignments;
	std::string curRead, line;
	curRead = "";
	// int readNumber = 0;

	readIndex = indexReads(readsFile);

	int poolSize = 10000;
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
        results[i] = myPool.push(processRead, curReadAlignments, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
        // std::cerr << "pushed job" << std::endl;
        jobsLoaded++;
	}

	// std::cerr << "end loop" << std::endl;

	// Load the remaining jobs as other jobs terminate
	int curJob = 0;
    std::pair<std::string, std::string> curRes;
    while(!f.eof() && jobsLoaded < jobsToProcess) {
    	// std::cerr << "got job" << std::endl;
    	// Get the job results
        curRes = results[curJob].get();
        std::cout << ">" << curRes.first << std::endl << curRes.second << std::endl;
        jobsCompleted++;
        
        // Load the next job
        curReadAlignments = getNextReadPile(f);
        while (curReadAlignments.size() == 0 and !f.eof()) {
        	curReadAlignments = getNextReadPile(f);
        }
        results[curJob] = myPool.push(processRead, curReadAlignments, readsDir, minSupport, windowSize, merSize, commonKMers, solidThresh, windowOverlap);
        jobsLoaded++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

	// std::cerr << "other loop" << std::endl;

	// Wait for the remaining jobs to terminate
	while(jobsCompleted < jobsLoaded) {
        // Get the job results
        // std::cerr << "we are looping" << std::endl;
        curRes = results[curJob].get();
        // std::cerr << "we looped" << std::endl;
        std::cout << ">" << curRes.first << std::endl << curRes.second << std::endl;
        jobsCompleted++;
        
        // Increment the current job nb, and loop if needed
        curJob++;
        if(curJob == poolSize) {
            curJob = 0;
        }
	}

}
