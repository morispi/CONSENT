#include "kMersProcessing.h"

std::map<std::string, int> getKMersPos(std::string sequence, unsigned merSize) {
	std::map<std::string, int> mers;

	for (unsigned i = 0; i < sequence.length() - merSize + 1; i++) {
		// mers.insert((sequence.substr(i, merSize)));
		if (mers[sequence.substr(i, merSize)] == 0) {
			mers[sequence.substr(i, merSize)] = i + 1;
		} else {
			mers[sequence.substr(i, merSize)] = -1;
		}
	}

	return mers;
}

std::map<std::string, int> getKMersCounts(std::vector<std::string> sequences, unsigned merSize) {
	std::map<std::string, int> merCounts;
	unsigned i;

	for (std::string s : sequences) {
		i = 0;
		while (s.length() >= merSize && i < s.length() - merSize + 1) {
			merCounts[(s.substr(i, merSize))]++;
			i++;
		}
	}

	return merCounts;
}