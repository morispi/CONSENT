#include "kMersProcessing.h"

std::map<std::string, std::vector<int>> getKMersPos(std::string sequence, unsigned merSize) {
	std::map<std::string, std::vector<int>> mers;

	for (unsigned i = 0; i < sequence.length() - merSize + 1; i++) {
			mers[sequence.substr(i, merSize)].push_back(i);
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