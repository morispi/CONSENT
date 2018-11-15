#include "kMersProcessing.h"

std::unordered_map<std::string, std::vector<unsigned>> getKMersPos(std::string sequence, unsigned merSize) {
	std::unordered_map<std::string, std::vector<unsigned>> mers;

	for (unsigned i = 0; i < sequence.length() - merSize + 1; i++) {
			mers[sequence.substr(i, merSize)].push_back(i);
	}

	return mers;
}

std::unordered_map<std::string, unsigned> getKMersCounts(std::vector<std::string>& sequences, unsigned merSize) {
	std::unordered_map<std::string, unsigned> merCounts;
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

std::unordered_map<std::string, std::vector<unsigned>> getKMersOccs(std::vector<std::string>& sequences, unsigned merSize) {
	std::unordered_map<std::string, std::vector<unsigned>> merOccs;

	for (unsigned i = 0; i < sequences.size(); i++) {
		for (unsigned j = 0; j < sequences[i].size() - merSize + 1; j++) {
			merOccs[sequences[i].substr(j, merSize)].push_back(i);
		}
	}

	return merOccs;
}