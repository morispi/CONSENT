#include <sstream>

struct Alignment {
	std::string qName;
	unsigned qLength;
	unsigned qStart;
	unsigned qEnd;
	bool strand;
	std::string tName;
	unsigned tLength;
	unsigned tStart;
	unsigned tEnd;
	std::string resMatches;
	std::string alBlockLen;
	std::string unordered_mapQual;

	bool operator<(const Alignment& a2) const {
		if (qName < a2.qName) {
			return true;
 		} else if (qName == a2.qName && qLength < a2.qLength) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart < a2.qStart) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd < a2.qEnd) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand < a2.strand) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName < a2.tName) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength < a2.tLength) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart < a2.tStart) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd < a2.tEnd) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd == a2.tEnd && resMatches < a2.resMatches) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd == a2.tEnd && resMatches == a2.resMatches && alBlockLen < a2.alBlockLen) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd == a2.tEnd && resMatches == a2.resMatches && alBlockLen == a2.alBlockLen && unordered_mapQual < a2.unordered_mapQual) {
 			return true;
 		} else {
 			return false;
 		}
	}

	Alignment() {
		
	}

	Alignment(std::string al) {
		std::string token;
		std::stringstream iss(al);
		getline(iss, token, '\t');
		qName = token;
		getline(iss, token, '\t');
		qLength = stoi(token);
		getline(iss, token, '\t');
		qStart = stoi(token);
		getline(iss, token, '\t');
		// Has to be -1, cause miniunordered_map flags as endPosition the nt following the last match
		qEnd = stoi(token) - 1;
		getline(iss, token, '\t');
		strand = token == "+" ? false : true;
		getline(iss, token, '\t');
		tName = token;
		getline(iss, token, '\t');
		tLength = stoi(token);
		getline(iss, token, '\t');
		tStart = stoi(token);
		getline(iss, token, '\t');
		// Has to be -1, cause miniunordered_map flags as endPosition the nt following the last match
		tEnd = stoi(token) - 1;
		getline(iss, token, '\t');
		resMatches = token;
		getline(iss, token, '\t');
		alBlockLen = token;
		getline(iss, token, '\t');
		unordered_mapQual = token;
	}
};