#include "utils.h"

bool fileExists(std::string fileName) {
	std::ifstream infile(fileName);
	return infile.good();
}

std::pair<int**, std::pair<int, int>> NeedlemanWunschLocalMatrix(std::string s1, std::string s2) {
	int **matrix = new int*[s1.length() + 1];
	for(int i = 0 ; i < s1.length() + 1; i++) {
    	matrix[i] = new int[s2.length() + 1];
    }

    for (int i = 0 ; i < s1.length() + 1 ; i++) {
    	matrix[i][0] = 0;
    }
	for (int j = 0 ; j < s2.length() + 1; j++) {
		matrix[0][j] = 0;
	}

	int maxI = 0;
	int maxJ = 0;
	int maximum = matrix[maxI][maxJ];
	for (int i = 1 ; i < s1.length() + 1; i++) {
		for (int j = 1 ; j < s2.length() + 1; j++) {
			int s = s1[i-1] == s2[j-1] ? 1 : -6;
			matrix[i][j] = std::max(std::max(0, matrix[i-1][j-1] + s), std::max(matrix[i-1][j] - 1, matrix[i][j-1] - 1));
			if (matrix[i][j] > maximum) {
				maximum = matrix[i][j];
				maxI = i;
				maxJ = j;
			}
		}
	}

	return std::make_pair(matrix, std::make_pair(maxI, maxJ));
}

std::pair<int, int> NeedlemanWunschLocalAlignments(std::string s1, std::string s2) {
	std::pair<int**, std::pair<int, int>> res = NeedlemanWunschLocalMatrix(s1, s2);
	int** matrix = res.first;
	std::pair<int, int> maxs = res.second;

	std::string al1, al2;
	// int i = s1.size();
	// int j = s2.size();
	int i = maxs.first;
	int j = maxs.second;

	int score, scoreDiag, scoreUp, scoreLeft, s;
	int editDistance = 0;

	while (i > 0 and j > 0 and matrix[i][j] != 0) {
		score = matrix[i][j];
		scoreDiag = matrix[i-1][j-1];
		scoreUp = matrix[i][j-1];
		scoreLeft = matrix[i-1][j];
		s = s1[i-1] == s2[j-1] ? 1 : -1;

		if (score == scoreDiag + s) {
			al1 = s1[i-1] + al1;
			al2 = s2[j-1] + al2;
			if (s1[i-1] != s2[j-1]) {
				editDistance++;
			}
			i--;
			j--;
		} else if (score == scoreLeft - 1) {
			al1 = s1[i-1] + al1;
			al2 = "-" + al2;
			i--;
			editDistance++;
		// score == scoreUp - 1
		} else {
			al1 = "-" + al1;
			al2 = s2[j-1] + al2;
			j--;
			editDistance++;
		}
	}

	while (i > 0 and matrix[i][j] != 0) {
		al1 = s1[i-1] + al1;
		al2 = "-" + al2;
		i--;
		editDistance++;
	}

	while (j > 0 and matrix[i][j] != 0) {
		al1 = "-" + al1;
		al2 = s2[j-1] + al2;
		j--;
		editDistance++;
	}

	for(int k = 0 ; k < s1.size() ; k++) {
		// std::cerr << i << std::endl;
    	delete [] matrix[k];
	}
	delete [] matrix;

	// // std::cerr << "editDistance : " << editDistance << std::endl;
	// std::cerr << "s1, from : " << i << " to " << maxs.first << std::endl;
	// std::cerr << "s2, from : " << j << " to " << maxs.second << std::endl;
	// return std::make_pair(al1, al2);

	return std::make_pair(j, maxs.second);
}

std::string reverseComplement(std::string seq) {
	std::string res = std::string(seq);
	for (int i = 0 ; i < seq.length() ; i++) {
		switch(seq[i]) {
			case 'A':
				res[seq.length() - i - 1] = 'T';
				break;
			case 'C':
				res[seq.length() - i - 1] = 'G';
				break;
			case 'G':
				res[seq.length() - i - 1] = 'C';
				break;
			case 'T':
				res[seq.length() - i - 1] = 'A';
				break;
			case 'a':
				res[seq.length() - i - 1] = 't';
				break;
			case 'c':
				res[seq.length() - i - 1] = 'g';
				break;
			case 'g':
				res[seq.length() - i - 1] = 'c';
				break;
			case 't':
				res[seq.length() - i - 1] = 'a';
				break;
		}
	}

	return res;
}