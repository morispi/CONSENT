#include "localAlignment.h"

std::pair<int**, std::pair<int, int>> NeedlemanWunschLocalMatrix(std::string s1, std::string s2) {
	int **matrix = new int*[s1.length() + 1];
	unsigned i, j;
	for(i = 0 ; i < s1.length() + 1; i++) {
    	matrix[i] = new int[s2.length() + 1];
    }

    for (i = 0 ; i < s1.length() + 1 ; i++) {
    	matrix[i][0] = 0;
    }
	for (j = 0 ; j < s2.length() + 1; j++) {
		matrix[0][j] = 0;
	}

	int maxI = 0;
	int maxJ = 0;
	int maximum = matrix[maxI][maxJ];
	int s;
	for (i = 1 ; i < s1.length() + 1; i++) {
		for (j = 1 ; j < s2.length() + 1; j++) {
			// s = s1[i-1] == s2[j-1] ? 5: -10;
			// matrix[i][j] = std::max(std::max(0, matrix[i-1][j-1] + s), std::max(matrix[i-1][j] - 5, matrix[i][j-1] - 5));
			s = s1[i-1] == s2[j-1] ? 1 : -5;
			matrix[i][j] = std::max(std::max(0, matrix[i-1][j-1] + s), std::max(matrix[i-1][j] - 1, matrix[i][j-1] - 1));
			if (matrix[i][j] >= maximum) {
				maximum = matrix[i][j];
				maxI = i;
				maxJ = j;
			}
		}
	}

	return std::make_pair(matrix, std::make_pair(maxI, maxJ));
}

std::pair<std::pair<int, int>, std::pair<int, int>> NeedlemanWunschLocalAlignments(std::string s1, std::string s2) {
	std::pair<int**, std::pair<int, int>> res = NeedlemanWunschLocalMatrix(s1, s2);
	int** matrix = res.first;
	std::pair<int, int> maxs = res.second;

	int i = maxs.first;
	int j = maxs.second;

	int score, scoreDiag, scoreLeft, s;

	while (i > 0 and j > 0 and matrix[i][j] != 0) {
		score = matrix[i][j];
		scoreDiag = matrix[i-1][j-1];
		scoreLeft = matrix[i-1][j];
		s = s1[i-1] == s2[j-1] ? 1 : -5;

		if (score == scoreDiag + s) {
			i--;
			j--;
		} else if (score == scoreLeft - 1) {
			i--;
		} else {
			j--;
		}
	}

	while (i > 0 and matrix[i][j] != 0) {
		i--;
	}

	while (j > 0 and matrix[i][j] != 0) {
		j--;
	}

	for(unsigned k = 0 ; k < s1.size() ; k++) {
    	delete [] matrix[k];
	}
	delete [] matrix;

	return std::make_pair(std::make_pair(i, maxs.first - 1), std::make_pair(j, maxs.second - 1));
}
