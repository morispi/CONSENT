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
	for (i = 1 ; i < s1.length() + 1; i++) {
		for (j = 1 ; j < s2.length() + 1; j++) {
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
	int i = maxs.first;
	int j = maxs.second;

	int score, scoreDiag, scoreLeft, s;
	int editDistance = 0;

	while (i > 0 and j > 0 and matrix[i][j] != 0) {
		score = matrix[i][j];
		scoreDiag = matrix[i-1][j-1];
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

	for(unsigned k = 0 ; k < s1.size() ; k++) {
    	delete [] matrix[k];
	}
	delete [] matrix;

	return std::make_pair(j, maxs.second);
}
