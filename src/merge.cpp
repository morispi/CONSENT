#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <string.h>

std::vector<std::string> getFilesNames(char** args) {
	std::vector<std::string> res;

	int i = 3;
	while (args[i] != NULL) {
		res.push_back(args[i]);
		i++;
	}

	return res;
}

std::string getHeader(std::string line) {
	std::stringstream iss(line);
	std::string token;
	getline(iss, token, '\t');
	return token;
}

int main (int argc, char* argv[]) {
	std::ofstream outFile;
	outFile.open(argv[1]);

	std::ifstream headers;
	headers.open(argv[2]);

	std::vector<std::string> filesNames = getFilesNames(argv);
	std::vector<std::ifstream*> files(filesNames.size());
	std::ifstream* f;
	for (unsigned long i = 0; i < filesNames.size(); i++) {
		f = new std::ifstream;
		files[i] = f;
		files[i]->open(filesNames[i]);
	}

	std::string header;
	std::string line;
	while(getline(headers, header)) {
		header = header.substr(1);
		for (unsigned i = 0; i < files.size(); i++) {
			while(getline(*files[i], line) && getHeader(line) == header) {
				outFile	<< line << std::endl;
			}
			if (!files[i]->eof()) {
				files[i]->seekg(-line.length()-1, files[i]->cur);
			}
		}
	}


	for (std::ifstream* f : files) {
		f->close();
	}
	outFile.close();
	headers.close();
}
