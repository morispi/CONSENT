#include <unistd.h>
#include<string>
#include <fstream>

std::pair<int**, std::pair<int, int>> NeedlemanWunschLocalMatrix(std::string s1, std::string s2);

std::pair<std::pair<int, int>, std::pair<int, int>> NeedlemanWunschLocalAlignments(std::string s1, std::string s2);