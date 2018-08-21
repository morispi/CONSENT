#include <unistd.h>
#include<string>
#include <fstream>

bool fileExists(std::string fileName);

std::pair<int**, std::pair<int, int>> NeedlemanWunschLocalMatrix(std::string s1, std::string s2);

std::pair<int, int> NeedlemanWunschLocalAlignments(std::string s1, std::string s2);

// std::string reverseComplement(std::string seq);

struct reverseComplement_t {
    
    reverseComplement_t();
   
    std::string operator()(std::string seq) const;

    private:
        char rev_comp_tab['T' + 1] = {0};
};

const static reverseComplement_t reverseComplement;