#include <string>

class rev_comp {
    public:
        static std::string run(std::string seq) {
        	rev_comp::build_instance();
        	
        	auto first = seq.begin(), last = seq.end();
                
        	while(true) {
        		if(first == last || first == --last) {
            	    if(seq.length() % 2) {
            		    *first = rev_comp::complement[*first];
                    }
            	    return seq;
                } else {
                	*first = rev_comp::complement[*first];
            	    *last = rev_comp::complement[*last];
            	    std::iter_swap(first, last);
            	    ++first;
                }
            }
        }

    protected:
        static rev_comp* _instance;
        static void build_instance() {
            if(_instance == nullptr) {
                _instance = new rev_comp();
            }
        }
        
        rev_comp() {
            this->complement['A'] = 'T';
            this->complement['T'] = 'A';
            this->complement['C'] = 'G';
            this->complement['G'] = 'C';
            this->complement['a'] = 't';
            this->complement['t'] = 'a';
            this->complement['c'] = 'g';
            this->complement['g'] = 'c';
        }
            
        ~rev_comp(){}
        
    private:
        static char complement['T'];
};

char rev_comp::complement['T'];
rev_comp* rev_comp::_instance = nullptr;