#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <vector>
#include <iterator>
#include <regex>
#include <stack> 


class Nussinov
{
    public:
        Nussinov(){} //default constructor
        void add_seqs(std::vector<std::string> strvec, std::vector<std::string> namevec) //Reads our vector of strings from main into our object
        {
            if (strvec.empty() || namevec.empty() || strvec.size() != namevec.size())
            {
                std::cout << "[ERROR]: Error reading in sequences from file. " << std::endl;
                exit(EXIT_FAILURE);
            }
            else
            {
                seqvec = strvec;
                svecsz = seqvec.size();
                namesvec = namevec;
            }
        }
        void runofanalysis(std::string seqstr) //This function will call other member functions in sequential order to perform the algorithm on a given string
        {
            set_RNA_seq(seqstr); //Set our sequence as per the argument string
            initialize_scoremat(); //Initialize scoring matrix to all 0s
            fillscmat(); //Fill the scoring matrix with corresponding values
            std::cout << std::endl;
            std::cout << "Now printing our scoring matrix for: " << namesvec[runiter] << std::endl;
            std::cout << std::endl;
            print_mat(); //print the scoring matrix
            //std::cout << std::endl;
            //traceback_mat(); [TO BE ADDED WHEN FUNCTIONAL]
            //std::cout << "Ok, here is what we got for our traceback: " << std::endl;
            //print_tbmat(); [TO BE ADDED WHEN FUNCTIONAL]
            //std::cout << std::endl;
            //make_structs() [TO BE ADDED WHEN FUNCTIONAL?]
            //std::cout << "Structures are as follows: " << std::endl;
            //print_structs(); [TO BE ADDED WHEN FUNCTIONAL]
        }
        void driver_func() //Driver function, calls the algorithm on all sequences in the input file. Will be used in final version barring unforseen circumstances
        {
            for (int i = 0; i < svecsz; i++)
            {
                runiter = i;
                runofanalysis(seqvec[i]);
                scoremat.clear(); //clear scoring matrix between runs
            }
        }
        void print_seqvec()
        {
            for (int i = 0; i < svecsz; i++)
            {
                std::cout << seqvec[i] << std::endl;
            }
        }
        bool is_complimentary(char a, char b) //checks our bases for complimentarity
        {
            return (a == 'A' && b == 'U') || (a == 'U' && b == 'A') || (a == 'C' && b == 'G') || (a == 'G' && b == 'C');
        }
        void test_function() //Test function for early submissions. Checks to see if our earliest iterations of progress work. 
        {
            runofanalysis(seqvec[0]);
        }
        int maxval(int ip, int jm, int ipjmd, int bifurcation) //Function returns our maximum value
        {
            std::vector<int> vec;
            vec.push_back(ip);
            vec.push_back(jm);
            vec.push_back(ipjmd);
            vec.push_back(bifurcation); 
            sort(vec.begin(), vec.end(), std::greater<int>()); //Sorts our vector by descending order, ensuring that the first value is either the largest one, or one of a few equal largest numbers
            return vec[0];
        }
        void set_RNA_seq(std::string str) //Adds the RNA sequence to our object in which we will utilize the Nussinov algorithm
        {
            inpseq = str;
            seqlen = inpseq.length();
        }
        void initialize_scoremat() //Initializes our scoring matrix with all 0s 
        {
            int init = 0;
            std::vector<int> emptyvec; //Made this to avoid seg faults primarily. 
            for (int i = 0; i < seqlen; i++)
            {
                scoremat.push_back(emptyvec);
                for (int j = 0; j < seqlen; j++)
                {
                    scoremat[i].push_back(init);
                }
            }
        }
        void print_vec(std::vector<int> vec) //Prints an individual vector
        {
            for (unsigned int i = 0; i < vec.size(); i++)
            {
                std::cout << vec[i];
            }
            std::cout << std::endl;
        }
        void print_mat() //Prints our entire scoring matrix
        {
            std::cout << "  " << inpseq << std::endl;
            for (int i = 0; i < seqlen; i++)
            {
                std::cout << inpseq[i] << ' ';
                print_vec(scoremat[i]);
            }
        }

        void traceback_mat(int i, int j) //recursively determines structure [IN PROGRESS]
        {
            std::stack<int> tbstack;
            tbstack.push(scoremat[0][seqlen - 1]);

        }

        int ret_bifurmax(int i, int j) //returns the maximum bifurcation value as an int
        {
            if (j - 2 > 0 && j - i > 3)
            {
                std::vector<int> vec;
                int k1, k2, kadd;
                for (int k = i + 1; k < j - 2; k++)
                {
                    k1 = scoremat[i][k];
                    k2 = scoremat[k + 1][j];
                    kadd = k1 + k2;
                    vec.push_back(kadd);
                }
                sort(vec.begin(), vec.end(), std::greater<int>());
                return vec[0];
            }
            else
            {
                return 0;
            }
        } 
        void fillfunc(int p, int j) //fills our scoring matrix diagonally
        {
            bool check;
            int diag, down, left, bifurcation; //diagonal, down, left
            for (int i = 0; i < p; i++)
            {
                check = is_complimentary(inpseq[i], inpseq[j]);
                diag = scoremat[i + 1][j - 1] + check; //bool cast as an int will return 0 if false and 1 if true, so we can do this to add score.
                down = scoremat[i + 1][j];
                left = scoremat[i][j - 1];
                bifurcation = ret_bifurmax(i, j);
                scoremat[i][j] = maxval(diag, down, left, bifurcation);
                j++;
            }
        }
        void fillscmat() //Driver function calls fillfunc until the matrix is filled. This ensures diagonal filling of the matrix
        {
            int j = 1;
            for (int i = seqlen - 1; i > 0; i--) //recursively calls our fillfunction
            {
                fillfunc(i, j);
                j++;
            }
            //traceback_mat(1, j);
        }

    private:
        int seqlen, svecsz, runiter; //sequence length, sequence vector size, run iterator
        std::vector<std::string> seqvec, namesvec; //vector holding our RNA strings for analysis
        std::string inpseq; //Input sequence
        std::vector<std::vector<int>> scoremat; //Scoring matrix, attempting with vector of vectors of ints
        std::vector<std::vector<std::pair<int, int>>> tbmat; //3D Traceback matrix (may need to be modified)
};


std::vector<std::string> make_seqs(std::vector<std::string> inpfil, std::vector<bool> bvec) //function returns a vector with only our sequence strings
{
    std::string seqline;
    std::vector<std::string> resvec;
    for (unsigned int i = 0; i < inpfil.size(); i++)
    {
        if (bvec[i] == true)
        {
            resvec.push_back(inpfil[i]);
        }
        else
        {
            continue;
        }
    }
    return resvec;
}

std::vector<std::string> makeseqnames(std::vector<std::string> inpfil, std::vector<bool> bvec)
{
    std::string seqline;
    std::vector<std::string> resvec;
    for (unsigned int i = 0; i < inpfil.size() - 1; i++)
    {
        if (bvec[i] == false && bvec[i + 1] == true)
        {
            resvec.push_back(inpfil[i]);
        }
        else
        {
            continue;
        }
    }
    return resvec;
}

bool isvalseq(std::string filelin)
{
    int len = filelin.length();
    for (int i = 0; i < len; i++)
    {
        if (filelin[i] != 'A' && filelin[i] != 'C' && filelin[i] != 'G' && filelin[i] != 'U')
        {
            return false;
        }
    }
    return true;
}


void printlvec(std::vector<std::string> linvec) //debugging function, shows what's actually input into our linevec object
{
    std::cout << "linvec.size() == " << linvec.size() << std::endl;
    for (unsigned int i = 0; i < linvec.size(); i++)
    {
        std::cout << linvec[i] << std::endl;
    }
}

void printboolvec(std::vector<bool> bvec)
{
    std::cout << "Printing our bool vector: " << std::endl;
    int sz = bvec.size();
    for (int i = 0; i < sz; i++)
    {
        std::cout << bvec[i] << " ";
    }
}

int main()
{
    std::regex r("\\S");
    std::string filenam, line; //Sequence and name of the file to be taken as a input
    std::vector<std::string> linevec, rnaseqs, seqnames; //Vector to store our lines from input file
    std::vector<bool> boolvec;
    bool linecheck;
    std::cout << "Hello! Please provide a FASTA format file with the RNA sequence(s) you would like a predicted secondary structure for: ";
    std::cin >> filenam;
    std::ifstream ifile(filenam);
    if (ifile.is_open())
    {
        while (getline(ifile, line))
        {
            std::string fileline = line;
            fileline = std::regex_replace(fileline, std::regex("\\r\\n|\\r|\\n"), "");
            linecheck = isvalseq(fileline);
            linevec.push_back(fileline);
            boolvec.push_back(linecheck);
        }
        ifile.close();
        printlvec(linevec);
        printboolvec(boolvec);
        rnaseqs = make_seqs(linevec, boolvec);
        seqnames = makeseqnames(linevec, boolvec);
        std::cout << "Ok, printing our used seq vec: " << std::endl; //remember to comment out in submission
        printlvec(rnaseqs);
        Nussinov nussi;
        nussi.add_seqs(rnaseqs, seqnames);
        //nussi.test_function();
        nussi.driver_func();
    }
    else
    {
        std::cout << "[ERROR] Error loading in input file. Terminating program: " << std::endl;
        exit(EXIT_FAILURE);
    }
    return 0;
}
