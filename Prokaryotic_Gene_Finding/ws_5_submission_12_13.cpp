#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <regex>
#include <unordered_map>
#include <map>
#include <cmath>

class Markov
{
public:
    Markov(){} //Default constructor
    void add_seqs(std::vector<std::string> svec1, std::vector<std::string> nvec1, std::vector<std::string> svec2, std::vector<std::string> nvec2, std::vector<std::string> tvec, std::vector<std::string> tesnam) //reads sequence vector from main into object
    {
        seqvecn = svec1;
        svecszn = svec1.size();
        namevecn = nvec1;
        seqveco = svec2;
        svecszo = seqveco.size();
        nameveco = nvec2;
        tseqs = tvec;
        tnames = tesnam;
    }

    void driver_func() //Driver function to cut down clutter in main()
    {
        bool check, chk, fincheck, chksm;
        build_charvec(); //for calculation purposes
        N_read_in_codons(); //build nmap
        O_read_in_codons(); //build omap
        n_vector_codon_occurences(); //put into a vector pre sorted as per std::map [NORF]
        n_build_vecovec(); //Put it all into a vector of vectors, sorted by codon transitioned from [NORF]
        n_percents(); //Calculate into percentages [NORF]
        o_vector_codon_occurences(); //put into a vector pre sorted as per std::map [ORF]
        o_build_vecovec(); //Put it all into a vector of vectors, sorted by codon transitioned from [ORF]
        o_percents(); //Calculate into percentages [ORF]
        acct_for_zeros_n(); //Account for transitions that never occur [NORF]
        acct_for_zeros_o(); //Account for transitions that never occur [ORF]
        set_zero_case(); //Set zero case bias doubles
        build_tscorevec(); //Build the scoring matrix
        read_sm_to_map(); //read our scoring matrix to a map object for ease of indexing
        score_known_seqs(); //score sequences from known NORFS and ORFS in their own vectors
        std::cout << "Output transition scores for known ORFS/NORFS to .csv format? [1 = yes, 0 = no]: ";
        std::cin >> chksm; 
        if (chksm == true)
        {
            output_scores_csv(); //output score vector to .csv format
        }
        else
        {
            std::cout << "Understandable, have a nice day!" << std::endl;
        }
        score_test_seqs(); //Score our test sequences
        std::cout << std::endl;
        std::cout << "Output final scoring results for your sequences? [1 = yes, 0 = no]: ";
        std::cin >> fincheck;
        if (fincheck == true)
        {
            output_final_results(); //output final results to a separate .txt file
        }
        else
        {
            std::cout << "Understandable, have a nice day! " << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Output raw NORF transition matrix? [1 = yes, 0 = no]: ";
        std::cin >> check;
        if (check == true)
        {
            output_sm_norf();
        }
        else
        {
            std::cout << "Understandable have a nice day! " << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Output raw ORF transition matrix? [1 = yes, 0 = no]: ";
        std::cin >> chk;
        if (chk == true)
        {
            output_sm_orf();
        }
        else
        {
            std::cout << "Understandable have a nice day! " << std::endl;
        }
    }

    double min_dub(std::vector<std::vector<std::pair<std::pair<std::string, std::string>, double>>> dvec)
    {
        double mindub = dvec[0][0].second;
        for (unsigned int i = 0; i < dvec.size(); i++)
        {
            for (unsigned int j = 0; j < dvec[i].size(); j++)
            {
                if (dvec[i][j].second < mindub && dvec[i][j].second != 0)
                {
                    mindub = dvec[i][j].second;
                }
            }
        }
        return mindub;
    }

    void set_zero_case() //sets minimum values for 0 case
    {
        zeron = min_dub(pctvecn);
        zeroo = min_dub(pctveco);
        zeron = zeron / 2.00;
        zeroo = zeroo / 2.00;
    }

    void build_charvec() //Builds a vector of characters storing all nucleotides from our input sequence vector, for determining total occurence of codons (NORF)
    {
        for (int i = 0; i < svecszn; i++)
        {
            for (unsigned int j = 0; j < seqvecn[i].length(); j++)
            {
                charvecn.push_back(seqvecn[i][j]);
            }
        }
        totnucn = charvecn.size();
        totcodon = (totnucn / 3) - 1; //Transitions
        for (int i = 0; i < svecszo; i++)
        {
            for (unsigned int j = 0; j < seqveco[i].length(); j++)
            {
                charveco.push_back(seqveco[i][j]);
            }
        }
        totnuco = charveco.size();
        totcodo = (totnuco / 3) - 1; //Transitions
    }

    void N_read_in_codons() //Reads in each codon in our input file into a hashmap, incrementing values as the number of times they appear (NORF)
    {
        for (int i = 0; i < svecszn; i++)
        {
            for (unsigned int j = 0; j < seqvecn[i].length() - 5; j += 3)
            {
                std::unordered_map<std::string, double> temp;
                std::string codon1, codon2; //Create codon string then add 3 nucleotides at a time to it, sliding by 3 each iteration, codon2 stores preceeding codon
                codon1 += seqvecn[i][j];
                codon1 += seqvecn[i][j + 1];
                codon1 += seqvecn[i][j + 2];
                codon2 += seqvecn[i][j + 3];
                codon2 += seqvecn[i][j + 4];
                codon2 += seqvecn[i][j + 5];
                std::pair<std::string, std::string> temppair(codon1, codon2);
                if (cmapn.find(temppair) == cmapn.end()) //Check to see if the codon key already exists. If not, add it with value of 1
                {
                    cmapn.insert(std::make_pair(temppair, 1));
                }
                else //If the key already exists, increment the value
                {
                    cmapn[temppair]++;
                }
            }
        }
        std::cout << std::endl;
    }

    double score_seq(std::string seq) //takes a sequence as an argument and returns a double of total score
    {
        double res = 0;
        for (unsigned int i = 0; i < seq.length() - 5; i += 3)
        {
            std::string codon1, codon2;
            codon1 += seq[i];
            codon1 += seq[i + 1];
            codon1 += seq[i + 2];
            codon2 += seq[i + 3];
            codon2 += seq[i + 4];
            codon2 += seq[i + 5];
            std::pair<std::string, std::string> temppair(codon1, codon2);
            double tdub = scoremap[temppair];
            res += tdub;
        }
        return res;
    }

    void score_known_seqs() //store sequence scores in vector of doubles
    {
        for (unsigned int i = 0; i < svecszn; i++)
        {
            double ndub = score_seq(seqvecn[i]);
            double szn = seqvecn[i].length();
            ndub = ndub / szn;
            knownnorf.push_back(ndub);
        }
        for (unsigned int j = 0; j < svecszo; j++)
        {
            double odub = score_seq(seqveco[j]);
            double sz = seqveco[j].length();
            odub = odub / sz;
            knownorf.push_back(odub);
        }
    }

    void output_final_results()
    {
        std::string fnm;
        std::cout << "Please choose a name for your output file (.txt format): ";
        std::cin >> fnm;
        std::ofstream opob(fnm, std::ofstream::out);
        opob << "Predictions for tested sequences: " << '\n';
        for (unsigned int i = 0; i < testscores.size(); i++)
        {
            if (testscores[i] >= -0.002000)
            {
                opob << tnames[i] << ": Predicted: Open Reading Frame, score = " << testscores[i] << '\n';
            }
            else
            {
                opob << tnames[i] << ": Predicted: Non-Open Reading Frame, score = " << testscores[i] << '\n';
            }
        }
        opob << "Thank you for using my Reading Frame Finder! [Cutoff = -0.002000] " << '\n';
        opob.close();
    }

    void O_read_in_codons() //Reads in each codon in our input file into a hashmap, incrementing values as the number of times they appear (ORF)
    {
        for (int i = 0; i < svecszo; i++)
        {
            for (unsigned int j = 0; j < seqveco[i].length() - 5; j += 3)
            {
                std::unordered_map<std::string, int> temp;
                std::string codon1, codon2; //Create codon string then add 3 nucleotides at a time to it, sliding by 3 each iteration, codon2 stores preceeding codon
                codon1 += seqveco[i][j];
                codon1 += seqveco[i][j + 1];
                codon1 += seqveco[i][j + 2];
                codon2 += seqveco[i][j + 3];
                codon2 += seqveco[i][j + 4];
                codon2 += seqveco[i][j + 5];
                std::pair<std::string, std::string> temppair(codon1, codon2);
                if (cmapo.find(temppair) == cmapo.end()) //Check to see if the codon key already exists. If not, add it with value of 1
                {
                    cmapo.insert(std::make_pair(temppair, 1));
                }
                else //If the key already exists, increment the value
                {
                    cmapo[temppair]++;
                }
            }
        }
    }

    void n_vector_codon_occurences() //Prints each found codon and the total number of times it appeared in our test set (NORF)
    {
        for (auto const &pair : cmapn)
        {
            codvecn.push_back(pair); //Add to vector. We can edit the name of the function and comment out the prints in later submissions
        }
    }

    void o_vector_codon_occurences() //Prints each found codon and the total number of times it appeared in our test set (NORF)
    {
        for (auto const &pair : cmapo)
        {
            codveco.push_back(pair); //Add to vector. We can edit the name of the function and comment out the prints in later submissions
        }
    }

    void output_sm_norf()
    {
        bool check;
        std::string onam;
        std::cout << "Hello! Output NORF Scoring matrix to .csv format? [1/0]";
        std::cin >> check;
        if (check == true)
        {
            std::cout << "Please choose a file name in .csv format: ";
            std::cin >> onam;
            std::ofstream oob(onam, std::ofstream::out);
            oob << " " << ',';
            for (unsigned int i = 0; i < pctvecn.size(); i++)
            {
                oob << pctvecn[i][0].first.first << ',';
            }
            oob << '\n';
            for (unsigned int j = 0; j < pctvecn.size(); j++)
            {
                oob << pctvecn[j][0].first.first << ',';
                for (unsigned int k = 0; k < pctvecn[j].size(); k++)
                {
                    oob << pctvecn[j][k].second << ',';
                }
                oob << '\n';
            }
            oob.close();
        }
        else
        {
            std::cout << "Ok, have a nice day then!" << std::endl;
        }
    }

    double ret_least()
    {
        if (zeron < zeroo)
        {
            return zeron;
        }
        else
        {
            return zeroo;
        }
        return zeroo;
    }

    void read_sm_to_map()
    {
        for (unsigned int i = 0; i < scorev.size(); i++)
        {
            for (unsigned int j = 0; j < scorev[i].size(); j++)
            {
                std::pair<std::string, std::string> codpair(scorev[i][j].first.first, scorev[i][j].first.second);
                scoremap.insert(std::make_pair(codpair, scorev[i][j].second));
            }
        }
    }

    void build_tscorevec() //build our scoring vector
    {
        for (unsigned int i = 0; i < pctvecn.size(); i++)
        {
            scorev.push_back(pctvecn[i]);
            for (unsigned int j = 0; j < pctvecn[i].size(); j++)
            {
                if (pctvecn[i][j].second == 0 && pctveco[i][j].second != 0)
                {
                    double sc = log2(zeroo);
                    scorev[i][j].second = sc;
                }
                else if (pctveco[i][j].second == 0 && pctvecn[i][j].second != 0)
                {
                    double scn = log2(zeron);
                    scorev[i][j].second == scn;
                }
                else if (pctveco[i][j].second == 0 && pctvecn[i][j].second == 0) //Ask about this case**
                {
                    double smallest = ret_least();
                    double scdub = log2(smallest);
                    scorev[i][j].second = smallest;
                }
                else
                {
                    double score = pctveco[i][j].second / pctvecn[i][j].second; //Divide our ORF score by our NORF score for the same transition
                    score = log2(score);                                        //Take the log2 of that score
                    scorev[i][j].second = score;
                }
            }
        }
    }

    void output_sm_orf()
    {
        bool check;
        std::string onam;
        std::cout << "Hello! Output ORF Scoring matrix to .csv format? [1/0]";
        std::cin >> check;
        if (check == true)
        {
            std::cout << "Please choose a file name in .csv format: ";
            std::cin >> onam;
            std::ofstream oob(onam, std::ofstream::out);
            oob << " " << ',';
            for (unsigned int i = 0; i < pctveco.size(); i++)
            {
                oob << pctveco[i][0].first.first << ',';
            }
            oob << '\n';
            for (unsigned int j = 0; j < pctveco.size(); j++)
            {
                oob << pctveco[j][0].first.first << ',';
                for (unsigned int k = 0; k < pctveco[j].size(); k++)
                {
                    oob << pctveco[j][k].second << ',';
                }
                oob << '\n';
            }
            oob.close();
        }
        else
        {
            std::cout << "Ok, have a nice day then!" << std::endl;
        }
    }

    void output_tscorem() //function outputs our transition scoring matrix
    {
        bool check;
        std::string onam;
        std::cout << "Hello! Output Transition Scoring matrix to .csv format? [1/0]";
        std::cin >> check;
        if (check == true)
        {
            std::cout << "Please choose a file name in .csv format: ";
            std::cin >> onam;
            std::ofstream oob(onam, std::ofstream::out);
            oob << " " << ',';
            for (unsigned int i = 0; i < scorev.size(); i++)
            {
                oob << scorev[i][0].first.first << ',';
            }
            oob << '\n';
            for (unsigned int j = 0; j < scorev.size(); j++)
            {
                oob << scorev[j][0].first.first << ',';
                for (unsigned int k = 0; k < scorev[j].size(); k++)
                {
                    oob << scorev[j][k].second << ',';
                }
                oob << '\n';
            }
            oob.close();
        }
        else
        {
            std::cout << "Ok, have a nice day then!" << std::endl;
        }
    }

    void n_build_vecovec()
    {
        std::vector<std::pair<std::pair<std::string, std::string>, double>> tmpvec;
        tmpvec.push_back(codvecn[0]);
        for (unsigned int i = 1; i < codvecn.size(); i++)
        {
            if (codvecn[i].first.first != codvecn[i - 1].first.first)
            {
                nwcodvec.push_back(tmpvec);
                tmpvec.clear();
                tmpvec.push_back(codvecn[i]);
            }
            else if (i + 1 == codvecn.size())
            {
                tmpvec.push_back(codvecn[i]);
                nwcodvec.push_back(tmpvec);
            }
            else
            {
                tmpvec.push_back(codvecn[i]);
            }
        }
    }

    std::vector<std::pair<std::pair<std::string, std::string>, double>> find_fullvec(std::vector<std::vector<std::pair<std::pair<std::string, std::string>, double>>> vec)
    {
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            if (vec[i].size() == vec.size())
            {
                return vec[i];
            }
        }
        return vec[0];
    }

    void score_test_seqs() //Score sequences from our test vector
    {
        for (unsigned int i = 0; i < tseqs.size(); i++)
        {
            double scr = score_seq(tseqs[i]);
            double szt = tseqs[i].length();
            scr = scr / szt;
            testscores.push_back(scr);
        }
    }

    void acct_for_zeros_n() //adds 0s for codon transitions not seen [NORF]
    {
        std::vector<std::pair<std::pair<std::string, std::string>, double>> compvec = find_fullvec(pctvecn); //finds first full vector
        for (unsigned int i = 0; i < pctvecn.size(); i++)
        {
            if (pctvecn[i].size() != compvec.size())
            {
                for (unsigned int j = 0; j < compvec.size(); j++)
                {
                    if (pctvecn[i][j].first.second != compvec[j].first.second)
                    {
                        std::pair<std::pair<std::string, std::string>, double> emptpair;
                        emptpair.first.first = pctvecn[i][0].first.first;
                        emptpair.first.second = compvec[j].first.second;
                        emptpair.second = 0;
                        auto itPos = pctvecn[i].begin() + j;
                        pctvecn[i].insert(itPos, emptpair);
                    }
                    else
                    {
                        continue;
                    }
                }
            }
            else
            {
                continue;
            }
        }
    }

    void acct_for_zeros_o() //adds 0s for codon transitions not seen [ORF]
    {
        std::vector<std::pair<std::pair<std::string, std::string>, double>> compvec = find_fullvec(pctveco); //finds first full vector
        for (unsigned int i = 0; i < pctveco.size(); i++)
        {
            if (pctveco[i].size() != compvec.size())
            {
                for (unsigned int j = 0; j < compvec.size(); j++)
                {
                    if (pctveco[i][j].first.second != compvec[j].first.second)
                    {
                        std::pair<std::pair<std::string, std::string>, double> emptpair;
                        emptpair.first.first = pctveco[i][0].first.first;
                        emptpair.first.second = compvec[j].first.second;
                        emptpair.second = 0;
                        auto itPos = pctveco[i].begin() + j;
                        pctveco[i].insert(itPos, emptpair);
                    }
                    else
                    {
                        continue;
                    }
                }
            }
            else
            {
                continue;
            }
        }
    }

    void output_scores_csv() //output scores to a csv file for distribution generation purposes
    {
        std::string filenam;
        std::cout << "Please provide a name for score file in .csv format: ";
        std::cin >> filenam;
        std::ofstream outp(filenam, std::ofstream::out);
        outp << "Known" << ',' << "Score" << ',';
        outp << '\n';
        for (unsigned int i = 0; i < knownnorf.size(); i++)
        {
            outp << "NORF" << ',' << knownnorf[i] << '\n';
        }
        for (unsigned int j = 0; j < knownorf.size(); j++)
        {
            outp << "ORF" << ',' << knownorf[j] << '\n';
        }
        outp.close();
    }

    double n_total(int it) //Return total transitions in NORF file for a given nucleotide
    {
        double sum = 0;
        for (unsigned int i = 0; i < nwcodvec[it].size(); i++)
        {
            sum += nwcodvec[it][i].second;
        }
        return sum;
    }

    double o_total(int it) //Return total transitions in ORF file for a given nucleotide
    {
        double sum = 0;
        for (unsigned int i = 0; i < ocodvec[it].size(); i++)
        {
            sum += ocodvec[it][i].second;
        }
        return sum;
    }

    void n_percents() //Store proportion of nucleotide transitions for each nucleotide as a double percentage of total transitions for that nucleotide (NORF)
    {
        for (unsigned int i = 0; i < nwcodvec.size(); i++)
        {
            pctvecn.emplace_back(nwcodvec[i]); //fix here (*)
            double todo = n_total(i);
            for (unsigned int j = 0; j < nwcodvec[i].size(); j++)
            {
                double pct = nwcodvec[i][j].second;
                pct = pct / todo;
                pctvecn[i][j].second = pct;
            }
        }
    }

    void o_percents() //Store proportion of nucleotide transitions for each nucleotide as a double percentage of total transitions for that nucleotide (ORF)
    {
        for (unsigned int i = 0; i < ocodvec.size(); i++)
        {
            pctveco.emplace_back(ocodvec[i]);
            double todo = o_total(i);
            for (unsigned int j = 0; j < ocodvec[i].size(); j++)
            {
                double pct = ocodvec[i][j].second;
                pct = pct / todo;
                pctveco[i][j].second = pct;
            }
        }
    }

    void o_build_vecovec()
    {
        std::vector<std::pair<std::pair<std::string, std::string>, double>> tmpvec;
        tmpvec.push_back(codveco[0]);
        for (unsigned int i = 1; i < codveco.size(); i++)
        {
            if (codveco[i].first.first != codveco[i - 1].first.first)
            {
                ocodvec.push_back(tmpvec);
                tmpvec.clear();
                tmpvec.push_back(codveco[i]);
            }
            else if (i + 1 == codveco.size())
            {
                tmpvec.push_back(codveco[i]);
                ocodvec.push_back(tmpvec);
            }
            else
            {
                tmpvec.push_back(codveco[i]);
            }
        }
    }

    int n_vec_location(std::string cdn)
    {
        for (unsigned int i = 0; i < pctvecn.size(); i++)
        {
            for (unsigned int j = 0; j < pctvecn[i].size(); j++)
            {
                if (pctvecn[i][0].first.first == cdn)
                {
                    return i;
                }
            }
        }
        return -1;
    }

    int o_vec_location(std::string cdn)
    {
        for (unsigned int i = 0; i < pctveco.size(); i++)
        {
            for (unsigned int j = 0; j < pctveco[i].size(); j++)
            {
                if (pctveco[i][0].first.first == cdn)
                {
                    return i;
                }
            }
        }
        return -1;
    }

    void print_percentages_n(int ind)
    {
        for (unsigned int i = 0; i < pctvecn[ind].size(); i++)
        {
            std::cout << "[" << pctvecn[ind][i].first.first << " -> " << pctvecn[ind][i].first.second << "  :  " << pctvecn[ind][i].second << "]" << std::endl;
        }
    }

    void print_percentages_o(int ind)
    {
        for (unsigned int i = 0; i < pctveco[ind].size(); i++)
        {
            std::cout << "[" << pctveco[ind][i].first.first << " -> " << pctveco[ind][i].first.second << "  :  " << pctveco[ind][i].second << "]" << std::endl;
        }
    }

    void print_norf_codons() //user input codon printing function for feedback
    {
        bool check, check2;
        std::string codn;
        std::cout << "Would you like to print transitions for a codon? (1/0) [NORF]: ";
        std::cin >> check;
        if (check == true)
        {
            std::cout << "Ok, what codon would you like to print?(Please print 3 uppercase valid nucelobases (A, T, G, or C)): ";
            std::cin >> codn;
            int index = n_vec_location(codn);
            if (index == -1)
            {
                std::cout << "[ERROR] Codon not found or invalid" << std::endl;
            }
            else
            {
                std::cout << "Transitions for: " << codn << std::endl;
                std::cout << std::endl;
                print_percentages_n(index);
            }
            std::cout << "Would you like to print transitions for another codon? ";
            std::cin >> check2;
            if (check2 == true)
            {
                print_norf_codons(); //call function again if we wanna do this again
            }
            else
            {
                std::cout << std::endl;
            }
        }
    }

    void print_orf_codons() //user input codon printing function for feedback
    {
        bool check, check2;
        std::string codn;
        std::cout << "Would you like to print transitions for a codon? (1/0) [ORF]: ";
        std::cin >> check;
        if (check == true)
        {
            std::cout << "Ok, what codon would you like to print?(Please print 3 uppercase valid nucelobases (A, T, G, or C)): ";
            std::cin >> codn;
            int index = o_vec_location(codn);
            if (index == -1)
            {
                std::cout << "[ERROR] Codon not found or invalid" << std::endl;
            }
            else
            {
                std::cout << "Transitions for: " << codn << std::endl;
                std::cout << std::endl;
                print_percentages_o(index);
            }
            std::cout << "Would you like to print transitions for another codon? ";
            std::cin >> check2;
            if (check2 == true)
            {
                print_orf_codons(); //call function again if we wanna do this again
            }
            else
            {
                std::cout << std::endl;
            }
        }
    }

private:
    int svecszn, totnucn, totcodon, svecszo, totnuco, totcodo;                                                 //Sequence vector size, total nucleotide count
    double tcdn, codcntn, tcdo, codcnto, zeron, zeroo;                                                         //Doubles for more accurate division/zero case
    std::vector<char> charvecn, charveco;                                                                      //Character vector holding every character in our sequence vector, for frequency purposes
    std::vector<std::string> seqvecn, namevecn, seqveco, nameveco, tseqs, tnames;                                             //Sequence vector and vector for names of sequences for NORF and ORF
    std::vector<std::pair<std::pair<std::string, std::string>, double>> codvecn, codveco;                      //Vector holding key value pairs for ease of data manipulation
    std::vector<std::vector<std::pair<std::pair<std::string, std::string>, double>>> pctvecn, pctveco, scorev; //Vector holding codon + proportion pairs, transition score vector
    std::map<std::pair<std::string, std::string>, double> cmapn, cmapo, scoremap;                              //maps to track transitions/store scores
    std::vector<std::vector<double>> codontabn, codontabo;                                                     //table for codons
    std::vector<double> knownnorf, knownorf, testscores;                                                                   // vector to store scores for known norf and orf values
    std::vector<std::vector<std::pair<std::pair<std::string, std::string>, double>>> nwcodvec, ocodvec;        //Separates by codon transition
};

std::vector<bool> make_boolvec(std::vector<std::string> inpvec) //Makes a vector of bools where true = a valid nucleotide sequence and false is everything else
{
    std::vector<bool> resvec;
    int sz = inpvec.size();
    for (int i = 0; i < sz; i++)
    {
        bool check = true;
        for (unsigned int j = 0; j < inpvec[i].length(); j++)
        {
            if (inpvec[i][j] != 'A' && inpvec[i][j] != 'C' && inpvec[i][j] != 'G' && inpvec[i][j] != 'T') //Edit made for U-T (reused from nussinov program)
            {
                check = false;
            }
        }
        resvec.push_back(check);
    }
    return resvec;
}

std::vector<std::string> make_seqs(std::vector<std::string> inpvec, std::vector<bool> bvec) //Creates a vector that stores valid nucleotide sequences (reused from Nussinov)
{
    std::vector<std::string> resvec;
    for (unsigned int i = 1; i < inpvec.size(); i++)
    {
        if (bvec[i] == true)
        {
            resvec.push_back(inpvec[i]);
        }
        else
        {
            continue;
        }
    }
    return resvec;
}

std::vector<std::string> make_names(std::vector<std::string> inpvec, std::vector<bool> bvec) //Creates vector to store sequence names (reused from nussinov)
{
    std::vector<std::string> resvec;
    for (unsigned int i = 0; i < inpvec.size() - 1; i++)
    {
        if (bvec[i] == false && bvec[i + 1] == true)
        {
            resvec.push_back(inpvec[i]);
        }
        else
        {
            continue;
        }
    }
    return resvec;
}

std::vector<std::string> make_strs(std::vector<std::string> inpvec) //Merges our sequence strings into a single string for each sequence (reused from nussinov)
{
    std::vector<std::string> resvec;
    std::string str;
    int sz = inpvec.size();
    for (int i = 0; i < sz; i++)
    {
        if (inpvec[i][0] == '>' && i == 0) //Add first line
        {
            resvec.push_back(inpvec[i]);
        }
        else if (inpvec[i][0] == '>' && i != 0) //Add any line with a carrot, but not before adding the cummulative str string
        {
            resvec.push_back(str);
            str.clear();
            resvec.push_back(inpvec[i]);
        }
        else if (i == sz - 1) //Add final string
        {
            str += inpvec[i];
            resvec.push_back(str);
        }
        else //Add consecutive sequence strings together
        {
            str += inpvec[i];
        }
    }
    return resvec;
}

std::string reverse_comp(std::string seqstr) //Function returns the reverse complement of an input string
{
    std::reverse(seqstr.begin(), seqstr.end());
    std::string res;
    for (unsigned int i = 0; i < seqstr.length(); i++)
    {
        if (seqstr[i] == 'A')
        {
            res += 'T';
        }
        else if (seqstr[i] == 'T')
        {
            res += 'A';
        }
        else if (seqstr[i] == 'C')
        {
            res += 'G';
        }
        else
        {
            res += 'C';
        }
    }
    return res;
}

std::vector<std::string> fix_cwise(std::vector<std::string> inpvec, std::vector<std::string> nvec) //Function to check for counterclockwise sequences and change to reverse complement
{
    std::vector<std::string> retvec;
    std::string substr = "Counterclockwise";
    for (unsigned int i = 0; i < inpvec.size(); i++)
    {
        if (nvec[i].find(substr) != std::string::npos)
        {
            std::string complem = reverse_comp(inpvec[i]);
            retvec.push_back(complem);
        }
        else
        {
            retvec.push_back(inpvec[i]);
        }
    }
    return retvec;
}

int main()
{
    std::vector<std::string> linevecn, strvecn, seqvecn1, seqvecn2, namevecn, lineveco, strveco, seqveco1, seqveco2, nameveco, testvec, testnames, tlvec, strtest;
    std::vector<bool> bvecn, bveco, bvect;
    std::string filenamen, linen, filenameo, lineo, testfi, tline;
    std::cout << "Please input test data set file in Fasta format [NORF]: ";
    std::cin >> filenamen;
    std::ifstream ifilen(filenamen);
    if (ifilen.is_open())
    {
        while (getline(ifilen, linen))
        {
            std::string fileline = linen;
            fileline = std::regex_replace(fileline, std::regex("\\r\\n|\\r|\\n"), "");
            linevecn.push_back(fileline);
        }
        ifilen.close();
    }
    else
    {
        std::cout << "[ERROR] Error opening test data set file [NORF]" << std::endl;
        exit(EXIT_FAILURE);
    }
    strvecn = make_strs(linevecn);
    bvecn = make_boolvec(strvecn);
    seqvecn1 = make_seqs(strvecn, bvecn);
    namevecn = make_names(strvecn, bvecn);
    seqvecn2 = fix_cwise(seqvecn1, namevecn);
    std::cout << "Please input test data set file in Fasta format [ORF]: ";
    std::cin >> filenameo;
    std::ifstream ifileo(filenameo);
    if (ifileo.is_open())
    {
        while (getline(ifileo, lineo))
        {
            std::string fileline = lineo;
            fileline = std::regex_replace(fileline, std::regex("\\r\\n|\\r|\\n"), "");
            lineveco.push_back(fileline);
        }
        ifileo.close();
    }
    else
    {
        std::cout << "[ERROR] Error opening test data set file [ORF]" << std::endl;
        exit(EXIT_FAILURE);
    }
    strveco = make_strs(lineveco);
    bveco = make_boolvec(strveco);
    nameveco = make_names(strveco, bveco);
    seqveco1 = make_seqs(strveco, bveco);
    seqveco2 = fix_cwise(seqveco1, nameveco);
    std::cout << "Ok, please provide a file to test sequences: ";
    std::cin >> testfi;
    std::ifstream ifilet(testfi);
    if (ifilet.is_open())
    {
        while (getline(ifilet, tline))
        {
            std::string fileline = tline;
            fileline = std::regex_replace(fileline, std::regex("\\r\\n|\\r|\\n"), "");
            tlvec.push_back(fileline);
        }
        ifilet.close();
    }
    else
    {
        std::cout << "[ERROR] Error opening file of sequences to test." << std::endl;
        exit(EXIT_FAILURE);
    }
    strtest = make_strs(tlvec);
    bvect = make_boolvec(strtest);
    testnames = make_names(strtest, bvect);
    testvec = make_seqs(strtest, bvect);
    Markov marky;
    marky.add_seqs(seqvecn2, namevecn, seqveco2, nameveco, testvec, testnames);
    marky.driver_func(); //like a "main()" function within our object, keeps clutter out of main()
    return 0;
}