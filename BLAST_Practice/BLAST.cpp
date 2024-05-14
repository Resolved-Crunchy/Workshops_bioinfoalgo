//This program was made as a dummy BLAST search tool for a class

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream> 
#include <regex>
#include <map>
#include <sys/time.h>

unsigned long msectime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000 + tv.tv_usec);
}

class database //Name could be workshopped, this will store our kmers in a hashmap and hopefully fascilitate quick lookup
{
    public:
        database(){} //default constructor
        unsigned long microsectime()
        {
            struct timeval tv;
            gettimeofday(&tv, NULL);
            return (tv.tv_sec * 1000000 + tv.tv_usec);
        }
        void setup_db(int size, std::vector<std::string> mainvec) //Takes user defined variables generated in main() and reads them into the database object
        {
            kmsz = size;
            dbvec = mainvec;
            dbveclen = dbvec.size();
            read_in_db();
        }
        void find_kmer() //Finds kmer keys in our query sequence and pushes their values (indexes) back into our vector of pairs of indexes (time intensive)
        {
            std::cout << "Beginning search of database..." << std::endl;
            begint = microsectime();
            for (int i = 0; i <= querysz - kmsz; i++)
            {
                std::pair<std::multimap<std::string, std::pair<int, int>>::iterator, std::multimap<std::string, std::pair<int, int>>::iterator> temp;
                std::string kmqry = kmer_query(i);
                temp = kmerdb.equal_range(kmqry);
                for (std::multimap<std::string, std::pair<int, int>>::iterator it = temp.first; it != temp.second; ++it)
                {
                    std::pair<int, int> indpair;
                    indpair.first = it->second.first;
                    indpair.second = it->second.second;
                    indexvec.push_back(indpair);
                }
                indlen = indexvec.size();
                sort_indvec();
            }
            endt = microsectime();
            time_spent = (double)(endt - begint);
            std::cout << "[Database Searched!] Time elapsed (microseconds): " << time_spent << std::endl;
        }
        void sort_indvec() //Here we implement a simple bubble sort on our index vector to sort by our second value in our vector of pairs
        {
            for (int i = 0; i < indlen - 1; i++)
            {
                for (int j = 0; j < indlen - i - 1; j++)
                {
                    std::pair<int, int> temp = indexvec[j];
                    if (temp.second > indexvec[j + 1].second)
                    {
                        indexvec[j] = indexvec[j + 1];
                        indexvec[j + 1] = temp;
                    }
                }
            }
        }
        std::string kmer_query(int index) //returns a string of kmer size from our query sequence to be searched against our hash table
        {
            std::string resstr;
            for (int i = index; i < (index + kmsz); i++)
            {
                resstr += query[i];
            }
            return resstr;
        }
        void add_query(std::string qry) //Adds query sequence into object
        {
            query = qry;
            querysz = query.length();
        }
        std::string ret_ext_kmer(int indx, int sz) //returns a kmer from our query sequence to be read into our hashmap
        {
            std::string retstr;
            for (int i = indx; i < indx + sz; i++)
            {
                retstr += query[i];
            }
            return retstr;
        }

        void read_in_db() //reads in our sequence database and extracts kmers into a multimap structure
        {
            std::cout << "Ok: Reading in your kmers, please be patient." << std::endl;
            begint = microsectime();
            for (int i = 0; i < dbveclen; i++)
            {
                int totalsz = dbvec[i].length();
                for (int j = 0; j <= (totalsz - kmsz); j++)
                {
                    std::pair<int, int> pospair;
                    pospair.first = j; 
                    pospair.second = i; //Store the index of our database vector as the second part of our pair, so we can retrieve the sequence itself
                    std::string kmerkey = ret_kmer(j, dbvec[i]); //utilize our ret_kmer function to return our kmer string
                    kmerdb.insert({kmerkey, pospair}); //Store our kmer string as a key and our position within the db string as an int "pos"
                }
            }
            endt = microsectime();
            time_spent = (double)(endt - begint);
            std::cout << "[SUCCESS] kmers added to the map! Time Spent (microseconds): " << time_spent << std::endl; //print statements help to see how long it takes
        }
        std::string ret_kmer(int indx, std::string dbstr) //returns a kmer from our database sequences to be read into our hashmap
        {
            std::string retstr;
            int seqlen = dbstr.length();
            for (int i = indx; i < indx + kmsz; i++)
            {
                retstr += dbstr[i];
            }
            return retstr;
        }
        void fill_query_index()
        {
            int len = extvec.size();
            for (int i = 0; i < len; i++)
            {
                std::regex rx(extvec[i]);
                for (auto it = std::sregex_iterator(query.begin(), query.end(), rx); it != std::sregex_iterator(); ++it)
                {
                    int sz = extvec[i].length();
                    int finind = (sz + it->position()) - 1;
                    qposvec.push_back(it->position());
                    qfinind.push_back(finind);
                }
            }
        }
        void fill_match_vec() //Fills a match vector with our string position, also extends kmers by full consecutive length
        {
            std::cout << "Filling match vector and extending kmers..." << std::endl;
            begint = microsectime();
            int count = 1;
            int fstind = indexvec[0].second; //first index of dbvec
            std::string fstkmer = ret_kmer(indexvec[0].first, dbvec[indexvec[0].second]); //first kmer, to be extended
            int firstpos = indexvec[0].first; //first position
            jvec.push_back(firstpos); //push back first position into jvec (may be used actually)
            std::pair<int, int> locpair(indexvec[0].first, indexvec[0].second); //location pair for matchstart (pos, dbind)
            matchstart.push_back(locpair); 
            kmervec.push_back(fstkmer);
            for (int i = 1; i < indlen; i++)
            {
                if (indexvec[i].second != indexvec[i - 1].second || indexvec[i].first != indexvec[i - 1].first + 1) //if dbvec seq != or position != 1 - pos
                {
                    std::pair<int, int> mpair;
                    mpair.first = fstind; //dbvec index
                    mpair.second = count; //count, counts kmers mashed together
                    matchstart.push_back(indexvec[i]); //if at terminator, push back indexvec values
                    matchvec.push_back(mpair); //push back dbvec index and # of consecutive kmers added
                    count = 1; //reset steps
                    fstind = indexvec[i].second;
                    std::string add_kmer = ret_kmer(indexvec[i].first, dbvec[fstind]);
                    firstpos = indexvec[i].first;
                    jvec.push_back(firstpos); //push back first positions to jvec
                    kmervec.push_back(add_kmer); //Add kmer, this loop is a starting loop
                }
                else if (i + 1 == indlen) //to terminate and get last valid value
                {
                    count += 1;
                    std::pair<int, int> mpair;
                    mpair.first = fstind;
                    mpair.second = count;
                    matchvec.push_back(mpair);
                }
                else //adds count to match
                {
                    count++;
                }
            }
            kmlen = kmervec.size();
            endt = microsectime();
            time_spent = (double)(endt - begint);
            std::cout << "[Success!] Kmers extended and matches stored! Time spent (microseconds): " << time_spent << std::endl;
        }

        std::pair<std::string, double> match_percent(std::string dbstr, std::string qmatch)
        {
            std::pair<std::string, double> retpair;
            int len = dbstr.length();
            double mmcount = 0;
            std::string retstr;
            for (int i = 0; i < len; i++)
            {
                if (dbstr[i] == qmatch[i])
                {
                    continue;
                }
                else
                {
                    mmcount += 1;
                    qmatch.insert(i, 1, '*');
                }
            }
            double qleng = qmatch.length();
            double mpct = mmcount / qleng;
            retpair.first = qmatch;
            retpair.second = mpct;
            return retpair;
        }
        void print_svec()
        {
            int len = kmervec.size();
            for (int i = 0; i < len; i++)
            {
                std::cout << "Seq: " << kmervec[i] << std::endl;
            }
        }
        std::string ext_kmer(std::string str, int extension, int ind, int dbind)
        {
            int klen = str.length();
            int start = (ind + klen); 
            std::string dbstr = dbvec[dbind];
            for (int i = start; i < (start + extension) - 1; i++)
            {
                str += dbstr[i];
            }
            int lastpos = (start + extension) - 2;
            lposvec.push_back(lastpos);
            return str;
        }
        std::string ret_db_section(int ind, int seqind, int fin)
        {
            std::string restr;
            std::string dbstr = dbvec[seqind];
            for (int i = ind; i <= fin; i++)
            {
                restr += dbstr[i];
            }
            return restr;
        }
        std::string ret_db_sectionc(int ind, int seqind, int fin)
        {
            std::string restr;
            std::string dbstr = dbvec[seqind];
            for (int i = ind; i <= fin; i++)
            {
                char c = dbstr[i] + 32;
                restr += c;
            }
            return restr;
        }
        void extend_kmers()
        {
            std::cout << "Extending kmers... " << std::endl;
            begint = microsectime();
            for (int i = 0; i < kmlen; i++)
            {
                std::string extended = ext_kmer(kmervec[i], matchvec[i].second, jvec[i], matchvec[i].first);
                extvec.push_back(extended);
            }
            endt = microsectime();
            time_spent = (double)(endt - begint);
            std::cout << "[Success!] Time elapsed (microseconds): " << time_spent << std::endl;
        }
        std::string ret_query_sect(int qpos1, int qpos2)
        {
            std::string restr;
            for (int i = qpos1; i < qpos2; i++)
            {
                restr += query[i];
            }
            return restr;
        }
        void final_kmer_ext()
        {
            std::cout << "Performing final kmer extension and consolidation... " << std::endl;
            begint = microsectime();
            int matches = 0;
            int len = extvec.size();
            bool prevflag = false;
            for (int i = 0; i < len; i++)
            {
                if (i + 1 != len && qposvec[i + 1] - qfinind[i] < 10 && matchstart[i].second == matchstart[i + 1].second && matchstart[i + 1].first - lposvec[i] == qposvec[i + 1] - qfinind[i])  
                {
                    cons += extvec[i];
                    qcons += extvec[i];
                    int mm = qposvec[i + 1] - qfinind[i]; //mismatches
                    std::string querysect = ret_query_sect(qfinind[i], qposvec[i + 1]);
                    std::string dbsect = ret_db_sectionc(lposvec[i], matchstart[i].second, matchstart[i + 1].first);
                    cons += dbsect;
                    qcons += querysect;
                    matches++;
                    prevflag = true;
                }
                else
                {
                    cons += extvec[i];
                    qcons += extvec[i];
                    consvec.push_back(cons);
                    finqvec.push_back(qcons);
                    cons.clear();
                    qcons.clear();
                    if (prevflag == true)
                    {
                        int start = matchstart[i - matches].first;
                        int seqind = matchstart[i].second;
                        int lpos = lposvec[i];
                        std::pair<int, std::pair<int, int>> posgroup;
                        std::pair<int, int> qpospair(qposvec[i - matches], qfinind[i]);
                        queryindex.push_back(qpospair);
                        posgroup.first = seqind;
                        posgroup.second.first = start;
                        posgroup.second.second = lpos;
                        db_indexes.push_back(posgroup);
                        matches = 0;
                        prevflag = false;
                    }
                    else
                    {
                        int start = matchstart[i].first;
                        int seqind = matchstart[i].second;
                        int lpos = lposvec[i];
                        std::pair<int, std::pair<int, int>> positions;
                        std::pair<int, int> qpospair(qposvec[i], qfinind[i]);
                        queryindex.push_back(qpospair);
                        positions.first = seqind;
                        positions.second.first = start;
                        positions.second.second = lpos;
                        db_indexes.push_back(positions);
                        matches = 0;
                        prevflag = false;
                    }
                }
            }
            endt = microsectime();
            time_spent = (double)(endt - begint);
            std::cout << "[SUCCESS!] Final consolidations perfrmed and data stored. Time elapsed (microseconds): " << time_spent << std::endl;
        }
        void print_ext()
        {
            int len = extvec.size();
            for (int i = 0; i < len; i++)
            {
                std::cout << "ext seq: " << extvec[i] << std::endl;
                std::cout << "Match.first: " << matchstart[i].first << "   Match.second: " << matchstart[i].second << std::endl;
            }
        }
        int mismatch_count(std::string str)
        {
            int count = 0;
            int len = str.length();
            for (int i = 0; i < len; i++)
            {
                if (str[i] == 'a' || str[i] == 'g' || str[i] == 't' || str[i] == 'c')
                {
                    count++;
                }
                else
                {
                    continue;
                }
            }
            return count;
        }
        double full_match(int matchsz)
        {
            double qz = querysz;
            double mz = matchsz;
            return (mz / qz) * 100;
        }
        double query_match(int matchsize, int subtsz)
        {
            double mz = matchsize;
            double sz = subtsz;
            return (sz / mz) * 100;
        }
        void results_print()
        {
            std::cout << std::endl;
            std::cout << "Printing results and query coverage: " << std::endl;
            begint = microsectime();
            std::cout << std::endl;
            std::cout << "Complete Query sequence: " << query << std::endl;
            std::cout << std::endl;
            std::cout << "[NOTE] Bases in lower case indicate a mismatch between database string and query string. " << std::endl;
            int len = consvec.size(); //choice to use consvec here is arbitrary, all vectors in question should share a size
            for (int i = 0; i < len; i++)
            {
                std::cout << std::endl;
                std::cout << "[Match found at DB seq_" << db_indexes[i].first + 1 << " From positions: " << db_indexes[i].second.first + 1 << " to " << db_indexes[i].second.second + 1 << "." << std::endl;
                std::string db_str = ret_db_section(db_indexes[i].second.first, db_indexes[i].first, db_indexes[i].second.second);
                std::cout << "Db Sequence:        " << db_str << std::endl;
                std::cout << std::endl;
                std::cout << "Consensus sequence: " << consvec[i] << std::endl;
                std::cout << std::endl;
                std::cout << "Query match from positions: " << queryindex[i].first + 1 << " to " << queryindex[i].second + 1 << std::endl;
                std::cout << "Query sequence:     " << finqvec[i] << std::endl;
                std::cout << std::endl;
                int matchsize = finqvec[i].length();
                int mm = mismatch_count(consvec[i]);
                int ident = consvec[i].size() - mm;
                int fmatchsize = matchsize - ident;
                double qmatch = query_match(matchsize, ident);
                double fmatch = full_match(ident);
                std::cout << qmatch << "% Identities (" << ident << "/" << consvec[i].size() << "),  " << fmatch << "% match to query overall.]" << std::endl;
                std::cout << std::endl;
            }
            endt = microsectime();
            time_spent = (double)(endt - begint);
            std::cout << "End of Results! Time elapsed (microsecond): " << time_spent << std::endl;
        }

    private:
        std::vector<std::string> dbvec, kmervec, extvec, resvec, consvec, finqvec; //Vector containing the strings from our input database file
        std::vector<std::pair<int, std::pair<int, int>>> db_indexes;
        std::vector<std::pair<int, int>> indexvec, matchvec, matchstart, gapvec, matchstart2, matchvec2, queryindex; //vector of pairs holding our indexes, vector of indexes + length of consecutive characters to be added.
        int kmsz, dbveclen, querysz, indlen, kmlen, mtqlen; //kmersize, size of our database vector of strings, length of our query sequence, and length of index vector
        std::multimap<std::string, std::pair<int, int>> kmerdb; //Hashmap utilizing our kmers as a string key and our position in the string as an int
        std::multimap<std::string, std::pair<int, int>> extkmer; //Multimap holding our extended kmer sequences
        std::string query, cons, qcons; //query sequence and consensus sequence
        std::vector<std::pair<std::string, int>> reskmer_sz; //vector containing extended kmers + size
        std::vector<int> jvec, lposvec, qposvec, qfinind; //vector containing our positions on our db strings
        long begint, endt;
        double time_spent;
};

std::vector<bool> ret_bvec(std::vector<char> buffvec)
{
    std::vector<bool> resvec;
    for (unsigned int i = 0; i < buffvec.size(); i++)
    {
        if (buffvec[i] == 'A' || buffvec[i] == 'T' || buffvec[i] == 'G' || buffvec[i] == 'C')
        {
            resvec.push_back(true);
        }
        else
        {
            resvec.push_back(false);
        }
    }
    return resvec;
}

std::vector<std::string> makedbvec(std::vector<char> buffvec, std::vector<bool> boolvec)
{
    std::string dbstr;
    std::vector<std::string> dbvec;
    for (unsigned int i = 0; i < buffvec.size(); i++)
    {
        if (boolvec[i] == true && i + 1 != buffvec.size() && boolvec[i + 1] == true)
        {
            dbstr += buffvec[i];
        }
        else if (boolvec[i] == true && i + 1 != buffvec.size() && boolvec[i + 1] == false)
        {
            dbstr += buffvec[i];
            dbvec.push_back(dbstr);
            dbstr.clear();
        }
        else if (i + 1 == buffvec.size())
        {
            dbstr += buffvec[i];
            dbvec.push_back(dbstr);
            break;
        }
        else
        {
            continue;
        }
    }
    return dbvec;
}

std::string pre_processquery(std::vector<char> buffvec) //This function seems to have fixed our string processing related issues
{
    std::string res;
    for (unsigned int i = 0; i < buffvec.size(); i++)
    {
        if (buffvec[i] == 'A' || buffvec[i] == 'T' || buffvec[i] == 'G' || buffvec[i] == 'C')
        {
            res += buffvec[i];
        }
        else
        {
            continue;
        }
    }
    return res;
}

int main()
{
    std::regex r("\\s"); //regular expression, eliminates whitespace
    std::regex seq("[ATCG]");
    long begintime, endtime, btime, e_time;
    double t_spent;
    int kmersz; //user defined kmer size
    std::cout << "Hello! What size would you like your kmers to be? ";
    std::cin >> kmersz;
    std::vector<std::string> seqdb; //vector to hold our database sequences in main
    std::vector<bool> bvec;
    std::string filestr, line, inpam, pamstr, query, qline, qnam, queryseq; //strings to read from files and to access files primarily
    std::cout << "Ok, please provide a database file in Fasta format: ";
    std::cin >> filestr;
    std::ifstream ifile(filestr);
    if (!ifile)
    {
        std::cout << "[ERROR] error opening sequence database file, terminating program." << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        btime = msectime();
        begintime = msectime();
        std::cout << "Ok, beginning to process from file: " << std::endl;
        std::istream_iterator<char> begin(ifile), end;
        std::vector<char> buffer(begin, end);
        bvec = ret_bvec(buffer);
        seqdb = makedbvec(buffer, bvec);
        e_time = msectime();
        t_spent = (double)(e_time - begintime);
        std::cout << "Database file read in. Time elapsed (microseconds): " << t_spent << std::endl;
        std::cout << "Ok, please provide a Query sequence file in Fasta format: ";
        std::cin >> query;
        std::cout << "Ok, beginning to read in your query now... " << std::endl;
        begintime = msectime();
        std::ifstream quefile(query);
        if (!quefile)
        {
            std::cout << "[ERROR], error with Query file" << std::endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            std::istream_iterator<char> begin(quefile), end;
            std::vector<char> buffer1(begin, end);
            std::string qry = pre_processquery(buffer1);
            endtime = msectime();
            t_spent = (double)(endtime - begintime);
            std::cout << "Query read in and processed! Time elapsed: " << t_spent << std::endl;
            database sdb;
            sdb.setup_db(kmersz, seqdb); 
            sdb.add_query(qry); 
            sdb.find_kmer(); //finds kmers in our query string that are present in our database and reads their indeces in the sequence database as a pair of ints into a vector 
            sdb.fill_match_vec();
            std::cout << "Ok, now printing our unextended kmers at their first found position: " << std::endl;
            sdb.print_svec();
            sdb.extend_kmers(); 
            std::cout << std::endl;
            std::cout << "Final search results: ";
            sdb.fill_query_index();
            sdb.final_kmer_ext();
            sdb.results_print();
            endtime = msectime();
            t_spent = (double)(endtime - btime);
            std::cout << "Total time elapsed (microseconds): " << t_spent << std::endl;
        }
    }
    return 0;
}