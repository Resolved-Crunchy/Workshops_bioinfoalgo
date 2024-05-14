//This is an alignment algorithm made for a class
//This algorithm should work for nucleic acid sequences in FASTA format, but will do nothing for Amino Acid sequences

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

int aff_gap_pen(int ingap, int exgap, int& gap_len){ //This function returns our affine gap penalty
    int affgap = ingap + (exgap * gap_len);
    return affgap;
}

int max_score(int up, int dia, int left, int& gap_len){ //This function returns the max between our upper, left, and diagonal scores and increments our gap length if a gap exists
    int max = 0;
    if (dia >= up && dia >= left){
        max = dia;
        gap_len = 0;
    }
    else if (up > dia && up >= left){
        max = up;
        gap_len += 1;
    }
    else if (left > dia && left >= up){
        max = left;
        gap_len += 1;
    }
    return max;
}

bool is_gap(int up, int dia, int left){ //This function tells us whether the score is a gap or a mismatch. If a gap, it returns true. Else it's false
    if (dia >= up && dia >= left){
        return false;
    }
    else
    {
        return true;
    }
}

bool g_type(int up, int left){ //This function sets where the gap is in our alignment. If true, it's on the top string
    if (up > left){
        return false;
    }
    else
    {
        return true;
    }
}

void align_tab(std::string seq1, std::string seq2, int mm_p, int gap_p, int ex_p, int match) //This function is supposed to perform the alignment then output 2 strings using std::cout to show alignment
{
    int i, j, k, l; //Iterators
    int u_sc, d_sc, l_sc; //Score values 
    int mat_x, mat_y; //Values for traversing our diagonal score matrix
    int sz1 = seq1.length(); //Initialize these szn ints to our input seqeunce lengths 
    int sz2 = seq2.length();
    std::string aligned_1 = ""; //I set these as empty strings because I couldn't get it to add anything otherwise
    std::string aligned_2 = "";
    char nu_tide1; //Stores our nucleotides at a given index
    char nu_tide2;
    int gap_len = 0;
    bool g_o_no, g_uol;
    char gapchar = '-';
    const int mat_mat[4][4] = {{match, mm_p, mm_p, mm_p},
                               {mm_p, match, mm_p, mm_p},
                               {mm_p, mm_p, match, mm_p},
                               {mm_p, mm_p, mm_p, match}}; //Scoring matrix for our diagonal scores. If mat_x == mat_y, we get our match score. Otherwise, we get a mismatch score
    int tab[sz1 + 1][sz2 + 1];
    tab[0][0] = 0; //Initialize our scoring table
    for (i = 1; i < sz1; i++)
    {
        for (j = 1; j <= sz2; j++){
            tab[0][i] = gap_p + (ex_p * (i - 1));
            tab[j][0] = gap_p + (ex_p * (j - 1));
        } 
    } //These for loops initialize the table by filling the first row and column with gaps
    for (k = 1; k <= sz1; k++){
        for (l = 1; l <= sz2; l++){
            nu_tide1 = seq1[k - 1];
            switch (nu_tide1)
            {
            case 'A':
                mat_x = 0;
                break;
            case 'T':
                mat_x = 1;
                break;
            case 'G':
                mat_x = 2;
                break;
            case 'C':
                mat_x = 3;
            }
            nu_tide2 = seq2[l - 1];
            switch (nu_tide2)
            {
            case 'A':
                mat_y = 0;
                break;
            case 'T':
                mat_y = 1;
                break;
            case 'G':
                mat_y = 2;
                break;
            case 'C':
                mat_y = 3;
            } //Switch statements employed here in order to return our score from our mat_matrix above
            u_sc = tab[k - 1][l] + aff_gap_pen(gap_p, ex_p, gap_len); //u_sc and l_sc are both gaps, so they are scored the same way save for their indeces
            d_sc = tab[k - 1][l - 1] + mat_mat[mat_x][mat_y]; //diagonal score can be either a match or a mismatch, so we score accordingly
            l_sc = tab[k][l - 1] + aff_gap_pen(gap_p, ex_p, gap_len);
            tab[k][l] = max_score(u_sc, d_sc, l_sc, gap_len); //Fill the spot we are at in with the max value of our above three integers
            g_o_no = is_gap(u_sc, d_sc, l_sc); //g_o_no = gap or not
            if (g_o_no == true){
                g_uol = g_type(u_sc, l_sc); //g_uol = upper or lower
                if (g_uol == true){
                    aligned_2 += seq2[l -1];
                    aligned_1 += gapchar;
                }
                else
                {
                    aligned_1 += seq1[k -1];
                    aligned_2 += gapchar;
                }
            }
            else
            {
                aligned_1 += seq1[k - 1];
                aligned_2 += seq2[l - 1];
            }
        }
    }
    std::cout << aligned_1 << std::endl;
    std::cout << aligned_2 << std::endl;
    return;
}

int main()
{
    int ingap, exgap, mm_p, mtch; //Initialize our initial gap penalty, our extenison penalty, our mismatch penalty, and our match score respectively
    bool pr_nuc_choice;
    std::string filestr, line;
    std::string seq_1, seq_2;
    std::cout << "Are we aligning a protein or nucleic acid sequence? (0 = Protein, 1 = Nucleic Acids) ";
    std::cin >> pr_nuc_choice; //Boolean functiont to select nucleic acids or proteins
    if (pr_nuc_choice == true){
        std::string snam_1, snam_2;
        std::cout << "Please provide a sequence file to read: (FASTA) ";
        std::cin >> filestr;
        std::ifstream ifile(filestr);
        if (ifile.is_open())
        {
            getline(ifile, line);
            snam_1 = line;
            snam_1.erase(snam_1.begin());
            getline(ifile, line);
            seq_1 = line;
            getline(ifile, line);
            snam_2 = line;
            snam_2.erase(snam_2.begin());
            getline(ifile, line);
            seq_2 = line;
        }
        std::cout << "Please input an initial gap penalty: ";
        std::cin >> ingap;
        std::cout << "Thank you, please input an extension penalty: ";
        std::cin >> exgap;
        std::cout << "Thank you, please input a mismatch penalty: ";
        std::cin >> mm_p;
        std::cout << "Thank you, almost done, please input a positive match score: ";
        std::cin >> mtch;
        std::cout << " Thank you, here are your sequences: " << std::endl;
        std::cout << snam_1 << ": " << seq_1 << std::endl;
        std::cout << snam_2 << ": " << seq_2 << std::endl;
        align_tab(seq_1, seq_2, mm_p, ingap, exgap, mtch);
        return 0;
    }
    else
    {
        std::cout << "This application does not support Amino Acid alignment, exiting now."
    }

    return 0;
}