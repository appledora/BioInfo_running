#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <time.h>

using namespace std;
#define MAX = 300;
int M[300][300];
char M_tb[300][300];
string A, B, A_al = "", B_al = "";
int lenA, lenB;

int align_nuc = 20; // amount of nucleotides per row
int a = 5;          // match
int b = -2;         // purine-purine / pyrimidine-pyrimidine
int c = -5;         // mismatch
int gap = 2;        // initial gap penalty. Gap penalty is lower than mismatch: two sequences from same species assumed.
float gap_ext = 1;  // bigger gap penalties for affine gap penalty

int gap_affinity(int gap, int gap_ext, int &length)
{
    int gap_aff = gap + (gap_ext * length);

    return gap_aff;
}

// Get maximal score and trace it
int max_score(int up, int diag, int left, char *ptr, int &length)
{
    int max = 0;

    if (diag >= up && diag >= left)
    {
        max = diag;
        *ptr = '\\';
        length = 0; // reset gap-length
    }
    else if (up > left)
    {
        max = up;
        *ptr = '|';
        length++; // add 1 gap-length
    }
    else
    {
        max = left;
        *ptr = '-';
        length++; // add 1 gap-length
    }

    return max;
}
// Initialise scoring matrix with first row and column
void init()
{
    M[0][0] = 0;      //scoring matrix
    M_tb[0][0] = 'n'; //direction matrix

    int i = 0, j = 0;

    for (j = 1; j <= lenA; j++)
    {
        M[0][j] = -(gap + (gap_ext * (j - 1))); // manually apply affine gap
        M_tb[0][j] = '-';
    }
    for (i = 1; i <= lenB; i++)
    {
        M[i][0] = -(gap + (gap_ext * (i - 1))); // manually apply affine gap
        M_tb[i][0] = '|';
    }
}

// Needleman and Wunsch algorithm
int alignment()
{
    int x = 0, y = 0;
    int scU, scD, scL; // scores for respectively cell above, diagonal and left
    char ptr, nuc;
    int i = 0, j = 0;
    int length = 0; // initial gap length

    // create substitution scoring matrix
    const int s[4][4] = {
        {a, b, c, c},
        {b, a, c, c},
        {c, c, a, b},
        {c, c, b, a}};

    for (i = 1; i <= lenB; i++)
    {
        for (j = 1; j <= lenA; j++)
        {
            nuc = A[j - 1];

            switch (nuc)
            {
            case 'C':
                x = 0;
                break;
            case 'T':
                x = 1;
                break;
            case 'A':
                x = 2;
                break;
            case 'G':
                x = 3;
            }

            nuc = B[i - 1];

            switch (nuc)
            {
            case 'C':
                y = 0;
                break;
            case 'T':
                y = 1;
                break;
            case 'A':
                y = 2;
                break;
            case 'G':
                y = 3;
            }

            scU = M[i - 1][j] - gap_affinity(gap, gap_ext, length); // get score if trace would go up
            scD = M[i - 1][j - 1] + s[x][y];                        // get score if trace would go diagonal
            scL = M[i][j - 1] - gap_affinity(gap, gap_ext, length); // get score if trace would go left

            M[i][j] = max_score(scU, scD, scL, &ptr, length); // get max score for current optimal global alignment

            M_tb[i][j] = ptr;
        }
    }
    i--;
    j--;

    while (i > 0 || j > 0)
    {
        switch (M_tb[i][j])
        {
        case '|':
            A_al += '-';
            B_al += B[i - 1];
            i--;
            break;

        case '\\':
            A_al += A[j - 1];
            B_al += B[i - 1];
            i--;
            j--;
            break;

        case '-':
            A_al += A[j - 1];
            B_al += '-';
            j--;
        }
    }

    reverse(A_al.begin(), A_al.end());
    reverse(B_al.begin(), B_al.end());

    return 0;
}

// Print the scoring matrix
void print_mtx()
{
    cout << "        ";
    for (int j = 0; j < lenA; j++)
    {
        cout << A[j] << "   ";
    }
    cout << "\n  ";

    for (int i = 0; i <= lenB; i++)
    {
        if (i > 0)
        {
            cout << B[i - 1] << " ";
        }
        for (int j = 0; j <= lenA; j++)
        {
            cout.width(3);
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// Print the traceback matrix
void print_tb()
{
    cout << "        ";
    for (int j = 0; j < lenA; j++)
    {
        cout << A[j] << "   ";
    }
    cout << "\n  ";

    for (int i = 0; i <= lenB; i++)
    {
        if (i > 0)
        {
            cout << B[i - 1] << " ";
        }
        for (int j = 0; j <= lenA; j++)
        {
            cout.width(3);
            cout << M_tb[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int NW(int gap, int gap_ext, int align_nuc)
{

    clock_t t;   // for timing execution
    t = clock(); // get time of start

    // Initialize traceback and F matrix (fill in first row and column)
    init();

    // Create alignment
    alignment();

    int score = M[lenB][lenA]; // get alignment score

    // print_mtx();
    // print_tb();

    cout << endl
         << "Alignments:" << endl;
    int start = 0;            // start of new line for printing alignments
    int cntr = 0;             // iterator for printing alignments
    int Al_n = A_al.length(); // length of alignment
    do
    {
        cout << start + 1 << " A: ";
        for (cntr = start; cntr < start + align_nuc; cntr++)
        {
            if (cntr < Al_n)
            {
                cout << A_al[cntr];
            }
            else
            {
                break;
            }
        }
        cout << " " << cntr << endl
             << start + 1 << " B: ";
        for (cntr = start; cntr < start + align_nuc; cntr++)
        {
            if (cntr < Al_n)
            {
                cout << B_al[cntr];
            }
            else
            {
                break;
            }
        }
        cout << " " << cntr << endl
             << endl;
        start += align_nuc;
    } while (start <= Al_n);

    // Show score and runtime
    t = clock() - t; // get time when finished
    cout << "Alignment score: " << score << endl;
    cout << "Alignment took " << t << " clicks (" << ((float)t) / CLOCKS_PER_SEC << "seconds).\n";

    return 0;
}
string reverse_complement(string T)
{
    reverse(T.begin(), T.end());
    for (int i = 0; i < T.length(); i++)
    {
        switch (T[i])
        {
        case 'A':
            T[i] = 'T';
            break;
        case 'T':
            T[i] = 'A';
            break;
        case 'G':
            T[i] = 'C';
            break;
        case 'C':
            T[i] = 'G';
            break;
        }
    }

    return T;
}
int main()
{
    A = "ATCAACATCATAGCCAGATGCCCAGAGATTAGAGCGCATGACAAGTAAAGGACGGTTGTCAGCGTCATAAGAGGTTTTACCTCCAAATGAAGAAATAACATCATGGTAACGCTGCATGAAGTAATCACGTTCTTGGTCAGTATGCAAATTAGCATAAGCAGCTTGCAGACCCATAATGTCAATAGATGTGGTAGAAGTCGTCATTTGGCGAGAAAGCTCAGTCTCAGGAGGAAGCGGAGCAGTCCAAATG";
    B = "AGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACGGAGCTTTCTCGCCACATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGCTGTAATT";

    B = reverse_complement(B);
    lenA = A.length();
    lenB = B.length();

    NW(gap, gap_ext, align_nuc);
    return 0;
}