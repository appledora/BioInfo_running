

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <bits/stdc++.h>

using namespace std;

int ind;
double penalty = -4;

double similarityScore(char a, char b);
double findMax(double array[], int length);

double similarityScore(char a, char b)
{
    double result;
    if (a == b)
    {
        result = 1;
    }
    else
    {
        result = penalty;
    }
    return result;
}

double findMax(double array[], int length)
{
    double max = array[0];
    ind = 0;

    for (int i = 0; i < length; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
            ind = i;
        }
    }
    return max;
}
string reverse_complement(string T)
{
    reverse(T.begin(), T.end());
    // cout << "reversed string \n"
    //      << T << endl;
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
    string seqA = "ATCAACATCATAGCCAGATGCCCAGAGATTAGAGCGCATGACAAGTAAAGGACGGTTGTCAGCGTCATAAGAGGTTTTACCTCCAAATGAAGAAATAACATCATGGTAACGCTGCATGAAGTAATCACGTTCTTGGTCAGTATGCAAATTAGCATAAGCAGCTTGCAGACCCATAATGTCAATAGATGTGGTAGAAGTCGTCATTTGGCGAGAAAGCTCAGTCTCAGGAGGAAGCGGAGCAGTCCAAATG"; // sequence A
    string seqB = "AGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACGGAGCTTTCTCGCCACATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGCTGTAATT";  // sequence B
    seqB = reverse_complement(seqB);
    // cout << "Sequence A" << endl;
    // cin >> seqA;
    // cout << "Sequence B" << endl;
    // cin >> seqB;
    cout << "You typed in " << endl
         << seqA << endl
         << seqB << endl;

    // initialize some variables
    int lengthSeqA = seqA.length();
    int lengthSeqB = seqB.length();

    // initialize matrix
    double matrix[lengthSeqA + 1][lengthSeqB + 1];
    for (int i = 0; i <= lengthSeqA; i++)
    {
        for (int j = 0; j <= lengthSeqB; j++)
        {
            matrix[i][j] = 0;
        }
    }

    double traceback[4];
    int I_i[lengthSeqA + 1][lengthSeqB + 1];
    int I_j[lengthSeqA + 1][lengthSeqB + 1];

    //start populating matrix
    for (int i = 1; i <= lengthSeqA; i++)
    {
        for (int j = 0; j <= lengthSeqB; j++)
        {
            //cout << i << " " << j << endl;
            traceback[0] = matrix[i - 1][j - 1] + similarityScore(seqA[i - 1], seqB[j - 1]);
            traceback[1] = matrix[i - 1][j] + penalty;
            traceback[2] = matrix[i][j - 1] + penalty;
            traceback[3] = 0;
            matrix[i][j] = findMax(traceback, 4);
            switch (ind)
            {
            case 0:
                I_i[i][j] = i - 1;
                I_j[i][j] = j - 1;
                break;
            case 1:
                I_i[i][j] = i - 1;
                I_j[i][j] = j;
                break;
            case 2:
                I_i[i][j] = i;
                I_j[i][j] = j - 1;
                break;
            case 3:
                I_i[i][j] = i;
                I_j[i][j] = j;
                break;
            }
        }
    }

    // // print the scoring matrix to console
    // for(int i=1;i<lengthSeqA;i++)
    // {
    // 	for(int j=1;j<lengthSeqB;j++)
    // 	{
    // 		cout << I_i[i][j] << " ";
    // 	}
    // 	cout << endl;
    // }
    // cout<<endl<<endl;
// for(int i=1;i<lengthSeqA;i++)
//     {
//     	for(int j=1;j<lengthSeqB;j++)
//     	{
//     		cout << I_j[i][j] << " ";
//     	}
//     	cout << endl;
//     }
    // find the max score in the matrix
    double matrix_max = 0;
    int i_max = 0, j_max = 0;
    for (int i = 1; i < lengthSeqA; i++)
    {
        for (int j = 1; j < lengthSeqB; j++)
        {
            if (matrix[i][j] > matrix_max)
            {
                matrix_max = matrix[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }

    cout << "Max score in the matrix is " << matrix_max << endl;

    // traceback

    int current_i = i_max  , current_j = j_max; // location of max_value
    int next_i = I_i[current_i][current_j];
    int next_j = I_j[current_i][current_j];
    int tick = 0;
    //char consensus_a[lengthSeqA + lengthSeqB + 2], consensus_b[lengthSeqA + lengthSeqB + 2];
    string consensus_a, consensus_b;

    while (((current_i != next_i) || (current_j != next_j)) && (next_j != 0) && (next_i != 0))
    {
        //cout << "current_i : " << current_i << "\ncurrent_j  : " << current_j << "\nnext_i  : " << next_i << "\nnext_j : " << next_j << endl;
        if (next_i == current_i)
        {
            consensus_a += '-'; // deletion in A
           // cout << "delA => " << consensus_a;
        }
        else
        {
            consensus_a += seqA[current_i - 1]; // match/mismatch in A
            //cout << "consA => " << consensus_a << "\n seqAchar : " << seqA[current_i - 1] << endl;
        }
        if (next_j == current_j)
        {
            consensus_b += '-'; // deletion in B
            //cout << "delB => " << consensus_b;
        }

        else
        {
            consensus_b += seqB[current_j - 1];
            //cout << "consB => " << consensus_b << "\n seqBchar : " << seqB[current_j - 1] << endl;
        } // match/mismatch in B

        current_i = next_i;
        current_j = next_j;
        next_i = I_i[current_i][current_j];
        next_j = I_j[current_i][current_j];
        tick++;
    }

    //print the consensus sequences
    cout << endl
         << " " << endl;
    cout << "Alignment:" << endl
         << endl;
    for (int i = 0; i < lengthSeqA; i++)
    {
        cout << seqA[i];
    };
    cout << "  and" << endl;
    for (int i = 0; i < lengthSeqB; i++)
    {
        cout << seqB[i];
    };
    cout << endl
         << consensus_a
         << endl
         << consensus_b << endl;

    // for (int i = consensus_a.size(); i >= 0; i--)
    //     cout << consensus_a[i];
    // cout << endl;
    // for (int j = consensus_b.size(); j >= 0; j--)
    //     cout << consensus_b[j];
    // cout << endl;

    return 0;
}
