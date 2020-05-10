#include <bits/stdc++.h>
using namespace std;
string x, y;
int scoreM[300][300];
char dirM[300][300];
int gap = -1;
string align_x = "";
string align_y = "";
int max_score = 0, max_i, max_j;
int align_score = 0;
int matchScore(char a, char b)
{
    if (a == b)
    {
        return 1;
    }
    else
        return 0;
}
void local_align(int lenX, int lenY)
{
    memset(scoreM, 0, sizeof(scoreM));
    memset(dirM, '.', sizeof(dirM));

    for (int i = 1; i <= lenY; i++)
    {
        for (int j = 1; j <= lenX; j++)
        {
            int s_diag = scoreM[i - 1][j - 1] + matchScore(y[i - 1], x[j - 1]);
            int s_left = scoreM[i][j - 1] + gap;
            int s_up = scoreM[i - 1][j] + gap;

            scoreM[i][j] = max(s_diag, max(s_left, max(s_up, 0)));
            if (s_diag >= s_up)
            {
                if (s_diag >= s_left)
                {
                    scoreM[i][j] = s_diag;
                    dirM[i][j] = 'd';
                }
                else
                {
                    scoreM[i][j] = s_left;
                    dirM[i][j] = 'l';
                }
            }
            else
            {
                if (s_up >= s_left)
                {
                    scoreM[i][j] = s_up;
                    dirM[i][j] = 'u';
                }
                else
                {
                    scoreM[i][j] = s_left;
                    dirM[i][j] = 'l';
                }
            }

            if (scoreM[i][j] > max_score)
            {
                max_score = scoreM[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
}
void traceback()
{
    {
        int i = max_i;
        int j = max_j;
        cout << i << " " << j << endl;
        while (true)
        {
            cout << dirM[i][j] << endl;
            cout << x[j - 1] << "-----" << y[i - 1] << endl
                 << endl;
            if (dirM[i][j] == '.')
            {
                break;
            }

            else if (dirM[i][j] == 'd')
            {
                align_x += x[j - 1];
                align_y += y[i - 1];
                if (x[j - 1] == y[i - 1])
                {
                    align_score += 1;
                }
                i--;
                j--;
            }

            else if (dirM[i][j] == 'l')
            {
                align_x += x[j - 1];
                align_y += '-';
                align_score += gap;
                cout<<"NOT DIAGONAL" <<endl;
                j--;
            }

            else if (dirM[i][j] == 'u')
            {
                align_x += '-';
                align_y += y[i - 1];
                align_score += gap;
                cout<<"NOT DIAGONAL" <<endl;
                i--;
            }
            cout << "score => " << align_score << endl;
        }

        reverse(align_x.begin(), align_x.end());
        reverse(align_y.begin(), align_y.end());

        cout << align_x << endl
             << align_y << endl;
    }
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

   x = "ATCAACATCATAGCCAGATGCCCAGAGATTAGAGCGCATGACAAGTAAAGGACGGTTGTCAGCGTCATAAGAGGTTTTACCTCCAAATGAAGAAATAACATCATGGTAACGCTGCATGAAGTAATCACGTTCTTGGTCAGTATGCAAATTAGCATAAGCAGCTTGCAGACCCATAATGTCAATAGATGTGGTAGAAGTCGTCATTTGGCGAGAAAGCTCAGTCTCAGGAGGAAGCGGAGCAGTCCAAATG";
    y = "AGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACGGAGCTTTCTCGCCACATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGCTGTAATT";
    y = reverse_complement(y);

    //x = "AACCCTA";
    //y = "AGCCTT";
    int lenX = x.size();
    int lenY = y.size();

    cout << lenX << " " << lenY << endl;
    local_align(lenX, lenY);

    // for (int i = 0; i <= lenY; i++)
    // {
    //     for (int j = 0; j <= lenX; j++)
    //     {
    //         cout << dirM[i][j] << " ";
    //     }

    //     cout << endl;
    // }
    // cout << endl;
    // for (int i = 0; i <= lenY; i++)
    // {
    //     for (int j = 0; j <= lenX; j++)
    //     {
    //         cout << scoreM[i][j] << " ";
    //     }

    //     cout << endl;
    // }

    traceback();
    cout << "Max Score : " << max_score << endl;
}