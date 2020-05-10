#include <bits/stdc++.h>
#include <map>
#include <time.h>

using namespace std;
string s, t, mergedRead, match, remaining;
//transition nodes
struct state
{
    int len, link;
    map<char, int> next;
};

struct result
{
    int bestpos;
    string match;
};

const int MAXLEN = 100000;
state st[MAXLEN * 2];
int sz = 0, last;

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

//initializing root node
void sa_init()
{
    st[0].len = 0;
    st[0].link = -1;
    //cout << "sz: " << sz << endl;
    sz++;
    last = 0;
}

//adding new character to the automata
void sa_extend(char c)
{
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    int p = last;
    while (p != -1 && !st[p].next.count(c))
    {
        st[p].next[c] = cur;
        p = st[p].link;
    }
    if (p == -1)
    {
        st[cur].link = 0;
    }
    else
    {
        int q = st[p].next[c];
        if (st[p].len + 1 == st[q].len)
        {
            st[cur].link = q;
        }
        else
        {
            int clone = sz++;
            st[clone].len = st[p].len + 1;
            st[clone].next = st[q].next;
            st[clone].link = st[q].link;
            while (p != -1 && st[p].next[c] == q)
            {
                st[p].next[c] = clone;
                p = st[p].link;
            }
            st[q].link = st[cur].link = clone;
        }
    }
    last = cur;
}

void collapse(int headpos_R1, int startpos1, int endpos2, int endpos1, int tailpos_R2, int lenS)
{
   // string qScore_r1 = "?????BBBDDDDDDDDFBFFCFEAFECHHHHGFHIHHEHHIHHHHCGH?DGHHHHEFEHH?CBEEHHHHIIIIIIIIIIIIIHIHHHHFGHHHHHHHHHHHHHHHHCCFFHHEEFFFFFFFDBBEDEEEEEFEEFCCBCEAEBBCEEEEFFFFEFFEECEECEEEFEFFFEEEFFEEACCCCECCEEFEEEAEEEFFCCEEEEFF?EAEDDDDEEE::?EECCEEFECEFEEFFE;>DDD?A:CAEEAE?";
    //string qScore_r2 = "<<???@@@BDDDDDDDFFFF7@FCHHHDHEFDFFFHFHFDG@?CGHDEFHHHG-EGFHFHGHEHHHHHHCFFHHHHF@EGBFHHHF?DHHFHHH@CFHHFFFFFFFFF:*=DDEDEEEEEEEE*3:@EEEEEEBBC;E)80AEE))4A?:A?CE:*8AEEEEEEE**?::C?CEA:?0:?:AEACAE?**00:?:*0:*0??AAEEEEEE########################################";

    string qScore_r1 = "~!@#$%^&*()_+:<>?";
    string qScore_r2 = "|}{+_&^%$#";

    reverse(qScore_r2.begin(), qScore_r2.end());

    string qScore_head_r1 = qScore_r1.substr(headpos_R1, startpos1);
    string qScore_head_r2 = qScore_r2.substr(0, startpos1);

    string qScore_tail_r1 = qScore_r1.substr(endpos2, lenS - endpos2);
    string qScore_tail_r2 = qScore_r2.substr(endpos1, tailpos_R2);

    mergedRead = s.substr(0, headpos_R1);
    for (int i = 0; i < qScore_head_r1.size(); i++)
    {
        if (t[i] != s[headpos_R1 + i])
        {
            if (qScore_head_r1[i] < qScore_head_r2[i])
            {
                mergedRead += t[i];
            }
            else
            {
                mergedRead += s[headpos_R1 + i];
            }
        }

        else
            mergedRead += t[i];
    }
    transform(match.begin(), match.end(), match.begin(), ::tolower);
    mergedRead += match;

    for (int i = 0; i < qScore_tail_r1.size(); i++)
    {
        //cout << qScore_tail_r1[i] << " --- " << qScore_tail_r2[i] << endl;
        if (t[endpos1 + i] != s[endpos2 + i])
        {
            if (qScore_tail_r1[i] < qScore_tail_r2[i])
            {
                mergedRead += t[endpos1 + i];
            }
            else
            {
                mergedRead += s[endpos2 + i];
            }
        }
        else
            mergedRead += t[endpos1 + i];
    }

    //cout << "mergedHead4 => " << mergedRead << endl;
    mergedRead += remaining;
    cout << "Final Read => " << mergedRead << endl
         << "Read Length : " << mergedRead.size() << endl;

    // cout << qScore_head_r1 << endl
    //      << qScore_head_r2 << endl;
    // cout << qScore_tail_r1 << endl
    //      << qScore_tail_r2 << endl;
}
//traversing the automata
result lcs(string S, string T)
{
    // for timing execution

    sa_init();
    for (int i = 0; i < S.size(); i++)
        sa_extend(S[i]);

    int v = 0, l = 0, best = 0, bestpos = 0;
    for (int i = 0; i < T.size(); i++)
    {
        while (v && !st[v].next.count(T[i]))
        {
            v = st[v].link;
            l = st[v].len;
        }
        if (st[v].next.count(T[i]))
        {
            v = st[v].next[T[i]];
            l++;
        }
        if (l > best)
        {
            best = l;
            bestpos = i;
        }
    }

    result res;
    res.bestpos = bestpos - best + 1;
    res.match = T.substr(bestpos - best + 1, best);

    return res;
}

int main()
{
    //s = "ATCAACATCATAGCCAGATGCCCAGAGATTAGAGCGCATGACAAGTAAAGGACGGTTGTCAGCGTCATAAGAGGTTTTACCTCCAAATGAAGAAATAACATCATGGTAACGCTGCATGAAGTAATCACGTTCTTGGTCAGTATGCAAATTAGCATAAGCAGCTTGCAGACCCATAATGTCAATAGATGTGGTAGAAGTCGTCATTTGGCGAGAAAGCTCAGTCTCAGGAGGAAGCGGAGCAGTCCAAATG";
    //t = "AGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACGGAGCTTTCTCGCCACATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGCTGTAATT";
    s = "abcdefghijklmux";
    t = "mnjkpqrstu";
    clock_t time;
    time = clock(); // get time of start
    //t = reverse_complement(t);
    result res1 = lcs(s, t);
    match = res1.match;
    cout << "INPUT STRINGS \n"
         << "Read1 : " << s << endl
         << "Read1 Length : " << s.size() << endl
         << endl;
    cout << "Read2 : " << t << endl
         << "Read2 Length : " << t.size() << endl
         << endl
         << "LCS : " << match << endl
         << "LCS Length : " << match.size()
         << endl;
    // clearing previous automata for 's' and preparing it for 't'
    for (int i = 0; i < sz; i++)
        st[i].next.clear();

    result res2 = lcs(t, s); // returns occurence index on r1

    int startpos1 = res1.bestpos; // the first occurence of the longest common substring in 't' / R2
    int startpos2 = res2.bestpos; // the first occurence of LCS in 's' /R1

    int headpos_R1 = startpos2 - startpos1; //start point of head substring in R1

    string head_r1 = s.substr(headpos_R1, startpos1);
    string head_r2 = t.substr(0, startpos1);

    int endpos1 = startpos1 + res1.match.size(); // where the match ends in string t
    int endpos2 = startpos2 + res2.match.size(); // do in string s

    int tailpos_R2 = s.size() - endpos2; // endposition of tail substring in R2

    string tail_r1 = s.substr(endpos2, s.size() - endpos2);
    string tail_r2 = t.substr(endpos1, tailpos_R2);
    remaining = t.substr((endpos1 + tailpos_R2), t.size() - endpos1 - tailpos_R2);
    //cout<< "Remaining end of r2 => " << remaining <<endl;

    //headpos_R1, startpos1, endpos2,endpos1,tailpos_R2, lenS

    collapse(headpos_R1, startpos1, endpos2, endpos1, tailpos_R2, s.size());
    time = clock() - time;
    cout << "\ntotal execution took " << time << "clicks (" << ((float)time) / CLOCKS_PER_SEC << "seconds).\n";
    // cout << head_r1 << endl
    //      << head_r2 << endl
    //      << "---------------------------" << endl
    //      << tail_r1 << endl
    //      << tail_r2 << endl;

    return 0;
}