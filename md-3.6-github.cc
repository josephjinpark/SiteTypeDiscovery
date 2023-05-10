// vi:foldmethod=marker

// CHANGE LOG
// ----------
//
// md-3.0
// - It is the integrated version of md-1.0 and md-2.0.
//   (cst1 command for md-1.0 and cst2 command for md-2.0)
// - add <in_dir_lam> and <in_dir_raw> arguments.
// - add assumption for chi-square tests on contingency tables:
//   all expected counts should be greater than 5.
//
// md-3.1
// - add targeting-only RESIDUES. (md-3.0: 22 residues, md-3.1: 32 residues)
// - rename variable 'bitcount' to 'ndc' (ndc: the number of dont-care residues)
// - detailed output format.
// - discard amo when it canot pass either expected counts filtering or relative
//   risk filtering in any cutoff.
//
// md-3.2
// - modify RESIDUES for easy understanding.
//
// md-3.3
// - put interaction residues into init_ir()
// - available motif length: 10, 12, 14, 16, 18mer
// - the number of contingency tables will be used: cont2, cont5
// - add filtering option 'minps': minimum of motif presence
// - add Fisher's exact test (command: fet1)
// 
// md-3.4
// - Solved 'O2 problem'
// - Add 'iac3' type residues A_, A1, A2, A3, A4, ..., C_, C1, C2, C3, C4, ...
// - Debugged (Segmentation fault in KISTI tachyon2)
//
// md-3.5
// - Simplified output format. It will reduce a size of output file.
// - Soft-coded variable 'ncont'(the number of contingency table).
//   Variable 'ncont' can be one of 0, 1, 2, ..., 9.
// - New argument 'psonly', which can be used when you want to see
//   motifs with positive enrichment p-values over all cutoffs.
//
// md-3.6
// - Revise output format. It will provide minscore as well as maxscore
// - Throw out the motifs with maxscore less than 0 (maxscore < 0)
//


#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "fisher2.h"

#define MOTIF      char*
#define MAX_SCORE  324





using namespace std;

/* struct */
struct Interaction {
    int  u;     // uamo index
    int* ipos;  // interacting positions on 3'UTR
};

struct Fcm {
    char   nmid[32];  // gene id
    char   maid[32];  // microarray id
    int    b;         // bin index
    char** lam;       // local au index
    vector<struct Interaction> vec_ia;
};

struct Chi2t {
    int b;               // bin index
    int n1, n2, n3, n4;  // counts in contingency table
    bool   rrsk;         // is relative risk greater than 1.0? (yes:1, no:0)
    double chi2;         // chi-square statistic
    double schi2;        // signed chi2
    // double pval;         // p-value from chi-square test
    // double logp;         // -log10(p-value)
    // double score;        // enrichment score
};

struct Fet {
    int b;
    int n1, n2, n3, n4;
    bool rrsk;
    double pval;
    double logp;
};


struct Ir {
    char main;  // main interaction residue
    char side;  // side interaction residue
};


typedef struct Ir* Uamo;


/* global variables */

int molen          = 0;  // motif length
int half_mask_size = 0;  // half mask size
int full_mask_size = 0;  // full mask size

vector<string>        vec_ir;
map<string,struct Ir> map_ir;

vector<int>* vec_halfmask;
int*         ndc;
vector<Uamo> vec_uamo;
Uamo tmp_uamo1;
Uamo tmp_uamo2;


vector<struct Fcm> vec_fc[10];
map<string,char**> map_lam;  // local au matrix



#ifdef DEBUG
// for time estimation
struct timeval tvs[2];
struct timeval tve[2];
#endif



// {{{ common

int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1) result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

vector<string> &split(const string &s, char delim, vector<string> &elems)
{
    // for split the string with given delimteres
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim))	{ elems.push_back(item); }
	return elems;
}

vector<string> split(const string &s, char delim)
{
	vector<string> elems;
	return split(s, delim, elems);
}


double clocksec() { return (double)clock() / CLOCKS_PER_SEC; }


// }}}
// {{{ init_map_ir

inline bool inc_check_by_mask(int target, int cand) { return (cand & target) == target; }

void init_map_ir(void)
{

    // LIST OF INTERACTION RESIDUES //

    vec_ir.clear();

    vec_ir.push_back("._"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {0, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");

    vec_ir.push_back("O_"); vec_ir.push_back("O1"); vec_ir.push_back("O2"); vec_ir.push_back("O3"); vec_ir.push_back("O4");  // {1, 0}
    vec_ir.push_back("O5"); vec_ir.push_back("O6"); vec_ir.push_back("O7"); vec_ir.push_back("O8"); vec_ir.push_back("O9");
    
    vec_ir.push_back("A_"); vec_ir.push_back("A1"); vec_ir.push_back("A2"); vec_ir.push_back("A3"); vec_ir.push_back("A4");  // {2, 0}
    vec_ir.push_back("A5"); vec_ir.push_back("A6"); vec_ir.push_back("A7"); vec_ir.push_back("A8"); vec_ir.push_back("A9");
    
    vec_ir.push_back("C_"); vec_ir.push_back("C1"); vec_ir.push_back("C2"); vec_ir.push_back("C3"); vec_ir.push_back("C4");  // {3, 0}
    vec_ir.push_back("C5"); vec_ir.push_back("C6"); vec_ir.push_back("C7"); vec_ir.push_back("C8"); vec_ir.push_back("C9");
    
    vec_ir.push_back("X_"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {4, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");

    vec_ir.push_back("D_"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {5, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    
    vec_ir.push_back("W_"); vec_ir.push_back("W1"); vec_ir.push_back("W2"); vec_ir.push_back("W3"); vec_ir.push_back("W4");  // {6, 0}
    vec_ir.push_back("W5"); vec_ir.push_back("W6"); vec_ir.push_back("W7"); vec_ir.push_back("W8"); vec_ir.push_back("W9");
    
    vec_ir.push_back("B1"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {7, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B2"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {8, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B3"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {9, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B4"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {10, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B5"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {11, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B6"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {12, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B7"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {13, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B8"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {14, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");
    vec_ir.push_back("B9"); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");  // {15, 0}
    vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  "); vec_ir.push_back("  ");

    for (int i=0; i<(int)vec_ir.size(); i++)
    {
        if ( vec_ir[i] == string("  ") ) continue;
        map_ir[vec_ir[i]].main = i/10;
        map_ir[vec_ir[i]].side = i%10;
    }

}

// }}}
// {{{ load_uamolist

void init_ndc(int molen, int full_mask_size)
{
    ndc = (int*)malloc(sizeof(int)*full_mask_size);
    memset(ndc, 0, sizeof(int)*full_mask_size);
    
    // m: mask index
    for (int m=0; m < full_mask_size; m++)
    {
        for (int i=0; i < molen; i++)
        {
            if ( (m >> i) & 1 ) ndc[m] += 1;
        }
    }
}

void init_vec_halfmask(int half_mask_size)
{
    vec_halfmask = (vector<int>*)malloc(sizeof(vector<int>)*half_mask_size);

	for (int i=0; i < half_mask_size; i++)
	{
        vector<int> vec_halfmask_i;
        vec_halfmask[i] = vec_halfmask_i;
		vec_halfmask[i].clear();

		for (int j=0; j < half_mask_size; j++)
		{
			if (inc_check_by_mask(i,j) == true)
			{
				vec_halfmask[i].push_back(j);
			}
		}
	}
}

void load_uamolist(char* in_file_uamo)
{
    #ifdef DEBUG
    gettimeofday(&tvs[0], NULL);
    #endif


	ifstream fin;
	string   line;
	
	fin.open(in_file_uamo, ios::in);
	if ( !fin.is_open() ) { fprintf(stderr, "error: cannot open file '%s'\n", in_file_uamo); exit(1); }

    for (int i=0; i<(int)vec_uamo.size(); i++) free(vec_uamo[i]);
    vec_uamo.clear();

	while ( getline(fin, line) )
	{
        if (molen == 0)
        {
            molen = strlen(line.c_str())/2;
            if (not (molen==10 or molen==12 or molen==14 or molen==16 or molen==18 or molen==20))
            {
                fprintf(stderr, "error: cannot support the motif length '%d'\n", molen);
                exit(1);
            }
            full_mask_size = ipow(2, molen);
            half_mask_size = ipow(2, molen/2);
            init_ndc(molen, full_mask_size);
            init_vec_halfmask(half_mask_size);
        }  // added in md-3.3
        

        Uamo uamo;
        uamo = (struct Ir*)malloc(sizeof(struct Ir)*molen);
        memset(uamo, 0, sizeof(struct Ir)*molen);

        for (int i=0; i<molen; i++) uamo[i] = map_ir[line.substr(i*2,2)];
		vec_uamo.push_back(uamo);
        
	}
	fin.close();
    
    
    tmp_uamo1 = (struct Ir*)malloc(sizeof(struct Ir)*molen);
    memset(tmp_uamo1, 0, sizeof(struct Ir)*molen);

    tmp_uamo2 = (struct Ir*)malloc(sizeof(struct Ir)*molen);
    memset(tmp_uamo2, 0, sizeof(struct Ir)*molen);


    #ifdef DEBUG
    gettimeofday(&tve[0], NULL); timersub(&tve[0], &tvs[0], &tve[0]);
    printf("load_uamolist      %4d.%06d sec\n", tve[0].tv_sec, tve[0].tv_usec);
    #endif

}

// }}}
// {{{ load_fc


void load_fc_lam(struct Fcm* fc, char* in_dir_lam)
{
    char in_file_lam[1024];
    ifstream fin;
    int utrlen = 0;

    if (map_lam.find(fc->nmid) != map_lam.end())
    {
        fc->lam = map_lam[fc->nmid]; //
    }
    else
    {
        sprintf(in_file_lam, "%s/%s.dat", in_dir_lam, fc->nmid);
        fin.open(in_file_lam, ios::in);
        if ( !fin.is_open() )
        {
            fc->lam = NULL; //
        }
        else
        {
            fin.seekg(0, ios_base::end);
            utrlen = sqrt(fin.tellg() / sizeof(char));
            fin.seekg(0, ios_base::beg);

            map_lam[fc->nmid] = (char**)malloc(sizeof(char*)*utrlen);
            for (int s=0; s<utrlen; s++)
            {
                map_lam[fc->nmid][s] = (char*)malloc(sizeof(char)*utrlen);
                fin.read(map_lam[fc->nmid][s], sizeof(char)*utrlen);
            }
            fc->lam = map_lam[fc->nmid]; //
        }
        fin.close();
    }
}


void load_fc_raw(struct Fcm* fc, char* in_dir_raw)
{
    ifstream fin;
    char in_file_raw[1024];

    sprintf(in_file_raw, "%s/%s/%s.dat", in_dir_raw, fc->maid, fc->nmid);
    fin.open(in_file_raw, ios::in);
    if ( !fin.is_open() ) { fprintf(stderr, "error: cannot open file '%s'\n", in_file_raw); exit(1); }
    fc->vec_ia.clear();
    while (true)
    {
        struct Interaction ia;
        ia.ipos = (int*)malloc(sizeof(int)*molen);
        fin.read((char*)(&ia.u), sizeof(int)      );
        fin.read((char*)ia.ipos, sizeof(int)*molen);
        if ( fin.eof() ) break;
        fc->vec_ia.push_back(ia);
    }
    fin.close();

}




void load_fc(char* in_file_fc, char* in_dir_lam, char* in_dir_raw)
{
    
    #ifdef DEBUG
    gettimeofday(&tvs[0], NULL); struct timeval tv1, tv2; timerclear(&tv1); timerclear(&tv2);
    #endif

	ifstream fin;
    string nmid;  // nmid: gene id
    string maid;  // maid: micro array id
	string line;
    vector<string> col;

    fin.open(in_file_fc, ios::in);
    if ( !fin.is_open() ) { fprintf(stderr, "error: cannot open file '%s'\n", in_file_fc); exit(1); }
    while ( getline(fin, line) )
    {

        // [0] init
        struct Fcm fc;
        col = split(line, '\t');
        strcpy(fc.nmid, col[0].c_str());  // fc.nmid
        strcpy(fc.maid, col[1].c_str());  // fc.maid
        fc.b   =   atoi(col[2].c_str());  // fc.b


        // [1] load lam(local au matrix) file
        #ifdef DEBUG
        gettimeofday(&tvs[1], NULL);
        #endif
        if ( in_dir_lam != NULL ) load_fc_lam(&fc, in_dir_lam);  // fc.lam
        #ifdef DEBUG
        gettimeofday(&tve[1], NULL); timersub(&tve[1], &tvs[1], &tve[1]); timeradd(&tv1, &tve[1], &tv1);
        #endif


        // [2] load raw file
        #ifdef DEBUG
        gettimeofday(&tvs[1], NULL);
        #endif
        load_fc_raw(&fc, in_dir_raw);  // fc.vec_ia
        #ifdef DEBUG
        gettimeofday(&tve[1], NULL); timersub(&tve[1], &tvs[1], &tve[1]); timeradd(&tv2, &tve[1], &tv2);
        #endif


        // [3] push_back
        vec_fc[fc.b].push_back(fc);
    }
    fin.close();

    #ifdef DEBUG
    gettimeofday(&tve[0], NULL); timersub(&tve[0], &tvs[0], &tve[0]);
    printf("load_fc            %4d.%06d sec\n", tve[0].tv_sec, tve[0].tv_usec);
    printf("- step1(load lam)  %4d.%06d sec\n", tv1.tv_sec, tv1.tv_usec);
    printf("- step2(load raw)  %4d.%06d sec\n", tv2.tv_sec, tv2.tv_usec);
    #endif

}


// }}}

// {{{ cst common

string masking(int u, int m)
{   
    int offset;
    int ir;
    string amo;

    for (int i=0; i<molen; i++)
    {
        offset = (molen-1) - i;

        if ( (m >> offset) & 1 )
        {
            ir = 0;
        }
        else
        {
            ir = vec_uamo[u][i].main * 10;

            if ( offset > 0 )
            {
                if ( (m >> (offset-1)) & 1 ) { ir += 0;                   }
                else                         { ir += vec_uamo[u][i].side; }
            }
            else                             { ir += vec_uamo[u][i].side; }
        }

        amo.append(vec_ir[ir]);
    }
    return amo;
}

double chi2stat(int n1, int n2, int n3, int n4, bool yc)
{
    int sum;
    int rsum[2];
    int csum[2];
    double e[4];
    double chi2;
    
    sum = n1 + n2 + n3 + n4;
    
    rsum[0] = n1 + n2;
    rsum[1] = n3 + n4;
    csum[0] = n1 + n3;
    csum[1] = n2 + n4;
    
    e[0] = double(rsum[0]) * double(csum[0]) / double(sum);
    e[1] = double(rsum[0]) * double(csum[1]) / double(sum);
    e[2] = double(rsum[1]) * double(csum[0]) / double(sum);
    e[3] = double(rsum[1]) * double(csum[1]) / double(sum);

    if (yc)
    {
        // with Yates correction
        chi2  = pow(fabs(n1-e[0])-0.5, 2) / e[0];
        chi2 += pow(fabs(n2-e[1])-0.5, 2) / e[1];
        chi2 += pow(fabs(n3-e[2])-0.5, 2) / e[2];
        chi2 += pow(fabs(n4-e[3])-0.5, 2) / e[3];
    }
    else
    {
        // without Yates correction
        chi2  = pow(fabs(n1-e[0]), 2) / e[0];
        chi2 += pow(fabs(n2-e[1]), 2) / e[1];
        chi2 += pow(fabs(n3-e[2]), 2) / e[2];
        chi2 += pow(fabs(n4-e[3]), 2) / e[3];
    }
	
	return chi2;
}


// }}}
// {{{ count_contbl


// {{{ 10mer-specific functions

int get_top_bit_10(int x) { return (x & 0x03E0) >> 5; }  // top    5 bit   1111100000  0x03E0
int get_bot_bit_10(int x) { return (x & 0x001F);      }  // bottom 5 bit   0000011111  0x001F
int merge_top_bot_10(int x, int y) { return (x << 5) | y; }

int get_maskd1_10(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;

    if ( uamo1[0].main != uamo2[0].main ) { mask=mask|(1<<9);                                   uamo1[0].side=0; uamo2[0].side=0; }
    if ( uamo1[1].main != uamo2[1].main ) { mask=mask|(1<<8); uamo1[0].side=0; uamo2[0].side=0; uamo1[1].side=0; uamo2[1].side=0; }
    if ( uamo1[2].main != uamo2[2].main ) { mask=mask|(1<<7); uamo1[1].side=0; uamo2[1].side=0; uamo1[2].side=0; uamo2[2].side=0; }
    if ( uamo1[3].main != uamo2[3].main ) { mask=mask|(1<<6); uamo1[2].side=0; uamo2[2].side=0; uamo1[3].side=0; uamo2[3].side=0; }
    if ( uamo1[4].main != uamo2[4].main ) { mask=mask|(1<<5); uamo1[3].side=0; uamo2[3].side=0; uamo1[4].side=0; uamo2[4].side=0; }
    if ( uamo1[5].main != uamo2[5].main ) { mask=mask|(1<<4); uamo1[4].side=0; uamo2[4].side=0; uamo1[5].side=0; uamo2[5].side=0; }
    if ( uamo1[6].main != uamo2[6].main ) { mask=mask|(1<<3); uamo1[5].side=0; uamo2[5].side=0; uamo1[6].side=0; uamo2[6].side=0; }
    if ( uamo1[7].main != uamo2[7].main ) { mask=mask|(1<<2); uamo1[6].side=0; uamo2[6].side=0; uamo1[7].side=0; uamo2[7].side=0; }
    if ( uamo1[8].main != uamo2[8].main ) { mask=mask|(1<<1); uamo1[7].side=0; uamo2[7].side=0; uamo1[8].side=0; uamo2[8].side=0; }
    if ( uamo1[9].main != uamo2[9].main ) { mask=mask|(1<<0); uamo1[8].side=0; uamo2[8].side=0; uamo1[9].side=0; uamo2[9].side=0; }

    return mask;
}

int get_maskd2_10(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;
    
    if ( uamo1[0].side != uamo2[0].side ) { mask=mask|(1<<9); }
    if ( uamo1[1].side != uamo2[1].side ) { mask=mask|(1<<8); }
    if ( uamo1[2].side != uamo2[2].side ) { mask=mask|(1<<7); }
    if ( uamo1[3].side != uamo2[3].side ) { mask=mask|(1<<6); }
    if ( uamo1[4].side != uamo2[4].side ) { mask=mask|(1<<5); }
    if ( uamo1[5].side != uamo2[5].side ) { mask=mask|(1<<4); }
    if ( uamo1[6].side != uamo2[6].side ) { mask=mask|(1<<3); }
    if ( uamo1[7].side != uamo2[7].side ) { mask=mask|(1<<2); }
    if ( uamo1[8].side != uamo2[8].side ) { mask=mask|(1<<1); }

    return mask;

}

// }}}
// {{{ 12mer-specific functions

int get_top_bit_12(int x) { return (x & 0x0FC0) >> 6; }  // top    6 bit   111111000000 0x0FC0
int get_bot_bit_12(int x) { return (x & 0x003F);      }  // bottom 6 bit   000000111111 0x003F
int merge_top_bot_12(int x, int y) { return (x << 6) | y; }

int get_maskd1_12(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;

    if ( uamo1[ 0].main != uamo2[ 0].main ) { mask=mask|(1<<11);                                     uamo1[ 0].side=0; uamo2[ 0].side=0; }
    if ( uamo1[ 1].main != uamo2[ 1].main ) { mask=mask|(1<<10); uamo1[ 0].side=0; uamo2[ 0].side=0; uamo1[ 1].side=0; uamo2[ 1].side=0; }
    if ( uamo1[ 2].main != uamo2[ 2].main ) { mask=mask|(1<< 9); uamo1[ 1].side=0; uamo2[ 1].side=0; uamo1[ 2].side=0; uamo2[ 2].side=0; }
    if ( uamo1[ 3].main != uamo2[ 3].main ) { mask=mask|(1<< 8); uamo1[ 2].side=0; uamo2[ 2].side=0; uamo1[ 3].side=0; uamo2[ 3].side=0; }
    if ( uamo1[ 4].main != uamo2[ 4].main ) { mask=mask|(1<< 7); uamo1[ 3].side=0; uamo2[ 3].side=0; uamo1[ 4].side=0; uamo2[ 4].side=0; }
    if ( uamo1[ 5].main != uamo2[ 5].main ) { mask=mask|(1<< 6); uamo1[ 4].side=0; uamo2[ 4].side=0; uamo1[ 5].side=0; uamo2[ 5].side=0; }
    if ( uamo1[ 6].main != uamo2[ 6].main ) { mask=mask|(1<< 5); uamo1[ 5].side=0; uamo2[ 5].side=0; uamo1[ 6].side=0; uamo2[ 6].side=0; }
    if ( uamo1[ 7].main != uamo2[ 7].main ) { mask=mask|(1<< 4); uamo1[ 6].side=0; uamo2[ 6].side=0; uamo1[ 7].side=0; uamo2[ 7].side=0; }
    if ( uamo1[ 8].main != uamo2[ 8].main ) { mask=mask|(1<< 3); uamo1[ 7].side=0; uamo2[ 7].side=0; uamo1[ 8].side=0; uamo2[ 8].side=0; }
    if ( uamo1[ 9].main != uamo2[ 9].main ) { mask=mask|(1<< 2); uamo1[ 8].side=0; uamo2[ 8].side=0; uamo1[ 9].side=0; uamo2[ 9].side=0; }
    if ( uamo1[10].main != uamo2[10].main ) { mask=mask|(1<< 1); uamo1[ 9].side=0; uamo2[ 9].side=0; uamo1[10].side=0; uamo2[10].side=0; }
    if ( uamo1[11].main != uamo2[11].main ) { mask=mask|(1<< 0); uamo1[10].side=0; uamo2[10].side=0; uamo1[11].side=0; uamo2[11].side=0; }

    return mask;
}

int get_maskd2_12(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;
    
    if ( uamo1[ 0].side != uamo2[ 0].side ) { mask=mask|(1<<11); }
    if ( uamo1[ 1].side != uamo2[ 1].side ) { mask=mask|(1<<10); }
    if ( uamo1[ 2].side != uamo2[ 2].side ) { mask=mask|(1<< 9); }
    if ( uamo1[ 3].side != uamo2[ 3].side ) { mask=mask|(1<< 8); }
    if ( uamo1[ 4].side != uamo2[ 4].side ) { mask=mask|(1<< 7); }
    if ( uamo1[ 5].side != uamo2[ 5].side ) { mask=mask|(1<< 6); }
    if ( uamo1[ 6].side != uamo2[ 6].side ) { mask=mask|(1<< 5); }
    if ( uamo1[ 7].side != uamo2[ 7].side ) { mask=mask|(1<< 4); }
    if ( uamo1[ 8].side != uamo2[ 8].side ) { mask=mask|(1<< 3); }
    if ( uamo1[ 9].side != uamo2[ 9].side ) { mask=mask|(1<< 2); }
    if ( uamo1[10].side != uamo2[10].side ) { mask=mask|(1<< 1); }

    return mask;

}

// }}}
// {{{ 14mer-specific functions

int get_top_bit_14(int x) { return (x & 0x3F80) >> 7; }  // top    7 bit   11111110000000 0x3F80
int get_bot_bit_14(int x) { return (x & 0x007F);      }  // bottom 7 bit   00000001111111 0x007F
int merge_top_bot_14(int x, int y) { return (x << 7) | y; }

int get_maskd1_14(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;

    if ( uamo1[ 0].main != uamo2[ 0].main ) { mask=mask|(1<<13);                                     uamo1[ 0].side=0; uamo2[ 0].side=0; }
    if ( uamo1[ 1].main != uamo2[ 1].main ) { mask=mask|(1<<12); uamo1[ 0].side=0; uamo2[ 0].side=0; uamo1[ 1].side=0; uamo2[ 1].side=0; }
    if ( uamo1[ 2].main != uamo2[ 2].main ) { mask=mask|(1<<11); uamo1[ 1].side=0; uamo2[ 1].side=0; uamo1[ 2].side=0; uamo2[ 2].side=0; }
    if ( uamo1[ 3].main != uamo2[ 3].main ) { mask=mask|(1<<10); uamo1[ 2].side=0; uamo2[ 2].side=0; uamo1[ 3].side=0; uamo2[ 3].side=0; }
    if ( uamo1[ 4].main != uamo2[ 4].main ) { mask=mask|(1<< 9); uamo1[ 3].side=0; uamo2[ 3].side=0; uamo1[ 4].side=0; uamo2[ 4].side=0; }
    if ( uamo1[ 5].main != uamo2[ 5].main ) { mask=mask|(1<< 8); uamo1[ 4].side=0; uamo2[ 4].side=0; uamo1[ 5].side=0; uamo2[ 5].side=0; }
    if ( uamo1[ 6].main != uamo2[ 6].main ) { mask=mask|(1<< 7); uamo1[ 5].side=0; uamo2[ 5].side=0; uamo1[ 6].side=0; uamo2[ 6].side=0; }
    if ( uamo1[ 7].main != uamo2[ 7].main ) { mask=mask|(1<< 6); uamo1[ 6].side=0; uamo2[ 6].side=0; uamo1[ 7].side=0; uamo2[ 7].side=0; }
    if ( uamo1[ 8].main != uamo2[ 8].main ) { mask=mask|(1<< 5); uamo1[ 7].side=0; uamo2[ 7].side=0; uamo1[ 8].side=0; uamo2[ 8].side=0; }
    if ( uamo1[ 9].main != uamo2[ 9].main ) { mask=mask|(1<< 4); uamo1[ 8].side=0; uamo2[ 8].side=0; uamo1[ 9].side=0; uamo2[ 9].side=0; }
    if ( uamo1[10].main != uamo2[10].main ) { mask=mask|(1<< 3); uamo1[ 9].side=0; uamo2[ 9].side=0; uamo1[10].side=0; uamo2[10].side=0; }
    if ( uamo1[11].main != uamo2[11].main ) { mask=mask|(1<< 2); uamo1[10].side=0; uamo2[10].side=0; uamo1[11].side=0; uamo2[11].side=0; }
    if ( uamo1[12].main != uamo2[12].main ) { mask=mask|(1<< 1); uamo1[11].side=0; uamo2[11].side=0; uamo1[12].side=0; uamo2[12].side=0; }
    if ( uamo1[13].main != uamo2[13].main ) { mask=mask|(1<< 0); uamo1[12].side=0; uamo2[12].side=0; uamo1[13].side=0; uamo2[13].side=0; }

    return mask;
}

int get_maskd2_14(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;
    
    if ( uamo1[ 0].side != uamo2[ 0].side ) { mask=mask|(1<<13); }
    if ( uamo1[ 1].side != uamo2[ 1].side ) { mask=mask|(1<<12); }
    if ( uamo1[ 2].side != uamo2[ 2].side ) { mask=mask|(1<<11); }
    if ( uamo1[ 3].side != uamo2[ 3].side ) { mask=mask|(1<<10); }
    if ( uamo1[ 4].side != uamo2[ 4].side ) { mask=mask|(1<< 9); }
    if ( uamo1[ 5].side != uamo2[ 5].side ) { mask=mask|(1<< 8); }
    if ( uamo1[ 6].side != uamo2[ 6].side ) { mask=mask|(1<< 7); }
    if ( uamo1[ 7].side != uamo2[ 7].side ) { mask=mask|(1<< 6); }
    if ( uamo1[ 8].side != uamo2[ 8].side ) { mask=mask|(1<< 5); }
    if ( uamo1[ 9].side != uamo2[ 9].side ) { mask=mask|(1<< 4); }
    if ( uamo1[10].side != uamo2[10].side ) { mask=mask|(1<< 3); }
    if ( uamo1[11].side != uamo2[11].side ) { mask=mask|(1<< 2); }
    if ( uamo1[12].side != uamo2[12].side ) { mask=mask|(1<< 1); }

    return mask;

}


// }}}
// {{{ 16mer-specific functions

int get_top_bit_16(int x) { return (x & 0xFF00) >> 8; }  // top    8 bit   1111111100000000   0xFF00
int get_bot_bit_16(int x) { return (x & 0x00FF);      }  // bottom 8 bit   0000000011111111   0x00FF
int merge_top_bot_16(int x, int y) { return (x << 8) | y; }

int get_maskd1_16(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;

    if ( uamo1[ 0].main != uamo2[ 0].main ) { mask=mask|(1<<15);                                     uamo1[ 0].side=0; uamo2[ 0].side=0; }
    if ( uamo1[ 1].main != uamo2[ 1].main ) { mask=mask|(1<<14); uamo1[ 0].side=0; uamo2[ 0].side=0; uamo1[ 1].side=0; uamo2[ 1].side=0; }
    if ( uamo1[ 2].main != uamo2[ 2].main ) { mask=mask|(1<<13); uamo1[ 1].side=0; uamo2[ 1].side=0; uamo1[ 2].side=0; uamo2[ 2].side=0; }
    if ( uamo1[ 3].main != uamo2[ 3].main ) { mask=mask|(1<<12); uamo1[ 2].side=0; uamo2[ 2].side=0; uamo1[ 3].side=0; uamo2[ 3].side=0; }
    if ( uamo1[ 4].main != uamo2[ 4].main ) { mask=mask|(1<<11); uamo1[ 3].side=0; uamo2[ 3].side=0; uamo1[ 4].side=0; uamo2[ 4].side=0; }
    if ( uamo1[ 5].main != uamo2[ 5].main ) { mask=mask|(1<<10); uamo1[ 4].side=0; uamo2[ 4].side=0; uamo1[ 5].side=0; uamo2[ 5].side=0; }
    if ( uamo1[ 6].main != uamo2[ 6].main ) { mask=mask|(1<< 9); uamo1[ 5].side=0; uamo2[ 5].side=0; uamo1[ 6].side=0; uamo2[ 6].side=0; }
    if ( uamo1[ 7].main != uamo2[ 7].main ) { mask=mask|(1<< 8); uamo1[ 6].side=0; uamo2[ 6].side=0; uamo1[ 7].side=0; uamo2[ 7].side=0; }
    if ( uamo1[ 8].main != uamo2[ 8].main ) { mask=mask|(1<< 7); uamo1[ 7].side=0; uamo2[ 7].side=0; uamo1[ 8].side=0; uamo2[ 8].side=0; }
    if ( uamo1[ 9].main != uamo2[ 9].main ) { mask=mask|(1<< 6); uamo1[ 8].side=0; uamo2[ 8].side=0; uamo1[ 9].side=0; uamo2[ 9].side=0; }
    if ( uamo1[10].main != uamo2[10].main ) { mask=mask|(1<< 5); uamo1[ 9].side=0; uamo2[ 9].side=0; uamo1[10].side=0; uamo2[10].side=0; }
    if ( uamo1[11].main != uamo2[11].main ) { mask=mask|(1<< 4); uamo1[10].side=0; uamo2[10].side=0; uamo1[11].side=0; uamo2[11].side=0; }
    if ( uamo1[12].main != uamo2[12].main ) { mask=mask|(1<< 3); uamo1[11].side=0; uamo2[11].side=0; uamo1[12].side=0; uamo2[12].side=0; }
    if ( uamo1[13].main != uamo2[13].main ) { mask=mask|(1<< 2); uamo1[12].side=0; uamo2[12].side=0; uamo1[13].side=0; uamo2[13].side=0; }
    if ( uamo1[14].main != uamo2[14].main ) { mask=mask|(1<< 1); uamo1[13].side=0; uamo2[13].side=0; uamo1[14].side=0; uamo2[14].side=0; }
    if ( uamo1[15].main != uamo2[15].main ) { mask=mask|(1<< 0); uamo1[14].side=0; uamo2[14].side=0; uamo1[15].side=0; uamo2[15].side=0; }

    return mask;
}

int get_maskd2_16(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;
    
    if ( uamo1[ 0].side != uamo2[ 0].side ) { mask=mask|(1<<15); }
    if ( uamo1[ 1].side != uamo2[ 1].side ) { mask=mask|(1<<14); }
    if ( uamo1[ 2].side != uamo2[ 2].side ) { mask=mask|(1<<13); }
    if ( uamo1[ 3].side != uamo2[ 3].side ) { mask=mask|(1<<12); }
    if ( uamo1[ 4].side != uamo2[ 4].side ) { mask=mask|(1<<11); }
    if ( uamo1[ 5].side != uamo2[ 5].side ) { mask=mask|(1<<10); }
    if ( uamo1[ 6].side != uamo2[ 6].side ) { mask=mask|(1<< 9); }
    if ( uamo1[ 7].side != uamo2[ 7].side ) { mask=mask|(1<< 8); }
    if ( uamo1[ 8].side != uamo2[ 8].side ) { mask=mask|(1<< 7); }
    if ( uamo1[ 9].side != uamo2[ 9].side ) { mask=mask|(1<< 6); }
    if ( uamo1[10].side != uamo2[10].side ) { mask=mask|(1<< 5); }
    if ( uamo1[11].side != uamo2[11].side ) { mask=mask|(1<< 4); }
    if ( uamo1[12].side != uamo2[12].side ) { mask=mask|(1<< 3); }
    if ( uamo1[13].side != uamo2[13].side ) { mask=mask|(1<< 2); }
    if ( uamo1[14].side != uamo2[14].side ) { mask=mask|(1<< 1); }

    return mask;

}

// }}}
// {{{ 18mer-specific functions

int get_top_bit_18(int x) { return (x & 0x03FE00) >> 9; }  // top    9 bit   111111111000000000   0x03FE00
int get_bot_bit_18(int x) { return (x & 0x0001FF);      }  // bottom 9 bit   000000000111111111   0x0001FF
int merge_top_bot_18(int x, int y) { return (x << 9) | y; }


int get_maskd1_18(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;

    if ( uamo1[ 0].main != uamo2[ 0].main ) { mask=mask|(1<<17);                                     uamo1[ 0].side=0; uamo2[ 0].side=0; }
    if ( uamo1[ 1].main != uamo2[ 1].main ) { mask=mask|(1<<16); uamo1[ 0].side=0; uamo2[ 0].side=0; uamo1[ 1].side=0; uamo2[ 1].side=0; }
    if ( uamo1[ 2].main != uamo2[ 2].main ) { mask=mask|(1<<15); uamo1[ 1].side=0; uamo2[ 1].side=0; uamo1[ 2].side=0; uamo2[ 2].side=0; }
    if ( uamo1[ 3].main != uamo2[ 3].main ) { mask=mask|(1<<14); uamo1[ 2].side=0; uamo2[ 2].side=0; uamo1[ 3].side=0; uamo2[ 3].side=0; }
    if ( uamo1[ 4].main != uamo2[ 4].main ) { mask=mask|(1<<13); uamo1[ 3].side=0; uamo2[ 3].side=0; uamo1[ 4].side=0; uamo2[ 4].side=0; }
    if ( uamo1[ 5].main != uamo2[ 5].main ) { mask=mask|(1<<12); uamo1[ 4].side=0; uamo2[ 4].side=0; uamo1[ 5].side=0; uamo2[ 5].side=0; }
    if ( uamo1[ 6].main != uamo2[ 6].main ) { mask=mask|(1<<11); uamo1[ 5].side=0; uamo2[ 5].side=0; uamo1[ 6].side=0; uamo2[ 6].side=0; }
    if ( uamo1[ 7].main != uamo2[ 7].main ) { mask=mask|(1<<10); uamo1[ 6].side=0; uamo2[ 6].side=0; uamo1[ 7].side=0; uamo2[ 7].side=0; }
    if ( uamo1[ 8].main != uamo2[ 8].main ) { mask=mask|(1<< 9); uamo1[ 7].side=0; uamo2[ 7].side=0; uamo1[ 8].side=0; uamo2[ 8].side=0; }
    if ( uamo1[ 9].main != uamo2[ 9].main ) { mask=mask|(1<< 8); uamo1[ 8].side=0; uamo2[ 8].side=0; uamo1[ 9].side=0; uamo2[ 9].side=0; }
    if ( uamo1[10].main != uamo2[10].main ) { mask=mask|(1<< 7); uamo1[ 9].side=0; uamo2[ 9].side=0; uamo1[10].side=0; uamo2[10].side=0; }
    if ( uamo1[11].main != uamo2[11].main ) { mask=mask|(1<< 6); uamo1[10].side=0; uamo2[10].side=0; uamo1[11].side=0; uamo2[11].side=0; }
    if ( uamo1[12].main != uamo2[12].main ) { mask=mask|(1<< 5); uamo1[11].side=0; uamo2[11].side=0; uamo1[12].side=0; uamo2[12].side=0; }
    if ( uamo1[13].main != uamo2[13].main ) { mask=mask|(1<< 4); uamo1[12].side=0; uamo2[12].side=0; uamo1[13].side=0; uamo2[13].side=0; }
    if ( uamo1[14].main != uamo2[14].main ) { mask=mask|(1<< 3); uamo1[13].side=0; uamo2[13].side=0; uamo1[14].side=0; uamo2[14].side=0; }
    if ( uamo1[15].main != uamo2[15].main ) { mask=mask|(1<< 2); uamo1[14].side=0; uamo2[14].side=0; uamo1[15].side=0; uamo2[15].side=0; }
    if ( uamo1[16].main != uamo2[16].main ) { mask=mask|(1<< 1); uamo1[15].side=0; uamo2[15].side=0; uamo1[16].side=0; uamo2[16].side=0; }
    if ( uamo1[17].main != uamo2[17].main ) { mask=mask|(1<< 0); uamo1[16].side=0; uamo2[16].side=0; uamo1[17].side=0; uamo2[17].side=0; }

    return mask;
}

int get_maskd2_18(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;
    
    if ( uamo1[ 0].side != uamo2[ 0].side ) { mask=mask|(1<<17); }
    if ( uamo1[ 1].side != uamo2[ 1].side ) { mask=mask|(1<<16); }
    if ( uamo1[ 2].side != uamo2[ 2].side ) { mask=mask|(1<<15); }
    if ( uamo1[ 3].side != uamo2[ 3].side ) { mask=mask|(1<<14); }
    if ( uamo1[ 4].side != uamo2[ 4].side ) { mask=mask|(1<<13); }
    if ( uamo1[ 5].side != uamo2[ 5].side ) { mask=mask|(1<<12); }
    if ( uamo1[ 6].side != uamo2[ 6].side ) { mask=mask|(1<<11); }
    if ( uamo1[ 7].side != uamo2[ 7].side ) { mask=mask|(1<<10); }
    if ( uamo1[ 8].side != uamo2[ 8].side ) { mask=mask|(1<< 9); }
    if ( uamo1[ 9].side != uamo2[ 9].side ) { mask=mask|(1<< 8); }
    if ( uamo1[10].side != uamo2[10].side ) { mask=mask|(1<< 7); }
    if ( uamo1[11].side != uamo2[11].side ) { mask=mask|(1<< 6); }
    if ( uamo1[12].side != uamo2[12].side ) { mask=mask|(1<< 5); }
    if ( uamo1[13].side != uamo2[13].side ) { mask=mask|(1<< 4); }
    if ( uamo1[14].side != uamo2[14].side ) { mask=mask|(1<< 3); }
    if ( uamo1[15].side != uamo2[15].side ) { mask=mask|(1<< 2); }
    if ( uamo1[16].side != uamo2[16].side ) { mask=mask|(1<< 1); }

    return mask;

}

// }}}
// {{{ 20mer-specific functions

int get_top_bit_20(int x) { return (x & 0x0FFC00) >> 10;   }  // top    10 bit   11111111110000000000   0x0FFC00
int get_bot_bit_20(int x) { return (x & 0x0003FF);         }  // bottom 10 bit   00000000001111111111   0x0003FF
int merge_top_bot_20(int x, int y) { return (x << 10) | y; }

int get_maskd1_20(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;

    if ( uamo1[ 0].main != uamo2[ 0].main ) { mask=mask|(1<<19);                                     uamo1[ 0].side=0; uamo2[ 0].side=0; }
    if ( uamo1[ 1].main != uamo2[ 1].main ) { mask=mask|(1<<18); uamo1[ 0].side=0; uamo2[ 0].side=0; uamo1[ 1].side=0; uamo2[ 1].side=0; }
    if ( uamo1[ 2].main != uamo2[ 2].main ) { mask=mask|(1<<17); uamo1[ 1].side=0; uamo2[ 1].side=0; uamo1[ 2].side=0; uamo2[ 2].side=0; }
    if ( uamo1[ 3].main != uamo2[ 3].main ) { mask=mask|(1<<16); uamo1[ 2].side=0; uamo2[ 2].side=0; uamo1[ 3].side=0; uamo2[ 3].side=0; }
    if ( uamo1[ 4].main != uamo2[ 4].main ) { mask=mask|(1<<15); uamo1[ 3].side=0; uamo2[ 3].side=0; uamo1[ 4].side=0; uamo2[ 4].side=0; }
    if ( uamo1[ 5].main != uamo2[ 5].main ) { mask=mask|(1<<14); uamo1[ 4].side=0; uamo2[ 4].side=0; uamo1[ 5].side=0; uamo2[ 5].side=0; }
    if ( uamo1[ 6].main != uamo2[ 6].main ) { mask=mask|(1<<13); uamo1[ 5].side=0; uamo2[ 5].side=0; uamo1[ 6].side=0; uamo2[ 6].side=0; }
    if ( uamo1[ 7].main != uamo2[ 7].main ) { mask=mask|(1<<12); uamo1[ 6].side=0; uamo2[ 6].side=0; uamo1[ 7].side=0; uamo2[ 7].side=0; }
    if ( uamo1[ 8].main != uamo2[ 8].main ) { mask=mask|(1<<11); uamo1[ 7].side=0; uamo2[ 7].side=0; uamo1[ 8].side=0; uamo2[ 8].side=0; }
    if ( uamo1[ 9].main != uamo2[ 9].main ) { mask=mask|(1<<10); uamo1[ 8].side=0; uamo2[ 8].side=0; uamo1[ 9].side=0; uamo2[ 9].side=0; }
    if ( uamo1[10].main != uamo2[10].main ) { mask=mask|(1<< 9); uamo1[ 9].side=0; uamo2[ 9].side=0; uamo1[10].side=0; uamo2[10].side=0; }
    if ( uamo1[11].main != uamo2[11].main ) { mask=mask|(1<< 8); uamo1[10].side=0; uamo2[10].side=0; uamo1[11].side=0; uamo2[11].side=0; }
    if ( uamo1[12].main != uamo2[12].main ) { mask=mask|(1<< 7); uamo1[11].side=0; uamo2[11].side=0; uamo1[12].side=0; uamo2[12].side=0; }
    if ( uamo1[13].main != uamo2[13].main ) { mask=mask|(1<< 6); uamo1[12].side=0; uamo2[12].side=0; uamo1[13].side=0; uamo2[13].side=0; }
    if ( uamo1[14].main != uamo2[14].main ) { mask=mask|(1<< 5); uamo1[13].side=0; uamo2[13].side=0; uamo1[14].side=0; uamo2[14].side=0; }
    if ( uamo1[15].main != uamo2[15].main ) { mask=mask|(1<< 4); uamo1[14].side=0; uamo2[14].side=0; uamo1[15].side=0; uamo2[15].side=0; }
    if ( uamo1[16].main != uamo2[16].main ) { mask=mask|(1<< 3); uamo1[15].side=0; uamo2[15].side=0; uamo1[16].side=0; uamo2[16].side=0; }
    if ( uamo1[17].main != uamo2[17].main ) { mask=mask|(1<< 2); uamo1[16].side=0; uamo2[16].side=0; uamo1[17].side=0; uamo2[17].side=0; }
    if ( uamo1[18].main != uamo2[18].main ) { mask=mask|(1<< 1); uamo1[17].side=0; uamo2[17].side=0; uamo1[18].side=0; uamo2[18].side=0; }
    if ( uamo1[19].main != uamo2[19].main ) { mask=mask|(1<< 0); uamo1[18].side=0; uamo2[18].side=0; uamo1[19].side=0; uamo2[19].side=0; }

    return mask;
}

int get_maskd2_20(Uamo uamo1, Uamo uamo2)
{
    int mask = 0;
    
    if ( uamo1[ 0].side != uamo2[ 0].side ) { mask=mask|(1<<19); }
    if ( uamo1[ 1].side != uamo2[ 1].side ) { mask=mask|(1<<18); }
    if ( uamo1[ 2].side != uamo2[ 2].side ) { mask=mask|(1<<17); }
    if ( uamo1[ 3].side != uamo2[ 3].side ) { mask=mask|(1<<16); }
    if ( uamo1[ 4].side != uamo2[ 4].side ) { mask=mask|(1<<15); }
    if ( uamo1[ 5].side != uamo2[ 5].side ) { mask=mask|(1<<14); }
    if ( uamo1[ 6].side != uamo2[ 6].side ) { mask=mask|(1<<13); }
    if ( uamo1[ 7].side != uamo2[ 7].side ) { mask=mask|(1<<12); }
    if ( uamo1[ 8].side != uamo2[ 8].side ) { mask=mask|(1<<11); }
    if ( uamo1[ 9].side != uamo2[ 9].side ) { mask=mask|(1<<10); }
    if ( uamo1[10].side != uamo2[10].side ) { mask=mask|(1<< 9); }
    if ( uamo1[11].side != uamo2[11].side ) { mask=mask|(1<< 8); }
    if ( uamo1[12].side != uamo2[12].side ) { mask=mask|(1<< 7); }
    if ( uamo1[13].side != uamo2[13].side ) { mask=mask|(1<< 6); }
    if ( uamo1[14].side != uamo2[14].side ) { mask=mask|(1<< 5); }
    if ( uamo1[15].side != uamo2[15].side ) { mask=mask|(1<< 4); }
    if ( uamo1[16].side != uamo2[16].side ) { mask=mask|(1<< 3); }
    if ( uamo1[17].side != uamo2[17].side ) { mask=mask|(1<< 2); }
    if ( uamo1[18].side != uamo2[18].side ) { mask=mask|(1<< 1); }

    return mask;

}


// }}}
// {{{ select_functions

int (*get_top_bit)(int);
int (*get_bot_bit)(int);
int (*merge_top_bot)(int, int);
int (*get_maskd1)(Uamo, Uamo);
int (*get_maskd2)(Uamo, Uamo);

void select_functions(int molen)
{
    if (molen == 10) {
        get_top_bit   = get_top_bit_10;
        get_bot_bit   = get_bot_bit_10;
        merge_top_bot = merge_top_bot_10;
        get_maskd1    = get_maskd1_10;
        get_maskd2    = get_maskd2_10;
    }
    else if (molen == 12)
    {
        get_top_bit   = get_top_bit_12;
        get_bot_bit   = get_bot_bit_12;
        merge_top_bot = merge_top_bot_12;
        get_maskd1    = get_maskd1_12;
        get_maskd2    = get_maskd2_12;
    }
    else if (molen == 14)
    {
        get_top_bit   = get_top_bit_14;
        get_bot_bit   = get_bot_bit_14;
        merge_top_bot = merge_top_bot_14;
        get_maskd1    = get_maskd1_14;
        get_maskd2    = get_maskd2_14;
    }
    else if (molen == 16)
    {
        get_top_bit   = get_top_bit_16;
        get_bot_bit   = get_bot_bit_16;
        merge_top_bot = merge_top_bot_16;
        get_maskd1    = get_maskd1_16;
        get_maskd2    = get_maskd2_16;
    }
    else if (molen == 18)
    {
        get_top_bit   = get_top_bit_18;
        get_bot_bit   = get_bot_bit_18;
        merge_top_bot = merge_top_bot_18;
        get_maskd1    = get_maskd1_18;
        get_maskd2    = get_maskd2_18;
    }
    else if (molen == 20)
    {
        get_top_bit   = get_top_bit_20;
        get_bot_bit   = get_bot_bit_20;
        merge_top_bot = merge_top_bot_20;
        get_maskd1    = get_maskd1_20;
        get_maskd2    = get_maskd2_20;
    }
    else
    {
        fprintf(stderr, "error: cannot support the motif length '%d'\n", molen);
        exit(1);
    }
}

// }}}



int check_side_mask(int maskm, int maskd2)
{

    int ret = 0;
    int offset;

    for (int i=0; i<molen-1; i++)
    {
        offset = (molen-1)-i;

        // cout << "i="<<i << "\t" << ((maskd2>>offset)&1) << "\t" << endl;

        if ( (maskd2 >> offset) & 1 )
        {

            // cout << "maskm(i):"<<((maskm>>offset)&1) << "\t" << "maskm(i+1):"<<((maskm>>(offset-1))&1) << endl;

            if ( (((maskm>>(offset))&1)==0) and (((maskm>>(offset-1))&1)==0) )
            {
                ret = 1;
                break;
            }
        }
    }

    return ret;
}



void (*count_contbl)(int, int**, int*);


void count_contbl1(int u, int** contbl_dm, int* contbl_dm_fc)
{
    int maskd1 = 0;   // difference main mask
    int maskd2 = 0;   // difference side mask
    int maskt  = 0;   // top    half mask
	int maskb  = 0;   // bottom half mask
	int maskm  = 0;   // mergred mask

    struct Fcm *pfcm;

    for (int b=0; b<10; b++)
    {

        (b == 0) ? memset(contbl_dm[b], 0, sizeof(int)*full_mask_size) : memcpy(contbl_dm[b], contbl_dm[b-1], sizeof(int)*full_mask_size);

        for (int i=0; i<(int)vec_fc[b].size(); i++)
        {
            pfcm = &vec_fc[b][i];
            memset(contbl_dm_fc, 0, sizeof(int)*full_mask_size);

            for (int j=0; j<(int)pfcm->vec_ia.size(); j++)
            {
                
                // cout << masking(u, 0) << "\t";
                // cout << "bin="<<b << "\t" << "fcm="<<i << "\t" << "iac="<<j << "\t" << "\t" << masking(pfcm->vec_ia[j].u, 0) << "\t";

                memcpy(tmp_uamo1, vec_uamo[u],                 sizeof(struct Ir)*molen);
                memcpy(tmp_uamo2, vec_uamo[pfcm->vec_ia[j].u], sizeof(struct Ir)*molen);

                maskd1 = (*get_maskd1)(tmp_uamo1, tmp_uamo2);
                maskd2 = (*get_maskd2)(tmp_uamo1, tmp_uamo2);

                // cout << "maskd1="<<maskd1 << "\t" << "maskd2="<<maskd2 << endl;

                maskt = (*get_top_bit)(maskd1);
                maskb = (*get_bot_bit)(maskd1);

                for (int m=0; m<(int)vec_halfmask[maskt].size(); m++) {
                for (int n=0; n<(int)vec_halfmask[maskb].size(); n++) {

                    maskm = (*merge_top_bot)(vec_halfmask[maskt][m], vec_halfmask[maskb][n]);
                    
                    // cout << "maskm="<<maskm << endl;

                    if (check_side_mask(maskm, maskd2) == 0)
                    {
                        contbl_dm_fc[maskm] = 1;
                    }

                }}

            }
            for (int m=0; m<full_mask_size; m++) { contbl_dm[b][m] += contbl_dm_fc[m]; }
        }
    }
}


void determine_iase(int* ipos, int mask, int* s, int* e)
{
    // determine interacting start,end region
    *s = -1; for (int i=0; i<molen; i++) { if ( (ipos[(molen-1)-i] != -1) and (((mask>>          i)&1)==0) ) { *s = ipos[(molen-1)-i]; break; } }
    *e = -1; for (int i=0; i<molen; i++) { if ( (ipos[          i] != -1) and (((mask>>(molen-1)-i)&1)==0) ) { *e = ipos[          i]; break; } }
}


void count_contbl2(int u, int** contbl_dm, int* contbl_dm_fc)
{
    int maskd1 = 0;   // difference main mask
    int maskd2 = 0;   // difference side mask
    int maskt = 0;   // top    half mask
	int maskb = 0;   // bottom half mask
	int maskm = 0;   // mergred mask
    int s, e;

    struct Fcm *pfcm;

    for (int b=0; b<10; b++)
    {
        (b == 0) ? memset(contbl_dm[b], 0, sizeof(int)*full_mask_size) : memcpy(contbl_dm[b], contbl_dm[b-1], sizeof(int)*full_mask_size);
        for (int i=0; i<(int)vec_fc[b].size(); i++)
        {
            pfcm = &vec_fc[b][i];
            if (pfcm->lam == NULL) continue;
            memset(contbl_dm_fc, 0, sizeof(int)*full_mask_size);

            for (int j=0; j<(int)pfcm->vec_ia.size(); j++)
            {
                memcpy(tmp_uamo1, vec_uamo[u],                 sizeof(struct Ir)*molen);
                memcpy(tmp_uamo2, vec_uamo[pfcm->vec_ia[j].u], sizeof(struct Ir)*molen);

                maskd1 = (*get_maskd1)(tmp_uamo1, tmp_uamo2);
                maskd2 = (*get_maskd2)(tmp_uamo1, tmp_uamo2);

                maskt = (*get_top_bit)(maskd1);
                maskb = (*get_bot_bit)(maskd1);

                for (int m=0; m<(int)vec_halfmask[maskt].size(); m++) {
                for (int n=0; n<(int)vec_halfmask[maskb].size(); n++) {
                    maskm = (*merge_top_bot)(vec_halfmask[maskt][m], vec_halfmask[maskb][n]);
                    if (check_side_mask(maskm, maskd2) == 0)
                    {
                        determine_iase(vec_fc[b][i].vec_ia[j].ipos, maskm, &s, &e);
                        if (s != -1 and e != -1 and s != e and vec_fc[b][i].lam[s][e]) contbl_dm_fc[maskm] = 1;
                    }
                }}
            }
            for (int m=0; m<full_mask_size; m++) { contbl_dm[b][m] += contbl_dm_fc[m]; }
        }
    }
}

// }}}
// {{{ assctest_cst

void assctest_cst(int uamo_start, int uamo_end, char* out_file_txt, int ps, int bf)
{

    #ifdef DEBUG
    gettimeofday(&tvs[0], NULL);
    #endif
    
    #ifdef DEBUG
    struct timeval tv1, tv2, tv3;
    timerclear(&tv1);
    timerclear(&tv2);
    timerclear(&tv3);
    gettimeofday(&tvs[1], NULL);
    #endif
    
    // [1] contbl_d, the number of down-regulated genes
    //     in 2x2 contingency table by cutoff. (=n1+n2)
    int* contbl_d = (int*)malloc(sizeof(int)*10);
    for (int b=0; b<10; b++)
    {
        contbl_d[b] = ( (b==0) ? 0 : contbl_d[b-1] ) + vec_fc[b].size();
        // printf("contbl_d[%d]=%d\n", b, contbl_d[b]);
	}

    // [2] contbl_dm, the number of down-regulated genes
    //     with the corresponding ambiguous motif. (=n1)
    int** contbl_dm = (int**)malloc(sizeof(int*)*10);
    for (int b=0; b<10; b++) { contbl_dm[b] = (int*)malloc(sizeof(int)*full_mask_size); }
    int* contbl_dm_fc = (int*)malloc(sizeof(int)*full_mask_size);


    // [3] write output
    
    bool is_passed;
    struct Chi2t chi2t;
    // struct Chi2t chi2tm;
    // struct Chi2t chi2tn;
    
    // char out_dir_txt[1024];
    // struct stat st = {0};
    // strcpy(out_dir_txt, out_file_txt);
    // dirname(out_dir_txt);
    // if (stat(out_dir_txt, &st) == -1) mkdir(out_dir_txt, 0775);
    
    #ifdef DEBUG
    gettimeofday(&tve[1], NULL);
    timersub(&tve[1], &tvs[1], &tv1);
    #endif
    
	ofstream fout;
	fout.open(out_file_txt, ios::out);  // ios::out=16, ios::binary=4
	for (int u=uamo_start; u<=uamo_end; u++)  // u: uamo index
	{
        
        // [4.1] contingency table
        
        #ifdef DEBUG
        gettimeofday(&tvs[1], NULL);
        #endif

        count_contbl(u, contbl_dm, contbl_dm_fc);
    
        #ifdef DEBUG
        gettimeofday(&tve[1], NULL);
        timersub(&tve[1], &tvs[1], &tve[1]);
        timeradd(&tv2, &tve[1], &tv2);
        #endif



        // [4.2] chi-square test

        #ifdef DEBUG
        gettimeofday(&tvs[1], NULL);
        #endif

		for (int m=0; m<full_mask_size; m++)  // m: mask index
		{

            //		D	!D
            //	 M	n1	n3
            //	!M	n2	n4
            
            is_passed = true;
            
            const char* separator = "";
            ostringstream oss;

            memset(&chi2t,  0, sizeof(struct Chi2t));
            // memset(&chi2tm, 0, sizeof(struct Chi2t)); chi2tm.schi2 = -1.0e+308;
            // memset(&chi2tn, 0, sizeof(struct Chi2t)); chi2tn.schi2 =  1.0e+308;

            if (not (contbl_dm[9][m] > ps)) continue;
            
            // bf = 111100011 = 483
            // bf = 110000011 = 387
            for (int b=0; b<9; b++)
            {
                if (((bf >> (8-b)) & 1) == 0) continue;

                chi2t.b  = b;
                chi2t.n1 = contbl_dm[b][m];
                chi2t.n2 = contbl_d[b]     - chi2t.n1;
                chi2t.n3 = contbl_dm[9][m] - chi2t.n1;
                chi2t.n4 = (contbl_d[9] - contbl_d[b]) - chi2t.n3;

                // cout << u << "\t" << m << "\t" << chi2t.b << "\t" << chi2t.n1 << "\t" << chi2t.n2 << "\t" << chi2t.n3 << "\t" << chi2t.n4 << endl;

                if (not (chi2t.n1>5 and chi2t.n2>5 and chi2t.n3>5 and chi2t.n4>5)) { is_passed = false; break; }

                chi2t.rrsk = (float(chi2t.n1)/float(chi2t.n1+chi2t.n2)) > (float(chi2t.n3)/float(chi2t.n3+chi2t.n4));
                // if (psonly == 1) { if (not (chi2t.rrsk == true)) { is_passed = false; break; } }
                
                chi2t.chi2  = chi2stat(chi2t.n1, chi2t.n2, chi2t.n3, chi2t.n4, false);
                // chi2t.pval  = pchisq(chi2t.chi2, 1, FALSE, FALSE);
                // chi2t.logp  = -log10l(chi2t.pval);
                // if (chi2t.logp > MAX_SCORE) chi2t.logp = MAX_SCORE;
                // chi2t.score = chi2t.rrsk ? chi2t.logp : -chi2t.logp;
                // if (chi2t.score >= mschi2t.score) memcpy(&mschi2t, &chi2t, sizeof(struct Chi2t));

                chi2t.schi2 = chi2t.rrsk ? chi2t.chi2 : -chi2t.chi2;
                // if (chi2t.schi2 > chi2tm.schi2) memcpy(&chi2tm, &chi2t, sizeof(struct Chi2t));
                // if (chi2t.schi2 < chi2tn.schi2) memcpy(&chi2tn, &chi2t, sizeof(struct Chi2t));
            
                oss << separator << chi2t.schi2;
                separator = ",";
            }

            if ( is_passed == true )
            {
                fout << u               << "\t";  // 0 uamo index
                fout << m               << "\t";  // 1 mask index
                fout << ndc[m]          << "\t";  // 2 the number of dont-care residues
                fout << masking(u,m)    << "\t";  // 3 amo string
                fout << contbl_dm[9][m] << "\t";  // 4 presence of amo (=n1+n3)
                fout << oss.str()       << endl;  // 5

                // fout << chi2tm.b        << "\t";  // 5 maximum
                // fout << chi2tm.schi2    << "\t";  // 6
                // fout << chi2tn.b        << "\t";  // 7 minimum
                // fout << chi2tn.schi2    << endl;  // 8

                // fout << mschi2t.chi2    << "\t";  // maximum chi-square statistic
                // fout << mschi2t.score   << "\t";  // maximum enrichment score
                // fout << mschi2t.b       << endl;  // cutoff index which gives maximum enrichment score
            }

		}

        #ifdef DEBUG
        gettimeofday(&tve[1], NULL);
        timersub(&tve[1], &tvs[1], &tve[1]);
        timeradd(&tv3, &tve[1], &tv3);
        #endif
        
	}

	fout.close();
    
    // release heap memory
    free(ndc);
    for (int i=0; i < half_mask_size; i++) vec_halfmask[i].clear();
    for (int i=0; i < (int)vec_uamo.size(); i++) free(vec_uamo[i]); vec_uamo.clear();
    free(contbl_d);
    for (int b=0; b<10; b++) free(contbl_dm[b]); free(contbl_dm);
    free(contbl_dm_fc);
    
    #ifdef DEBUG
    gettimeofday(&tve[0], NULL);
    timersub(&tve[0], &tvs[0], &tve[0]);
    printf("assctest_chi2      %4d.%06d sec\n", tve[0].tv_sec, tve[0].tv_usec);
    printf("- step1(init)      %4d.%06d sec\n", tv1.tv_sec, tv1.tv_usec);
    printf("- step2(contbl)    %4d.%06d sec\n", tv2.tv_sec, tv2.tv_usec);
    printf("- step3(chi2t)     %4d.%06d sec\n", tv3.tv_sec, tv3.tv_usec);
    #endif
}


// }}}
// {{{ cst1, cst2


int usage_cst(int mode)
{
    fprintf(stderr, "\n"); if (mode == 1) {
    fprintf(stderr, "usage: md cst1 <in_file_uamo> <uamo_start> <uamo_end> <in_file_fc> <in_dir_raw> <out_file_txt> <ps> <bf>\n"); } else {
    fprintf(stderr, "usage: md cst2 <in_file_uamo> <uamo_start> <uamo_end> <in_file_fc> <in_dir_lam> <in_dir_raw> <out_file_txt> <ps> <bf>\n"); }
    fprintf(stderr, "\n");
    fprintf(stderr, "inputs: in_file_uamo  text file containing non-redundant unambiguous motifs\n");
    fprintf(stderr, "        uamo_start    index of the first unambiguous motif to be tested (0-based)\n");
    fprintf(stderr, "        uamo_end      index of the last unambiguous motif to be tested (0-based)\n");
    fprintf(stderr, "        in_file_fc    sorted list of genes by fold change value\n"); if (mode == 2) {
    fprintf(stderr, "        in_dir_lam    directory which contains local AU matrix binary files\n"); }
    fprintf(stderr, "        in_dir_raw    directory which contains mRNA:miRNA interactions (unambiguous motif index, interacting 3'utr positions)\n");
    fprintf(stderr, "        out_file_txt  text file where output result will be writeen\n");
    fprintf(stderr, "        ps            minimum of motif presence\n");
    fprintf(stderr, "        bf            flag for bins to be tested\n");
    fprintf(stderr, "\n");
    return 1;
}


int cst1(int argc, char *argv[])
{
    if (argc < 9) return usage_cst(1);
    
	char* in_file_uamo =      argv[1];
	int   uamo_start   = atoi(argv[2]);
    int   uamo_end     = atoi(argv[3]);
    char* in_file_fc   =      argv[4];
    char* in_dir_raw   =      argv[5];
    char* out_file_txt =      argv[6];
    int   ps           = atoi(argv[7]);
    int   bf           = atoi(argv[8]);

    init_map_ir();
    load_uamolist(in_file_uamo);
    load_fc(in_file_fc, NULL, in_dir_raw);
    
    select_functions(molen);
    count_contbl = count_contbl1;
    assctest_cst(uamo_start, uamo_end, out_file_txt, ps, bf);

    return 0;
}


int cst2(int argc, char *argv[])
{
    if (argc < 10) return usage_cst(2);
    
	char* in_file_uamo =      argv[1];
	int   uamo_start   = atoi(argv[2]);
    int   uamo_end     = atoi(argv[3]);
    char* in_file_fc   =      argv[4];
    char* in_dir_lam   =      argv[5];
    char* in_dir_raw   =      argv[6];
    char* out_file_txt =      argv[7];
    int   ps           = atoi(argv[8]);
    int   bf           = atoi(argv[9]);

    init_map_ir();
    load_uamolist(in_file_uamo);
    load_fc(in_file_fc, in_dir_lam, in_dir_raw);

    select_functions(molen);
    count_contbl = count_contbl2;
    assctest_cst(uamo_start, uamo_end, out_file_txt, ps, bf);

    return 0;
}


// }}}



int x2p(int argc, char *argv[])
{
    if (argc < 3) return 1;
    
	char* in_file_txt  = argv[1];
    char* out_file_txt = argv[2];

    // ------ //

    ifstream fin;
    ofstream fout;

    string line;
    vector<string> col;

    fin.open(in_file_txt, ios::in);
	fout.open(out_file_txt, ios::out);

    double chi2m, pvalm;
    double chi2n, pvaln;

    if ( !fin.is_open() ) { fprintf(stderr, "error: cannot open file '%s'\n", in_file_txt); exit(1); }
    while ( getline(fin, line) )
    {
        col = split(line, '\t');

        chi2m = fabs(atof(col[6].c_str()));
        pvalm = pchisq(chi2m, 1, FALSE, FALSE);
        
        chi2n = fabs(atof(col[8].c_str()));
        pvaln = pchisq(chi2n, 1, FALSE, FALSE);

        fout << pvalm << "\t" << pvaln << endl;

    }

    fin.close();
    fout.close();



    return 0;
}


// {{{ assctest_fet

/*

double fexact_pval(int n1, int n2, int n3, int n4)
{
    int nrow = 2;
    int ncol = 2;
    double table[4] = {n1, n2, n3, n4};
    double expect = -1.0;
    double percnt = 100.0;
    double emin = 0.0;
    double prt = 0.0;
    double pre = 1.0;
    int workspace = 3000000;

    fexact(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre, &workspace);

    return pre;

    // ACATTCC  258 773 804 8480    2.003776e-47    miR-206_7mer-m8;miR-613_7mer-m8;miR-1_7mer-m8;
}


void assctest_fet(int uamo_start, int uamo_end, int minps, int ncont, char* out_file_txt)
{
    
    #ifdef DEBUG
    gettimeofday(&tvs[0], NULL);
    #endif
    
    #ifdef DEBUG
    struct timeval tv1, tv2, tv3;
    timerclear(&tv1);
    timerclear(&tv2);
    timerclear(&tv3);
    gettimeofday(&tvs[1], NULL);
    #endif
    
    // [1] contbl_d, the number of down-regulated genes
    //     in 2x2 contingency table by cutoff. (=n1+n2)
    int* contbl_d = (int*)malloc(sizeof(int)*10);
    for (int b=0; b<10; b++)
    {
        contbl_d[b] = ( (b==0) ? 0 : contbl_d[b-1] ) + vec_fc[b].size();
        // printf("contbl_d[%d]=%d\n", b, contbl_d[b]);
	}

    // [2] contbl_dm, the number of down-regulated genes
    //     with the corresponding ambiguous motif. (=n1)
    int** contbl_dm = (int**)malloc(sizeof(int*)*10);
    for (int b=0; b<10; b++) { contbl_dm[b] = (int*)malloc(sizeof(int)*full_mask_size); }
    int* contbl_dm_fc = (int*)malloc(sizeof(int)*full_mask_size);


    // [3] write output
    bool is_passed;
    char buf[65536];

    string out_fet;
    struct Fet fet;
    struct Fet msfet;
    
    // char out_dir_txt[1024];
    // struct stat st = {0};
    // strcpy(out_dir_txt, out_file_txt);
    // dirname(out_dir_txt);
    // if (stat(out_dir_txt, &st) == -1) mkdir(out_dir_txt, 0775);
    
    #ifdef DEBUG
    gettimeofday(&tve[1], NULL);
    timersub(&tve[1], &tvs[1], &tv1);
    #endif
    
	ofstream fout;
	fout.open(out_file_txt, ios::out);  // ios::out=16, ios::binary=4
	for (int u=uamo_start; u<=uamo_end; u++)  // u: uamo index
	{
        
        // [4.1] contingency table
        
        #ifdef DEBUG
        gettimeofday(&tvs[1], NULL);
        #endif

        count_contbl(u, contbl_dm, contbl_dm_fc);
    
        #ifdef DEBUG
        gettimeofday(&tve[1], NULL);
        timersub(&tve[1], &tvs[1], &tve[1]);
        timeradd(&tv2, &tve[1], &tv2);
        #endif



        // [4.2] fisher's exact test

        #ifdef DEBUG
        gettimeofday(&tvs[1], NULL);
        #endif

		for (int m=0; m<full_mask_size; m++)  // m: mask index
		{

            //		D	!D
            //	M	n1	n3
            //	!M	n2	n4

            is_passed  = true;
            msfet.logp = 0.0;
            out_fet.clear();

            if (not (contbl_dm[9][m] >= minps)) continue;

            for (int b=0; b<ncont; b++)
            {
                fet.b  = b;
                fet.n1 = contbl_dm[b][m];
                fet.n2 = contbl_d[b]     - fet.n1;
                fet.n3 = contbl_dm[9][m] - fet.n1;
                fet.n4 = (contbl_d[9] - contbl_d[b]) - fet.n3;
                if (not (fet.n1>5 and fet.n2>5 and fet.n3>5 and fet.n4>5)) { is_passed = false; break; }

                fet.rrsk = (float(fet.n1)/float(fet.n1+fet.n2)) > (float(fet.n3)/float(fet.n3+fet.n4));
                if (not (fet.rrsk == true)) { is_passed = false; break; }

                fet.pval = fexact_pval(fet.n1, fet.n2, fet.n3, fet.n4);
                fet.logp = -log10l(fet.pval);
                if (fet.logp > MAX_SCORE) fet.logp = MAX_SCORE;
                if (fet.logp >= msfet.logp) memcpy(&msfet, &fet, sizeof(struct Fet));

                sprintf(buf, "%d,%d,%d,%d,%d,%d,na,%f;", fet.b, fet.n1, fet.n2, fet.n3, fet.n4, fet.rrsk, fet.logp);
                out_fet.append(buf);
            }

            if ( is_passed == true )
            {
                fout << u               << "\t";  // uamo index
                fout << m               << "\t";  // mask index
                fout << ndc[m]          << "\t";  // the number of dont-care residues
                fout << masking(u,m)    << "\t";  // amo string
                fout << contbl_dm[9][m] << "\t";  // presence of amo (=n1+n3)
                fout << "na"            << "\t";  //
                fout << msfet.logp      << "\t";  // maximum -log10(p-value)
                fout << out_fet         << endl;  // chi-square test details
            }
		}

        #ifdef DEBUG
        gettimeofday(&tve[1], NULL);
        timersub(&tve[1], &tvs[1], &tve[1]);
        timeradd(&tv3, &tve[1], &tv3);
        #endif
        
	}

	fout.close();
    
    // release heap memory
    free(ndc);
    for (int i=0; i < half_mask_size; i++) vec_halfmask[i].clear();
    for (int i=0; i < (int)vec_uamo.size(); i++) free(vec_uamo[i]); vec_uamo.clear();
    free(contbl_d);
    for (int b=0; b<10; b++) free(contbl_dm[b]); free(contbl_dm);
    free(contbl_dm_fc);
    
    #ifdef DEBUG
    gettimeofday(&tve[0], NULL);
    timersub(&tve[0], &tvs[0], &tve[0]);
    printf("assctest_fet       %4d.%06d sec\n", tve[0].tv_sec, tve[0].tv_usec);
    printf("- step1(init)      %4d.%06d sec\n", tv1.tv_sec, tv1.tv_usec);
    printf("- step2(contbl)    %4d.%06d sec\n", tv2.tv_sec, tv2.tv_usec);
    printf("- step3(fet)       %4d.%06d sec\n", tv3.tv_sec, tv3.tv_usec);
    #endif
}

*/

// }}}
// {{{ fet1: Fisher's exact test

/*

int usage_fet1(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: md fet1 <in_file_uamo> <uamo_start> <uamo_end> <in_file_fc> <in_dir_raw> <minps> <ncont> <out_file_txt>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "inputs: in_file_uamo  text file containing non-redundant unambiguous motifs\n");
    fprintf(stderr, "        uamo_start    index of the first unambiguous motif to be tested (0-based)\n");
    fprintf(stderr, "        uamo_end      index of the last unambiguous motif to be tested (0-based)\n");
    fprintf(stderr, "        in_file_fc    sorted list of genes by fold change value\n");
    fprintf(stderr, "        in_dir_raw    directory which contains mRNA:miRNA interactions (unambiguous motif index, interacting 3'utr positions)\n");
    fprintf(stderr, "        minps         minimum of motif presence\n");
    fprintf(stderr, "        ncont         how many contingency tables are used? (2:10%,20%, 5:10%,20%,30%,40%,50%)\n");
    fprintf(stderr, "        out_file_txt  text file where output result will be writeen\n");
    fprintf(stderr, "\n");
    return 1;
}


int fet1(int argc, char *argv[])
{
    if (argc < 9) return usage_fet1();
    
	char* in_file_uamo =      argv[1];
	int   uamo_start   = atoi(argv[2]);
    int   uamo_end     = atoi(argv[3]);
    char* in_file_fc   =      argv[4];
    char* in_dir_raw   =      argv[5];
    int   minps        = atoi(argv[6]);
    int   ncont        = atoi(argv[7]);
    char* out_file_txt =      argv[8];
    

    init_map_ir();
    load_uamolist(in_file_uamo);
    load_fc(in_file_fc, NULL, in_dir_raw);
    
    select_functions(molen);
    count_contbl = count_contbl1;
    assctest_fet(uamo_start, uamo_end, minps, ncont, out_file_txt);

    return 0;
}

*/

// }}}

// {{{ rst: rank sum test

/*

int rst(int argc, char *argv[])
{
    return 0;
}

*/

// }}}




int usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "program: motif discovery (for targeting project)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "contact: Sukjun Kim (dandyrilla@naver.com)\n");
    fprintf(stderr, "         The Baek Research Group of Computational Biology, SNU\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "version: 3.6\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: md <command> [inputs]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "command: cst1  chi-square test\n");
    fprintf(stderr, "         cst2  chi-square test with local AU considered\n");
    fprintf(stderr, "         x2p   convert chi-squre into p-value\n");
    //fprintf(stderr, "         fet1  Fisher's exact test\n");
    //fprintf(stderr, "         rst   Wilcoxon's rank sum test\n");
    fprintf(stderr, "\n");
    return 1;
}



int main(int argc, char *argv[])
{
    if (argc < 2) return usage();
    if      (strcmp(argv[1], "cst1") == 0) return cst1(argc-1, argv+1);
    else if (strcmp(argv[1], "cst2") == 0) return cst2(argc-1, argv+1);
    else if (strcmp(argv[1], "x2p") == 0)  return x2p(argc-1, argv+1);
    //else if (strcmp(argv[1], "fet1") == 0) return fet1(argc-1, argv+1);
    //else if (strcmp(argv[1], "rst" ) == 0) return rst(argc-1, argv+1);
    else { fprintf(stderr, "unrecognized command '%s'\n", argv[1]); return 1; }
	return 0;
}



