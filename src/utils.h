#ifndef UTILS
#define UTILS

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <string>
#include "typedefs.h"

using namespace std;

int vecIntSum(const vector <int>& x);
double mult_llh_no_constant(const vector <int>& n, const vector <double>& pi);
int numGenotypes (const int nA, const int ploidy);
void printVersion();
double probSubErr(const char& phred, const double& eps_S);
double ph33toProb(const char& phred);
bool compareRefPos(baseInfo snp1, baseInfo snp2);
void readPosFile(const string& posInput, unordered_map <string, locusInfo>& posMap, vector <string>& locusNames);
void readPosFile2(const string& posInput,
				unordered_map <string, locusInfo>& posMap,  // edited by this function instead of returned
				vector <string>& locusNames);
void readRefSeqs(const string& refInput, const vector <string>& locusNames, unordered_map <string, locusInfo>& posMap);
void buildAlleles(const unordered_map <string, locusInfo>& posMap, const vector <string>& locusNames,
			unordered_map <string, alleleTable>& locusAlleleTables);
void splitString(const string& inString, const char delim, vector <string>& split);
double logSumExp(const vector <double>& x);

#endif // UTILS
