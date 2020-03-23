#ifndef CALCULATE_LLH
#define CALCULATE_LLH

#include <string>
#include <regex>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"

using namespace std;

void calculate_llh(const string bamInput,
					const vector <string>& locusNames,
					const unordered_map <string, locusInfo>& posMap,
					unordered_map <string, genoTable>& locusGenoTables,
					const double eps_I,
					const double eps_S);
void calculate_llh_logSumExp(const string bamInput,
					const vector <string>& locusNames,
					const unordered_map <string, locusInfo>& posMap,
					unordered_map <string, genoTable>& locusGenoTables,
					const double eps_I,
					const double eps_S);

#endif // CALCULATE_LLH
