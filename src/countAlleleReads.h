#ifndef COUNTALLELEREADS
#define COUNTALLELEREADS

#include <string>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"

using namespace std;

void countAlleleReads(const string bamInput,
					const unordered_map <string, locusInfo>& posMap,
					indAlleleCounts& indAlleleTable);
void buildA(const int nA,
						const int ploidy,
						vector <vector <int>>& allGenos);
void calcPosterior(const indAlleleCounts& batchGenosInd,
				const int ploidy,
				const vector <string>& locusNames,
				const double eps,
				const double c,
				const double depth,
				const bool reportAll,
				vector <string>& batchGenos_toWrite);

#endif // COUNTALLELEREADS
