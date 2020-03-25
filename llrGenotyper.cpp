/*
 simple genotyper for output of microTyper
 g++ llrGenotyper.cpp -c -I /usr/local/include/bamtools -L /usr/local/lib -lbamtools -lz -o llrGenotyper.o
 g++ llrGenotyper.o utils.o -o genoCaller
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <unordered_map>
#include "utils.h"

using namespace std;

void addGenotype(const vector <string>& split, genoTable& gT, const bool priorUse,
				const unordered_map <string, unordered_map<string, double>>& priorMap){
	genotype g;
	allele a;
	a.name = split[2];
	g.A1 = a;
	a.name = split[3];
	g.A2 = a;
	g.llh = stod(split[4]);
	g.A1_perfect = stoi(split[5]);
	g.A2_perfect = stoi(split[6]);

	// incorporating prior by just multiplying llh by prior,
	//   and then callGenotype will calculate appropriately
	if(priorUse){
		try {
			g.llh += priorMap.at(gT.locName).at(g.A1.name + "/" + g.A2.name);
		}
		catch (const out_of_range& err){
			cerr << "Error: could not find prior value for locus " <<
				gT.locName << " genotype " << g.A1.name + "/" + g.A2.name << endl;
			exit(EXIT_FAILURE);
		}
	}
	// add to genotype table
	gT.gTable.push_back(g);
}

// call genotypes with a flat prior
void callGenotype(const genoTable& gT, ofstream& genoOutFile, const double& c, const double& minPerfect){
	int maxl = 0;
	vector <double> rSum (gT.gTable.size());
	rSum[0] = gT.gTable[0].llh;
	for(int i = 1, max = gT.gTable.size(); i < max; i++){
		if (gT.gTable[i].llh > gT.gTable[maxl].llh) maxl = i;
		rSum[i] = gT.gTable[i].llh;
	}
	double p_max = gT.gTable[maxl].llh - logSumExp(rSum);
	// filtering by probability and number of perfect reads for the two alleles in the genotype
	if(p_max >= c && (gT.gTable[maxl].A1_perfect + gT.gTable[maxl].A2_perfect >= minPerfect)){
		genoOutFile << gT.indName << "\t" << gT.locName << "\t" << gT.gTable[maxl].A1.name << "\t" <<
			gT.gTable[maxl].A2.name << "\t" << exp(p_max) << "\t" << gT.gTable[maxl].A1_perfect <<
			"\t" << gT.gTable[maxl].A2_perfect << "\n";
	} else {
		// leaving genotype field blank for missing data
		genoOutFile << gT.indName << "\t" << gT.locName << "\t" << "\t" << "\t" <<
		exp(p_max) << "\t" << gT.gTable[maxl].A1_perfect << "\t" << gT.gTable[maxl].A2_perfect << "\n";
	}
}

int main(int argc, char* argv[]){

	string mhgenosInput; // -f
	string outputName ("mh_genotypes.txt"); // -o
	string priorFile; // -p
	double c = log(.95); // -c probability threshold to accept a genotype: currently assumes a constant prior on genotypes
	double minPerfect = 0; // -m minimum number of perfect reads for the alleles in that genotype

	// get input options and check
	std::string x;
	for (int i = 1; i < argc; i++) {
		x = string(argv[i]);
		if (x == "-f") {
			mhgenosInput = string(argv[++i]);
		} else if (x == "-o") {
			outputName = string(argv[++i]);
		} else if (x == "-c") {
			c = log(atof(argv[++i]));
		} else if (x == "-m") {
			minPerfect = atof(argv[++i]);
		} else if (x == "-p") {
			priorFile = string(argv[++i]);
		} else {
			cerr << "Error: Option " << x << " not recognized." << endl;
			exit(EXIT_FAILURE);
		}
	}

	// input error checking
	if (mhgenosInput.length() == 0) {
		cerr << "Error: No mhgenos file input." << endl;
		exit(EXIT_FAILURE);
	}
	if (outputName.length() == 0) {
		cerr << "Error: No reference file input." << endl;
		exit(EXIT_FAILURE);
	}

	// determine if a prior was input
	bool priorUse = false;
	if(priorFile.length() > 0) priorUse = true;

	// load prior if applicable
	// prior input file should have a header line and
	//	be tab delimited with columns of:
	//  locus \t allele1 \t allele2 \t priorProbability
	unordered_map <string, unordered_map <string, double> > priorMap;
	if(priorUse){
		vector <string> locusNames;
		ifstream pFile (priorFile);
		if(!pFile.is_open()){
			cerr << "Error: Prior file specified (-p) could not be opened." << endl;
			exit(EXIT_FAILURE);
		}
		string line;
		vector <string> split;
		getline(pFile, line); // skip header
		while (getline(pFile, line)){
			splitString(line, '\t', split);
			// if wrong length, skip
			if(split.size() != 4) {
				split.clear();
				continue;
			}
			if(split[3].length() == 0 || stod(split[3]) < 0) {
				cerr << "Error: missing or negative value in prior file." << endl;
				exit(EXIT_FAILURE);
			}
			if(stod(split[3]) == 0) {
				cerr << "Error: priors cannot be zero(but they can be very small)." << endl;
				exit(EXIT_FAILURE);
			}
			if(priorMap.count(split[0]) == 1){ // if in map already
					priorMap[split[0]][split[1] + "/" + split[2]] = stod(split[3]);
			} else {
				// add to map
				unordered_map <string, double> tempMap;
				tempMap[split[1] + "/" + split[2]] = stod(split[3]);
				priorMap[split[0]] = tempMap;
				// add to locus names vector
				locusNames.push_back(split[0]);
			}
			split.clear();
		}
		pFile.close();
		// normalize priors and convert to log
		for(int i = 0, m = locusNames.size(); i < m; i++){
			unordered_map <string, double> * tM;
			tM = &priorMap[locusNames[i]];
			double rsum = 0;
			for(auto it = (*tM).begin(); it != (*tM).end(); it++) rsum += it->second;
			for(auto it = (*tM).begin(); it != (*tM).end(); it++) it->second = log(it->second / rsum);
		}
	}

	ifstream llhFile (mhgenosInput);
	if(!llhFile.is_open()){
		cerr << "Error: Log-likelihood file specified (-f) could not be opened." << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	ofstream genoOutFile (outputName, ofstream::trunc);
	genoOutFile << "Indiv\tLocus\tAllele1\tAllele2\tPr_geno\tA1_perfect\tA2_perfect\n";
	string curInd ("");
	string curLoc ("");
	genoTable gT;
	getline(llhFile, line); // skip header
	// have to "bookend" and treat first and last outside of loop
	vector <string> split;
	getline(llhFile, line); // first entry
	splitString(line, '\t', split);
	gT.indName = split[0];
	gT.locName = split[1];

	addGenotype(split, gT, priorUse, priorMap);
	split.clear();

	while (getline(llhFile, line)){
		splitString(line, '\t', split);
		// if same, add to list
		if(split[0] == gT.indName && split[1] == gT.locName){
			addGenotype(split, gT, priorUse, priorMap);
		} else {
			// call genotype
			callGenotype(gT, genoOutFile, c, minPerfect);

			// start new genoTable
			gT.gTable.clear();
			gT.indName = split[0];
			gT.locName = split[1];
			addGenotype(split, gT, priorUse, priorMap);
		}
		split.clear();
	}
	// call last genotype
	callGenotype(gT, genoOutFile, c, minPerfect);

	llhFile.close();
	genoOutFile.close();

}
