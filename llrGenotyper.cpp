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
#include "utils.h"

using namespace std;

void addGenotype(const vector <string>& split, genoTable& gT){
	genotype g;
	allele a;
	a.name = split[2];
	g.A1 = a;
	a.name = split[3];
	g.A2 = a;
	g.llh = stod(split[4]);
	g.A1_perfect = stoi(split[5]);
	g.A2_perfect = stoi(split[6]);
	gT.gTable.push_back(g);
}

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
		} else {
			cerr << "Error: Option " << x << " not recognized." << endl;
			return 1;
		}
	}

	// input error checking
	if (mhgenosInput.length() == 0) {
		cerr << "Error: No mhgenos file input." << endl;
		return 1;
	}
	if (outputName.length() == 0) {
		cerr << "Error: No reference file input." << endl;
		return 1;
	}

	ifstream llhFile (mhgenosInput);
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

	addGenotype(split, gT);
	split.clear();

	while (getline(llhFile, line)){
		splitString(line, '\t', split);
		// if same, add to list
		if(split[0] == gT.indName && split[1] == gT.locName){
			addGenotype(split, gT);
		} else {
			// call genotype
			callGenotype(gT, genoOutFile, c, minPerfect);

			// start new genoTable
			gT.gTable.clear();
			gT.indName = split[0];
			gT.locName = split[1];
			addGenotype(split, gT);
		}
		split.clear();
	}
	// call last genotype
	callGenotype(gT, genoOutFile, c, minPerfect);

	llhFile.close();
	genoOutFile.close();

}
