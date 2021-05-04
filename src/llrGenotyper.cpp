/*
 simple genotyper for output of microTyper
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
				const unordered_map <string, unordered_map<string, double>>& priorMap,
				const bool bool_counts){
	genotype g;
	allele a;
	a.name = split[2];
	g.A1 = a;
	a.name = split[3];
	g.A2 = a;
	if(g.A1.name == g.A2.name){
		g.hom = true;
	} else {
		g.hom = false;
	}
	if(bool_counts){
		g.llh = 0;
	} else {
		g.llh = stod(split[4]);
	}
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
void callGenotype(genoTable& gT, ofstream& genoOutFile, const double& c, const double& minPerfect,
				const bool bool_counts, const double eps){
	// if using counts, calculate llh for each genotype here
	if(bool_counts){
		vector <string> allele_names;
		vector <int> allele_counts;
		// get counts of each allele
		for(int i = 0, max = gT.gTable.size(); i < max; i++){
			if(gT.gTable[i].hom){
				allele_names.push_back(gT.gTable[i].A1.name);
				allele_counts.push_back(gT.gTable[i].A1_perfect);
			}
		}
		// template vector of probability for multinomial with error rates, so just have to update one or two values
		vector <double> template_pi (allele_counts.size(), (eps / (allele_counts.size() - 1)));
		// calculate llh for each genotype
		for(int i = 0, max = gT.gTable.size(); i < max; i++){
			vector <double> pi (template_pi);
			if(gT.gTable[i].hom){
				// update one value of pi
				for(int j = 0, m = pi.size(); j < m; j++){
					if(allele_names[j] == gT.gTable[i].A1.name){
						pi[j] = 1 - eps;
						break;
					}
				}
			} else {
				// update two values of pi
				int place = 0;
				for(int j = 0, m = pi.size(); j < m; j++){
					if(allele_names[j] == gT.gTable[i].A1.name){
						pi[j] = (.5 * (1 - eps)) + (.5 * eps / (allele_names.size() - 1));
						// add perfect allele counts to output
						gT.gTable[i].A1_perfect = allele_counts[j];
						place++;
					} else if(allele_names[j] == gT.gTable[i].A2.name){
						pi[j] = (.5 * (1 - eps)) + (.5 * eps / (allele_names.size() - 1));
						// add perfect allele counts to output
						gT.gTable[i].A2_perfect = allele_counts[j];
						place++;
					}
					if(place == 2) break;
				}
			}
			// calc llh
			gT.gTable[i].llh += mult_llh_no_constant(allele_counts, pi);
		}
	}

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

	cerr << "Just a friendly notice that this function is deprecated. Only use it if you specifically prefer it to the newer option of mtype2." << endl;

	string mhgenosInput; // -f
	string outputName ("mh_genotypes.txt"); // -o
	string priorFile; // -p
	double c = log(.95); // -c probability threshold to accept a genotype: currently assumes a constant prior on genotypes
	double minPerfect = 0; // -m minimum number of perfect reads for the alleles in that genotype
	double eps = .01;
	bool bool_counts = false;

	ios_base::sync_with_stdio(false);
	cin.tie(NULL);

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
		} else if (x == "--version") {
			 printVersion();
		} else if (x == "--count") {
			 bool_counts = true;
		} else if (x == "-eps") {
			 eps = stod(argv[++i]);
		} else {
			cerr << "Error: Option " << x << " not recognized." << endl;
			exit(EXIT_FAILURE);
		}
	}

	if(!bool_counts){
		cerr << "Just a friendly warning: you are not using the --count option. " <<
		"The --count option is typically recommended. Calling genotypes without" <<
		" --count may lead to a higher rate of calling heterozygous genotypes in " <<
		"the presence of low levels of index hopping, contamination etc. Basically, " <<
		"if the only source of error is sequencing error and is more or less independent by position in a read " <<
		"you are fine. Otherwise, use the --count option." << endl;
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

	bool use_stdin = false;
	if(mhgenosInput == "-") use_stdin = true;
	istream * llhIn;
	ifstream llhFile;
	if(!use_stdin){
		llhFile.open(mhgenosInput);
		if(!llhFile.is_open()){
			cerr << "Error: Log-likelihood file specified (-f) could not be opened." << endl;
			exit(EXIT_FAILURE);
		}
		llhIn = &llhFile;
	} else {
		llhIn = &cin;
	}


	string line;
	ofstream genoOutFile (outputName, ofstream::trunc);
	genoOutFile << "Indiv\tLocus\tAllele1\tAllele2\tPr_geno\tA1_perfect\tA2_perfect\n";
	string curInd ("");
	string curLoc ("");
	genoTable gT;
	getline((*llhIn), line); // skip header
	// have to "bookend" and treat first and last outside of loop
	vector <string> split;
	// first entry
	getline((*llhIn), line);

	splitString(line, '\t', split);
	gT.indName = split[0];
	gT.locName = split[1];

	addGenotype(split, gT, priorUse, priorMap, bool_counts);
	split.clear();

	while (getline((*llhIn), line)){
		splitString(line, '\t', split);
		// if same, add to list
		if(split[0] == gT.indName && split[1] == gT.locName){
			addGenotype(split, gT, priorUse, priorMap, bool_counts);
		} else {
			// call genotype
			callGenotype(gT, genoOutFile, c, minPerfect, bool_counts, eps);

			// start new genoTable
			gT.gTable.clear();
			gT.indName = split[0];
			gT.locName = split[1];
			addGenotype(split, gT, priorUse, priorMap, bool_counts);
		}
		split.clear();
	}
	// call last genotype
	callGenotype(gT, genoOutFile, c, minPerfect, bool_counts, eps);

	if (llhFile.is_open()) llhFile.close();
	genoOutFile.close();

}
