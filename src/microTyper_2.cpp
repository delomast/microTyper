/*
 microTyper 2
 using read counts for valid alleles
 only considering allele within an individual
 TA Delomas Spring 2021
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"
#include "countAlleleReads.h"

using namespace std;

int main(int argc, char* argv[]){

	vector <string> bamInput; // -f
	string outputName ("microtyper_llh.mhgenos"); // -o
	string refInput; // -r fasta format
	string posInput; // -p
	double eps = .01; // this is the error rate for sequencing an incorrect allele -eps
	int batchSize = 100; // number of individuals to process and write at once -b
	unsigned int numThreads = 1; // number of threads to use -t
	int ploidy = 2; // -ploidy
	bool justAlleleCounts = false; // --justCount to just get allele counts
	double c = log(.99); // min post prob to call genotype -c
	double depth = 10; // min depth (total depth = sum of all alleles) to call genotype -d
	bool reportAll = false; // true to report all genotypes and posterior probabilities (ignores c and depth). --all

	ios_base::sync_with_stdio(false);

	// get input options and check
	string x;
	for (int i = 1; i < argc; i++) {
		x = string(argv[i]);
		if (x == "-f") {
			i++;
			for(; i < argc; i++){
				if(argv[i][0] == '-') {
					i--;
					break;
				}
				bamInput.push_back(string(argv[i]));
			}
		} else if (x == "-o") {
			outputName = string(argv[++i]);
		} else if (x == "-r") {
			 refInput = string(argv[++i]);
		} else if (x == "-p") {
			 posInput = string(argv[++i]);
		} else if (x == "-eps") {
			 eps = atof(argv[++i]);
		} else if (x == "-b") {
			 batchSize = stoi(argv[++i]);
		} else if (x == "-t") {
			 numThreads = stoi(argv[++i]);
		} else if (x == "-ploidy") {
			 ploidy = stoi(argv[++i]);
		} else if (x == "-c") {
			c = log(atof(argv[++i]));
		} else if (x == "-d") {
			depth = atof(argv[++i]);
		} else if (x == "--version") {
			 printVersion();
		} else if (x == "--justCount") {
			 justAlleleCounts = true;
		} else if (x == "--all") {
			 reportAll = true;
		} else {
			cerr << "Error: Option " << x << " not recognized." << endl;
			return 1;
		}
	}

	// optionally output results to std_out
	bool use_stdout = false;
	if (outputName == "-") use_stdout = true;

	// input error checking
	if (bamInput.size() == 0) {
		cerr << "Error: No bam file input." << endl;
		return 1;
	}
	if (refInput.length() == 0) {
		cerr << "Error: No reference file input." << endl;
		return 1;
	}
	if (posInput.length() == 0) {
		cerr << "Error: No position file input." << endl;
		return 1;
	}
	if(batchSize < 1){
		cerr << "Error: batch size must be 1 or greater." << endl;
		return 1;
	}
	if(eps <= 0 || eps >= 1){
		cerr << "Error: eps must be between 0 and 1." << endl;
		return 1;
	}
	if(ploidy <= 0){
		cerr << "Error: ploidy must be greater than 0." << endl;
		return 1;
	}

	// read in position input file
	unordered_map <string, locusInfo> posMap;
	vector <string> locusNames; // list of locus names - to make iterating through loci simple
	readPosFile2(posInput, posMap, locusNames);

	// read in reference sequences
	// adding reference value to struct baseInfo for each locus: for S, the ref base; D, the value 'I'; I, the value 'D'
	// and checking that all loci in position file have reference sequences
	readRefSeqs(refInput, locusNames, posMap);


	// Now, add each locus to the allele count table
	indAlleleCounts oneIndAlleleTable;
	for(int i = 0, max1 = locusNames.size(); i < max1; i++){
		unordered_map <string, int> tempMap;
		oneIndAlleleTable.locusMap[locusNames[i]] = tempMap;
	}

	// opening output file and writing header line
	ofstream llhOut;
	if(use_stdout) {
		if(justAlleleCounts){
			cout << "Indiv" << "\t" << "Locus" << "\t" << "Allele" <<
				"\t" << "Count" << endl;
		} else {
			cout << "Indiv" << "\t" << "Locus";
			for(int i = 1; i <= ploidy; i++) cout << "\t" << "Allele" << i;
			for(int i = 1; i <= ploidy; i++) cout << "\t" << "Allele" << i <<"_count";
			cout << "\t" << "p" << endl;
		}

	} else {
		llhOut.open(outputName, ofstream::trunc);
		if(!llhOut.is_open()){
			cerr << "Could not open " << outputName << "to write output." << endl;
			exit(EXIT_FAILURE);
		}
		if(justAlleleCounts){
			llhOut << "Indiv" << "\t" << "Locus" << "\t" << "Allele" <<
				"\t" << "Count" << endl;
		} else {
			llhOut << "Indiv" << "\t" << "Locus";
			for(int i = 1; i <= ploidy; i++) llhOut << "\t" << "Allele" << i;
			for(int i = 1; i <= ploidy; i++) llhOut << "\t" << "Allele" << i <<"_count";
			llhOut << "\t" << "p" << endl;
		}
	}
	// for each individual
	int indivPos = 0;
	int nIndiv = bamInput.size();
	while(indivPos < nIndiv){
		int curBatchSize = min(nIndiv - indivPos, batchSize);
		vector <indAlleleCounts> batchGenos (curBatchSize);
		for(int i = 0; i < curBatchSize; i++) batchGenos[i] = oneIndAlleleTable;
		// start parallel
		#pragma omp parallel for shared(batchGenos, posMap) num_threads(numThreads)
		for(int b = 0; b < curBatchSize; b++){
			// count reads
			countAlleleReads(bamInput[indivPos + b], posMap, batchGenos[b]);
		}
		// end parallel


		// if allele counts only requested, output those and move to the next batch
		// "Indiv" << "\t" << "Locus" << "\t" << "Allele" << "\t" << "Count" << endl;
		if(justAlleleCounts){
			for(int b = 0; b < curBatchSize; b++){ // for each ind
				for(int i = 0, max1 = locusNames.size(); i < max1; i++){ // for each locus
					unordered_map <string, int> * tempMap;
					tempMap = &batchGenos[b].locusMap[locusNames[i]];
					if((*tempMap).empty()){
						if(use_stdout){
							cout << batchGenos[b].indName << "\t" << locusNames[i] << "\t" << "noReads" << "\t" << "noReads" << endl;
						} else {
							llhOut << batchGenos[b].indName << "\t" << locusNames[i] << "\t" << "noReads" << "\t" << "noReads" << endl;
						}
					} else { // loop over all alleles found in that individual
						for(auto it = (*tempMap).cbegin(); it != (*tempMap).cend(); ++it){
							if(use_stdout){
								cout << batchGenos[b].indName << "\t" << locusNames[i] << "\t" << it->first << "\t" << it->second << endl;
							} else {
								llhOut << batchGenos[b].indName << "\t" << locusNames[i] << "\t" << it->first << "\t" << it->second << endl;
							}
						}
					}
				}
			}
		} else { // otherwise, another parallel for loop to calculate llh values
			// calculate llh
			vector <vector <string>> batchGenos_toWrite (curBatchSize);
			#pragma omp parallel for shared(batchGenos_toWrite, locusNames) num_threads(numThreads)
			for(int b = 0; b < curBatchSize; b++){
				// calculate posterior and return lines to write
				calcPosterior(batchGenos[b], ploidy, locusNames, eps, c, depth, reportAll, batchGenos_toWrite[b]);
			}
			// end parallel

			//write genotypes
			for(int b = 0; b < curBatchSize; b++){
				for(int l = 0, m = batchGenos_toWrite[b].size(); l < m; l++){
					if(use_stdout){
						cout << batchGenos_toWrite[b][l] << endl;
					} else {
						llhOut << batchGenos_toWrite[b][l] << endl;
					}
				}
			}
		}
		indivPos += curBatchSize;
	}

	// closing output file
	if(llhOut.is_open()) llhOut.close();

	return 0;
}
