/*
 microTyper
 TA Delomas Spring 2020
*/
// g++ microhap_dev.cpp -I /usr/local/include/bamtools -L /usr/local/lib -lbamtools -lz
// g++ microTyper_codeBlocks.cpp -I /usr/local/include/bamtools -L /usr/local/lib -lbamtools -lz -fopenmp -o mtype

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <regex>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"
#include "calculate_llh.h"

using namespace std;

int main(int argc, char* argv[]){

	vector <string> bamInput; // -f
	string outputName ("microtyper_llh.mhgenos"); // -o
	string refInput; // -r fasta format
	string posInput; // -p
	double eps_S = .01; // this is the additional error rate for sequencing substitutions (in addition to phred) -eps_S
	double eps_I = .01; // this is the error rate for sequencing indels -eps_I
	int batchSize = 100; // number of individuals to process and write at once -b
	unsigned int numThreads = 1; // number of threads to use -t
	bool call_with_lse = false; // to use log-sum-exp --manySNPs

	// get input options and check
	std::string x;
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
		} else if (x == "-eps_S") {
			 eps_S = atof(argv[++i]);
		} else if (x == "-eps_I") {
			 eps_I = atof(argv[++i]);
		} else if (x == "-b") {
			 batchSize = stoi(argv[++i]);
		} else if (x == "-t") {
			 numThreads = stoi(argv[++i]);
		} else if (x == "--manySNPs") {
			 call_with_lse = true;
		} else {
			cerr << "Error: Option " << x << " not recognized." << endl;
			return 1;
		}
	}

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
	if(eps_S < 0 || eps_S > 1){
		cerr << "Error: eps_S must be between 0 and 1 inclusive." << endl;
		return 1;
	}
	if(eps_I < 0 || eps_I > 1){
		cerr << "Error: eps_I must be between 0 and 1 inclusive." << endl;
		return 1;
	}

	// read in position input file
	unordered_map <string, locusInfo> posMap;
	vector <string> locusNames; // list of locus names - to make iterating through loci simple
	readPosFile(posInput, posMap, locusNames);



	// read in reference sequences
	// adding reference value to struct baseInfo for each locus: for S, the ref base; D, the value 'I'; I, the value 'D'
	// and checking that all loci in position file have reference sequences
	readRefSeqs(refInput, locusNames, posMap);


	// Now, for each locus, build reference table of all possible alleles and of all possible genotypes
	unordered_map <string, alleleTable> locusAlleleTables;
	buildAlleles(posMap, locusNames, locusAlleleTables);


	// now make a (diploid) genotype table for each locus
	unordered_map <string, genoTable> locusGenoTables;
	for(int i = 0, max = locusNames.size(); i < max; i++){
		genoTable tGeno;
		tGeno.locName = locusNames[i];
		alleleTable * tAT;
		tAT = &locusAlleleTables[locusNames[i]];
		for(int j = 0, max2 = (*tAT).alleleList.size(); j < max2; j++){
			for(int k = j; k < max2; k++){
				genotype tG;
				tG.A1 = (*tAT).alleleList[j];
				tG.A2 = (*tAT).alleleList[k];
				if (k == j){
					tG.hom = true;
				} else {
					tG.hom = false;
				}
				tG.llh = 0;
				tG.A1_perfect = 0;
				tG.A2_perfect = 0;
				tGeno.gTable.push_back(tG);
			}
		}
		locusGenoTables[locusNames[i]] = tGeno;
	}

	// opening output file and writing header line
	ofstream llhOut (outputName, ofstream::trunc);
	llhOut << "Indiv" << "\t" << "Locus" << "\t" << "Allele1" <<
				"\t" << "Allele2" << "\t" << "LLH" << "\t" <<
				"A1_perfect_count" << "\t" << "A2_perfect_count" << "\n";

	// for each individual
	int indivPos = 0;
	int nIndiv = bamInput.size();
	while(indivPos < nIndiv){
		int curBatchSize = min(nIndiv - indivPos, batchSize);
		vector <unordered_map <string, genoTable> > batchGenos (curBatchSize);
		for(int i = 0; i < curBatchSize; i++)batchGenos[i] = locusGenoTables;
		// start parallel
		#pragma omp parallel for shared(batchGenos, posMap) num_threads(numThreads)
		for(int b = 0; b < curBatchSize; b++){
			// calculate llh across loci and save to batchGenos
			if(call_with_lse){
				calculate_llh_logSumExp(bamInput[indivPos + b], locusNames, posMap, batchGenos[b],
					eps_I, eps_S);
			} else {
				calculate_llh(bamInput[indivPos + b], locusNames, posMap, batchGenos[b],
					eps_I, eps_S);
			}
		}
		// end parallel

		// write out genotypes
		for(int b = 0; b < curBatchSize; b++){
			genoTable * tGT;
			for(int i = 0, max = locusNames.size(); i < max; i++){ 	// for each locus
				tGT = &batchGenos[b][locusNames[i]];
				// output likelihoods and genotypes
				for(int j = 0, max2 = (*tGT).gTable.size(); j < max2; j++){
					// tab delimited
					// columns are Indiv	Locus	Allele1	Allele2	LLH	A1_perfect	A2_perfect
					llhOut << (*tGT).indName << "\t" << (*tGT).locName << "\t" << (*tGT).gTable[j].A1.name <<
						"\t" << (*tGT).gTable[j].A2.name << "\t" << (*tGT).gTable[j].llh << "\t" <<
						(*tGT).gTable[j].A1_perfect << "\t" << (*tGT).gTable[j].A2_perfect << "\n";
				}
			}
		}

		indivPos += curBatchSize;
	}



	// closing output file
	llhOut.close();

	return 0;
}


























