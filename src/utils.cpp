// utility functions for microtyper

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <string>
#include "typedefs.h"

using namespace std;

// print version number
void printVersion(){
	cout << "Version 1.0" << endl;
	exit(EXIT_SUCCESS);
}

// calculate probability of wrong base from PHRED score
double ph33toProb(const char& phred){
	double Q = static_cast<int>(phred) - 33;
	if(Q < 0){
		cerr << "Error: Unexpected PHRED score. Is it +33 or +64? Only +33 is supported." << endl;
		exit(EXIT_FAILURE);
	}
	return pow(10, (Q / -10));
}

// calculate probability of substitution error from eps_S and phred score
double probSubErr(const char& phred, const double& eps_S){
	double e2 = ph33toProb(phred);
	return (eps_S * (1 - e2)) + ((2/3) * eps_S * e2) + ((1 - eps_S) * e2);
}

// comparison for baseInfo objects based on reference position
// used to sort SNPs after reading in the position file so that
// haplotypes are in order of reference position
bool compareRefPos(baseInfo snp1, baseInfo snp2){
	return snp1.refPos < snp2.refPos;
}

// splits a string based on a 1-char delimiter
void splitString(const string& inString, const char delim, vector <string>& split){
	string tCell ("");
	for(int i = 0, max = inString.length(); i < max; i++){
		if (inString[i] == delim){
			split.push_back(tCell);
			tCell = "";
		} else if (inString[i] == '\n' || inString[i] == '\r') {
			// to make multi-system compatible, just skipping all end of line characters
			// this also means have to add a push back once loop is completed
			continue;
		} else {
			tCell += inString[i];
		}
	}
	split.push_back(tCell); // add last entry
}

// log - sum - exp function
// for a set of values x, returns
// log(sum(e^x1 + e^x2 + ...))
double logSumExp(const vector <double>& x){
	if(x.size() < 1){
		cerr << "Internal error: vector of length 0 to logSumExp." << endl;
		exit(EXIT_FAILURE);
	}
	// find max value
	double maxV = x[0];
	for(int i = 1, max = x.size(); i < max; i++) if (x[i] > maxV) maxV = x[i];

	// calculate and return
	double sum = 0;
	for(int i = 0, max = x.size(); i < max; i++) sum += exp(x[i] - maxV);
	return maxV + log(sum);
}


// read in position file
// create map of loci and information as well as vector of locus names
void readPosFile(const string& posInput,
				unordered_map <string, locusInfo>& posMap,  // edited by this function instead of returned
				vector <string>& locusNames){ // edited by this function instead of returned
	ifstream posFile (posInput);
	if(!posFile.is_open()){
		cerr << "Error: variant file specified (-p) could not be opened." << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	getline(posFile, line); // skip header
	while (getline(posFile, line)){
		if(line.substr(0,1) == "\t" || line.length() == 0) continue; // skip blank lines
		int pos = 0;
		vector <string> cell; //Locus	RefPos	Type	ValidAlt
		splitString(line, '\t', cell);
		if(cell.size() != 4){
			cerr << "Error: Line in the position file with wrong number of fields: " << endl;
			cerr << line << endl;
			exit(EXIT_FAILURE);
		}
		// splitting validAlt
		vector <char> validAlt;
		bool nextIsAlt = true;
		for(int i = 0, max = cell[3].size(); i < max; i++){
			if(nextIsAlt){
				if(cell[3][i] == ','){
					cerr << "Error: The validAlt field for a SNP in " << cell[0] << " is not formatted appropriately." <<
						" Likely error is either a comma at the beginning or multiple commas in a row." << endl;
					exit(EXIT_FAILURE);
				}
				validAlt.push_back(cell[3][i]);
				nextIsAlt = false;
			} else {
				if(cell[3][i] != ','){
					cerr << "Error: The validAlt field for a SNP in " << cell[0] << " is not formatted appropriately." <<
						" Likely error is a alternate allele that is more than one character in length." << endl;
					exit(EXIT_FAILURE);
				}
				nextIsAlt = true;
			}
		}
		if(posMap.count(cell[0]) == 1){ // already exists in map
			baseInfo tempBase;
			tempBase.refPos = stoi(cell[1]) - 1;
			if(cell[2] == "S"){
				tempBase.type = 0;
			} else if(cell[2] == "D"){
				tempBase.type = 1;
			} else if(cell[2] == "I"){
				tempBase.type = 2;
			} else {
				cerr << "Error: Unrecognized type for locus " << cell[0] << endl;
				exit(EXIT_FAILURE);
			}
			tempBase.validAlt = validAlt;

			posMap[cell[0]].snps.push_back(tempBase);
		} else {
			locusInfo tempLoc;
			tempLoc.name = cell[0];
			locusNames.push_back(cell[0]); // add to list of locus names
			baseInfo tempBase;
			tempBase.refPos = stoi(cell[1]) - 1;
			if(cell[2] == "S"){
				tempBase.type = 0;
			} else if(cell[2] == "D"){
				tempBase.type = 1;
			} else if(cell[2] == "I"){
				tempBase.type = 2;
			} else {
				cerr << "Error: Unrecognized type for locus " << cell[0] << endl;
				exit(EXIT_FAILURE);
			}
			tempBase.validAlt = validAlt;

			tempLoc.snps.push_back(tempBase);
			posMap[cell[0]] = tempLoc;
		}
	}
	posFile.close();

	// put snps in order
	locusInfo * tlocInfo;
	for(int i=0, max = locusNames.size(); i < max; i++){
		tlocInfo = &posMap[locusNames[i]];
		sort((*tlocInfo).snps.begin(), (*tlocInfo).snps.end(), compareRefPos);
	}
}



// read in reference sequences and info add to posMap
void readRefSeqs(const string& refInput,
				const vector <string>& locusNames,
				unordered_map <string, locusInfo>& posMap){
	ifstream refFile (refInput);
	if(!refFile.is_open()){
		cerr << "Error: reference file specified (-r) could not be opened." << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	string curSeqName;
	string curSeq;
	bool start = true;
	while(getline(refFile, line)){
		if(line.substr(0,1) == ">"){
			if (start){
				start = false;
			} else {
				if(posMap.count(curSeqName) == 1){
					posMap[curSeqName].refSequence = curSeq; // assign reference sequence to map
				} else {
					cerr << "Warning: " << curSeqName << " not in position file. It will not be genotyped." << endl;
				}
			}
			curSeqName = line.substr(1); // name of ref sequence
			curSeq = "";
		} else {
			curSeq = curSeq + line;
		}
	}
	refFile.close();
	// add the last value, if applicable
	if(posMap.count(curSeqName) == 1){
		posMap[curSeqName].refSequence = curSeq; // assign reference sequence to map
	} else {
		cerr << "Warning: " << curSeqName << " not in position file. It will not be genotyped." << endl;
	}

	// adding reference value to struct baseInfo for each locus: for S, the ref base; D, the value 'I'; I, the value 'D'
	// and checking that all loci in position file have reference sequences
	for(int i = 0, max = locusNames.size(); i < max; i++){ // for each locus
		for(int j = 0, max2 = posMap[locusNames[i]].snps.size(); j < max2; j++){ // for each snp
			if(posMap[locusNames[i]].refSequence.length() == 0){
				cerr << "Error: " << locusNames[i] << " reference sequence either not found or has 0 length.";
				exit(EXIT_FAILURE);
			}
			if(posMap[locusNames[i]].refSequence.length() < posMap[locusNames[i]].snps[j].refPos + 1){
				cerr << "Error: " << locusNames[i] << " reference sequence is not long enough to contain all of its SNPs.";
				exit(EXIT_FAILURE);
			}
			switch(posMap[locusNames[i]].snps[j].type){
				case 0: // S
					posMap[locusNames[i]].snps[j].refValue =
						posMap[locusNames[i]].refSequence[posMap[locusNames[i]].snps[j].refPos];
					break;
				case 1: // D
					posMap[locusNames[i]].snps[j].refValue = 'I';
					break;
				case 2: // I
					posMap[locusNames[i]].snps[j].refValue = 'D';
					break;
				default:
					// should not get here b/c checked when read in
					cerr << "Error: Unrecognized type for locus " << posMap[locusNames[i]].name << endl;
					exit(EXIT_FAILURE);
			}
		}
	}
}


// build vectors of all possibly alleles for each locus
void buildAlleles(const unordered_map <string, locusInfo>& posMap, const vector <string>& locusNames,
			unordered_map <string, alleleTable>& locusAlleleTables){
	for(int i = 0, max = locusNames.size(); i < max; i++){ // for each locus
		alleleTable tempAT;
		tempAT.name = locusNames[i];

		// this is to generate a vector of all possible alleles
		const locusInfo * tempLI;
		tempLI = &posMap.at(locusNames[i]);
		vector <vector <char> > tempVec;
		for(int j = 0, max2 = (*tempLI).snps.size(); j < max2; j++){ // for each SNP in the locus

			// first check that all reference and alt alleles are different
			vector <char> allSNPalleles;
			allSNPalleles = (*tempLI).snps[j].validAlt;
			allSNPalleles.push_back((*tempLI).snps[j].refValue);
			for(int k = 0, m = allSNPalleles.size(); k < m - 1; k++){
				for(int l = k + 1; l < m; l++){
					if(allSNPalleles[k] == allSNPalleles[l]){
						cerr << "Error: two SNP alleles with the same symbol in locus " << locusNames[i] << endl;
						cerr << "The position of the SNP is " << (*tempLI).snps[j].refPos << endl;
						cerr << "The reference allele is " << (*tempLI).snps[j].refValue << endl;
						cerr << "The alternative alleles are : ";
						for(int n = 0, nmax = (*tempLI).snps[j].validAlt.size(); n < nmax; n++) cerr << (*tempLI).snps[j].validAlt[n] << " ";
						cerr << endl;
						exit(EXIT_FAILURE);
					}
				}
			}

			if(j == 0){ // start vector of alleles
				vector <char> temp2(1);
				temp2[0] = (*tempLI).snps[j].refValue;
				tempVec.push_back(temp2);
				for(int k = 0, max3 = (*tempLI).snps[j].validAlt.size(); k < max3; k++){
					temp2[0] = (*tempLI).snps[j].validAlt[k];
					tempVec.push_back(temp2);
				}
				continue;
			}
			vector <vector <char> > newTempVec;
			for(int k = 0, max3 = tempVec.size(); k < max3; k++){ // for each allele already in the list
				vector <char> newAllele (tempVec[k]);
				newAllele.push_back((*tempLI).snps[j].refValue); // add the reference allele
				newTempVec.push_back(newAllele);
				int lastPos = newAllele.size() - 1;
				for(int l = 0, max4 = (*tempLI).snps[j].validAlt.size(); l < max4; l++){
					newAllele[lastPos] = (*tempLI).snps[j].validAlt[l];
					newTempVec.push_back(newAllele);
				}
			}
			tempVec = newTempVec; // replace old list with new;
		}
		// now make the alleles and compile into the alleleTable
		for(int j = 0, max2 = tempVec.size(); j < max2; j++){
			allele tempAllele;
			tempAllele.name = "";
			for(int k = 0, max3 = tempVec[j].size(); k < max3; k++) tempAllele.name += tempVec[j][k];
			tempAllele.snpAlleles = tempVec[j];
			tempAT.alleleList.push_back(tempAllele);
		}
		locusAlleleTables[locusNames[i]] = tempAT; // and save allele table to unordered map
	}


}
















