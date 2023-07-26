// functions for microtyper 2.0
// genotyping with better efficiency and ability to scale to large loci with many SNPs

#include <string>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"

// counts reads of alleles present in an individual
// only the alleles present in the individual are considered
// list of alleles is built on the fly

void countAlleleReads(const string bamInput,
					const unordered_map <string, locusInfo>& posMap,
					indAlleleCounts& indAlleleTable){
	// add individual name
	indAlleleTable.indName = bamInput;
	// open bam file
	BamTools::BamReader reader;
	if (!reader.Open(bamInput)) {
		cerr << "Error: could not open " << bamInput << endl;
		return;
	}

	// make a vector to translate between reference ID and reference name
	const BamTools::RefVector references = reader.GetReferenceData();
	vector <string> rL;
	for(int i = 0, max1 = references.size(); i < max1; i++) rL.push_back(references[i].RefName);

	// check that which reference sequences are in the posMap
	// position is reference ID and value is true it should be skipped (nothing in position file)
	vector <bool> skipLocus (rL.size(), true);
	for(int i = 0, max1 = rL.size(); i < max1; i++){
		if(posMap.count(rL[i]) == 1) skipLocus[i] = false;
	}

	// for each read
	BamTools::BamAlignment al;
	string rName;
	unordered_map <string, int> * alleleCounts;
	const locusInfo * lInfo;
	while (reader.GetNextAlignment(al)) {

		// check if paired and throw error if it is
		if(al.IsPaired()){
			cerr << "Error: paired read found in " << bamInput <<
			  ". Only single end reads are currently supported." << endl;
			exit(EXIT_FAILURE);
		}

		// skip if unmapped
		if(!al.IsMapped()) continue;

		// if different from last, get the appropriate info
		if(rName != rL[al.RefID]){
			if(skipLocus[al.RefID]) continue; // if refID is not in the position file, skip this alignment
			rName = rL[al.RefID];
			alleleCounts = &indAlleleTable.locusMap[rName];
			lInfo = &posMap.at(rName); // using .at() to pass as constant and avoid compiler error
		}

		// now identify allele
		unsigned int snpsTested = 0; // count of snps covered by read, to make sure none skipped for perfect counts
        int refCur = al.Position; // next position in reference sequence
        int queCur = 0;	// next position in query sequence (read)
        int dist;
        string curAllele ("");
        int k = 0; // current position in the CIGAR string
        int cigarLength = al.CigarData.size();
		for(int j = 0, max1 = (*lInfo).snps.size(); j < max1; j++){ // for each SNP
			dist = (*lInfo).snps[j].refPos - refCur;
			// progress through CIGAR until reach the next relevant position
			while (k < cigarLength){
				if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
					if(dist - static_cast<int>(al.CigarData[k].Length) < 0) break; // SNP is in this stretch or was skipped
					queCur += al.CigarData[k].Length;
					refCur += al.CigarData[k].Length;
				} else if (al.CigarData[k].Type == 'I'){
					if((*lInfo).snps[j].type == 2 && dist == 0) break; // type is insertion and next position is location of insertion
					queCur += al.CigarData[k].Length;
				}else if (al.CigarData[k].Type == 'S'){
					queCur += al.CigarData[k].Length;
				} else if (al.CigarData[k].Type == 'D' || al.CigarData[k].Type == 'N'){
					if(dist - static_cast<int>(al.CigarData[k].Length) < 0) break; // SNP is in this stretch or was skipped
					refCur += al.CigarData[k].Length;
				} else if (al.CigarData[k].Type == 'H' || al.CigarData[k].Type == 'P'){
					// do nothing
				} else {
					cerr << "Error: unrecognized CIGAR character in " << bamInput << endl;
					exit(EXIT_FAILURE);
				}
				dist = (*lInfo).snps[j].refPos - refCur;
				k++; // haven't reached SNP yet, progress to next CIGAR character
			}
			if(k >= cigarLength || dist < 0) break; // a SNP was not present in the read (dist should never be less than zero, though)

			switch((*lInfo).snps[j].type){
				case 0: // S
					if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
						if(al.QueryBases[queCur + dist] == 'N') break; // skip if base not called
						snpsTested++;
						curAllele += al.QueryBases[queCur + dist];
					}
					break;
				case 1: // D
					if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
						snpsTested++;
						curAllele += 'I';
					} else if (al.CigarData[k].Type == 'D' || al.CigarData[k].Type == 'N'){
						snpsTested++;
						curAllele += 'D';
					}
					break;
				case 2: // I
					if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
						snpsTested++;
						curAllele += 'D';
					} else if (al.CigarData[k].Type == 'I'){
						snpsTested++;
						curAllele += 'I';
					}
					break;
				default:
					// should not get here b/c checked when read in
					cerr << "Error: Unrecognized type for locus " << (*lInfo).name << endl;
			}
		}
		if(snpsTested == (*lInfo).snps.size()){
			// either increase the count or add to table (new allele found)
			if((*alleleCounts).count(curAllele) == 1){
				(*alleleCounts)[curAllele] += 1;
			} else {
				// validate allele
				bool validAllele = true;
				for(int i = 0, max1 = (*lInfo).snps.size(); i < max1; i++){
					bool validSNP = false;
					if(curAllele[i] == (*lInfo).snps[i].refValue) continue;
					for(int j = 0, max2 = (*lInfo).snps[i].validAlt.size(); j < max2; j++){
						if((*lInfo).snps[i].validAlt[j] == curAllele[i]){
							validSNP = true;
							break;
						}
					}
					if(!validSNP){
						validAllele = false;
						break;
					}
				}
				// add to allele count table
				if(validAllele) (*alleleCounts)[curAllele] = 1;
			}
		}
	}
	reader.Close();
}

// recursively build up genotypes
void buildA(const int nA,
						const int ploidy,
						vector <vector <int>>& allGenos){
	if(allGenos.size() == 0){
		vector <int> tempVec;
		allGenos.push_back(tempVec);
	}
	if(allGenos[0].size() < ploidy - 1){
		buildA(nA, ploidy - 1, allGenos);
	}
	vector <vector <int>> newGenos;
	newGenos.reserve(allGenos.size());
	if(allGenos[0].size() == 0){
		for(int i = 0; i < nA; i++){
			vector <int> tempVec (1, i);
			newGenos.push_back(tempVec);
		}
	} else {
		for(int i = 0; i < nA; i++){
			for(int j = 0, m = allGenos.size(); j < m; j++){
				if(allGenos[j].back() > i) continue;
				vector <int> tempVec (allGenos[j]);
				tempVec.push_back(i);
				newGenos.push_back(tempVec);
			}
		}
	}
	allGenos = newGenos;
}

// calculates LLR, then posterior with uniform prior, and returns lines to write as output
// returning lines to write rather than data structure b/c makes outputting
// much cleaner (no need to loop over many things b/c ploidy may be variable)
void calcPosterior(const indAlleleCounts& batchGenosInd,
				const int ploidy,
				const vector <string>& locusNames,
				const double eps,
				const double c,
				const double depth,
				const bool reportAll,
				vector <string>& batchGenos_toWrite){

	// reserve some capacity to limit the number of reallocations
	int capReq = 0;
	for(int i = 0, max1 = locusNames.size(); i < max1; i++){
		int nA = batchGenosInd.locusMap.at(locusNames[i]).size();
		if(nA == 0){
			capReq += 1; // missing genotype;
		} else {
			capReq += numGenotypes (nA, ploidy);
		}
	}
	batchGenos_toWrite.reserve(capReq);

	// build partial line for missing genotype to reuse as needed
	string missLine ("\t");
	for(int i = 0; i < ploidy; i++) missLine += "\t\t";

	for(int i = 0, max1 = locusNames.size(); i < max1; i++){ // for each locus
		const unordered_map <string, int> * tempMap;
		tempMap = &batchGenosInd.locusMap.at(locusNames[i]);
		if((*tempMap).empty()){ // determine if empty or not
			batchGenos_toWrite.push_back(batchGenosInd.indName + "\t" + locusNames[i] + missLine);
		} else { // loop over all genotypes considering alleles found in that individual
			vector <string> alleles;
			vector <int> counts;
			alleles.reserve((*tempMap).size());
			counts.reserve((*tempMap).size());
			// get alleles and counts in convenient fashion
			for(auto it = (*tempMap).cbegin(); it != (*tempMap).cend(); ++it){
				alleles.push_back(it->first);
				counts.push_back(it->second);
			}
			// compare to minimum depth if necessary
			if(!reportAll && vecIntSum(counts) < depth){
				batchGenos_toWrite.push_back(batchGenosInd.indName + "\t" + locusNames[i] + missLine);
				continue; // output missing genotype and move to next locus
			}
			// for all genotypes
			vector <vector <int>> allGenos; // vector <int> gives indices of alleles in the genotype
			int nA = alleles.size();
			buildA(nA, ploidy, allGenos);

			// calc LLH
			vector <double> LLH (allGenos.size());
			for(int j = 0, max2 = LLH.size(); j < max2; j++){
				vector <double> pi (counts.size(), 0);
				// build pi vector
				double aDose = 1.0 / ploidy;
				for(int k = 0, max3 = allGenos[j].size(); k < max3; k++){
					pi[allGenos[j][k]] += aDose * (1 - eps);
					for(int i2 = 0, max4 = pi.size(); i2 < max4; i2++){
						if(i2 == allGenos[j][k]) continue;
						pi[i2] += aDose * eps / (max4 - 1);
					}
				}
				LLH[j] = mult_llh_no_constant(counts, pi);
			}

			// scale to create posterior and output
			double denom = logSumExp(LLH);
			for(int j = 0, max2 = LLH.size(); j < max2; j++){
				LLH[j] = LLH[j] - denom; // log of posterior with flat prior
				if(reportAll){
					//and output
					string tempLine ("");
					for(int k = 0; k < ploidy; k++) tempLine += "\t" + alleles[allGenos[j][k]];
					for(int k = 0; k < ploidy; k++) tempLine += "\t" + to_string(counts[allGenos[j][k]]);
					tempLine += "\t" + to_string(exp(LLH[j]));
					batchGenos_toWrite.push_back(batchGenosInd.indName + "\t" + locusNames[i] + tempLine);
				}
			}
			if(!reportAll){
				// depth already checked above, just checking posterior prob
				int maxPostPos = 0;
				double maxPost = LLH[0];
				for(int j = 1, max2 = LLH.size(); j < max2; j++){
					if(LLH[j] > maxPost){
						maxPostPos = j;
						maxPost = LLH[j];
					}
				}
				if(maxPost > c){
					string tempLine ("");
					for(int k = 0; k < ploidy; k++) tempLine += "\t" + alleles[allGenos[maxPostPos][k]];
					for(int k = 0; k < ploidy; k++) tempLine += "\t" + to_string(counts[allGenos[maxPostPos][k]]);
					tempLine += "\t" + to_string(exp(LLH[maxPostPos]));
					batchGenos_toWrite.push_back(batchGenosInd.indName + "\t" + locusNames[i] + tempLine);
				} else {
					batchGenos_toWrite.push_back(batchGenosInd.indName + "\t" + locusNames[i] + missLine);
				}
			}
		}
	}
}
