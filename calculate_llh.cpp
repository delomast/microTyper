// calculate log-likelihood of reads given genotypes
// for microtyper

#include <string>
#include <regex>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"

void calculate_llh(const string bamInput,
					const vector <string>& locusNames,
					const unordered_map <string, locusInfo>& posMap,
					unordered_map <string, genoTable>& locusGenoTables,
					const double eps_I,
					const double eps_S){
	// add individual name to geno tables
	regex bamEnding ("\\.bam$");
	string indName = regex_replace(bamInput, bamEnding, "");
	for(int i = 0, max = locusNames.size(); i < max; i++) locusGenoTables[locusNames[i]].indName = indName;
	BamTools::BamReader reader;
	if (!reader.Open(bamInput)) {
		cerr << "Error: could not open " << bamInput << endl;
		return;
	}

	const BamTools::RefVector references = reader.GetReferenceData();
	// make a vector to translate between reference ID and reference name
	vector <string> rL;
	for(int i = 0, max = references.size(); i < max; i++) rL.push_back(references[i].RefName);

	// check that which reference sequences are in the locusAlleleTables / posMap
	// position is reference ID and value is true it should be skipped (nothing in position file)
	vector <bool> skipLocus (rL.size(), true);
	for(int i = 0, max = rL.size(); i < max; i++){
		if(posMap.count(rL[i]) == 1) skipLocus[i] = false;
	}

	// for each read
	BamTools::BamAlignment al;
	string rName;
	genoTable * gT;
	const locusInfo * lInfo;
	regex testRegex;
	bool regexBool;
	while (reader.GetNextAlignment(al)) {

		// check if paired and throw error if it is
		if(al.IsPaired()){
			cerr << "Error: paired read found in " << bamInput <<
			  ". Only single end reads are currently supported." << endl;
			exit(EXIT_FAILURE);
		}

		// check if aligned to reverse strand and throw error if it did
		if(al.IsReverseStrand()){
			cerr << "Error: read mapped to reverse strand found in " << bamInput <<
			  ". Only alignments to the forward strand are currently supported." << endl;
			exit(EXIT_FAILURE);
		}

		// skip if unmapped
		if(!al.IsMapped()) continue;

		// if different from last, get the appropriate info
		if(rName != rL[al.RefID]){
			if(skipLocus[al.RefID]) continue; // if refID is not in the position file, skip this alignment
			rName = rL[al.RefID];
			gT = &locusGenoTables[rName];
			lInfo = &posMap.at(rName);
			regexBool = (*lInfo).regExBool;
			if(regexBool) testRegex.assign((*lInfo).regEx);
		}

		// skip alignment if regex is present and it does not match
		if(regexBool){
			if(!regex_search(al.QueryBases, testRegex)) continue;
		}

		// add to likelihood of each genotype
		for(int i = 0, max = (*gT).gTable.size(); i < max; i++){ // for each genotype
			double A1_lh = 1;
			double A2_lh = 1;
			int A1_pCounts = 0;
			int A2_pCounts = 0;
			bool A1_p = true;
			bool A2_p = true;
			unsigned int snpsTested = 0; // count of snps covered by read, to make sure none skipped for perfect counts
			for(int j = 0, max2 = (*gT).gTable[i].A1.snpAlleles.size(); j < max2; j++){ // for each SNP
				int refCur = al.Position; // next position in reference sequence
				int queCur = 0;	// next position in query sequence (read)
				int dist;
				switch((*lInfo).snps[j].type){
					case 0: // S
						for(int k = 0, max3 = al.CigarData.size(); k < max3; k++){ // progress through CIGAR until reach the relevant position
							dist = (*lInfo).snps[j].refPos - refCur;
							if(dist < 0) break; // alignment skipped position
							if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
								if(dist - static_cast<int>(al.CigarData[k].Length) < 0){
									// snp is in this stretch
									if(al.QueryBases[queCur + dist] == 'N') break; // skip if base not called
									snpsTested++;
									if (al.QueryBases[queCur + dist] == (*gT).gTable[i].A1.snpAlleles[j]){
										// read matches the allele
										A1_lh *= 1 - probSubErr(al.Qualities[queCur + dist], eps_S);
									} else {
										A1_lh *= probSubErr(al.Qualities[queCur + dist], eps_S) / 3;
										A1_p = false;
									}
									if ((*gT).gTable[i].hom) break;
									if (al.QueryBases[queCur + dist] == (*gT).gTable[i].A2.snpAlleles[j]){
										// read matches the allele
										A2_lh *= 1 - probSubErr(al.Qualities[queCur + dist], eps_S);
									} else {
										A2_lh *= probSubErr(al.Qualities[queCur + dist], eps_S) / 3;
										A2_p = false;
									}
									break;
								} else {
									queCur += al.CigarData[k].Length;
									refCur += al.CigarData[k].Length;
								}
							} else if (al.CigarData[k].Type == 'I'){
								queCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'D'){
								refCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'N'){
								refCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'S'){
								queCur += al.CigarData[k].Length;
							}
						}
						break;
					case 1: // D
						for(int k = 0, max3 = al.CigarData.size(); k < max3; k++){
							dist = (*lInfo).snps[j].refPos - refCur;
							if(dist < 0) break; // alignment skipped position
							if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
								if(dist - static_cast<int>(al.CigarData[k].Length) < 0){
									snpsTested++;
									// snp is in this stretch and it is NOT a deletion allele
									if ('I' == (*gT).gTable[i].A1.snpAlleles[j]){
										// read matches the allele
										A1_lh *= 1 - eps_I;
									} else {
										A1_lh *= eps_I;
										A1_p = false;
									}
									if ((*gT).gTable[i].hom) break;
									if ('I' == (*gT).gTable[i].A2.snpAlleles[j]){
										A2_lh *= 1 - eps_I;
									} else {
										A2_lh *= eps_I;
										A2_p = false;
									}
									break;
								} else {
									queCur += al.CigarData[k].Length;
									refCur += al.CigarData[k].Length;
								}
							} else if (al.CigarData[k].Type == 'I'){
								queCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'D'){
								// if this covers current SNP it is a deletion allele
								if(dist - static_cast<int>(al.CigarData[k].Length) < 0){
									snpsTested++;
									if ('I' != (*gT).gTable[i].A1.snpAlleles[j]){
										A1_lh *= 1 - eps_I;
									} else {
										A1_lh *= eps_I;
										A1_p = false;
									}
									if ((*gT).gTable[i].hom) break;
									if ('I' != (*gT).gTable[i].A2.snpAlleles[j]){
										A2_lh *= 1 - eps_I;
									} else {
										A2_lh *= eps_I;
										A2_p = false;
									}
									break;
								} else {
									refCur += al.CigarData[k].Length;
								}

							} else if (al.CigarData[k].Type == 'N'){
								refCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'S'){
								queCur += al.CigarData[k].Length;
							}
						}
						break;
					case 2: // I
						for(int k = 0, max3 = al.CigarData.size(); k < max3; k++){
							dist = (*lInfo).snps[j].refPos - refCur;
							if(dist < 0) break; // alignment skipped position
							if(al.CigarData[k].Type == 'M' || al.CigarData[k].Type == '=' || al.CigarData[k].Type == 'X'){
								if(dist - static_cast<int>(al.CigarData[k].Length) < 0){
									// snp is in this stretch and it is NOT an insertion allele
									snpsTested++;
									if ('D' == (*gT).gTable[i].A1.snpAlleles[j]){
										A1_lh *= 1 - eps_I;
									} else {
										A1_lh *= eps_I;
										A1_p = false;
									}
									if ((*gT).gTable[i].hom) break;
									if ('D' == (*gT).gTable[i].A2.snpAlleles[j]){
										A2_lh *= 1 - eps_I;
									} else {
										A2_lh *= eps_I;
										A2_p = false;
									}
									break;
								} else {
									queCur += al.CigarData[k].Length;
									refCur += al.CigarData[k].Length;
								}
							} else if (al.CigarData[k].Type == 'I'){
								if(dist - static_cast<int>(al.CigarData[k].Length) < 0){
									snpsTested++;
									if ('D' != (*gT).gTable[i].A1.snpAlleles[j]){
										A1_lh *= 1 - eps_I;
									} else {
										A1_lh *= eps_I;
										A1_p = false;
									}
									if ((*gT).gTable[i].hom) break;
									if ('D' != (*gT).gTable[i].A2.snpAlleles[j]){
										A2_lh *= 1 - eps_I;
									} else {
										A2_lh *= eps_I;
										A2_p = false;
									}
									break;
								} else {
									queCur += al.CigarData[k].Length;
								}
							} else if (al.CigarData[k].Type == 'D'){
									refCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'N'){
								refCur += al.CigarData[k].Length;
							} else if (al.CigarData[k].Type == 'S'){
								queCur += al.CigarData[k].Length;
							}
						}
						break;
					default:
						// should not get here b/c checked when read in
						cerr << "Error: Unrecognized type for locus " << (*lInfo).name << endl;
				}
			}
			if ((*gT).gTable[i].hom){
				(*gT).gTable[i].llh += log(A1_lh);
				// only saving perfect counts in A1_perfect for homozygous loci
				if(A1_p && snpsTested == (*gT).gTable[i].A1.snpAlleles.size()) (*gT).gTable[i].A1_perfect++;
			} else {
				(*gT).gTable[i].llh += log(0.5 * (A1_lh + A2_lh));
				if(snpsTested == (*gT).gTable[i].A1.snpAlleles.size()){
					if(A1_p) (*gT).gTable[i].A1_perfect++;
					if(A2_p) (*gT).gTable[i].A2_perfect++;
				}
			}
		}
	}
	reader.Close();

}
