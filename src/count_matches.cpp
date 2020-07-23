// counts reads that perfectly match each allele
// to be used by alternate genotyping methods than the original
// llr method

#include <string>
#include <vector>
#include <unordered_map>
#include <math.h>
#include "api/BamReader.h"
#include "typedefs.h"
#include "utils.h"

void count_matches(const string bamInput,
					const vector <string>& locusNames,
					const unordered_map <string, locusInfo>& posMap,
					unordered_map <string, genoTable>& locusGenoTables){
	// add individual name to geno tables
	for(int i = 0, max = locusNames.size(); i < max; i++) locusGenoTables[locusNames[i]].indName = bamInput;
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
		}

		// add to likelihood of each genotype
		for(int i = 0, max = (*gT).gTable.size(); i < max; i++){ // for each genotype
			// skipping if a het locus
			//  it would be better to just loop over alleles, but... this is much easier to add in given what is already written
			if(!(*gT).gTable[i].hom) continue;
			int A1_pCounts = 0;
			bool A1_p = true;
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
										break;
									} else {
										A1_p = false;
										break;
									}
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
										break;
									} else {
										A1_p = false;
										break;
									}
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
										break;
									} else {
										A1_p = false;
										break;
									}
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
										break;
									} else {
										A1_p = false;
										break;
									}
								} else {
									queCur += al.CigarData[k].Length;
									refCur += al.CigarData[k].Length;
								}
							} else if (al.CigarData[k].Type == 'I'){
								if(dist - static_cast<int>(al.CigarData[k].Length) < 0){
									snpsTested++;
									if ('D' != (*gT).gTable[i].A1.snpAlleles[j]){
										break;
									} else {
										A1_p = false;
										break;
									}
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
			if(A1_p && snpsTested == (*gT).gTable[i].A1.snpAlleles.size()) (*gT).gTable[i].A1_perfect++;
		}
	}
	reader.Close();

}
