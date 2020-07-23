// declaration of structs used in microtyper

#ifndef TYPEDEFS
#define TYPEDEFS

#include <vector>
#include <string>

using namespace std;

struct baseInfo {
	int refPos; // position in reference with 0-offset
	int type; // S:0, D:1, I:2
	vector <char> validAlt;
	char refValue; // for S, the ref base; D, the value 'I'; I, the value 'D'
};

struct locusInfo{
	string name;
	string refSequence;
	vector <baseInfo> snps;
};

struct allele{
	string name; // name of allele - concatenation of snpAlleles, in order
	vector <char> snpAlleles; // vector of the alleles of the snps, in order
};

struct alleleTable{
	string name; // name of the locus
	vector <allele> alleleList;
};

struct genotype{
	allele A1;
	allele A2;
	bool hom; // true for homozygous, false for heterozygous
	double llh;
	int A1_perfect; // count of reads matching A1 at all SNPs
	int A2_perfect;
};

struct genoTable{
	string locName; // name of the locus
	string indName; // name of the individual
	vector <genotype> gTable;
	int numAlleles; // number of alleles at this locus
};




#endif // TYPEDEFS
