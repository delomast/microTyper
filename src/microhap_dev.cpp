// microhap_dev.cpp : Defines the entry point for the application.
//
//g++ microhap_dev.cpp -I /usr/local/include/bamtools -L /usr/local/lib -lbamtools -lz

#include <iostream>
#include <string>
#include <vector>
#include "api/BamMultiReader.h"

using namespace std;

int main(int argc, char* argv[]){

	string bamTestFile;

	// get input options and check
	std::string x;
	for (int i = 1; i < argc; i++) {
		x = string(argv[i]);
		if (x == "-f") {
			bamTestFile = string(argv[++i]);
		}
	}

	if (bamTestFile.length() == 0) {
		cerr << "No bam file input for testing." << endl;
		return 1;
	}

	BamTools::BamMultiReader reader;
	if (!reader.OpenFile(bamTestFile)) {
		cerr << "could not open BAM" << endl;
		return 1;
	}

//	cout << "file opened" << endl; // testing

	const BamTools::SamHeader header = reader.GetHeader();
	const BamTools::RefVector references = reader.GetReferenceData();

	BamTools::BamAlignment al;
	for (int i = 0; i < 10; i++) {
		//reader.GetNextAlignmentCore(al);
		reader.GetNextAlignment(al);
		cout << al.AlignedBases << endl;
	}

	reader.Close();

	return 0;
}
