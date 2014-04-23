#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <vector>
#include <hash_map>
#include <map>
#include <set>

#include "buildIndex_info.h"
#include "build2ndLevelIndex_info.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 2)
	{
		cout << "Executable <InputIndexPrefix> <Output2ndLevelIndexPrefix> " << endl;
		exit(0);
	}

	string inputIndexPrefix = argv[1];
	string Output2ndLevelIndexPrefix = argv[2];

	cout << "inputIndexPrefix: " << inputIndexPrefix << endl;
	cout << "Output2ndLevelIndexPrefix: " << Output2ndLevelIndexPrefix << endl; 

	inputIndexPrefix.append("/");

	string chrom_bit_file = inputIndexPrefix; chrom_bit_file.append("_chrom"); ifstream chrom_bit_file_ifs(chrom_bit_file.c_str(),ios::binary);

	string parameter_file = inputIndexPrefix; parameter_file.append("_parameter"); ifstream parameter_file_ifs(parameter_file.c_str(),ios::binary);

	SecondLevelIndex_Info* secondLevelIndexInfo = new SecondLevelIndex_Info(parameter_file_ifs);
	
	cout << "start to load whole genome" << endl;
	char *chrom;

	chrom = (char*)malloc((secondLevelIndexInfo->indexSize) * sizeof(char));
	chrom_bit_file_ifs.read((char*)chrom, (secondLevelIndexInfo->indexSize) * sizeof(char)); 

	string genomeString = chrom;

	free(chrom);
	cout << endl << "start to divide Genome to chromosomes" << endl;
	secondLevelIndexInfo->Genome2ChromString(genomeString);
	cout << "finish dividing Genome to chromosomes" << endl;

	cout << endl <<  "start to create subfolders" << endl;
	secondLevelIndexInfo->creatAllSubFolder(Output2ndLevelIndexPrefix);
	cout << "finish creating subfolders" << endl;

	cout << endl <<  "start to generate 2nd level index chrom file" << endl;
	secondLevelIndexInfo->generate2ndLevelIndexChrom(Output2ndLevelIndexPrefix);
	cout << "finish generating 2nd level index chrom file" << endl;

	cout << endl <<  "start to generate 2nd level index with original size" << endl;
	secondLevelIndexInfo->generate2ndLevelIndexOriginalSize(Output2ndLevelIndexPrefix);
	cout << "finish generating 2nd level index with original size" << endl;

	cout << endl << "start to generate 2nd level index with Compressed size" << endl;
	secondLevelIndexInfo->generate2ndLevelIndexCompressedSize(Output2ndLevelIndexPrefix);
	cout << "finish generating 2nd level index with Compressed size" << endl;

	cout << endl << "finish all 2nd level index" << endl;

	chrom_bit_file_ifs.close();
	parameter_file_ifs.close();

	delete(secondLevelIndexInfo);

	//free(chrom);
	return 0;
}