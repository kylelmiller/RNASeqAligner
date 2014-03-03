#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <iterator>

using namespace std;

void sepratefusion(string input, size_t range)
{
	ifstream ifs(input.c_str());

	string interchr = input; interchr.append(".interchr"); ofstream ofs_interchr(interchr.c_str());

	string longrange = input; longrange.append(".longrange"); ofstream ofs_longrange(longrange.c_str());

	string inversion = input; inversion.append(".inversion"); ofstream ofs_inversion(inversion.c_str());

	string del_dup = input; del_dup.append(".del_dup"); ofstream ofs_del_dup(del_dup.c_str());

	size_t interchrcount = 0, longrangecount = 0, inversioncount = 0, del_dupcount = 0;

	if (ifs.is_open())
	{
		string line;

		//getline(ifs,line);

		
		while (getline(ifs,line))
		{
			if (line == "" || line.find("track") != string::npos)
				continue;

			char skip[1000], gene1[10000], gene2[10000], chrname[1000];

			size_t start, end;

			char strand1, strand2;

			//chr11~chr1	101985125	50732261	JUNC_606	3	+-	255,0,0	2	35,25,93,318,	0,51252900,	1.098612	6	GTAG	0	0	0.000000	
			//25	25	0	3	0	3	2	1	0	3	0	10	101985090	50732236	50731945,317M|	101985018,108M|	0	0	0.295873	0.653333	
			//317	317	108	108	doner_exact_matched	not_matched	CTTCTGGTCAGAGATACTTCTTAAA	gTTTTCAGTTGACAGCCAGAGTGCT	0	from_fusion	intergenic		YAP1,	-,	

			sscanf(line.c_str(), "%s\t%llu\t%llu\t%s\t%s\t%c%c", chrname, &start, &end, skip, skip, &strand1, &strand2);

			string chrstr = chrname;

			string chr1 = chrstr.substr(0, chrstr.find("~"));

			string chr2 = chrstr.substr(chrstr.find("~") + 1, chrstr.length() - 1);

			size_t intron;

			if (start > end)
				intron = start - end;
			else 
				intron = end - start;

			if (chr1 != chr2)
			{
				ofs_interchr << line << endl;
				++interchrcount;
			}
			else if (intron > range)
			{
				ofs_longrange << line << endl;
				++longrangecount;
			}
			else if (strand1 != strand2)
			{
				ofs_inversion << line << endl;
				++inversioncount;
			}
			else
			{
				ofs_del_dup << line << endl;
				++del_dupcount;
			}
		}
	}

	cout << "interchrcount:\t"<<interchrcount<<endl;
	cout << "inversioncount:\t"<<inversioncount<<endl;
	cout << "longrangecount:\t"<<longrangecount<<endl;
	cout << "del_dupcount:\t"<<del_dupcount<<endl;
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "infile unspliced_reads unique_spliced_reads multiple_spliced_reads" <<endl;
		exit(0);
	}
	
	string input = argv[1];

	size_t range = atoi(argv[2]);

	sepratefusion(input, range);
}