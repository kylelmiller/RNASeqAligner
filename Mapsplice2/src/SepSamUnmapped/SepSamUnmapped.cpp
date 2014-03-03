#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
#include <string>
using namespace std;

int 
main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "in_sam unspliced spliced del ins clip unmapped head"<<endl;
		exit(0);
	}
	
	cout << argv[1] << endl;

	ofstream ofs_normal(argv[1]);

	ofstream* ofs_fusion;//;(argv[2]);

	if (string(argv[1]) == string(argv[2]))
		ofs_fusion = &ofs_normal;
	else
	{
		ofs_fusion = new ofstream;

		ofs_fusion->open(argv[2]);
	}

	ofstream ofs_unmapped(argv[3]);

	ofstream ofs_head(argv[4]);
 
	string line;
	int normal = 0, fusion = 0, unmapped = 0;

	for (size_t i = 5; i < argc; ++i)
	{
		ifstream input_fs(argv[i]);

		if( !input_fs ) 
		{
			fprintf(stderr,"error: open extended file error\n");exit(1);
		}

		while(getline(input_fs, line))
		{
			if (line == "" || line[0] == '@')
			{
				ofs_head <<line<< endl;
				continue;
			}

			//count++;

			char tagname[1000], chrom[100], mapped[100];
			int strand, offset, something, pid;

			//HWI-EAS217:4:100:1000:1884#0/1_6458116	16	chr9	35092770	0	51M226N49M	*	0	0	CAGTCCAGAGGAGGCGGGGCGCGGAGCGCGGCCAGAAGCCAGTAGAGAGCCCCTCAGCAAAAGGGCCCCAGTGCCCCGCGCCGCGC
			//GCGCCAGCATTTCC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
			sscanf(line.c_str(), "%s\t%d\t%s\t%d\t%d\t%s", tagname, &strand, chrom, &offset, &something, mapped);

			string mappedstr = mapped;

			if (mappedstr.find("*") != string::npos)
			{
				ofs_unmapped << line << endl; ++unmapped;
			}
			else if (line.find("ZF:Z:FUS") != string::npos)
			{
				*ofs_fusion << line << endl; ++fusion;
			}
			else
			{
				ofs_normal << line << endl; ++normal;
			}
		}
	}

	cout << "unmapped " <<unmapped<<endl;
	cout << "normal " <<normal<<endl;
	cout << "fusion " <<fusion<<endl;
}
