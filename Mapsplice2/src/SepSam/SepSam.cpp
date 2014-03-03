#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <string.h>

using namespace std;

int 
main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "in_sam unspliced spliced del ins clip unmapped head"<<endl;
		exit(0);
	}
	ifstream input_fs(argv[1]);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	cout << argv[1] << endl;

	ofstream ofs_unspliced(argv[2]);

	ofstream ofs_spliced(argv[3]);

	ofstream ofs_small_del(argv[4]);

	ofstream ofs_small_ins(argv[5]);

	ofstream ofs_clipped(argv[6]);

	ofstream ofs_unmapped(argv[7]);

	ofstream ofs_head(argv[8]);
 
	string line;
	int unspliced = 0, spliced = 0, small_del = 0, small_ins = 0, unmapped = 0, clipped = 0;
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
		else if (mappedstr.find("S") != string::npos)
		{
			ofs_clipped << line << endl; ++clipped;
		}
		else if (mappedstr.find("D") != string::npos)
		{
			ofs_small_del << line << endl; ++small_del;
		}
		else if(mappedstr.find("I") != string::npos)
		{
			ofs_small_ins << line << endl; ++small_ins;
		}
		else if(mappedstr.find("N") != string::npos)
		{
			ofs_spliced << line << endl; ++spliced;
		}
		else if (mappedstr.find("M") != string::npos)
		{
			ofs_unspliced << line << endl; ++unspliced;
		}
		else 
			cout << line << endl;
	}

	cout << "unmapped " <<unmapped<<endl;
	cout << "small_del " <<small_del<<endl;
	cout << "small_ins " <<small_ins<<endl;
	cout << "clipped " <<clipped<<endl;
	cout << "spliced " <<spliced<<endl;
	cout << "unspliced " <<unspliced<<endl;
}
