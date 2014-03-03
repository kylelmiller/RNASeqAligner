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

	string del = input; del.append(".del"); ofstream ofs_del(del.c_str());

	string dup = input; dup.append(".dup"); ofstream ofs_dup(dup.c_str());

	string header = input; header.append(".header"); ofstream ofs_header(header.c_str());

	size_t interchrcount = 0, longrangecount = 0, inversioncount = 0, del_count = 0, dup_count = 0;

	if (ifs.is_open())
	{
		string line;

		getline(ifs,line);

		ofs_header << line << endl;

		while (getline(ifs,line))
		{
			if (line == "" || line.find("track") != string::npos)
				continue;

			if (line.find("inter_chr") != string::npos)
			{
				ofs_interchr << line << endl;
				++interchrcount;
			}
			else if (line.find("inversion") != string::npos)
			{
				ofs_inversion << line << endl;
				++inversioncount;
			}
			else if (line.find("long_range") != string::npos)
			{
				ofs_longrange << line << endl;
				++longrangecount;
			}
			else if (line.find("deletion") != string::npos)
			{
				ofs_del << line << endl;
				++del_count;
			}
			else if (line.find("tandem_dup") != string::npos)
			{
				ofs_dup << line << endl;
				++dup_count;
			}
		}
	}

	cout << "interchrcount:\t"<<interchrcount<<endl;
	cout << "inversioncount:\t"<<inversioncount<<endl;
	cout << "longrangecount:\t"<<longrangecount<<endl;
	cout << "del_count:\t"<<del_count<<endl;
	cout << "dup_count:\t"<<dup_count<<endl;
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