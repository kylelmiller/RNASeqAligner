#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <set>
#include <map>
using namespace std;


int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr,"error: too few arguments\n");exit(1);
	}

	ifstream input_reads(argv[1]);
	
	ofstream output_reads(argv[2]);

	map<string, int> found_reads;

	size_t count = 0;

	while(!input_reads.eof())
	{
		string readline;

		getline(input_reads, readline);

		if (readline.empty())
			continue;
		char readname[1000];

		sscanf(readline.c_str(), "%s", readname);

		string readnamestr = readname;

		string base_read_name = readnamestr.substr(0, readnamestr.length() - 2);

		if (readnamestr[readnamestr.length() - 1] == '1')
		{
			map<string, int>::iterator found_iter = found_reads.find(base_read_name);

			if (found_iter == found_reads.end())
			{
				found_reads.insert(make_pair(base_read_name, 1));
			}
			else
			{
				if (found_iter->second == 1)
				{
					cout << "found 1 before? " << endl;
					cout << found_iter->first<< endl;
				}
				else if (found_iter->second == 2)
				{
					found_iter->second = 3;
				}
				else
				{
					cout << "other end found before"<<endl;
					cout << found_iter->first << '\t' << found_iter->second << endl;
				}
			}
		}
		else if (readnamestr[readnamestr.length() - 1] == '2')
		{
			map<string, int>::iterator found_iter = found_reads.find(base_read_name);

			if (found_iter == found_reads.end())
			{
				if (base_read_name == "UNC11-SN627_80:8:2208:20539:200800")
					cout << readline << endl;
				found_reads.insert(make_pair(base_read_name, 2));
			}
			else
			{
				if (found_iter->second == 2)
				{
					cout << "found 2 before? " << endl;
					cout << found_iter->first<< endl;
					cout << readline<<endl;
				}
				else if (found_iter->second == 1)
				{
					found_iter->second = 3;
				}
				else
				{
					cout << "other end found before"<<endl;
					cout << found_iter->first << '\t' << found_iter->second << endl;
				}
			}
		}

		++count;
	}

	for (size_t i = 3; i < argc; ++i)
	{
		ifstream input_reads2(argv[i]);

		while(!input_reads2.eof())
		{
			string readline;

			string line2,line3,line4;

			getline(input_reads2, readline);

			getline(input_reads2, line2);

			getline(input_reads2, line3);

			getline(input_reads2, line4);

			if (readline.empty())
				continue;

			string base_read_name = readline.substr(1, readline.length() - 3);

			//cout << base_read_name << endl;

			int endid = 0;

			map<string, int>::iterator found_iter = found_reads.find(base_read_name);

			if (found_iter != found_reads.end())
			{
				output_reads << readline << endl;

				output_reads << line2 << endl;

				output_reads << line3 << endl;

				output_reads << line4 << endl;
			}
			//else if (found_iter != found_reads.end() )
			//{
			//	//cout << "reaad appeared " << endl;
			//	//cout << readline << endl;
			//}
		}
	}

	return 0;
}