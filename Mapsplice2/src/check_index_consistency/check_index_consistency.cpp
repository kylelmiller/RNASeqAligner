#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
using namespace std;


int main(int argc, const char** argv)
{
	if (argc < 3)
	{
		cout << "too few arguments"<<endl;
		exit(0);
	}
	int return_flag = 0;
	ifstream input_fs1(argv[1]);
	ifstream input_fs2(argv[2]);
	map<string, pair<int, int> > ref_names;
	string line;
	getline(input_fs1, line);
	getline(input_fs1, line);
	getline(input_fs1, line);
	while(getline(input_fs1, line))
	{
		char name_tmp[1000];
		int seq_len = -1;
		sscanf(line.c_str(), "%*s\t%[^\t]\t%d", name_tmp, &seq_len);
		string cur_name = name_tmp;
		if(cur_name.find(' ') != string::npos)
		{
			cerr << "Error: Reference name in Bowtie Index contains space:" << endl;
			cerr << "'" << cur_name << "' contains space" << endl;
			return_flag = 1;
			return return_flag;
		}
		if(ref_names.find(cur_name) == ref_names.end())
			ref_names.insert(make_pair(cur_name, make_pair(seq_len, 1)));
		else
		{
			cerr << "Error: duplicate reference name in Bowtie Index:" << endl;
			cerr << "'" <<line << "' detected multiple times" << endl;
			return_flag = 1;
		}
	}
	while(getline(input_fs2, line))
	{
		char name_tmp[1000];
		int seq_len = -1;
		sscanf(line.c_str(), "%[^\t]\t%d", name_tmp, &seq_len);
		string cur_name = name_tmp;
		if(cur_name.find(' ') != string::npos)
		{
			cerr << "Error: Reference name in Reference Sequence contains space:" << endl;
			cerr << "'" << cur_name << "' contains space" << endl;
			return_flag = 1;
			return return_flag;
		}
		if(ref_names.find(cur_name) != ref_names.end())
		{
			if(ref_names[cur_name].first != seq_len)
			{
			  cerr << "Error: Bowtie Index not consistent with Reference Sequence" << endl;
				cerr << "'" << cur_name << "' has length " << ref_names[cur_name].first << " in Bowtie Index, but has length " << seq_len << " in Reference Sequence" << endl;
			  return_flag = 1;
			}
			ref_names[cur_name].second ++;
		}
		else
		{
			cerr << "Error: Bowtie Index not consistent with Reference Sequence" << endl;
			cerr << "'" << cur_name << "' does not exist in Bowtie Index" << endl;
			return_flag = 1;
		}
	}
	map<string, pair<int, int> >::iterator it;
	for(it = ref_names.begin(); it != ref_names.end(); it++)
	{
		if(it->second.second == 1)
		{
			cerr << "Error: Bowtie Index not consistent with Reference Sequence" << endl;
			cerr << "'" << it->first << "' does not exist in Reference Sequence" << endl;
			return_flag = 1;
		}
		else if(it->second.second > 2)
		{
			cerr << "Error: duplicate reference name in Reference Sequence" << endl;
			cerr << "'" << it->first  << "' detected multiple times" << endl;
			return_flag = 1;
		}
	}
	input_fs1.close();
	input_fs2.close();
	return return_flag;
}