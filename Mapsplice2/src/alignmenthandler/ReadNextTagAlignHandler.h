#ifndef READNEXTTAGALIGNHANDLER_H
#define READNEXTTAGALIGNHANDLER_H

#include "sharedlib.h"

template<typename T>
class ReadNextTagAlignHandler {

public:
	typedef T value_type;
	typedef size_t size_type;
	typedef T& reference;
	typedef const T& const_reference;

private:
	ifstream ifs;
	string prev_tag_name;
	string prev_line;
	size_t prev_read_id;
	size_t prev_end_id;
	string tag_name;	
	int min_ins;

public:

	ReadNextTagAlignHandler()
	{
	}

	ReadNextTagAlignHandler(const char* file_name, int min_insert) : prev_read_id(0), min_ins(min_insert)
	{
		ifs.open(file_name);

		if (!ifs.is_open())
		{
			cerr<<"can't open file: "<< file_name << endl;
			exit(1);
		}
	}

	bool Init(const char* file_name, int min_insert)
	{
		Clear();

		min_ins = min_insert;

		ifs.open(file_name);

		if (!ifs.is_open())
		{
			cerr<<"can't open file: "<< file_name << endl;
			exit(1);
		}

		return true;
	}

	bool Clear()
	{
		ifs.close();

		prev_tag_name.clear();

		prev_line.clear();

		tag_name.clear();

		prev_read_id = 0;

		min_ins = 0;

		return true;
	}

	size_t ReadNextTagAlign(vector<value_type>& stored_tag_aligns, size_t& read_id)
	{
		stored_tag_aligns.clear();

		tag_name = "";

		if (!prev_tag_name.empty())
		{
			tag_name = prev_tag_name;

			read_id = prev_read_id;

			stored_tag_aligns.push_back(value_type(prev_line, min_ins));
		}

		while (!ifs.eof() )
		{
			string line;

			getline(ifs, line);

			if (line.empty() || line[0] == '@')
				continue;

			char tag_name_chr[1000];

			sscanf(line.c_str(), "%s", tag_name_chr);

			size_t cur_read_id; 

			sscanf(tag_name_chr, "%llu", &cur_read_id);

			if (prev_tag_name.empty())
			{
				prev_tag_name = tag_name_chr;

				tag_name = prev_tag_name;

				read_id = cur_read_id;

				prev_read_id = cur_read_id;

				stored_tag_aligns.push_back(value_type(line, min_ins));
			}
			else if (prev_tag_name == tag_name_chr)
			{
				stored_tag_aligns.push_back(value_type(line, min_ins));
			}
			else
			{
				prev_tag_name = tag_name_chr;

				prev_line = line;

				prev_read_id = cur_read_id;

				return stored_tag_aligns.size();
			}
		}

		prev_tag_name.clear();

		prev_line.clear();

		return stored_tag_aligns.size();
	}

	size_t ReadNextTagAlignPE(pair<vector<value_type>, vector<value_type> >& stored_tag_aligns_pe, size_t& read_id)
	{
		stored_tag_aligns_pe.first.clear();

		stored_tag_aligns_pe.second.clear();

		tag_name = "";

		if (!prev_tag_name.empty())
		{
			tag_name = prev_tag_name;

			size_t pairid = prev_end_id;

			//sscanf(tag_name.c_str(), "%llu", &read_id);

			if (pairid % 2 == 0)
				stored_tag_aligns_pe.second.push_back(value_type(prev_line, min_ins));
			else
				stored_tag_aligns_pe.first.push_back(value_type(prev_line, min_ins));
		}

		while (!ifs.eof() )
		{
			string line;

			getline(ifs, line);

			if (line.empty() || line[0] == '@')
				continue;

			char tag_name_chr[1000];

			sscanf(line.c_str(), "%s", tag_name_chr);

			size_t cur_read_id; 

			sscanf(tag_name_chr, "%llu", &cur_read_id);

			string tag_name_chr_str = tag_name_chr;

			int pairid = atoi(tag_name_chr_str.substr(tag_name_chr_str.length() - 1, 1).c_str());

			size_t st_idx = tag_name_chr_str.find("~");			

			if (st_idx != string::npos)
				tag_name_chr_str = tag_name_chr_str.substr(st_idx + 1, tag_name_chr_str.length() - 3 - st_idx);	
			else
				tag_name_chr_str = tag_name_chr_str.substr(0, tag_name_chr_str.length() - 1);	

			if (prev_tag_name.empty())
			{
				prev_tag_name = tag_name_chr_str;

				tag_name = prev_tag_name;

				read_id = cur_read_id;

				prev_read_id = cur_read_id;

				prev_end_id =  pairid;

				if (pairid % 2 == 0)
					stored_tag_aligns_pe.second.push_back(value_type(line, min_ins));
				else
					stored_tag_aligns_pe.first.push_back(value_type(line, min_ins));
			}
			else if (tag_name_chr_str == prev_tag_name )
			{
				if (pairid % 2 == 0)
					stored_tag_aligns_pe.second.push_back(value_type(line, min_ins));
				else
					stored_tag_aligns_pe.first.push_back(value_type(line, min_ins));
			}
			else
			{
				prev_tag_name = tag_name_chr_str;

				prev_line = line;

				prev_read_id = cur_read_id;

				prev_end_id =  pairid;

				return stored_tag_aligns_pe.first.size() + stored_tag_aligns_pe.second.size();
			}
		}

		prev_tag_name.clear();

		prev_line.clear();

		return stored_tag_aligns_pe.first.size() + stored_tag_aligns_pe.second.size();
	}
};

#endif