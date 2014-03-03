#include "ReadNextTagAlignHandler.h"

template<class T>
ReadNextTagAlignHandler<T>::ReadNextTagAlignHandler(const char* file_name) : prev_read_id(0)
{
	ifs.open(file_name);

	if (!ifs.is_open())
	{
		cerr<<"can't open file: "<< file_name << endl;
		exit(1);
	}
}

template<class T> 
bool 
ReadNextTagAlignHandler<T>::Init(const char* file_name)
{
	Clear();

	ifs.open(file_name);

	if (!ifs.is_open())
	{
		cerr<<"can't open file: "<< file_name << endl;
		exit(1);
	}
}

template<class T> 
bool 
ReadNextTagAlignHandler<T>::Clear()
{
	ifs.close();

	prev_tag_name.clear();

	prev_line.clear();

	tag_name.clear();

	prev_read_id = 0;

	return true;
}

template<class T> 
size_t
ReadNextTagAlignHandler<T>::ReadNextTagAlign(vector<value_type>& stored_tag_aligns, size_t& read_id)
{
	stored_tag_aligns.clear();

	tag_name = "";

	if (!prev_tag_name.empty())
	{
		tag_name = prev_tag_name;

		read_id = prev_read_id;

		stored_tag_aligns.push_back(value_type(prev_line));
	}

	while (!ifs.eof() )
	{
		string line;

		getline(ifs, line);

		if (line.empty())
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

			stored_tag_aligns.push_back(value_type(line));
		}
		else if (prev_read_id == cur_read_id)
		{
			stored_tag_aligns.push_back(value_type(line));
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

template<class T> 
size_t 
ReadNextTagAlignHandler<T>::ReadNextTagAlignPE(pair<vector<value_type>, vector<value_type> >& stored_tag_aligns_pe, size_t& read_id)
{
	stored_tag_aligns_pe.first.clear();

	stored_tag_aligns_pe.second.clear();

	tag_name = "";

	if (!prev_tag_name.empty())
	{
		tag_name = prev_tag_name;

		read_id = prev_read_id;

		//sscanf(tag_name.c_str(), "%llu", &read_id);

		if (read_id % 2 == 0)
			stored_tag_aligns_pe.second.push_back(value_type(prev_line));
		else
			stored_tag_aligns_pe.first.push_back(value_type(prev_line));
	}

	while (!ifs.eof() )
	{
		string line;

		getline(ifs, line);

		if (line.empty())
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

			if (cur_read_id % 2 == 0)
				stored_tag_aligns_pe.second.push_back(value_type(line));
			else
				stored_tag_aligns_pe.first.push_back(value_type(line));
		}
		else if (cur_read_id == prev_read_id || (cur_read_id == prev_read_id + 1 && cur_read_id % 2 == 0) )
		{
			if (cur_read_id % 2 == 0)
				stored_tag_aligns_pe.second.push_back(value_type(line));
			else
				stored_tag_aligns_pe.first.push_back(value_type(line));
		}
		else
		{
			prev_tag_name = tag_name_chr;

			prev_line = line;

			prev_read_id = cur_read_id;

			return stored_tag_aligns_pe.first.size() + stored_tag_aligns_pe.second.size();
		}
	}

	prev_tag_name.clear();

	prev_line.clear();

	return stored_tag_aligns_pe.first.size() + stored_tag_aligns_pe.second.size();
}