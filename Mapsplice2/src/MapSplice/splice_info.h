#ifndef SPLICE_INFO_H
#define SPLICE_INFO_H

#include <string>
#include "mismatch.h"
using namespace std;

class Bwtmap_Info;

class Jump_Code
{
public:
	int len;
	string type;

	Jump_Code(int _len, string _type)
	{
		len=_len;
		type=_type;
	}
	
	~Jump_Code()
	{}
	
	string toString()
	{
		return int_to_str(len) + type;
	}
	
	string toString(int add_value)
	{
		return int_to_str(len + add_value) + type;
	}
	
	string toStringConvertDeletion(int max_del_len)
	{
		if(type == "N" && len <= max_del_len)
		{
			return int_to_str(len) + "D";
		}
		return int_to_str(len) + type;
	}
};

class Splice_Info
{
public:
	Splice_Info()
	{ 
		mapped_len = -1;	
	}

	Splice_Info(const Splice_Info& copy_info)
	{
		start_pos = copy_info.start_pos;
		end_pos = copy_info.end_pos;
		start_contig = copy_info.start_contig;
		end_contig = copy_info.end_contig;
		chrom = copy_info.chrom;
		strand = copy_info.strand;
		buffer_len = copy_info.buffer_len;
		mapped_len = copy_info.mapped_len;
		start_seg_no = copy_info.start_seg_no;
		end_seg_no = copy_info.end_seg_no;
		for(size_t i = 0; i< copy_info.jump_code.size(); i++)
		{
			jump_code.push_back(copy_info.jump_code[i]);
		}
		for(size_t i = 0; i< copy_info.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(copy_info.junc_flank_seq[i]);
		}
	}

	void copy(const Splice_Info& copy_info)
	{
		start_pos = copy_info.start_pos;
		end_pos = copy_info.end_pos;
		start_contig = copy_info.start_contig;
		end_contig = copy_info.end_contig;
		chrom = copy_info.chrom;
		strand = copy_info.strand;
		buffer_len = copy_info.buffer_len;
		mapped_len = copy_info.mapped_len;
		start_seg_no = copy_info.start_seg_no;
		end_seg_no = copy_info.end_seg_no;
		for(size_t i = 0; i< copy_info.jump_code.size(); i++)
		{
			jump_code.push_back(copy_info.jump_code[i]);
		}
		for(size_t i = 0; i< copy_info.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(copy_info.junc_flank_seq[i]);
		}
	}
	
	void Get_Extra_Info(string& chrom_seq)
	{
		mapped_len = 0;
		int	end_tmp = start_pos - 1;
		for(size_t i = 0; i < jump_code.size(); i++)
		{
			if(jump_code[i].type == "M")
			{
				end_tmp += jump_code[i].len;
				mapped_len += jump_code[i].len;
			}
			else if(jump_code[i].type == "N")
			{
				int junc_start = end_tmp + 1;
				end_tmp += jump_code[i].len;
				int junc_end =  end_tmp;
				string flank_seq = chrom_seq.substr(junc_start, 2) + chrom_seq.substr(junc_end - 1, 2);
				junc_flank_seq.push_back(flank_seq);
			}
			else   // "I"
			{
				mapped_len += jump_code[i].len;
			}
		}	
	}
	
	void set_value(string line, bool pairend)
	{
	}

	void set_value_single(string line)
	{
	}

	void set_value_pairend(string line)
	{
	}

	~Splice_Info()
	{
	}

	int GetMapLen()
	{
	  if(mapped_len == -1)
	  {
	  	mapped_len = 0;
	  	for(size_t i = 0; i < jump_code.size(); i++)
	  	{
	  		if(jump_code[i].type == "M" || jump_code[i].type == "I")
	  			mapped_len += jump_code[i].len;
	  	}
		}
	  return mapped_len;
	}

	int Count_Gap()
	{
		int num_gap = 0;
		for(size_t i = 0; i < jump_code.size(); i++)
	  {
	  	if(jump_code[i].type != "M")
	  		num_gap ++;
	  }
	  return num_gap;
	}

	vector<Jump_Code> jump_code;
	vector<string> junc_flank_seq;
	string chrom;
	string strand;
	int start_pos;
	int end_pos;
	int buffer_len;
	int mapped_len;
	int start_seg_no;
	int end_seg_no;
	size_t start_contig;
	size_t end_contig;
};


class Fusion_Splice
{
public:
	Splice_Info first_splice;
	Splice_Info second_splice;
	string flank_seq;

	Fusion_Splice()
	{}
	
	Fusion_Splice (Splice_Info& my_first_splice, Splice_Info& my_second_splice, string& my_flank_seq)
	{
	 	first_splice = my_first_splice;
		second_splice = my_second_splice;
		flank_seq = my_flank_seq; 
	}
	
	Fusion_Splice (const Fusion_Splice& fusion_splice_info)
	{
		first_splice = fusion_splice_info.first_splice;
		second_splice = fusion_splice_info.second_splice;
		flank_seq = fusion_splice_info.flank_seq;
	}
	
	~Fusion_Splice ()
	{}
};


#endif

