#ifndef SPLICE_INFO_H
#define SPLICE_INFO_H

#include <string>
#include <string.h>

#include "utilities.h"

using namespace std;

class Jump_Code
{
public:
	int len; // length fo current feature
	string type; //feature type, M/N/I/D/S

	Jump_Code(int _len, string _type)
	{
		len=_len;
		type=_type;
	}
	
	~Jump_Code()
	{}

	string toString()
	{
		return Utilities::int_to_str(len) + type;
	}
};

class Splice_Info
{
public:

	/*
	 * Member Variables
	 */
	vector<Jump_Code> jump_code;
	vector<Jump_Code> final_jump_code;
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

	Splice_Info()
	{ 
		mapped_len = -1;	
	}

	Splice_Info(const Splice_Info& other)
	{
		copy(other);
	}

	void copy(const Splice_Info& other)
	{
		start_pos = other.start_pos;
		end_pos = other.end_pos;
		start_contig = other.start_contig;
		end_contig = other.end_contig;
		chrom = other.chrom;
		strand = other.strand;
		buffer_len = other.buffer_len;
		mapped_len = other.mapped_len;
		start_seg_no = other.start_seg_no;
		end_seg_no = other.end_seg_no;
		for(size_t i = 0; i< other.jump_code.size(); i++)
		{
			jump_code.push_back(other.jump_code[i]);
		}
		for(size_t i = 0; i< other.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(other.junc_flank_seq[i]);
		}
	}

	void appendJumpCode(Splice_Info* anotherSpliceInfo)
	{
		for(int tmp = 0; tmp < anotherSpliceInfo->jump_code.size(); tmp++)
			jump_code.push_back(anotherSpliceInfo->jump_code[tmp]);
	}

	void getFinalJumpCode()
	{
		int jumpCodeSize = jump_code.size();
		if (jumpCodeSize == 0)
		{
			return;
		}
		else if(jumpCodeSize == 1)
		{
			if((jump_code[0].len <= 0)||(jump_code[0].type != "M"))
				return;
			else
			{	
				final_jump_code.push_back(jump_code[0]);
				return;
			}
		}
		else
		{	
			Jump_Code currentJumpCode = jump_code[0];
			Jump_Code nextJumpCode = jump_code[1];
			for(int tmp = 0; tmp < jump_code.size() - 1; tmp++)
			{
				int tmp2 = tmp+1;
				nextJumpCode = jump_code[tmp2];
				if(currentJumpCode.type != nextJumpCode.type)
				{	
					final_jump_code.push_back(currentJumpCode);
					currentJumpCode = nextJumpCode;
				}
				else
				{
					currentJumpCode.len = currentJumpCode.len + nextJumpCode.len;
				}
			}
			final_jump_code.push_back(currentJumpCode);
		}
	}

	int getEndBaseMapPos_jump_code(int startMapPos)
	{
		int tmpMapPos = 0;
		int jumpCodeSize = jump_code.size();
		for(int tmp = 0; tmp < jumpCodeSize; tmp++)
		{
			if((jump_code[tmp].type == "S")||(jump_code[tmp].type == "I"))
			{
				continue;
			}
			else if((jump_code[tmp].type == "M")||(jump_code[tmp].type == "N")||(jump_code[tmp].type == "D"))
			{
				tmpMapPos += jump_code[tmp].len;
			}
			else
			{
				cout << "invalid jumpCode type, error in splice_info.h: getEndBaseMapPos_jump_code !" << endl;
			}
		}
		int endBaseMapPos = startMapPos + tmpMapPos - 1;
		return endBaseMapPos;
	}

	bool allFinalJumpCodeValid()
	{
		for(int i = 0; i<final_jump_code.size(); i++)
			if(final_jump_code[i].len <= 0)
				return false;

		return true;
	}
};

#endif

