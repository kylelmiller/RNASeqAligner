#ifndef SPLICE_INFO_H
#define SPLICE_INFO_H

#include <string>
#include <string.h>
#include <stdexcept>

#include "utilities.h"

using namespace std;

class Jump_Code
{
public:
	int _length; // length fo current feature
	string _type; //feature type, M/N/I/D/S

	Jump_Code()
	{
		_length = 0;
		_type = "";
	}

	Jump_Code(int length, string type)
	{
		_length=length;
		_type=type;
	}
	
	~Jump_Code()
	{}

	string toString()
	{
		return Utilities::int_to_str(_length) + _type;
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

	// merges up all of the jump codes into the
	// cigar string which will be written to the sam file
	void getFinalJumpCode()
	{
		switch(jump_code.size())
		{
			case 0:
				throw invalid_argument("Bad Jump Code!");
				break;
			case 1:
				if(jump_code.back()._length <= 0 || jump_code.back()._type != "M")
					throw invalid_argument("Bad Jump Code!");
				else
					final_jump_code.push_back(jump_code.back());
				break;
			default:
				Jump_Code currentJumpCode;
				Jump_Code nextJumpCode;
				for(int i=1; i<jump_code.size(); i++)
				{
					currentJumpCode = jump_code[i-1];
					nextJumpCode = jump_code[i];

					// if the codes are different (say S and M) then
					// just add the current jump code and move on
					if(currentJumpCode._type != nextJumpCode._type)
						final_jump_code.push_back(currentJumpCode);
					else  // Different, then add lengths and keep type
						nextJumpCode._length += currentJumpCode._length;
				}
				final_jump_code.push_back(nextJumpCode);
				break;
		}
	}

	int getEndBaseMapPos_jump_code(int startMapPos)
	{
		int tmpMapPos = 0;
		int jumpCodeSize = jump_code.size();
		for(int tmp = 0; tmp < jumpCodeSize; tmp++)
		{
			if((jump_code[tmp]._type == "S")||(jump_code[tmp]._type == "I"))
			{
				continue;
			}
			else if((jump_code[tmp]._type == "M")||(jump_code[tmp]._type == "N")||(jump_code[tmp]._type == "D"))
			{
				tmpMapPos += jump_code[tmp]._length;
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
			if(final_jump_code[i]._length <= 0)
				return false;

		return true;
	}
};

#endif

