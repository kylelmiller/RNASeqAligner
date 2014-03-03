#include <algorithm>
#include <math.h>
#include <string>
#include <stdio.h>
#include "sam_info.h"

class Local_Align
{

public:
	Local_Align()
	{
		q = 40;
		r = 5;
		k = 10;
		d = q + k * r;
		match = 20;
		mismatch = -40;
		canonical_score = 10;
		non_canonical_score = -40;
		S = 0;
		ED = 0;
		ID = 0;
		I = 0;
	}

	~Local_Align()
	{
		reset();
	}

	void set(string& text, string& pattern, const char* start_seq, const char* end_seq, bool start_penalty, bool end_penalty, bool _debug)
	{
		reset();
		text_ptr = &text;
		pattern_ptr = &pattern;
		text_size = (int)(text.length());
		pattern_size = (int)(pattern.length());
		S = new int[(text_size + 1) * (pattern_size + 1)];
		ED = new int[(text_size + 1) * (pattern_size + 1)];
		ID = new int[(text_size + 1) * (pattern_size + 1)];
		I = new int[(text_size + 1) * (pattern_size + 1)];
		canonical_start_seq = start_seq;
		canonical_end_seq = end_seq;
		start_intron_penalty = start_penalty;
		end_intron_penalty = end_penalty;
		debug = _debug;
	}

	void reset()
	{
		delete[] S;
		S = 0;
		delete[] ED;
		ED = 0;
		delete[] ID;
		ED = 0;
		delete[] I;
		I = 0;
	}

	int donor_score(int pos)
	{
		string flank_seq = text_ptr->substr(pos - 1, 2);
		if(flank_seq == canonical_start_seq)
			return canonical_score;
		else
			return non_canonical_score;
	}

	int acceptor_score(int pos)
	{
		if(pos < 2)
			return -100;
		string flank_seq = text_ptr->substr(pos - 2, 2);
		if(flank_seq == canonical_end_seq)
			return canonical_score;
		else
			return non_canonical_score;
	}

	int substitution_score(char a, char b)
	{
		if(a == b)
			return match;
		else
			return mismatch;
	}

	void build_matrices()
	{
		if(debug)
			cout << "building matrices" << endl;
		for(int i = 0; i <= text_size; i ++)
		{
			for(int j = 0; j <= pattern_size; j++)
			{
				if(i == 0 && j == 0)
				{
					S[i * (pattern_size + 1) + j] = 0;
					ED[i * (pattern_size + 1) + j] = (start_intron_penalty ? -q : 0);
					ID[i * (pattern_size + 1) + j] = (start_intron_penalty ? (donor_score(1) - d) : 0);
					I[i * (pattern_size + 1) + j] = -q;
				}
				else if(i == 0)
				{
					I[i * (pattern_size + 1) + j] = I[i * (pattern_size + 1) + j - 1] - r;
					S[i * (pattern_size + 1) + j] = I[i * (pattern_size + 1) + j];
					ED[i * (pattern_size + 1) + j] = S[i * (pattern_size + 1) + j] - ((j == pattern_size && !end_intron_penalty) ? 0 : q);
					ID[i * (pattern_size + 1) + j] = S[i * (pattern_size + 1) + j] + ((j == pattern_size && !end_intron_penalty) ? 0 : (donor_score(1) - d));
				}
				else if(j == 0)
				{
					ED[i * (pattern_size + 1) + j] = ED[(i - 1) * (pattern_size + 1) + j] - (start_intron_penalty ? r : 0);
					ID[i * (pattern_size + 1) + j] = ID[(i - 1) * (pattern_size + 1) + j];   
					S[i * (pattern_size + 1) + j] = max(ED[i * (pattern_size + 1) + j], ID[i * (pattern_size + 1) + j] + (start_intron_penalty ? acceptor_score(i) : 0));
					I[i * (pattern_size + 1) + j] = S[i * (pattern_size + 1) + j] - q;
				}
				else
				{
					I[i * (pattern_size + 1) + j] = max(I[i * (pattern_size + 1) + j - 1] - r, S[i * (pattern_size + 1) + j - 1] - q - r);	
					ED[i * (pattern_size + 1) + j] = max(ED[(i - 1) * (pattern_size + 1) + j] - ((j == pattern_size && !end_intron_penalty) ? 0 : r), S[(i - 1) * (pattern_size + 1) + j] -((j == pattern_size && !end_intron_penalty) ? 0 : (q + r)));
					ID[i * (pattern_size + 1) + j] = max(ID[(i - 1) * (pattern_size + 1) + j], S[(i - 1) * (pattern_size + 1) + j] + ((j == pattern_size && !end_intron_penalty) ? 0 : (donor_score(i) - d)));
					int max_tmp1 = max(S[(i - 1) * (pattern_size + 1) + j - 1] + substitution_score((*text_ptr)[i - 1], (*pattern_ptr)[j - 1]), ED[i * (pattern_size + 1) + j]);
					int max_tmp2 = max(I[i * (pattern_size + 1) + j], ID[i * (pattern_size + 1) + j] + ((j == pattern_size && !end_intron_penalty) ? 0 : acceptor_score(i)));
					S[i * (pattern_size + 1) + j] = max(max_tmp1, max_tmp2);
				}
			}
		}
		if(debug)
			cout << "building matrices done" << endl;
	}

	int trace_back(Sam_Info& my_sam, int& num_mismatch, int start_pos)
	{
		if(debug)
			cout << "trace back" << endl;
		int* current_matrix = S;
		int i = text_size, j = pattern_size;
		string text_align;
		string pattern_align;
		while(true)
		{
			if(current_matrix == S)
			{
				if(i == 0 && j == 0)
					break;
				else if(i > 0 && S[i * (pattern_size + 1) + j] == ED[i * (pattern_size + 1) + j])
					current_matrix = ED;
				else if(j > 0 && S[i * (pattern_size + 1) + j] == I[i * (pattern_size + 1) + j])
					current_matrix = I;
				else if(i > 0 && S[i * (pattern_size + 1) + j] == ID[i * (pattern_size + 1) + j] + ((( j == 0 && !start_intron_penalty) || (j == pattern_size && !end_intron_penalty)) ? 0 : acceptor_score(i)) )
					current_matrix = ID;
				else   // mismatch
				{
					text_align.insert(0, (*text_ptr), i - 1, 1);
					pattern_align.insert(0, (*pattern_ptr), j - 1, 1);
					i--;
					j--;
				}
			}
			else if(current_matrix == ID)
			{
				text_align.insert(0, (*text_ptr), i - 1, 1);   //deletion
				pattern_align.insert(0, "-");
				if(i > 1 && ID[i * (pattern_size + 1) + j] == ID[(i - 1) * (pattern_size + 1) + j])
					i--;
				else
				{
					current_matrix = S;
					i--;
				}
			}
			else if(current_matrix == ED)
			{
				text_align.insert(0, (*text_ptr), i - 1, 1);   //deletion
				pattern_align.insert(0, "-");
				if(i > 1 && ED[i * (pattern_size + 1) + j] == ED[(i - 1) * (pattern_size + 1) + j] - ( ((j == 0 && !start_intron_penalty) || (j == pattern_size && !end_intron_penalty)) ? 0 : r))
					i--;
				else
				{
					current_matrix = S;
					i--;
				}
			}
			else if(current_matrix == I)
			{
				text_align.insert(0, "-");     // insertion
				pattern_align.insert(0, (*pattern_ptr), j - 1, 1);
				if(j > 1 && I[i * (pattern_size + 1) + j] == I[i * (pattern_size + 1) + j - 1] - r)
					j--;
				else
				{
					current_matrix = S;
					j--;
				}
			}
		}
		if(debug)
			cout << "trace back done" << endl;
		get_cigar_code(text_align, pattern_align, my_sam, num_mismatch, start_pos);
		return S[(text_size + 1) * (pattern_size + 1) - 1];
	}

	void get_cigar_code(string& text_align, string& pattern_align, Sam_Info& my_sam, int& num_mismatch, int start_pos)
	{
		if(debug)
		{
			cout << "get jump code" << endl;
			cout << "text align is " << text_align << endl;
			cout << "pattern align is " << pattern_align << endl;
		}
		my_sam.start_pos = -1;
		int span_len = 0;
		string span_type = "";
		num_mismatch = 0;
		for(size_t i = 0; i < text_align.length(); i++)
		{
			if(text_align[i] != '-' && pattern_align[i] == '-')
			{
				if(span_type != "")
				{
					if(span_type != "N")
					{
						Jump_Code new_jumpcode(span_len, span_type);
						my_sam.jump_code.push_back(new_jumpcode);
						span_len = 1;
						span_type = "N";
					}
					else
						span_len ++;
				}
				else if(my_sam.start_pos != -1)
				{
					span_type = "N";
					span_len  = 1;
				}
			}
			else
			{
				if(my_sam.start_pos == -1)
				{
					if(start_intron_penalty && i != 0)
					{
						my_sam.start_pos = start_pos;
						Jump_Code new_jumpcode1(0, "M");
						Jump_Code new_jumpcode2((int)i, "N");
						my_sam.jump_code.push_back(new_jumpcode1);
						my_sam.jump_code.push_back(new_jumpcode2);
					}
					else
						my_sam.start_pos = (int)i + start_pos;
				}
				if(text_align[i] == '-' && pattern_align[i] != '-')
				{
					if(span_type != "")
					{
						if(span_type != "I")
						{
							Jump_Code new_jumpcode2(span_len, span_type);
							my_sam.jump_code.push_back(new_jumpcode2);
							span_len = 1;
							span_type = "I";
						}
						else
							span_len ++;
					}
					else
					{
						span_len = 1;
						span_type = "I";
					}
				}
				else if(text_align[i] != '-' && pattern_align[i] != '-')
				{
					if(text_align[i] != pattern_align[i])
					{
						num_mismatch ++;
					}
					if(span_type != "")
					{
						if(span_type != "M")
						{
							Jump_Code new_jumpcode(span_len, span_type);
							my_sam.jump_code.push_back(new_jumpcode);
							span_len = 1;
							span_type = "M";
						}
						else
							span_len ++;
					}
					else
					{
						span_len = 1;
						span_type = "M";
					}
				}
			}
		}
		if(span_type != "N")
		{
			Jump_Code new_jumpcode1(span_len, span_type);
			my_sam.jump_code.push_back(new_jumpcode1);
		}
		else if(end_intron_penalty)
		{
			Jump_Code new_jumpcode1(span_len, span_type);
			Jump_Code new_jumpcode2(0, "M");
			my_sam.jump_code.push_back(new_jumpcode1);
			my_sam.jump_code.push_back(new_jumpcode2);
		}
		if(debug)
			cout << "get jump code done" << endl;
	}
	
	int align(string& text, string& pattern, const char* start_seq, const char* end_seq, bool start_penalty, bool end_penalty, Sam_Info& my_sam, int& num_mismatch, int start_pos, bool _debug)
	{
		set(text, pattern, start_seq, end_seq, start_penalty, end_penalty, _debug);
		build_matrices();
		return trace_back(my_sam, num_mismatch, start_pos);
	}

private:
	int* S;
	int* ED;
	int* ID;
	int* I;
	string* pattern_ptr;
	string* text_ptr;
	int text_size;
	int pattern_size;
	string canonical_start_seq;
	string canonical_end_seq;
	int q;     // deletion gap open penalty
	int r;     // deletion gap extension penalty
	int k;     // intron / deletion cutoff
	int d;     // intron gap penalty
	int match;  // match bonus
	int mismatch; // mismatch penalty
	int canonical_score;  //canonical bonus	
	int non_canonical_score; // non-canonical penalty
	bool start_intron_penalty;
	bool end_intron_penalty;
	bool debug;
};








