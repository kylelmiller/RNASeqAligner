#ifndef __GENOME_SCAN_H_INCLUDED__
#define __GENOME_SCAN_H_INCLUDED__

#include <vector>
#include <string>

#include "masks.h"
#include "wordPair.h"

using namespace std;

class GenomeScan{
public:

	bool Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, string& flank_seq, size_t& msimatch_num);

	bool Double_anchored_score_least_mis(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, size_t& mismatch_num);

	bool Double_anchored_score_ins(string tobe_fixed_str, string chrom_seq, size_t max_mismatch, size_t& prefix_length, size_t& comb_bits, size_t& num_mismatch);

private:
	/*
	 * Constants
	 */
	static const size_t bit_GT_upper = 3;
	static const size_t bit_GT_lower = 1;

	static const size_t bit_TG_upper = 3;
	static const size_t bit_TG_lower = 2;

	static const size_t bit_AG_upper = 1;
	static const size_t bit_AG_lower = 0;

	static const size_t bit_GA_upper = 2;
	static const size_t bit_GA_lower = 0;

	static const size_t bit_GC_upper = 2;
	static const size_t bit_GC_lower = 1;

	static const size_t bit_CG_upper = 1;
	static const size_t bit_CG_lower = 2;

	static const size_t bit_AT_upper = 1;
	static const size_t bit_AT_lower = 1;

	static const size_t bit_TA_upper = 2;
	static const size_t bit_TA_lower = 2;

	static const size_t bit_AC_upper = 0;
	static const size_t bit_AC_lower = 1;

	static const size_t bit_CA_upper = 0;
	static const size_t bit_CA_lower = 2;

	static const size_t bit_CT_upper = 1;
	static const size_t bit_CT_lower = 3;

	static const size_t bit_TC_upper = 2;
	static const size_t bit_TC_lower = 3;

	static const size_t bit_GTAG_upper = bit_GT_upper << 2 | bit_AG_upper;
	static const size_t bit_GTAG_lower = bit_GT_lower << 2 | bit_AG_lower;
	static const size_t bit_GTAG = bit_GTAG_upper << 4 | bit_GTAG_lower;

	static const size_t bit_GATG_upper = bit_GA_upper << 2 | bit_TG_upper;
	static const size_t bit_GATG_lower = bit_GA_lower << 2 | bit_TG_lower;
	static const size_t bit_GATG = bit_GATG_upper << 4 | bit_GATG_lower;

	static const size_t bit_GCAG_upper = bit_GC_upper << 2 | bit_AG_upper;
	static const size_t bit_GCAG_lower = bit_GC_lower << 2 | bit_AG_lower;
	static const size_t bit_GCAG = bit_GCAG_upper << 4 | bit_GCAG_lower;

	static const size_t bit_GACG_upper = bit_GA_upper << 2 | bit_CG_upper;
	static const size_t bit_GACG_lower = bit_GA_lower << 2 | bit_CG_lower;
	static const size_t bit_GACG = bit_GACG_upper << 4 | bit_GACG_lower;

	static const size_t bit_ATAC_upper = bit_AT_upper << 2 | bit_AC_upper;
	static const size_t bit_ATAC_lower = bit_AT_lower << 2 | bit_AC_lower;
	static const size_t bit_ATAC = bit_ATAC_upper << 4 | bit_ATAC_lower;

	static const size_t bit_CATA_upper = bit_CA_upper << 2 | bit_TA_upper;
	static const size_t bit_CATA_lower = bit_CA_lower << 2 | bit_TA_lower;
	static const size_t bit_CATA = bit_CATA_upper << 4 | bit_CATA_lower;

	static const size_t bit_CTAC_upper = bit_CT_upper << 2 | bit_AC_upper;
	static const size_t bit_CTAC_lower = bit_CT_lower << 2 | bit_AC_lower;
	static const size_t bit_CTAC = bit_CTAC_upper << 4 | bit_CTAC_lower;

	static const size_t bit_CATC_upper = bit_CA_upper << 2 | bit_TC_upper;
	static const size_t bit_CATC_lower = bit_CA_lower << 2 | bit_TC_lower;
	static const size_t bit_CATC = bit_CATC_upper << 4 | bit_CATC_lower;

	static const size_t bit_CTGC_upper = bit_CT_upper << 2 | bit_GC_upper;
	static const size_t bit_CTGC_lower = bit_CT_lower << 2 | bit_GC_lower;
	static const size_t bit_CTGC = bit_CTGC_upper << 4 | bit_CTGC_lower;

	static const size_t bit_CGTC_upper = bit_CG_upper << 2 | bit_TC_upper;
	static const size_t bit_CGTC_lower = bit_CG_lower << 2 | bit_TC_lower;
	static const size_t bit_CGTC = bit_CGTC_upper << 4 | bit_CGTC_lower;

	static const size_t bit_GTAT_upper = bit_GT_upper << 2 | bit_AT_upper;
	static const size_t bit_GTAT_lower = bit_GT_lower << 2 | bit_AT_lower;
	static const size_t bit_GTAT = bit_GTAT_upper << 4 | bit_GTAT_lower;

	static const size_t bit_TATG_upper = bit_TA_upper << 2 | bit_TG_upper;
	static const size_t bit_TATG_lower = bit_TA_lower << 2 | bit_TG_lower;
	static const size_t bit_TATG = bit_TATG_upper << 4 | bit_TATG_lower;

	/*
	 * Member Variables
	 */
	WordPair m_five_prim_suffix, m_three_prim_prefix;

	size_t m_matched_flank, m_matched_bads;

	/*
	 * Methods
	 */
	bool FixHoleCheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	bool FixHoleCheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width);

	size_t Fixhole_score_selective_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	size_t Fixhole_score_selective_insert_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr);

	string FlankString(size_t matched_flank, bool bad);

};

inline bool
GenomeScan::FixHoleCheckFirstTime(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	mins = s;
	loc = score_buf_width;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	m_matched_bads = combine_bads;

	if (combine_bads == 0)
	{
		combine_words = ((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask);

		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 6;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				m_matched_flank = combine_words;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 2;
				m_matched_flank = combine_words;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			{
				m_matched_flank = combine_words;
			}
			break;
		}
	}

	return false;
}

inline bool
GenomeScan::FixHoleCheckBeforeMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;
	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_ATAC:
			{
				prim = 1;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				prim = 3;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				prim = 4;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				prim = 2;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			if (s < mins)
			{
				mins = s;
				loc = score_buf_width - i;
				m_matched_bads = combine_bads;
				m_matched_flank = combine_words;
			}
			break;
		}
	}
	else
	{
		if (s < mins)
		{
			mins = s;
			loc = score_buf_width - i;
			m_matched_bads = combine_bads;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;
		switch (combine_words)
		{
		case bit_ATAC:
			{
				if ( s < mins)
				{
					prim = 1;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_CTAC:
			{
				//if ( s < mins)
				//{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				//}
				if (mins == 0)
					return true;
			}
			break;
		case bit_CTGC:
			{
				if ( s < mins)
				{
					prim = 3;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GCAG:
			{
				if ( s < mins)
				{
					prim = 4;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
			break;
		case bit_GTAT:
			{
				if ( s < mins)
				{
					prim = 2;
					mins = s;
					loc = score_buf_width - i;
					m_matched_flank = combine_words;
					m_matched_bads = combine_bads;
				}
				//if (mins == 0)
				//	return true;
			}
			break;
		default:
			break;
		}
	}
	return false;
}

inline bool
GenomeScan::FixHoleCheckAfterGTAGMatch(const size_t& five_prim_mask, const size_t& three_prim_mask, const size_t& i, size_t& prim, size_t& mins, size_t& loc, const size_t& s, const size_t& score_buf_width)
{
	size_t combine_words, combine_bads;

	combine_bads = (m_five_prim_suffix.bads & five_prim_mask) | (m_three_prim_prefix.bads & three_prim_mask);

	if (combine_bads == 0)
	{
		combine_words = (((m_five_prim_suffix.upper & five_prim_mask) | (m_three_prim_prefix.upper & three_prim_mask)) << 4
			| (m_five_prim_suffix.lower & five_prim_mask) | (m_three_prim_prefix.lower & three_prim_mask)) >> i;

		switch (combine_words)
		{
		case bit_CTAC:
			{
				prim = 6;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		case bit_GTAG:
			{
				prim = 5;
				mins = s;
				loc = score_buf_width - i;
				m_matched_flank = combine_words;
				m_matched_bads = combine_bads;
				if (mins == 0)
					return true;
			}
		}
	}
	return false;
}

inline size_t
GenomeScan::Fixhole_score_selective_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl;

	size_t five_prim_mask = LAST_THIRD_FOUTH;
	size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;
	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) |
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
		return mins;

	five_prim_mask <<= 1;
	three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off
	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on
	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off

	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;
	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;

		if (s <= left_mismatch)
		{
			if (!prim)
			{
				if (FixHoleCheckBeforeMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (prim < 5)
			{
				if (FixHoleCheckAfterMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}
			else if (s < mins)
			{
				if (FixHoleCheckAfterGTAGMatch(five_prim_mask, three_prim_mask, i, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
					return mins;
			}

		}
		selector1 <<= 1;
		selector2 <<= 1;
		five_prim_mask <<= 1;
		three_prim_mask <<= 1;
	}

	return mins;
}

inline size_t
GenomeScan::Fixhole_score_selective_insert_var_mask(const WordPair & read_word_dup, const WordPair & comb_chrom_seq, size_t& loc, size_t& prim, const size_t& left_mismatch, size_t& rbits, const Masks* mask_ptr)
{
	//get the score with LPrefix + Rbuf + Rsuffix
	size_t s = read_word_dup.score(comb_chrom_seq, mask_ptr->prefix_seg_bits_on);

	//cout << "prefix score: "<<s << endl;

	//size_t five_prim_mask = LAST_THIRD_FOUTH;
	//size_t three_prim_mask = LAST_TWO_BIT;

	size_t mins = left_mismatch + 1;

	prim = 0;

	//generate the bits where 1 indicates the mismatches.
	register size_t bits = ((comb_chrom_seq.upper ^ read_word_dup.upper) |
		(comb_chrom_seq.lower ^ read_word_dup.lower) | comb_chrom_seq.bads | read_word_dup.bads) & mask_ptr->comb_seg_bits_on;

	rbits = bits;

	if (s <= left_mismatch)// && FixHoleCheckFirstTime(five_prim_mask, three_prim_mask, 0, prim, mins, loc, s, mask_ptr->score_seg_buf_width))
	{
		mins = s;

		loc = mask_ptr->score_seg_buf_width;

		if (!mins)
			return mins;
	}

	//five_prim_mask <<= 1;
	//three_prim_mask <<= 1;

	//the following are two pointers each point to the bit to be turned on and turned off

	                                              // - m_max_mismatch
	size_t selector1 = mask_ptr->comb_seg_first_selector_rt;//LEAST_SIG_BIT << (m_seed_width - m_max_mismatch); //to be turned on

	                                             // + m_max_mismatch
	size_t selector2 = mask_ptr->comb_seg_first_selector_lt;//to be turned off

	//check out the accumulation score of mismatches when the selectors are moving...
	//size_t score_buff_width = m_//m_masks.score_buff_width;

	                                     // + m_max_mismatch
	for (size_t i = 1; i <= mask_ptr->score_seg_buf_width/*score_buff_width*/; ++i)
	{
		s += (selector1 & bits ) ? 1 : 0;
		s -= (selector2 & bits ) ? 1 : 0;

		if (mins > s)
		{
			mins = s;

			loc = mask_ptr->score_seg_buf_width - i;
			if (!mins)
				return mins;
		}

		selector1 <<= 1;
		selector2 <<= 1;
	}

	return mins;
}

inline string
GenomeScan::FlankString(size_t matched_flank, bool bad)
{
	string flankstr = "";
	if (bad)
		flankstr = "BADS";
	else
	{
		for (size_t i = 7; i >= 4; --i)
		{
			if (matched_flank & (LEAST_SIG_BIT << i))
			{
				if (matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "T";
				else
					flankstr += "G";
			}
			else
			{
				if (matched_flank & (LEAST_SIG_BIT << (i - 4)))
					flankstr += "C";
				else
					flankstr += "A";
			}
		}
	}

	return flankstr;
}

bool
GenomeScan::Double_anchored_score(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, string& flank_seq, size_t& msimatch_num)
{
	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);

	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;

	//if(do_noncanonical)

	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);

	//else

	score = Fixhole_score_selective_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);

	flank_seq = FlankString(m_matched_flank, m_matched_bads > 0);

	if(!do_noncanonical)
	{
		if(flank_seq != "GTAG" && flank_seq != "CTAC" && flank_seq != "ATAC" && flank_seq != "GTAT" && flank_seq != "CTGC" && flank_seq != "GCAG")
			return false;
	}

	if(score <= max_mismatch)
	{
		msimatch_num = score;
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}
		return true;
	}
	else
		return false;
}

bool      // select longer prefix
GenomeScan::Double_anchored_score_least_mis(string tobe_fixed_str, string doner_str, string acceptor_str, size_t& prefix_length, size_t max_mismatch, size_t& comb_bits, bool do_noncanonical, size_t& mismatch_num)
{

	size_t tobe_fixed_len = tobe_fixed_str.length();

	Masks mask(tobe_fixed_len);

	WordPair to_be_fixed_wp(tobe_fixed_str);

	WordPair pre_wp;

	to_be_fixed_wp.duplicate_self(tobe_fixed_len, pre_wp);

	m_five_prim_suffix = WordPair(doner_str);

	m_five_prim_suffix.left_shift(2);

	m_three_prim_prefix = WordPair(acceptor_str);

	WordPair prefix_wp(doner_str);

	prefix_wp.left_shift(tobe_fixed_len - 2);

	WordPair comb_chrom_seq;

	prefix_wp.ps_combine(mask.prefix_seg_bits_on, mask.suffix_seg_bits_on, mask.comb_seg_bits_on, m_three_prim_prefix, comb_chrom_seq);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = max_mismatch;

	size_t score;

	//if(do_noncanonical)

	//	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);

	//else

	score = Fixhole_score_selective_insert_var_mask(pre_wp, comb_chrom_seq, max_loc, prim, left_mismatches, rbits, &mask);

	if(!do_noncanonical)
	{
		string flank_string = FlankString(m_matched_flank, m_matched_bads > 0);
		if(flank_string != "GTAG" && flank_string != "CTAC" && flank_string != "ATAC" && flank_string != "GTAT" && flank_string != "CTGC" && flank_string != "GCAG")
			return false;
	}

	if(score <= max_mismatch)
	{
		prefix_length = max_loc;
		if (score)
		{
			mismatch_num = score;
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			//string comb_chrom_str1 = doner_str.substr(0, max_loc) + acceptor_str.substr(acceptor_str.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}
		return true;
	}
	else
		return false;
}

bool
GenomeScan::Double_anchored_score_ins(string tobe_fixed_str, string chrom_seq, size_t max_mismatch, size_t& prefix_length, size_t& comb_bits, size_t& num_mismatch)
{
	size_t tobe_fixed_len = chrom_seq.length();

	Masks mask(tobe_fixed_len);

	string prefix_read = tobe_fixed_str.substr(0, tobe_fixed_len);

	string suffix_read = tobe_fixed_str.substr(tobe_fixed_str.length() - tobe_fixed_len, tobe_fixed_len);

	string comb_read = prefix_read + suffix_read;

	WordPair comb_read_wp(comb_read);

	WordPair dup_chrom_wp, chrom_wp(chrom_seq);

	chrom_wp.duplicate_self(tobe_fixed_len, dup_chrom_wp);

	size_t max_loc = 0, prim = 0;

	m_matched_flank = 0;

	m_matched_bads = 0;

	size_t rbits;

	size_t left_mismatches = 2;

	size_t score = Fixhole_score_selective_insert_var_mask(comb_read_wp, dup_chrom_wp, max_loc, prim, left_mismatches, rbits, &mask);

	vector<pair<size_t, pair<char, char> > > mis_pos;

	if(/*max_loc != 0 && max_loc != chrom_seq.length() && */score <= max_mismatch)
	{
		//cout << "score: " << score << endl;
		(num_mismatch) = score;
		prefix_length = max_loc;
		if (score)
		{
			size_t seg1_suffix_len = tobe_fixed_len - max_loc;

			string mapped_read_str = comb_read.substr(0, max_loc) + comb_read.substr(comb_read.length() - seg1_suffix_len, seg1_suffix_len);

			size_t seg1_mask_prefix = mask.suffix_seg_bits_on >> seg1_suffix_len << seg1_suffix_len;

			size_t seg1_mask_suffix = mask.suffix_seg_bits_on >> max_loc;

			comb_bits = ((rbits >> tobe_fixed_len) & seg1_mask_prefix) + (rbits & seg1_mask_suffix);
		}

		return true;
	}
	else
		return false;
}

#endif
