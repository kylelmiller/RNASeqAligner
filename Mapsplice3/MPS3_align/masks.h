#ifndef __MASKS_H_INCLUDED__
#define __MASKS_H_INCLUDED__

#include "constantDefinitions.h"

using namespace std;

struct Masks{
	/*
	 * Member Variables
	 */
	size_t comp_buff_width;
	size_t comb_seg_bits_on;
	size_t suffix_seg_bits_on;
	size_t prefix_seg_bits_on;
	size_t comb_seg_first_selector_rt;
	size_t comb_seg_first_selector_lt;
	size_t score_seg_buf_width;

	/*
	 * Constructor
	 */
	Masks (size_t seg_len)
	{
		//duplicate reads
		comp_buff_width = 2 * seg_len;
		comb_seg_bits_on = (ALL_BITS_ON >> (BITS_SUPPORTED - comp_buff_width));
		suffix_seg_bits_on = (ALL_BITS_ON >> (BITS_SUPPORTED - seg_len));
		prefix_seg_bits_on = suffix_seg_bits_on << seg_len;
		comb_seg_first_selector_rt = LEAST_SIG_BIT;
		comb_seg_first_selector_lt = LEAST_SIG_BIT << seg_len;
		score_seg_buf_width = seg_len;
	}
};

#endif
