#ifndef SAMREC_H
#define SAMREC_H


#include "sharedlib.h"
#include "JunctionSeed.h"




class SamRec {

public:
	string tag_name;
	unsigned short strand_t;
	string chrom_name;
	size_t start, end;
	unsigned short confid;
	string splice_way;
	string ori_splice_way;
	string mapped_seq;
	string clippled_way;
	unsigned short mis_match;

	set<FILTERED_TYPE> filter_type;

	vector<JunctionSeed*> corresponding_juncs;

	vector<pair<size_t, int> > spliceway_vec;

	bool wrong_format;

	size_t tagidx;

	double filter_score;

	double junc_anchor_len;

	double ave_intron_len;

	double ave_junc_mis;

	double pair_rate;

	vector<string> junc_id;

	int best;

	string alters;

	string qual_str;

	bool isunique, isexonic, isspliced, issmallins, issmalldel, iscanonical, issemicanonical, isnoncanoical, isunmapped, isclipped;

	size_t canon_count, noncanon_count;

	double canon_rate;

	int matched_id;

	string mate_match;

	size_t mate_offset;

	long mate_diff;

	string mate_match2;

	size_t mate_offset2;

	long mate_diff2;

	string cur_line;

	//bool is_insert;

	size_t intron_size;

	size_t mappedlen;

	size_t mappedlen1;

	size_t mappedlen2;

	unsigned short min_anchor;

	unsigned short max_anchor;

	bool is_fusion;

	bool is_swapped;

	string chrom_name2;

	unsigned short strand_t2;

	char strand1, strand2;

	size_t start2, end2;

	string splice_way2;

	string ori_splice_way2;

	vector<pair<size_t, int> > spliceway_vec2;

	string syn_seq;

	size_t fusion_prefix_st, fusion_prefix_end, fusion_suffix_st, fusion_suffix_end;

	size_t fusion_prefix_len, fusion_suffix_len;

	vector<SpliceWay> left_splice_ways, right_splice_ways;

	vector<pair<size_t, char> > mis_info;

	string mis_info_str;

	size_t m_contig_len;

	//bool is_spliced;

	SamRec() {}

	SamRec(const string& line, int min_ins);

	SamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, 
		const string& alt, const string& qualstr, const string& matematch, size_t mateoffest, int matediff, string line);

	string tostring(size_t total, size_t idx) const;

	//
	string tostandfusion() const;
	void set_unmapped();

	bool modify_jumpcode_by_filtered_junc();

	bool clip_by_small_anchor(bool add_S/*JunctionSeed**/);

	bool clip_by_end_mismatch(bool add_S, size_t max_insert_len/*JunctionSeed**/);
	
	void setfusionbit();

};



class PairedSamRec {

public:

	PairedSamRec(long dist1, long dist2, unsigned short tm, size_t is, size_t maplen, const string& chrom_name1, const string& chrom_name2, 
		size_t strand1, size_t strand2, SamRec* sam_rec1, SamRec* sam_rec2, size_t prefix_len, size_t suffix_len, vector<SpliceWay> l_splice_ways, vector<SpliceWay> r_splice_ways) : 
		mate_dist1(dist1), mate_dist2(dist2), intron_size(is), mappedlen(maplen), total_mismatch(tm), total_filter_score(sam_rec1->filter_score + sam_rec2->filter_score), 
		total_anchor_len(sam_rec1->junc_anchor_len + sam_rec2->junc_anchor_len), total_ave_mismatch(sam_rec1->ave_junc_mis + sam_rec2->ave_junc_mis),
		total_pairing_rate(sam_rec1->pair_rate + sam_rec1->pair_rate), contiglen(sam_rec1->m_contig_len + sam_rec1->m_contig_len), prefixlen(prefix_len), suffixlen(suffix_len),
		left_splice_ways(l_splice_ways), right_splice_ways(r_splice_ways)
	{
		min_anchor = sam_rec1->min_anchor < sam_rec2->min_anchor ? sam_rec1->min_anchor : sam_rec2->min_anchor;

		if (abs(mate_dist1) < abs(mate_dist2))
		{
			mate_dist = abs(mate_dist1);

			real_mate_dist = mate_dist1;

			inner_dist = mate_dist1;

			outter_dist = mate_dist2;
		}
		else
		{
			mate_dist = abs(mate_dist2);

			real_mate_dist = mate_dist2;

			outter_dist = mate_dist1;

			inner_dist = mate_dist2;
		}

		is_same_chrom = (chrom_name1 == chrom_name2);

		is_diff_strand = ((strand1 & IS_REVERSE) != (strand2 & IS_REVERSE));

		paired_sam_rec = make_pair(sam_rec1, sam_rec2);
	}

	long mate_dist1;
	long mate_dist2;
	long inner_dist;
	long outter_dist;
	size_t mate_dist;
	long real_mate_dist;
	size_t intron_size;
	size_t min_anchor;
	size_t mappedlen;
	size_t contiglen;
	size_t prefixlen, suffixlen;
	unsigned short total_mismatch;
	double total_filter_score;
	double total_anchor_len;
	double total_ave_mismatch;
	double total_pairing_rate;
	bool is_same_chrom;
	bool is_diff_strand;
	pair<SamRec*, SamRec*> paired_sam_rec;
	vector<SpliceWay> left_splice_ways;
	vector<SpliceWay> right_splice_ways;
};

#endif