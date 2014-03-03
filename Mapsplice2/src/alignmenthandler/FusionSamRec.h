#ifndef FusionSamRec_H
#define FusionSamRec_H


#include "sharedlib.h"
#include "JunctionSeed.h"
#include "SamRec.h"

class FusionSamRec : public SamRec {

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

	double ave_intron_len;

	double ave_junc_mis;

	vector<string> junc_id;

	int best;

	string alters;

	string qual_str;

	bool isunique, isexonic, isspliced, issmallins, issmalldel, iscanonical, issemicanonical, isnoncanoical, isunmapped, isclipped;

	size_t canon_count, noncanon_count;

	double canon_rate;

	int matched_id;

	char mate_match;

	size_t mate_offset;

	int mate_diff;

	string cur_line;

	//bool is_insert;

	size_t intron_size;

	size_t mappedlen;

	unsigned short min_anchor;

	unsigned short max_anchor;

	string chrom_name2;

	char strand1, strand2;

	size_t start2, end2;

	string splice_way2;

	string ori_splice_way2;

	vector<pair<size_t, int> > spliceway_vec2;

	//bool is_spliced;

	FusionSamRec() {}

	FusionSamRec(const string& line, int min_ins);

	FusionSamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& alt, const string& qualstr, 
		char matematch, size_t mateoffest, int matediff, string line);

	string tostring(size_t total, size_t idx) const;

	void set_unmapped();

	bool modify_jumpcode_by_filtered_junc();

	bool clip_by_small_anchor(bool add_S/*JunctionSeed**/);

};

class PairedFusionSamRec {

public:

	PairedFusionSamRec(size_t dist1, size_t dist2, unsigned short tm, size_t is, size_t maplen, FusionSamRec* sam_rec1, FusionSamRec* sam_rec2) : mate_dist1(dist1), mate_dist2(dist2),
		intron_size(is), mappedlen(maplen), total_mismatch(tm) 
	{
		if (mate_dist1 < mate_dist2)
			mate_dist = mate_dist1;
		else
			mate_dist = mate_dist2;

		is_same_chrom = (sam_rec1->chrom_name == sam_rec2->chrom_name);

		is_diff_strand = (sam_rec1->strand_t == sam_rec2->strand_t);

		paired_sam_rec = make_pair(sam_rec1, sam_rec2);
	}

	size_t mate_dist1;
	size_t mate_dist2;
	size_t mate_dist;
	size_t intron_size;
	size_t mappedlen;
	unsigned short total_mismatch;
	bool is_same_chrom;
	bool is_diff_strand;
	pair<FusionSamRec*, FusionSamRec*> paired_sam_rec;
};

#endif