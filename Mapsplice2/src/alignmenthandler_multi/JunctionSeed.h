#ifndef JUNCTIONSEED_H
#define JUNCTIONSEED_H

#include "sharedlib.h"
#include "SpliceWay.h"
#include "SamRec.h"
//chr1    11420   133803  JUNC_1  3       -       11420   133803  255,0,0 2       10,66,  0,122450,       0       6       CTAC    0.99965 0.44492 2       2       2

class JunctionSeed{

public:
	string m_juncname;
	unsigned int m_hits;
	char m_strand;
	char m_strand1, m_strand2;
	unsigned short m_kinds;

	string m_chrom;
	string m_chrom2;
	size_t m_start;
	size_t m_end;

	size_t m_max_prefix_len;
	size_t m_max_suffix_len;

	size_t m_max_fusion_prefix_len;
	size_t m_max_fusion_suffix_len;

	size_t m_start_blockoffset;
	size_t m_end_blockoffset;

	size_t m_encompass_reads_count;

	double m_entropy;
	unsigned short m_flankcase;
	string m_flankstring;
	double m_intronlen;
	double m_lpq;
	double m_il_score;
	unsigned short m_min_mismatch, m_max_mismatch, m_sum_mismatch;
	double m_ave_mismatch;

	unsigned int m_left_exon, m_right_exon;

	vector<unsigned short> m_prefix_count;

	unsigned int m_positive_count, m_negative_count;

	bool m_three_prime_known, m_five_prime_known, m_pair_known;

	map<size_t, size_t> m_three_prime_known_id, m_five_prime_known_id;

	//map<size_t, int> m_mapped_idx;

	map<string, int> m_ins;

	unsigned short m_min_anchor_difference;

	unsigned m_max_min_prefix, m_max_min_suffix;

	unsigned int m_unique_count, m_multi_count;

	unsigned int m_paired_mutiple_count, m_paired_unique_count;

	unsigned int m_encompass_unique_count, m_encompass_multi_count;

	unsigned int m_paired_count, m_single_count;

	unsigned int m_left_paired_count, m_right_paired_count;

	static size_t m_junc_count;

	size_t m_junc_id;

	vector<pair<SamRec, SamRec> > doner_side_spanning_pairs;

	vector<pair<SamRec, SamRec> > accetpr_side_spanning_pairs;

	vector<SamRec > single_spanning;

	vector<size_t > m_fusion_encompassing_reads_doner;

	vector<size_t > m_fusion_encompassing_reads_acceptor;

	FILTERED_TYPE m_filtered_type;

	vector<SpliceWayTrue> left_splice_ways;

	vector<SpliceWayTrue> right_splice_ways;

	vector<pair<size_t, size_t> > left_exons;

	vector<pair<size_t, size_t> > right_exons;

	vector<vector<int> > left_paths;

	vector<vector<int> > right_paths;

	bool m_is_fusion;

	JunctionSeed();

	JunctionSeed(int st, int ed);

	JunctionSeed(const string& juncname, unsigned int hits, char strand, unsigned short kinds, size_t max_prefix_len, size_t max_suffix_len, 
		size_t start_blockoffset, size_t end_blockoffset,
		double entropy, unsigned short flankcase, const string& flankstring, double intronlen, double lpq, unsigned short min_mismatch,
		unsigned short max_mismatch, double ave_mismatch, size_t start, size_t end, const string& chrom, unsigned int unique_count, unsigned int multi_count,
		unsigned int paired_count, unsigned int left_paired_count, unsigned int right_paired_count, unsigned int paired_mutiple_count, unsigned int paired_unique_count,
		unsigned int single_count, unsigned short min_anchor_difference);

	JunctionSeed(size_t loc, size_t suffix_len, size_t fusion_prefix_len, size_t fusion_suffix_len, size_t rw, size_t tagidx, unsigned short mis, size_t strand, size_t strand2, size_t start, size_t end, const string& chrom, const string& chrom2,
		size_t sam_count, const string& mate_match, int mate_diff, const vector<SpliceWay>& l_splice_ways, const vector<SpliceWay>& r_splice_ways, bool is_fusion, SamRec* samrecprt, string insert = "");

	bool inc_hits(size_t idx, size_t suffix_len, size_t fusion_prefix_len, size_t fusion_suffix_len, size_t tagidx, unsigned short mis, size_t strand, size_t sam_count, const string& mate_match, int mate_diff, 
		const vector<SpliceWay>& l_splice_ways, const vector<SpliceWay>& r_splice_ways, SamRec* samrecprt, string insert = "");

	void reset_splice_ways(vector<SamRec*>& fusion_encompassing_reads_doner_filtered,		
						   vector<SamRec*>& fusion_encompassing_reads_acceptor_filtered,
						   vector<pair<SamRec*, SamRec*> >& doner_side_spanning_pairs_filtered,
						   vector<pair<SamRec*, SamRec*> >& accetpr_side_spanning_pairs_filtered,
						   vector<SamRec >& single_spanning);

	void reset_splice_ways(vector<SamRec >& m_fusion_encompassing_reads);

	void clear_splice_ways();

	size_t spliceways2exons_doner();

	size_t spliceways2exons_acceptor();

	size_t hits2exons(vector<pair<size_t, size_t> >& exons, vector<bool>& exon_regions, size_t left_most);

	size_t construct_graph(vector<pair<size_t, size_t> >& exons, vector<SpliceWayTrue>& splice_ways, vector<vector<int> >& graph, char strand);

	void generate_fusion_struct(vector<SamRec >& m_fusion_encompassing_reads);

	bool FindRegion(pair<size_t, size_t>& cur_region, vector<pair<size_t, size_t> >& sorted_regions, size_t& find_region_idx);

	void set_coverage();

	void set_entropy();

	void set_flankstring(/*const string& flankstring*/);

	void set_flankstring(const string& flankstring);

	void set_pq_score(size_t intron_len, size_t chrom_size);

	void set_il_score(size_t junc_st, size_t junc_end);

	void set_ave_mis();

	void set_block_offset(size_t start, size_t end);

	bool clear();

	string to_normal_junction(size_t junc_id, bool isdel=false);

	string to_normal_junction_bed(size_t junc_id);

	string to_insert_junction(size_t junc_id);

};

struct FusionJuncRegion {
	size_t m_doner_st, m_doner_end, m_acceptor_st, m_acceptor_end;

	JunctionSeed* m_junc_seed_ptr;

	FusionJuncRegion(JunctionSeed* junc_seed_ptr);

	FusionJuncRegion();

	FusionJuncRegion(size_t doner_st, size_t doner_end, size_t acceptor_st, size_t acceptor_end);

	void set(size_t doner_st, size_t doner_end, size_t acceptor_st, size_t acceptor_end);
};


#endif