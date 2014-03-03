#ifndef JUNCTIONHANDLER_H
#define JUNCTIONHANDLER_H


#include "sharedlib.h"
#include "SamRec.h"
#include "JunctionSeed.h"

enum Jump_Type {M, N};

const size_t NUM_BP_PER_LINE = 70;

class Jump_Code
{
public:
	size_t len;
	Jump_Type type;

	Jump_Code();

	Jump_Code(size_t _len, Jump_Type _type);
};

class Junc_Seq
{
public:
	size_t start;
	vector<Jump_Code> jump_code;
	vector< pair<size_t, size_t> > sequence;
	char strand;
	unsigned short flank_case;
	size_t junc_id;
	size_t total_len;

	Junc_Seq();

	Junc_Seq(const Junc_Seq& my_junc_seq);

	Junc_Seq(char _strand, 	unsigned short _flank_case, size_t _junc_id, Jump_Code& _new_code);

	Junc_Seq(char _strand, 	unsigned short _flank_case, size_t _junc_id, size_t _start);

	~Junc_Seq();
};

struct JuncSeedSp{

public:
	size_t m_max_prefix_len;
	size_t m_max_suffix_len;

	size_t m_fusion_prefix_len;
	size_t m_fusion_suffix_len;

	size_t m_start;
	size_t m_end;	

	string m_ins_str;

	SamRec* m_sam_rec_ptr;

	bool m_is_fusion;

	JuncSeedSp(size_t max_prefix_len, size_t max_suffix_len, size_t start, size_t end, const string& ins_str, SamRec* sam_rec_ptr, bool is_fusion);

	void set(size_t max_prefix_len, size_t max_suffix_len, size_t start, size_t end, const string& ins_str, SamRec* sam_rec_ptr, bool is_fusion);

};

typedef hash_map<size_t, JunctionSeed> JUNC_HASH_COMB;

typedef hash_map<string, JUNC_HASH_COMB> CHROM_JUNC_HASH_COMB;

bool comp_junc_sort(const JunctionSeed* lhs, const JunctionSeed* rhs);

bool comp_fusion_junc_sort(const FusionJuncRegion& lhs, const FusionJuncRegion& rhs);

bool comp_junc_sort_by_end(const JunctionSeed* lhs, const JunctionSeed* rhs);

class JunctionHandler {

private:

	vector<string> m_junc_files;

	string m_head_line;

	vector<string> m_alignment_files;

	size_t m_max_read_width;

	int m_insert_len;

	int m_delete_len;

	string m_chrom_dir;

	size_t m_canonical, m_semi_canonical, m_non_canonical, m_small_del, m_small_ins, m_fusion, m_fusion_canon, m_fusion_semicanon, m_fusion_noncanon;

	size_t m_min_anchor;

	size_t m_min_mismatch;

	size_t m_min_junc_anchor;

	size_t m_min_fusion_coverage;

	map<size_t, vector<Junc_Seq> > m_head_seq;

	map<size_t, vector<Junc_Seq> > m_tail_seq;

	static map<string, string> m_loaded_chromosomes;

	hash_map<string, size_t> m_chrom_sizes;

	size_t m_do_filter;

public:

	size_t m_total_junc;

	size_t m_total_encompass;

	map<string, vector<JunctionSeed*> > m_junc_sort;

	map<string, vector<FusionJuncRegion> > m_fuson_junc_sort;

	vector<SamRec > m_fusion_encompassing_reads;

	CHROM_JUNC_HASH_COMB m_junc_hash;

	JunctionHandler();

	bool Init(const vector<string>& alignment_files, size_t max_read_width, int insert_len, int delete_len, const string& chrom_dir, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor, size_t min_fusion_coverage, int do_filter);

	bool Clear();

	bool ClearJunction();

	size_t ReadJunction(vector<string>& junc_files);

	void UpdateJuncTable(const string& chromname, size_t combined_offset, size_t prefixlen, size_t suffixlen, size_t fusion_prefixlen, size_t fusion_suffixlen, size_t start, size_t end, const string& ins_str, SamRec& samrec, size_t sam_count, bool is_fusion = false);

	void SamRec2Junc(SamRec& samrec, size_t sam_count, hash_map<string, hash_map<size_t, JuncSeedSp> >& cur_read_junc);

	void SamRecVec2Junc(vector<SamRec*>& samrecs);

	void ReadSam();

	void LoadFlankString();

	void SortJunc();

	void SortJuncByEnd();

	void LoadJuncToSortVec();

	void LoadUnfilteredJuncToSortVec();

	void LoadFusionJuncToSortVec();

	void SortFusionJunc();

	void WriteJunction(string junction_file, string junction_ins_file, string junction_del_file, string junction_fusion_file, string filtered_junc);

	void WriteFusionJunctionWithEncompassingAlignments(string junction_fusion_file, JunctionHandler* filtered_junction_ptr);

	void GenerateFusionStruct();

	void CollectStats();

	void ClearStats();

	void WriteStats(string stat_file, string headline);

	void MarkFiltered(bool paired, JunctionHandler& junc_db);

	//vector<JunctionSeed*>& sorted_juncs = (junc_handler->m_junc_sort.find(expressed_regions_chrom_iter->first))->second;

	//vector<pair<size_t, size_t> >& sorted_regions = expressed_regions_chrom_iter->second;

	//DisjointSet& cur_disjointSet = (m_disjointset_map.find(expressed_regions_chrom_iter->first))->second;


	void FindCorrespondingJunction(SamRec* sam_rec);

	void FindMinPreviousExon(vector<JunctionSeed*>& junc_list, size_t curr_junc);
		
	void FindMinNextExon(vector<JunctionSeed*>& junc_list, size_t curr_junc);
	
	void GetMinimumExon();

	size_t CheckJuncSeqNum(map<size_t, vector<Junc_Seq> >::iterator& _it_head, map<size_t, vector<Junc_Seq> >::iterator& _it_tail);

	void OutPutJuncSeq(size_t junc_id, string chrom, string& chromseq, size_t max_seq_thresh, ofstream& out_fs);

	void AddJunSeq(map<size_t, vector<Junc_Seq> >& junc_seq_list, Junc_Seq& my_junc_seq, size_t max_seq_thresh);

	void GetJuncHead(Junc_Seq& my_junc_seq, vector<JunctionSeed*>& junction_set, size_t curr_junc, size_t min_anchor, size_t max_anchor, size_t max_seq_thresh);

	void GetJuncTail(Junc_Seq& my_junc_seq, vector<JunctionSeed*>& junction_set, size_t curr_junc, size_t min_anchor, size_t max_anchor, size_t chrom_len, size_t max_seq_thresh);

	void GenerateJunctionSequence(size_t min_anchor, size_t max_anchor, size_t max_seq_thresh, string output_file);

	map<string, vector<JunctionSeed*> >& GetJuncSortList();

	void ConfirmJunction(JunctionHandler& database_handler, size_t range, bool clear_internal_result);

	void ReadChromSize(string chrom_size_file);

	void Load2SimpJunc(hash_map<string, vector<pair<int, int> > >& m_original_junc_locs);

	int FindFragmentLengthForward(vector<vector<int> >& paths, vector<pair<size_t, size_t> >& exons, size_t end_point, hash_set<size_t>& stored_fragment_lengths);

	int FindFragmentLengthReverse(vector<vector<int> >& paths, vector<pair<size_t, size_t> >& exons, size_t end_point, hash_set<size_t>& stored_fragment_lengths);

	double PrintExonExpression(vector<vector<int> >& paths, vector<pair<size_t, size_t> >& exons, vector<vector<int> >& exon_express, ofstream* ofs_ptr, double& ks_score);

	void GenerateMeanAndSD(vector<int> expressed_regions, double& mean, double& chisqure);

	void KSTest(vector<int> expressed_regions, double& max_distance, ofstream* ofs);

};




//struct JuncForSort{
//	size_t juncst;
//	size_t juncend;
//	unsigned short hits;
//	unsigned short kinds;
//	string blocks;
//	string blocksoffset;
//	double rank;
//	double lpq;
//	string flankstring;
//	unsigned short flankcase;
//	double intronlen;
//	string juncname;
//
//	unsigned short min_mismatch, max_mismatch;
//
//	double ave_mismatch;
//
//	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
//		const double& rk, const double& l, const string& jn) : juncst(jst), juncend(jend), hits(hts), 
//		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l), juncname(jn) {}
//
//	JuncForSort(const int& jst, const int& jend, const int& hts, const int& kds, const string& blks, const string& blksoft,
//		const double& rk, const double& l, const string& fs, const int fc, const double i, const string& jn, unsigned short min_mis, unsigned short max_mis, double ave_mis) : juncst(jst), juncend(jend), hits(hts), 
//		kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), lpq(l), flankstring(fs), flankcase(fc), intronlen(i), juncname(jn), min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis) {}
//};

#endif