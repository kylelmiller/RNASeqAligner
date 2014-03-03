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

	size_t m_canonical, m_semi_canonical, m_non_canonical, m_small_del, m_small_ins;

	size_t m_min_anchor;

	size_t m_min_mismatch;

	size_t m_min_junc_anchor;

	map<size_t, vector<Junc_Seq> > m_head_seq;

	map<size_t, vector<Junc_Seq> > m_tail_seq;

public:

	map<string, vector<JunctionSeed*> > m_junc_sort;

	map<string, vector<FusionJuncRegion> > m_fuson_junc_sort;

	CHROM_JUNC_HASH_COMB m_junc_hash;

	JunctionHandler();

	bool Init(const vector<string>& alignment_files, size_t max_read_width, int insert_len, int delete_len, const string& chrom_dir, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor);

	bool Clear();

	size_t ReadJunction(vector<string>& junc_files);

	void UpdateJuncTable(const string& chromname, size_t combined_offset, size_t prefixlen, size_t suffixlen, size_t start, size_t end, const string& ins_str, const SamRec& samrec, size_t sam_count, bool is_fusion = false);

	void SamRec2Junc(SamRec& samrec, size_t sam_count, hash_map<string, hash_map<size_t, JuncSeedSp> >& cur_read_junc);

	void SamRecVec2Junc(vector<SamRec*>& samrecs);

	void ReadSam();

	void LoadFlankString();

	void SortJunc();

	void SortJuncByEnd();

	void LoadJuncToSortVec();

	void LoadFusionJuncToSortVec();

	void SortFusionJunc();

	void WriteJunction(string junction_file, string junction_ins_file, string junction_del_file, string junction_fusion_file, string filtered_junc);

	void WriteFusionJunctionWithEncompassingAlignments(string junction_fusion_file);

	void GenerateFusionStruct();

	void CollectStats();

	void ClearStats();

	void WriteStats(string stat_file, string headline);

	void MarkFiltered(bool paired);

	void FindCorrespondingJunction(SamRec* sam_rec);

	void FindMinPreviousExon(vector<JunctionSeed*>& junc_list, size_t curr_junc);
		
	void FindMinNextExon(vector<JunctionSeed*>& junc_list, size_t curr_junc);
	
	void GetMinimumExon();

	size_t CheckJuncSeqNum(map<size_t, vector<Junc_Seq> >::iterator& _it_head, map<size_t, vector<Junc_Seq> >::iterator& _it_tail);

	void OutPutJuncSeq(string chrom, string& chromseq, size_t max_seq_thresh, ofstream& out_fs);

	void AddJunSeq(map<size_t, vector<Junc_Seq> >& junc_seq_list, Junc_Seq& my_junc_seq, size_t max_seq_thresh);

	void GetJuncHead(Junc_Seq& my_junc_seq, vector<JunctionSeed*>& junction_set, size_t curr_junc, size_t min_anchor, size_t max_anchor, size_t max_seq_thresh);

	void GetJuncTail(Junc_Seq& my_junc_seq, vector<JunctionSeed*>& junction_set, size_t curr_junc, size_t min_anchor, size_t max_anchor, size_t chrom_len, size_t max_seq_thresh);

	void GenerateJunctionSequence(size_t min_anchor, size_t max_anchor, size_t max_seq_thresh, string output_file);

	map<string, vector<JunctionSeed*> >& GetJuncSortList();

	void ConfirmJunction(JunctionHandler& database_handler, size_t range, bool clear_internal_result);

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