#ifndef ALIGNMENTHANDLER_H
#define ALIGNMENTHANDLER_H


#include "sharedlib.h"
#include "SamRec.h"
//#include "ReadNextTagAlignHandler.h"
#include "JunctionHandler.h"
#include "disjointset.h"
#include "UnionExpressedRegions.h"

class AlignmentHandler {
private:
	
	
	
	
	vector<string> m_sam_files;
	
	string m_filtered_alignment_file;
	ofstream m_ofs_filtered_alignment;
	ofstream m_ofs_fusion_std;
	ofstream m_ofs_fusion_paired;
	ofstream m_ofs_single;
	ofstream m_ofs_paired;
	ofstream m_ofs_to_mapper;
	ofstream m_ofs_bothunspliced;
	ofstream m_ofs_onespliced;

	ofstream m_ofs_unspliced_fusion_paired;
	ofstream m_ofs_unspliced_single;
	ofstream m_ofs_unspliced_paired;

	ofstream m_ofs_unmapped;

	//ReadNextTagAlignHandler<SamRec> m_sam_file_handler;
	//ReadNextTagAlignHandler<FusionSamRec> m_fusion_sam_file_handler;
	
	CHROM_JUNC_HASH_COMB* m_junc_hash_ptr;
	pair<vector<SamRec>, vector<SamRec> > m_sam_rec_pe;
	pair<vector<SamRec*>, vector<SamRec*> > m_sam_rec_pe_ptr;
	//vector<PairedSamRec> m_fusion_paired_reads_ptr;

	string m_junction_file;
	size_t m_max_read_width;
	static string m_chrom_dir;
	
	int m_min_ins;
	
	int m_max_del;
	size_t m_min_anchor;
	size_t m_min_mismatch;
	size_t m_min_junc_anchor;
	size_t m_min_fusion_coverage;

	ofstream large_mismatch_ofs, short_50_75_ofs;

	ifstream m_ifs_alignments;

	//stats variables
	//size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_semi_canonical, m_non_canonical, m_paired, m_single, m_unmapped, m_clipped;

	//parallel buffer
	static size_t m_cur_read_alignments_index;// =  0;

	static size_t m_cur_line_index;

	static vector<string> m_stored_lines;

	static vector<SamRec> m_stored_sam_recs;

	static vector<pair<vector<SamRec>::iterator, vector<SamRec>::iterator> > m_same_reads_alignments;

	static vector<pair<vector<SamRec*>, vector<SamRec*> > > m_same_reads_alignments_paired;

	static hash_map<string, vector<pair<int, int> > > m_normal_paired_locs;

	static hash_map<string, vector<pair<int, int> > > m_original_junc_locs;

	static vector<bool> m_same_reads_alignments_paired_is_unspliced;

	static int m_min_ins_static;

	static int m_do_filter;

	static size_t m_max_pair_dist;

	static JunctionHandler m_junction_handler;

	static JunctionHandler m_junction_handler_intermediate;

	static JunctionHandler m_junction_handler_filtered;

	static JunctionHandler m_junction_handler_annotated;

	static JunctionHandler* m_junction_handler_prev_ptr;

	static double m_entrpy_weight;

	static double m_pqlen_weight;

	static double m_ave_mis_weight;

	static bool m_add_S;

	static size_t m_max_hits;

	static size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_semi_canonical, m_non_canonical, m_paired, m_single, m_fusion_paired, m_unmapped, m_clipped, m_fusion;

	//static size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_canoical, m_semi_canonical, m_non_canonical;

	//spliced, unspliced, insertion, deletion vs paired, single, fusion_paired

	static size_t m_spliced_paired, m_spliced_single, m_spliced_fusion_paired;
	
	static size_t m_unspliced_paired, m_unspliced_single, m_unspliced_fusion_paired;

	static size_t m_insertion_paired, m_insertion_single, m_insertion_fusion_paired;

	static size_t m_deletion_paired, m_deletion_single, m_deletion_fusion_paired;

	static size_t m_canoical_paired, m_canoical_single, m_canoical_fusion_paired;

	static size_t m_semi_canonical_paired, m_semi_canonical_single, m_semi_canonical_fusion_paired;

	static size_t m_non_canonical_paired, m_non_canonical_single, m_non_canonical_fusion_paired;

	static size_t m_clipped_paired, m_clipped_single, m_clipped_fusion_paired;

	//fusion vs canonical, semicanonical noncanical

	static size_t m_fusion_canonical, m_fusion_semi_canonical, m_fusion_non_canonical;

	static size_t m_fusion_multiple, m_fusion_unique;

	///////////

	static size_t m_both_unspliced_multiple, m_both_unspliced_unique, m_both_unspliced_paired, m_both_unspliced_single, m_both_unspliced_fusion_paired;

	//paired, single and fusion_paired vs multiple and unique

	static size_t m_paired_multiple, m_paired_unique, m_single_multiple, m_single_unique, m_fusion_paired_multiple, m_fusion_paired_unique;

	static size_t m_both_unspliced_paired_multiple, m_both_unspliced_paired_unique;
	
	static size_t m_both_unspliced_single_multiple, m_both_unspliced_single_unique;

	static size_t m_both_unspliced_fusion_paired_multiple, m_both_unspliced_fusion_paired_unique;

	///////////Final stats///////////////////////

	//static size_t m_final_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_semi_canonical, m_non_canonical, m_paired, m_single, m_fusion_paired, m_unmapped, m_clipped, m_fusion;

	//static size_t m_both_unspliced_multiple, m_both_unspliced_unique, m_both_unspliced_paired, m_both_unspliced_single, m_both_unspliced_fusion_paired;

	//static size_t m_both_unspliced_paired_multiple, m_both_unspliced_paired_unique;
	//
	//static size_t m_both_unspliced_single_multiple, m_both_unspliced_single_unique;

	//static size_t m_both_unspliced_fusion_paired_multiple, m_both_unspliced_fusion_paired_unique;

	/////////////////////

	static string m_prev_tag;

	static bool m_is_paired;
	/////
	static hash_map<string, vector<bool> > m_mapped_pos;

	static hash_map<string, vector<pair<size_t, size_t> > > m_expressed_regions;

	char* m_chrom_size_file;

	static hash_map<string, DisjointSet> m_disjointset_map;

	static hash_map<string, hash_map<size_t, UnionExpressedRegions > > m_unioned_expressed_regions_map;

	

private:
	void MarkCanonNoncanonByReads(vector<SamRec>& read_sam, JunctionHandler* junc_hash_ptr);
	void FilterByMinAnchor(vector<SamRec>& read_sam);
	void FilterByMapLen(vector<SamRec>& read_sam);
	void FilterCanonNonCanonByReads(vector<SamRec>& read_sam, vector<SamRec* >& filtered_read_sam_ptr);
	bool EstablishPairing(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr, vector<PairedSamRec>& fusion_paired_reads_ptr);
	void FilterPairedSamRec(vector<PairedSamRec>& pairedsamrec_vec, pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr);
	void FilterSingleMulti(vector<SamRec*>& sam_rec_ptr);
	void FilterByFilteredJunction(vector<SamRec*>& sam_rec_ptr, JunctionHandler* junc_handler);
	void RemoveDup(vector<SamRec*>& sam_rec_ptr);
	void FilterByMisMatch(vector<SamRec*>& sam_rec_ptr);
	void AllocateHitsMemory();
	static void SamRec2Hits(SamRec* sam_rec_ptr);
	static void SetHits(vector<SamRec*>& samrecs);
	void Hits2ExpressedRegions();
	void WriteAlignment(vector<SamRec*>& sam_rec_ptr);
	void WriteAlignment(vector<SamRec*>& sam_rec_ptr, ofstream& cur_ofs);
	void WriteAlignmentUnspliced(vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator pair_samvec_iter, bool is_unspliced);
	void SetBitInfo(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr);

	void SetFusionBitInfo(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr);

	void CollectStats(const vector<SamRec*>& sam_rec_ptr);
	void CollectStats(const vector<SamRec>& sam_rec_ptr);

	void UnionSets(JunctionHandler* junc_handler, hash_map<string, vector<pair<size_t, size_t> > >* expressed_regions);

	void FilterFusionByIsolatedExons(JunctionHandler* junc_handler);

	void FindFusionJuncRegionVec(JunctionHandler* junc_handler, vector<PairedSamRec>& fusion_paired_reads_ptr);

	void FindFusionJuncRegion(FusionJuncRegion& cur_fusion_region, PairedSamRec& cur_paired_sam_rec, vector<FusionJuncRegion>& sorted_fusion_regions);

	static bool FindRegion(pair<size_t, size_t>& cur_region, vector<pair<size_t, size_t> >& sorted_regions, size_t& find_region_idx);

	void FindSamVecContigLen(vector<SamRec*>& sam_rec_ptr);

	bool FindSamContigLen(vector<SamRec*>& sam_rec_ptr);

	void LoadSamVec2SamPtr(vector<SamRec>& sam_rec, vector<SamRec*>& sam_rec_ptr);
	
	void SetSamRecUnpaired(pair<vector<SamRec>, vector<SamRec> >& m_sam_rec_pe);

	void EstimateFragmentLength(JunctionHandler* junc_handler);	

	void ProcessPEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessPERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessPEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessSERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessSEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessSEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessFile(string alignment_file, JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered, int fileter_by_junc);

	static void* ParseLines(void * str);

	static void* ProcessPEReadMulti(void * str);

	static void* ProcessSEReadMulti(void * str);

	static void* ProcessPEReadFilterJuncMulti(void * str);

	static void* ProcessSEReadFilterJuncMulti(void * str);

	static void* ProcessPEReadSelectBestMulti(void * str);

	static void* ProcessSEReadSelectBestMulti(void * str);

	static bool EstablishPairingStatic(void * str);

	static bool EstablishFusionPairingStatic(void * str);

	static void* FilterPairedSamRecStatic(void * str);

	static void* MarkCanonNoncanonByReadsStatic(void * str);

	static void* FilterByFilteredJunctionStatic(void * str);

	static void* RemoveDupStatic(void * str);

	static void* FindFusionJuncRegionVecStatic(void * str);

	static void* FindFusionJuncRegionStatic(void * str);

	static void* FilterSingleMultiStatic(void * str);

	static void* RemoveFusionUnPairedStatic(void * str);

	static void* FilterByMisMatchStatic(void * str);

	static void* CollectStatsStatic(void * str);

	static void* CollectStatsStatic(void * str, bool both_unspliced);

	static void* SetFusionBitInfoStatic(void * str);

	static void* SetBitInfoStatic(void * str);

	static void* GenerateAlignment(void * str, bool fusion_std = true);

	static void SetPairedType(vector<SamRec*>& sam_rec_ptr, PAIRED_TYPE paired_type);

	static void SetSamRecUnpaired(pair<vector<SamRec*>, vector<SamRec*> >& m_sam_rec_pe_ptr);

	static bool MergeStdFusion(vector<SamRec>::iterator sam_rec_ptr1, vector<SamRec>::iterator sam_rec_ptr2);

	

public:
	AlignmentHandler(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, string filtered_alignment_file,
		size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor, char* chr_sz_file,
		size_t min_fusion_coverage,
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	bool Init(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, 
		string filtered_alignment_file, size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor,
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	bool Clear();

	bool ClearStats();

	void FilterAlignment();

	void FilterAlignmentFiltered();

	void WriteStats(string stat_file, string headline);

	void WriteFinalStats(string stat_file, string headline);

	
};

bool comp_mate_dist(size_t lhs, size_t rhs);

bool comp_intron_dist(size_t lhs, size_t rhs);

bool comp_filter_score_range(double lhs, double rhs);

bool comp_dup(const SamRec* lhs, const SamRec* rhs);

bool check_overlap(vector<pair<size_t, int> >& spliceway_vec_1, vector<pair<size_t, int> >& spliceway_vec_2);




#endif
