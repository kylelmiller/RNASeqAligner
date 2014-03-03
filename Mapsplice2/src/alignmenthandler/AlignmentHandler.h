#ifndef ALIGNMENTHANDLER_H
#define ALIGNMENTHANDLER_H


#include "sharedlib.h"
#include "SamRec.h"
#include "ReadNextTagAlignHandler.h"
#include "JunctionHandler.h"
#include "FusionSamRec.h"
#include "disjointset.h"
#include "UnionExpressedRegions.h"

class AlignmentHandler {
private:
	bool m_is_paired;
	bool m_add_S;
	size_t m_max_pair_dist;
	size_t m_max_hits;
	vector<string> m_sam_files;
	double m_entrpy_weight;
	double m_pqlen_weight;
	double m_ave_mis_weight;
	string m_filtered_alignment_file;
	ofstream m_ofs_filtered_alignment;
	ofstream m_ofs_fusion_std;
	ofstream m_ofs_fusion_paired;
	ofstream m_ofs_single;
	ofstream m_ofs_paired;
	ofstream m_ofs_to_mapper;
	ReadNextTagAlignHandler<SamRec> m_sam_file_handler;
	ReadNextTagAlignHandler<FusionSamRec> m_fusion_sam_file_handler;
	JunctionHandler m_junction_handler;
	JunctionHandler m_junction_handler_filtered;
	CHROM_JUNC_HASH_COMB* m_junc_hash_ptr;
	pair<vector<SamRec>, vector<SamRec> > m_sam_rec_pe;
	pair<vector<SamRec*>, vector<SamRec*> > m_sam_rec_pe_ptr;
	//vector<PairedSamRec> m_fusion_paired_reads_ptr;

	string m_junction_file;
	size_t m_max_read_width;
	string m_chrom_dir;
	int m_do_filter;
	int m_min_ins;
	int m_max_del;
	size_t m_min_anchor;
	size_t m_min_mismatch;
	size_t m_min_junc_anchor;

	ofstream large_mismatch_ofs, short_50_75_ofs;

	//stats variables
	size_t m_spliced, m_unspliced, m_insertion, m_deletion, m_unique, m_multiple, m_canoical, m_semi_canonical, m_non_canonical, m_paired, m_single, m_unmapped, m_clipped;

	hash_map<string, vector<bool> > m_mapped_pos;

	hash_map<string, vector<pair<size_t, size_t> > > m_expressed_regions;

	char* m_chrom_size_file;

	hash_map<string, DisjointSet> m_disjointset_map;

	hash_map<string, hash_map<size_t, UnionExpressedRegions > > m_unioned_expressed_regions_map;

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
	void SamRec2Hits(SamRec* sam_rec_ptr);
	void SetHits(vector<SamRec*>& samrecs);
	void Hits2ExpressedRegions();
	void WriteAlignment(vector<SamRec*>& sam_rec_ptr);
	void WriteAlignment(vector<SamRec*>& sam_rec_ptr, ofstream& cur_ofs);
	void SetBitInfo(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr);

	void SetFusionBitInfo(pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr);

	void CollectStats(const vector<SamRec*>& sam_rec_ptr);
	void CollectStats(const vector<SamRec>& sam_rec_ptr);

	void UnionSets(JunctionHandler* junc_handler, hash_map<string, vector<pair<size_t, size_t> > >* expressed_regions);

	void FindFusionJuncRegionVec(JunctionHandler* junc_handler, vector<PairedSamRec>& fusion_paired_reads_ptr);

	void FindFusionJuncRegion(FusionJuncRegion& cur_fusion_region, PairedSamRec& cur_paired_sam_rec, vector<FusionJuncRegion>& sorted_fusion_regions);

	bool FindRegion(pair<size_t, size_t>& cur_region, vector<pair<size_t, size_t> >& sorted_regions, size_t& find_region_idx);

	void FindSamVecContigLen(vector<SamRec*>& sam_rec_ptr);

	bool FindSamContigLen(vector<SamRec*>& sam_rec_ptr);

	void LoadSamVec2SamPtr(vector<SamRec>& sam_rec, vector<SamRec*>& sam_rec_ptr);
	
	void SetSamRecUnpaired(pair<vector<SamRec>, vector<SamRec> >& m_sam_rec_pe);

	void ProcessPEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessPERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessPEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessSERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessSEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered);

	void ProcessFile(string alignment_file, JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered, int fileter_by_junc);

public:
	AlignmentHandler(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, string filtered_alignment_file,
		size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor, char* chr_sz_file,
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	bool Init(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, 
		string filtered_alignment_file, size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor,
		double entrpy_weight = 0.097718, double pqlen_weight = 0.66478, double ave_mis_weight = -0.21077);

	bool Clear();

	bool ClearStats();

	void FilterAlignment();

	void WriteStats(string stat_file, string headline);
};


bool comp_min_anchor(const SamRec& lhs, const SamRec& rhs);

bool comp_sam_rec_maplen(const SamRec& lhs, const SamRec& rhs);

bool comp_mate_dist(size_t lhs, size_t rhs);

bool comp_intron_dist(size_t lhs, size_t rhs);

bool comp_dist(const PairedSamRec& lhs, const PairedSamRec& rhs);

bool better_than(const PairedSamRec& lhs, const PairedSamRec& rhs);

bool comp_mis(const SamRec* lhs, const SamRec* rhs);

bool comp_dup(const SamRec* lhs, const SamRec* rhs);

bool comp_filterscore(const SamRec* lhs, const SamRec* rhs);

bool check_overlap(vector<pair<size_t, int> >& spliceway_vec_1, vector<pair<size_t, int> >& spliceway_vec_2);




#endif
