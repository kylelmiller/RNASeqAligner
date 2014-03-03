#ifndef MAPSPLICE_H
#define MAPSPLICE_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
#include "ebwt.h"
#include "hit.h"
#include "pat.h"
#include "bwtmap_info.h"
#include "read_block.h"
#include "DoubleAnchorScore.h"
#include "sam_info.h"
#include "quality_score.h"
#include "cluster.h"
#include "bm.h"
#include "sbndm.h"
#include "stats.h"
#include "local_align.h"


const int MAX_JUMP = 3;      // maximum missed segments between two segments
const int HEAD_SEG_NO = 1;   // first segment number
const int CONTIG_NO_RELATIONSHIP = 0;
const int CONTIG_RELATIONSHIP = 1;
const int CONTIG_CONTINUE_SEARCH = -1;

const int FIX_NO_RELATIONSHIP = 0;
const int FIX_TOO_CLOSE = 1;
const int FIX_DOUBLE_ANCHOR = 2;
const int FIX_TOO_FAR = 3;
const int FIX_INSERTION = 4;

#define SET_PARAM_PTR(p) \
	assert(!empty(p->bufa().patFw)); \
	String<Dna5>& patFw  = p->bufa().patFw;  \
	patFw.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patRc  = p->bufa().patRc;  \
	patRc.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& qual = p->bufa().qual; \
	qual.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& qualRev = p->bufa().qualRev; \
	qualRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patFwRev  = p->bufa().patFwRev;  \
	patFwRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<Dna5>& patRcRev  = p->bufa().patRcRev;  \
	patRcRev.data_begin += 0; /* suppress "unused" compiler warning */ \
	String<char>& name   = p->bufa().name;   \
	name.data_begin += 0; /* suppress "unused" compiler warning */ \
	uint32_t      patid  = p->patid();       \
	params.setPatId(patid);

inline bool bwtmap_compare(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)   //sort compare function
{
	if(seg1.pair_no != seg2.pair_no)
		return seg1.pair_no < seg2.pair_no;
	if(seg1.chrom != seg2.chrom)
		return seg1.chrom < seg2.chrom;
	if(seg1.strand != seg2.strand)
		return seg1.strand < seg2.strand;
	if (seg1.start != seg2.start)
		return seg1.start < seg2.start;
	if(seg1.end != seg2.end)
		return seg1.end <seg2.end;
	else
		return seg1.start_seg_no < seg2.start_seg_no;
}

inline bool sort_sam_by_pos(const Sam_Info& sam1, const Sam_Info& sam2)
{
	if(sam1.chrom != sam2.chrom)
		return sam1.chrom < sam2.chrom;
	if(sam1.start_pos != sam2.start_pos)
		return sam1.start_pos < sam2.start_pos;
	else
		return sam1.end_pos < sam2.end_pos;
}

/*inline bool comp_dist(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs)
{
	if (lhs.mappedlen != rhs.mappedlen)
		return lhs.mappedlen > rhs.mappedlen;
	if (lhs.total_mismatch + lhs.insert_size != rhs.total_mismatch + rhs.insert_size)
		return lhs.total_mismatch + lhs.insert_size < rhs.total_mismatch + rhs.insert_size; 
	if(lhs.insert_size != rhs.insert_size)
		return lhs.insert_size < rhs.insert_size;
	else
		return lhs.mate_dist < rhs.mate_dist;
}*/

inline bool comp_dist_maplen(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs)
{
		return lhs.mappedlen > rhs.mappedlen;
}

inline bool comp_dist_mismatch(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs)
{
		return lhs.total_mismatch + lhs.insert_size < rhs.total_mismatch + rhs.insert_size; 
}

inline bool comp_dist_insert(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs)
{
		return lhs.insert_size < rhs.insert_size;
}

inline bool comp_dist_matedist(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs)
{
		return lhs.mate_dist_diff < rhs.mate_dist_diff;
}

inline bool comp_dist_intron(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs)
{
	return lhs.intron_size < rhs.intron_size;
}

void sort_pair_alignment(vector<Paired_Sam_Info>::iterator _it_begin, vector<Paired_Sam_Info>::iterator _it_end, vector<Paired_Sam_Info>::iterator& it_reserve, int mate_diff, int intron_diff)
{
	vector<Paired_Sam_Info>::iterator it_begin = _it_begin;
	vector<Paired_Sam_Info>::iterator it_end = _it_end;
	vector<Paired_Sam_Info>::iterator it_index;
		
	sort(it_begin, it_end, comp_dist_maplen);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if(it_index->mappedlen < it_begin->mappedlen)
			break;
	}
	it_end = it_index;
	
	sort(it_begin, it_end, comp_dist_mismatch);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if(it_index->total_mismatch + it_index->insert_size > it_begin->total_mismatch + it_begin->insert_size)
			break;
	}
	it_end = it_index;	
	
	sort(it_begin, it_end, comp_dist_insert);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if(it_index->insert_size > it_begin->insert_size)
			break;
	}
	it_end = it_index;	
	
	sort(it_begin, it_end, comp_dist_matedist);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if(it_index->mate_dist_diff - it_begin->mate_dist_diff > mate_diff)
			break;
	}
	it_end = it_index;	
	
	sort(it_begin, it_end, comp_dist_intron);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if((it_begin->intron_size == 0 && it_index->intron_size != 0) || (it_index->intron_size - it_begin->intron_size > intron_diff))
			break;
	}
	it_reserve = it_index;	
}

inline bool comp_dist_maplen_s(const Sam_Info* lhs, const Sam_Info* rhs)
{
		return lhs->mapped_len > rhs->mapped_len;
}

inline bool comp_dist_mismatch_s(const Sam_Info* lhs, const Sam_Info* rhs)
{
		return lhs->mismatchs.size() + lhs->insert_len + lhs->num_junc_penalty < rhs->mismatchs.size() + rhs->insert_len + rhs->num_junc_penalty;
}

inline bool comp_dist_insert_s(const Sam_Info* lhs, const Sam_Info* rhs)
{
		return lhs->insert_len < rhs->insert_len;
}

inline bool comp_dist_intron_s(const Sam_Info* lhs, const Sam_Info* rhs)
{
		return lhs->intron_len < rhs->intron_len;
}

void sort_single_alignment(vector<Sam_Info*>::iterator _it_begin, vector<Sam_Info*>::iterator _it_end, vector<Sam_Info*>::iterator& it_reserve, int intron_diff)
{
	vector<Sam_Info*>::iterator it_begin = _it_begin;
	vector<Sam_Info*>::iterator it_end = _it_end;
	vector<Sam_Info*>::iterator it_index;
		
	sort(it_begin, it_end, comp_dist_maplen_s);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if((*it_index)->mapped_len < (*it_begin)->mapped_len)
			break;
	}
	it_end = it_index;
	
	sort(it_begin, it_end, comp_dist_mismatch_s);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if((*it_index)->mismatchs.size() + (*it_index)->insert_len + (*it_index)->num_junc_penalty > (*it_begin)->mismatchs.size() + (*it_begin)->insert_len + (*it_begin)->num_junc_penalty )
			break;
	}
	it_end = it_index;	
	
	sort(it_begin, it_end, comp_dist_insert_s);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if((*it_index)->insert_len > (*it_begin)->insert_len)
			break;
	}
	it_end = it_index;	
	
	sort(it_begin, it_end, comp_dist_intron_s);
	for(it_index = it_begin + 1; it_index != it_end; it_index++)
	{
		if(((*it_begin)->intron_len == 0 && (*it_index)->intron_len != 0) || ((*it_index)->intron_len - (*it_begin)->intron_len > intron_diff))
			break;
	}
	it_reserve = it_index;	
}

/*inline int equal_dist(const Paired_Sam_Info& lhs, const Paired_Sam_Info& rhs, int mate_dist_diff, int intron_diff)
{
	if (lhs.mappedlen > rhs.mappedlen)
		return -1;
	else if (lhs.total_mismatch + lhs.insert_size < rhs.total_mismatch + rhs.insert_size)
		return -1; 
	else if (lhs.insert_size < rhs.insert_size)	
		return -1;
	else if ((int)(rhs.mate_dist) - (int)(lhs.mate_dist) > mate_dist_diff)
		return -1;
	else if((int)(rhs.intron_size) - (int)(lhs.intron_size) > intron_diff)
		return 1;
	else	
		return 0;	
}*/

inline	bool in_range(string src_chrom, int src_start, int src_end, string dst_chrom, int dst_start, int dst_end)
{
	if(src_chrom == dst_chrom && src_start >= dst_start && src_end <= dst_end)
		return true;
	else
		return false;
}

class in_dist
{
public:
	in_dist(Bwtmap_Info* _my_bwt, int _max_dist): my_bwt(_my_bwt), max_dist(_max_dist){}

	~in_dist() 
	{
		fusion_region.clear();
	}

	bool operator() (const Sam_Info& my_sam) 
	{
		if(my_bwt->chrom == my_sam.chrom && (abs(my_bwt->start - my_sam.end_pos) < max_dist || abs(my_bwt->end - my_sam.start_pos) < max_dist))
			return true;
		return false;
	}

private:
	Bwtmap_Info* my_bwt;
	vector<Cluster*> fusion_region;
	int max_dist;
};

inline bool fusion_cluster_cmp(const Cluster& fusion_cluster1, const Cluster& fusion_cluster2)
{
	if(fusion_cluster1.start1 != fusion_cluster2.start1)
		return fusion_cluster1.start1 < fusion_cluster2.start1;
	else
		return fusion_cluster1.end1 < fusion_cluster2.end1;
}

class MapSplice
{
private: 
	vector<Bwtmap_Info>* bwt_vector;   // mapped segments vector 
	vector<Bwtmap_Info>* bwt_vector_normal;
	vector<Bwtmap_Info> bwt_vector_repeats;
	Read_Block* r_block_pe1;               // read sequence
	Read_Block* r_block_pe2;
	map<string, size_t>* chrom_map;    // reference sequence name and coresponding length
	map<string, string>* reference_sequence; // reference sequence
	GreedyDFSRangeSource* bt_ptr;      // bowtie mapping class
	PatternSourcePerThread* patsrc;    // bowtie read source
	HitSinkPerThread* sink;            // bowtie hitsink
	Ebwt< String<Dna, Alloc<> > >* original_ebwt;
	Ebwt< String<Dna, Alloc<> > >* original_ebwtBw;
	Ebwt< String<Dna, Alloc<> > >* juncdb_ebwt;
	Ebwt< String<Dna, Alloc<> > >* juncdb_ebwtBw;
	Ebwt< String<Dna, Alloc<> > >* fusiondb_ebwt;
	Ebwt< String<Dna, Alloc<> > >* fusiondb_ebwtBw;
	EbwtSearchParams<String<Dna> >* params_ptr;  // bowtie mapping parameter
	ReadBuf* full_readbufa;             // readbuf of current read
	ReadBuf* full_readbufb;
	GenomeScan genome_scan;           // double anchor fixing class
	pthread_mutex_t* outfs_lock;
	pthread_mutex_t* fusion_fs_lock;
	map<string, vector<Cluster> >* fusion_cluster;
	SBNDM_BACKWARD anchor_search_head;
	SBNDM_FORWARD anchor_search_tail;
	Local_Align local_align;

	int max_read_length;    // maximum length of read
	int segment_length;     // length per segment
	int new_seg_thresh;
	int min_intron_length;  // minimum intron length allowed
	int max_intron_length_double_anchor;  // maximum intron length allowed
	int max_intron_length_single_anchor;  // maximum intron length allowed
	int trunc_len;
	int extend_len;
	int max_ins;
	int max_ins_mismatch;
	int max_del;
	int anchor_length;
	size_t max_double_splice_mismatch;  // maximum mismatches allowed in seed
	size_t max_single_splice_mismatch;
	size_t max_append_mismatch;	
	int min_mapped_len;              // minimum consequtive segments allowed for report
	bool fa;
	bool pair_end;           // whether read is paired
	bool juncdb;
	bool fusiondb;
	bool debug;
	bool output_unmapped;
	bool output_unmapped_pe;
	bool double_anchor_noncanonical;
	bool single_anchor_noncanonical;
	bool fusion_double_anchor_noncanonical;
	bool fusion_single_anchor_noncanonical;
	uint32_t seed;               // bowtie generate random seed

	bool nofw;                   // no forward mapping for bowtie
	bool norc;                   // no reverse mapping for bowtie
	int pe1_segment_number;     // number of segment for pe1 read
	int pe2_segment_number;			// number of segment for pe2 read
	int pe1_read_length;
	int pe2_read_length;
	bool pe1_unspliced;
	bool pe1_spliced;
	bool pe1_aligned;
	bool pe1_fusion;
	bool pe2_unspliced;
	bool pe2_spliced;
	bool pe2_aligned;
	bool pe2_fusion;
	bool* pe1_no_repeats;
	bool* pe2_no_repeats; 
	vector<bool>* pe1_segment_no_repeats;
	vector<bool>* pe2_segment_no_repeats; 
	size_t num_pe1_alignment;
	size_t num_pe2_alignment;
	bool read_paired;
	bool fusion;
	uint32_t max_anchor_hits;
	uint32_t max_remap_hits;
	uint32_t max_repeat_hits;
	uint32_t max_repeat_all_hits;
	uint32_t max_segment_hits;
	size_t anchor_number;
	size_t remap_number;
	bool optimize_for_repeats;
	size_t maximum_alignments;
	int maximum_alignments_output;
	bool do_pairing;
	bool splice_only;
	bool try_hard_ins;
	bool do_local_align;
	int local_align_max_mis;
	int local_align_max_gap;
	int local_align_max_dist;
	int try_hard_max_ins;
	int try_hard_max_mis;
	int pairing_mate_dist_diff;
	int pairing_intron_diff;
	vector<Sam_Info> sam_vec;
	vector<Sam_Info> fusion_candidate;
	vector<Paired_Sam_Info> paired_sam_vec;
	vector<Sam_Info*> pe1_sam;
	vector<Sam_Info*> pe2_sam;
	string global_qual_string;    // global quality string for Fasta read
	ofstream* sam_fs;  
	ofstream* fusion_fs;       
	ofstream* unmapped_fs1;
	ofstream* unmapped_fs2; 
	Mapping_Stats* mapping_stats;
	string quality_scale;
	int mate_distance_training_set;
	int mate_distance_training_count;
	vector<int> mate_distance_training_array;
	int adaptive_mate_pair_dist;
	int max_pair_dist;
	int min_fusion_distance;
	int max_semicanonical_intron_length;
	int min_semicanonical_anchor_length;


public:
	MapSplice()
	{}

	MapSplice(int seg_len, int _new_seg_thresh, int min_thr, int max_thr_d, int max_thr_s, int ext_len, int max_insertion, int max_deletion, bool fasta, 
		bool pairend,  uint32_t _seed, 
		size_t _max_double_splice_mismatch, size_t _max_single_splice_mismatch, size_t _max_append_mismatch, int _min_mapped_len, int _anchor_len, bool _juncdb, bool _fusiondb,
		bool _debug, uint32_t _max_segment_hits, uint32_t _max_anchor_hits, uint32_t _max_remap_hits, uint32_t _max_repeat_hits, 
		uint32_t _max_repeat_all_hits, size_t _maximum_alignments, bool _do_pairing, bool _splice_only, int _max_pair_dist, bool _optimize_for_repeats, bool _output_unmapped, bool _output_unmapped_pe, 
		bool _fusion, bool _double_anchor_noncanonical, bool _single_anchor_noncanonical, bool _fusion_double_anchor_noncanonical, bool _fusion_single_anchor_noncanonical, bool _do_local_align, 
		int _local_align_max_mis, int _local_align_max_gap, int _local_align_max_dist, string _quality_scale, int _mate_distance_training_set, int _min_fusion_distance)
	{
		segment_length = seg_len;
		new_seg_thresh = _new_seg_thresh;
		trunc_len = seg_len / 2;
		extend_len = ext_len;
		min_intron_length = min_thr;
		max_intron_length_double_anchor = max_thr_d;
		max_intron_length_single_anchor = max_thr_s;
		max_ins = max_insertion;
		max_del = max_deletion;
		pair_end = pairend;
		fa = fasta;
		seed = _seed;
		anchor_length = _anchor_len;
		max_double_splice_mismatch = _max_double_splice_mismatch;
		max_single_splice_mismatch = _max_single_splice_mismatch;
		max_append_mismatch =  _max_append_mismatch;
		min_mapped_len = _min_mapped_len;
		pair_end = pairend;
		juncdb = _juncdb;
		fusiondb = _fusiondb;
		debug = _debug;
		max_anchor_hits= _max_anchor_hits;
		max_remap_hits = _max_remap_hits;
		double_anchor_noncanonical = _double_anchor_noncanonical;
		single_anchor_noncanonical = _single_anchor_noncanonical;
		max_repeat_hits = _max_repeat_hits;
		optimize_for_repeats = _optimize_for_repeats;
		do_pairing = pair_end && _do_pairing;
		splice_only = _splice_only;
		max_pair_dist = _max_pair_dist + max_intron_length_double_anchor;
		max_repeat_all_hits = _max_repeat_all_hits;
		output_unmapped = _output_unmapped;
		output_unmapped_pe = _output_unmapped_pe;
		fusion = _fusion;
		max_segment_hits = _max_segment_hits;
		maximum_alignments = _maximum_alignments;
		max_ins_mismatch = 1;
		max_read_length = 2000;
		global_qual_string.assign(max_read_length, 'I');
		fusion_double_anchor_noncanonical = _fusion_double_anchor_noncanonical;
		fusion_single_anchor_noncanonical = _fusion_single_anchor_noncanonical;
		pairing_mate_dist_diff = 500;
		pairing_intron_diff = 5000;
		try_hard_ins = true;
		try_hard_max_ins = 1;
	  try_hard_max_mis = 0;
	  do_local_align = _do_local_align;
	  local_align_max_mis = _local_align_max_mis;
	  local_align_max_gap = _local_align_max_gap;
	  local_align_max_dist = _local_align_max_dist;
	  quality_scale = _quality_scale;
	  mate_distance_training_set = _mate_distance_training_set;
	  mate_distance_training_count = 0;
	  adaptive_mate_pair_dist = 0;
	  min_fusion_distance = _min_fusion_distance;
	  min_semicanonical_anchor_length = 15;
	  if(double_anchor_noncanonical)
	    max_semicanonical_intron_length = max_intron_length_double_anchor;
	  else
	  	max_semicanonical_intron_length = min(50000, max_intron_length_double_anchor);
	  maximum_alignments_output = 10;
	}

	~MapSplice()
	{
	}

	void set_parameter(vector<Bwtmap_Info>* bwt_vec, Read_Block* r_b_pe1, Read_Block* r_b_pe2, map<string, size_t>* _chrom_map, GreedyDFSRangeSource* _bt, PatternSourcePerThread* _patsrc, HitSinkPerThread* _sink, Ebwt< String<Dna, Alloc<> > >* _original_ebwt, Ebwt< String<Dna, Alloc<> > >* _original_ebwtBw, Ebwt< String<Dna, Alloc<> > >* _juncdb_ebwt, Ebwt< String<Dna, Alloc<> > >* _juncdb_ebwtBw, Ebwt< String<Dna, Alloc<> > >* _fusiondb_ebwt, Ebwt< String<Dna, Alloc<> > >* _fusiondb_ebwtBw, ReadBuf* _full_readbufa, ReadBuf* _full_readbufb, EbwtSearchParams<String<Dna> >* _params, map<string, string>* _reference_sequence, ofstream* _sam_fs, ofstream* _unmapped_fs1, ofstream* _unmapped_fs2, pthread_mutex_t* _outfs_lock, bool* _pe1_no_repeats, bool* _pe2_no_repeats, vector<bool>* _pe1_segment_no_repeats, vector<bool>* _pe2_segment_no_repeats, map<string, vector<Cluster> >* _fusion_cluster, ofstream* _fusion_fs, pthread_mutex_t* _fusion_fs_lock, Mapping_Stats* _mapping_stats)
	{
		bwt_vector_normal = bwt_vec;
		set_bwt_vector_normal();
		r_block_pe1 = r_b_pe1;
		r_block_pe2 = r_b_pe2;
		chrom_map = _chrom_map;
		bt_ptr = _bt;
		patsrc = _patsrc;
		sink = _sink;
		original_ebwt = _original_ebwt;
		original_ebwtBw = _original_ebwtBw;
		juncdb_ebwt = _juncdb_ebwt;
		juncdb_ebwtBw = _juncdb_ebwtBw;
		fusiondb_ebwt = _fusiondb_ebwt;
		fusiondb_ebwtBw = _fusiondb_ebwtBw;
		full_readbufa = _full_readbufa;
		full_readbufb = _full_readbufb;
		params_ptr = _params;
		reference_sequence = _reference_sequence;
		sam_fs = _sam_fs;
		outfs_lock = _outfs_lock;
		pe1_no_repeats = _pe1_no_repeats;
		pe2_no_repeats = _pe2_no_repeats;
		pe1_segment_no_repeats = _pe1_segment_no_repeats;
		pe2_segment_no_repeats = _pe2_segment_no_repeats;
		fusion_cluster = _fusion_cluster;
		fusion_fs = _fusion_fs;
		unmapped_fs1 = _unmapped_fs1;
		unmapped_fs2 = _unmapped_fs2;
		fusion_fs_lock = _fusion_fs_lock;
		mapping_stats = _mapping_stats;
	}

	void set_bwt_vector_normal()
	{
		bwt_vector = bwt_vector_normal;
	}

	void set_bwt_vector_repeats()
	{
		bwt_vector = &bwt_vector_repeats;
	}

	int check_contig(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)  //check 2 bwtmap segment are connected
	{
		if(seg1.pair_no != seg2.pair_no || seg1.chrom!=seg2.chrom || seg1.strand!=seg2.strand || seg2.start - seg1.end > 1)
			return CONTIG_NO_RELATIONSHIP;
		else if(seg1.end + 1 == seg2.start) //position related
		{
			if( (seg1.strand=="+" && seg2.start_seg_no == seg1.end_seg_no + 1) || (seg1.strand=="-"&&seg2.start_seg_no == seg1.end_seg_no - 1))     // strand +
			{
				return CONTIG_RELATIONSHIP;        // case 1: contig
			} 
		}
		else if (seg2.start - seg1.end < 1)  //need to continue search
			return CONTIG_CONTINUE_SEARCH;
		return CONTIG_NO_RELATIONSHIP;                 // case 0: no relationship
	}

	/*void combine_contig()      // combined connected segments
	{
	for (vector<Bwtmap_Info>::iterator it=(*bwt_vector).begin(); it + 1 != (*bwt_vector).end(); )
	{
	for(vector<Bwtmap_Info>::iterator it2= it + 1; ; it2++)
	{
	if(it2 == (*bwt_vector).end())
	{
	it++;
	break;
	}
	int contig_status= check_contig(*it, *(it2));
	if (contig_status == CONTIG_RELATIONSHIP)   //contig
	{
	(*it).end_seg_no = (*it2).end_seg_no;
	(*it).end = (*it2).end;
	for(size_t i = 0; i < (*it2).mismatches.size(); i++)
	{
	(*it).mismatches.push_back((*it2).mismatches[i]);
	}
	(*bwt_vector).erase(it2);
	break;
	}
	else if(contig_status == CONTIG_NO_RELATIONSHIP)
	{
	it++;
	break;
	}
	}
	}
	}*/

	void combine_contig()      // combined connected segments
	{
		for (vector<Bwtmap_Info>::iterator it=(*bwt_vector).begin(); it + 1 != (*bwt_vector).end(); )
		{
			for(vector<Bwtmap_Info>::iterator it2= it + 1; ; it2++)
			{
				if(it2 == (*bwt_vector).end())
				{
					it++;
					break;
				}
				int contig_status= check_contig(*it, *(it2));
				if (contig_status == CONTIG_RELATIONSHIP)   //contig
				{
					int segment_number = (it->pair_no == 1 ? pe1_segment_number : pe2_segment_number);
					Read_Block* r_block = (it->pair_no == 1 ? r_block_pe1 : r_block_pe2);
					if((it->strand == "+" && it->start_seg_no == HEAD_SEG_NO && it->end_seg_no == HEAD_SEG_NO && it->num_mis > 0)
						|| (it->strand == "-" && it->start_seg_no == segment_number && it->end_seg_no == segment_number && it->num_mis > 0))
					{
						Splice_Info new_head_splice;
						it->to_splice_info(new_head_splice, r_block->get_seg_len(it->start_seg_no));
						it->copy(*it2);
						it->splice_head.push_back(new_head_splice);
					}
					else if((it2->strand == "+" && it2->start_seg_no == segment_number && it2->end_seg_no == segment_number && it2->num_mis > 0)
						|| (it2->strand == "-" && it2->start_seg_no == HEAD_SEG_NO && it2->end_seg_no == HEAD_SEG_NO && it2->num_mis > 0))
					{
						Splice_Info new_tail_splice;
						it2->to_splice_info(new_tail_splice, r_block->get_seg_len(it2->start_seg_no));
						it->splice_tail.push_back(new_tail_splice);
					}
					else
					{
						(*it).end_seg_no = (*it2).end_seg_no;
						(*it).end = (*it2).end;
					}
					(*bwt_vector).erase(it2);
					break;
				}
				else if(contig_status == CONTIG_NO_RELATIONSHIP)
				{
					it++;
					break;
				}
			}
		}
	}

	int check_relation(const Bwtmap_Info& seg1, const Bwtmap_Info& seg2)
	{
		if(seg1.pair_no == seg2.pair_no && seg1.chrom == seg2.chrom && seg1.strand == seg2.strand)
		{
			int seg_gap = (abs(seg2.start_seg_no - seg1.end_seg_no) - 1) * segment_length;
			int ins_gap = seg_gap - max_ins;
			int del_gap = seg_gap + max_del;
			int min_gap = seg_gap + min_intron_length;
			int max_gap = seg_gap + max_intron_length_double_anchor;
			int dist = seg2.start - seg1.end - 1;
			if(dist < ins_gap)
				return FIX_TOO_CLOSE;
			else if(dist < seg_gap)
				return FIX_INSERTION; 								// insertion
			else if(dist <= del_gap)              
				return FIX_DOUBLE_ANCHOR;     	      // append and deletion
			else if(dist < min_gap)
				return FIX_TOO_CLOSE;  								// too close
			else if(dist <= max_gap)
				return FIX_DOUBLE_ANCHOR;        			// in range
			else 
				return FIX_TOO_FAR;        					  // too far
		}
		else
			return FIX_NO_RELATIONSHIP;            				//no relationship	
	}
	
	bool try_hard_insertion(string& pending_seq, string& pending_chrom_seq, size_t& prefix_length, size_t& insert_length, size_t& mismatch_bits, bool head)
	{
		bool insertion_fixed = false;
		size_t current_best_mismatch = pending_seq.length();
		int i;
		if(pending_seq.length() == pending_chrom_seq.length())
			i = 1;
		else
			i = 0;
		for(; (i + (int)(pending_seq.length() - pending_chrom_seq.length())) <= try_hard_max_ins; i++)
		{
			string chrom_seq;
			if(head)
				chrom_seq = pending_chrom_seq.substr(i);
			else
				chrom_seq = pending_chrom_seq.substr(0, pending_chrom_seq.length() - i);
			size_t this_mismatch;
			size_t this_prefix_length = 0;
			size_t this_mismatch_bits = 0;
			if(genome_scan.Double_anchored_score_ins(pending_seq, chrom_seq, try_hard_max_mis, this_prefix_length, this_mismatch_bits, this_mismatch))
			{
				if(head && this_prefix_length == 0)
					continue;
				else if(!head && this_prefix_length == chrom_seq.length())
					continue;
				else if(this_mismatch < current_best_mismatch)
				{
					current_best_mismatch = this_mismatch;
					prefix_length = this_prefix_length;
					mismatch_bits = this_mismatch_bits;
					insert_length = i + pending_seq.length() - pending_chrom_seq.length();
					insertion_fixed = true;
					if(current_best_mismatch == 0)
						break;
				}
			}
		}
		return insertion_fixed;
	}

	void generate_sam_head(size_t current_contig_index)  
	{                                             //has fix head
		for(size_t i=0;i<(*bwt_vector)[current_contig_index].splice_head.size();i++)
		{
			Sam_Info new_sam_info;
			new_sam_info.start_contig_index=current_contig_index;
			new_sam_info.chrom=(*bwt_vector)[current_contig_index].chrom;
			new_sam_info.pair_no = (*bwt_vector)[current_contig_index].pair_no;
			new_sam_info.strand = (*bwt_vector)[current_contig_index].strand;
			new_sam_info.start_seg_no = (*bwt_vector)[current_contig_index].splice_head[i].start_seg_no;
			if((*bwt_vector)[current_contig_index].splice_head[i].jump_code[0].len!=0)   // no 0M case
			{
				new_sam_info.start_pos=(*bwt_vector)[current_contig_index].splice_head[i].start_pos;
				for(size_t j=0;j<(*bwt_vector)[current_contig_index].splice_head[i].jump_code.size();j++)
				{
					new_sam_info.jump_code.push_back((*bwt_vector)[current_contig_index].splice_head[i].jump_code[j]);
				}
			}
			generate_sam_internal(new_sam_info,current_contig_index, false, 0, false);
		}
		if((*bwt_vector)[current_contig_index].splice_head.size() == 0 || (*bwt_vector)[current_contig_index].incomplete_head)  //no complete head fixed
		{
			Sam_Info new_sam_info;
			new_sam_info.start_contig_index=current_contig_index;
			new_sam_info.start_seg_no=(*bwt_vector)[current_contig_index].start_seg_no;
			new_sam_info.chrom=(*bwt_vector)[current_contig_index].chrom;
			new_sam_info.pair_no = (*bwt_vector)[current_contig_index].pair_no;
			new_sam_info.strand = (*bwt_vector)[current_contig_index].strand;
			new_sam_info.imcomplete_single_anchor = (*bwt_vector)[current_contig_index].incomplete_head;
			generate_sam_internal(new_sam_info,current_contig_index,false,0,false);
		} 
	}

	void generate_sam_internal(Sam_Info &my_sam_info, size_t current_contig_index, bool skip_edge_seg, int overlap_len, bool _has_special_case)
	{
		bool has_special_case = _has_special_case;
		Read_Block* r_block = ((*bwt_vector)[current_contig_index].pair_no == 1) ? r_block_pe1 : r_block_pe2;
		if(my_sam_info.start_pos == -1) // started without fix head
			my_sam_info.start_pos = (*bwt_vector)[current_contig_index].start;
		int match_len= r_block->get_seg_len((*bwt_vector)[current_contig_index].start_seg_no, (*bwt_vector)[current_contig_index].end_seg_no);
		if(my_sam_info.jump_code.size()==0)
		{
			Jump_Code new_jumpcode(match_len, "M");
			my_sam_info.jump_code.push_back(new_jumpcode);
		}
		else
		{
			match_len -= overlap_len;
			if(skip_edge_seg)
				match_len -= r_block->get_seg_len((*bwt_vector)[current_contig_index].start_seg_no);	
			if(my_sam_info.jump_code[my_sam_info.jump_code.size()-1].type == "M")	
				my_sam_info.jump_code[my_sam_info.jump_code.size()-1].len+=match_len;
			else
			{
				Jump_Code new_jumpcode(match_len, "M");
				my_sam_info.jump_code.push_back(new_jumpcode);
			}
		}
		if((*bwt_vector)[current_contig_index].splice_internal.size()==0)   //no more internal splice, generate tail and finish
		{
			generate_sam_tail(my_sam_info, current_contig_index); 
		}
		else      //more internal splice
		{
			for(size_t i = 0; i < (*bwt_vector)[current_contig_index].splice_internal.size(); i++)
			{
				bool next_skip_edge_seg=false;
				int next_overlap_len=0;
				int seg_gap=abs((*bwt_vector)[current_contig_index].end_seg_no - (*bwt_vector)[(*bwt_vector)[current_contig_index].splice_internal[i].end_contig].start_seg_no);
				if(seg_gap==2)  //gap hole
				{
					next_overlap_len=extend_len;
				}
				else if(seg_gap==1)  // continuous hole
				{
					next_overlap_len=trunc_len;
				}
				else if(seg_gap==0)
				{
					if(has_special_case)
						continue;
					next_overlap_len=extend_len;
					next_skip_edge_seg=true;
					has_special_case = true;
				}
				Sam_Info new_sam_info(my_sam_info);
				new_sam_info.jump_code[my_sam_info.jump_code.size()-1].len-=next_overlap_len;
				if(next_skip_edge_seg) // special case
					new_sam_info.jump_code[my_sam_info.jump_code.size()-1].len -= r_block->get_seg_len((*bwt_vector)[current_contig_index].end_seg_no);
				if((*bwt_vector)[current_contig_index].splice_internal[i].jump_code[0].type == "M")
					new_sam_info.jump_code[my_sam_info.jump_code.size()-1].len += (*bwt_vector)[current_contig_index].splice_internal[i].jump_code[0].len;
				else
					new_sam_info.jump_code.push_back((*bwt_vector)[current_contig_index].splice_internal[i].jump_code[0]);
				for(size_t j=1;j<(*bwt_vector)[current_contig_index].splice_internal[i].jump_code.size();j++)
				{
					new_sam_info.jump_code.push_back((*bwt_vector)[current_contig_index].splice_internal[i].jump_code[j]);
				}
				generate_sam_internal(new_sam_info,(*bwt_vector)[current_contig_index].splice_internal[i].end_contig,next_skip_edge_seg,next_overlap_len, has_special_case);
			}
		}
	}

	void generate_sam_tail(Sam_Info &my_sam_info, size_t current_contig_index)
	{            // has tail
		for(size_t i=0;i<(*bwt_vector)[current_contig_index].splice_tail.size();i++)
		{
			Sam_Info new_sam_info(my_sam_info);
			new_sam_info.end_seg_no = (*bwt_vector)[current_contig_index].splice_tail[i].end_seg_no;
			if((*bwt_vector)[current_contig_index].splice_tail[i].jump_code[0].type == "M")
				new_sam_info.jump_code[my_sam_info.jump_code.size()-1].len += (*bwt_vector)[current_contig_index].splice_tail[i].jump_code[0].len;
			else
				new_sam_info.jump_code.push_back((*bwt_vector)[current_contig_index].splice_tail[i].jump_code[0]);
			for(size_t j=1;j<(*bwt_vector)[current_contig_index].splice_tail[i].jump_code.size();j++) 
			{
				if(j+1<(*bwt_vector)[current_contig_index].splice_tail[i].jump_code.size() && (*bwt_vector)[current_contig_index].splice_tail[i].jump_code[j+1].type == "M" &&(*bwt_vector)[current_contig_index].splice_tail[i].jump_code[j+1].len==0)
				{ //0M case
					if((*bwt_vector)[current_contig_index].splice_tail[i].jump_code[j].type =="I")					
						new_sam_info.jump_code.push_back((*bwt_vector)[current_contig_index].splice_tail[i].jump_code[j]);
					break;
				}
				new_sam_info.jump_code.push_back((*bwt_vector)[current_contig_index].splice_tail[i].jump_code[j]);
			}
			print_saminfo(new_sam_info);
		}
		if((*bwt_vector)[current_contig_index].splice_tail.size() == 0 || (*bwt_vector)[current_contig_index].incomplete_tail)   // no complete tail fixed
		{
			my_sam_info.end_seg_no=(*bwt_vector)[current_contig_index].end_seg_no;
			my_sam_info.imcomplete_single_anchor = (*bwt_vector)[current_contig_index].incomplete_tail;
			print_saminfo(my_sam_info);
		}  
	}

	bool check_duplicate_jumpcode(vector<Jump_Code> &jump1, vector<Jump_Code> &jump2) 
	{
		if(jump1.size() != jump2.size())
			return false;
		for(size_t i=0;i<jump1.size();i++)
		{
			if(jump1[i].len !=jump2[i].len || jump1[i].type != jump2[i].type)
				return false;
		}
		return true;
	}

	bool check_duplicate_sam(vector<Sam_Info> &my_vector, Sam_Info &my_sam_info)    // check duplicate alignments when report
	{
		for(size_t i=0; i<my_vector.size();i++)
		{
			if(my_sam_info.start_pos == my_vector[i].start_pos && my_sam_info.chrom == my_vector[i].chrom && my_sam_info.pair_no == my_vector[i].pair_no && my_sam_info.strand == my_vector[i].strand && check_duplicate_jumpcode(my_sam_info.jump_code, my_vector[i].jump_code))
				return true;
		}
		return false;
	}

	bool check_duplicate_splice(vector<Splice_Info> &my_vector, Splice_Info &my_info)    // check duplicate alignments when report
	{
		for(size_t i=0; i<my_vector.size();i++)
		{
			if(check_duplicate_jumpcode(my_info.jump_code, my_vector[i].jump_code))
				return true;
		}
		return false;
	}

	void print_saminfo(Sam_Info &my_sam_info)
	{
		for(size_t i = 0; i < my_sam_info.jump_code.size(); i++)
		{
			if(my_sam_info.jump_code[i].len <= 0)
				return;
		}
		Read_Block* r_block = (my_sam_info.pair_no == 1 ? r_block_pe1 : r_block_pe2);
		int read_length = my_sam_info.pair_no == 1 ? pe1_read_length : pe2_read_length;
		my_sam_info.Get_Extra_Info((*reference_sequence)[my_sam_info.chrom], max_del, r_block);
		if(my_sam_info.mapped_len == read_length || (min_mapped_len > 0 && my_sam_info.mapped_len >= min_mapped_len))   // only output above threshold sam
		{
			size_t num_current_alignment;
			if(my_sam_info.pair_no == 1)
			{
				pe1_aligned = true;
				num_current_alignment = num_pe1_alignment;
				num_pe1_alignment ++;
			}
			else
			{
				pe2_aligned = true;
				num_current_alignment = num_pe2_alignment;
				num_pe2_alignment ++;
			}
			if(!(my_sam_info.imcomplete_single_anchor) && num_current_alignment < maximum_alignments && !check_duplicate_sam(sam_vec, my_sam_info))
			{
				sam_vec.push_back(my_sam_info);
			}
		}
		if(fusion &&  my_sam_info.mapped_len < read_length && !check_duplicate_sam(fusion_candidate, my_sam_info))
			fusion_candidate.push_back(my_sam_info);
	}

	void print_sam_info_to_file(Sam_Info &my_sam_info)    // report single read alignment
	{
		if(my_sam_info.intron_len > 0)
			(my_sam_info.pair_no == 1 ? pe1_spliced : pe2_spliced) = true;
		else
		  (my_sam_info.pair_no == 1 ? pe1_unspliced : pe2_unspliced) = true;
		if(my_sam_info.printed || (splice_only && my_sam_info.intron_len == 0))
			return;
		Read_Block* r_block = (my_sam_info.pair_no == 1 ? r_block_pe1 : r_block_pe2);
		int segment_number = (my_sam_info.pair_no == 1 ? pe1_segment_number : pe2_segment_number);
		int read_length = (my_sam_info.pair_no == 1 ? pe1_read_length : pe2_read_length);
		string quality_string;                                          // get quality string        
		if(debug)
			cout << "print sam info, start segment no is " << my_sam_info.start_seg_no << ", end segment no is " << my_sam_info.end_seg_no << endl;
		if(!fa)
		{
			if(my_sam_info.strand == "+")
			{
				for(int i = my_sam_info.start_seg_no; i <= (int)(my_sam_info.end_seg_no); i++)
				{
					quality_string.append(r_block->get_seg_qual(i));
				}
			}
			else    
			{
				for(int i = my_sam_info.start_seg_no; i >= (int)(my_sam_info.end_seg_no); i--)
				{
					quality_string.append(r_block->get_revcom_seg_qual(i));
				}
			}
		}
		else
			quality_string = global_qual_string.substr(0, r_block->get_seg_len((int)(my_sam_info.start_seg_no), (int)(my_sam_info.end_seg_no)));
		int quality_score = 255;     // default quality score for Fasta read
		if(!fa && quality_string != global_qual_string.substr(0,quality_string.size()) )
			quality_score = GetQualityScore(my_sam_info.mismatchs, quality_string, quality_scale);   // compute quality score
		size_t sam_flag;
		if(my_sam_info.strand == "+")
			sam_flag = 0;
		else
			sam_flag = 16;
		(*sam_fs) << r_block->read_id << "\t" << sam_flag << "\t" << my_sam_info.chrom << "\t" << (my_sam_info.start_pos + 1) << "\t"
			<< quality_score << "\t";
		if(my_sam_info.strand == "+" && my_sam_info.start_seg_no != HEAD_SEG_NO)
		{
			(*sam_fs) << r_block->get_seg_len(HEAD_SEG_NO, my_sam_info.start_seg_no - 1) << "S";
		}
		if(my_sam_info.strand == "-" && my_sam_info.start_seg_no != segment_number)
		{
			(*sam_fs) << r_block->get_seg_len(my_sam_info.start_seg_no + 1, segment_number) << "S";
		}
		for(size_t i = 0; i < my_sam_info.jump_code.size(); i++)    // output jump code
		{
			(*sam_fs) << my_sam_info.jump_code[i].len;
			if(my_sam_info.jump_code[i].len <= max_del && my_sam_info.jump_code[i].type == "N")
				(*sam_fs) << "D";
			else
				(*sam_fs) << my_sam_info.jump_code[i].type;
		}
		if(my_sam_info.strand == "+" && my_sam_info.end_seg_no != segment_number)
		{
			(*sam_fs) << r_block->get_seg_len(my_sam_info.end_seg_no + 1, segment_number) << "S";
		}
		if(my_sam_info.strand == "-" && my_sam_info.end_seg_no != HEAD_SEG_NO)
		{
			(*sam_fs) << r_block->get_seg_len(HEAD_SEG_NO, my_sam_info.end_seg_no - 1) << "S";
		}
		(*sam_fs) << "\t*\t0\t0\t";               
		if(my_sam_info.strand == "+")     //output sequence
		{
			for(int i = HEAD_SEG_NO; i <= segment_number; i++)
			{
				(*sam_fs) << r_block->get_seg_seq(i);
			}
		}
		else                   //out put - sequence
		{
			for(int i = segment_number; i >= HEAD_SEG_NO; i--)
			{
				(*sam_fs) << r_block->get_revcom_seg_seq(i);
			}
		}
		(*sam_fs) << "\t";
		if(fa)
			(*sam_fs) << global_qual_string.substr(0, read_length);
		else
		{
			if(my_sam_info.strand == "+")
			{
				for(int i = HEAD_SEG_NO; i <= segment_number; i++)
				{
					(*sam_fs) << r_block->get_seg_qual(i);
				}
			}
			else    
			{
				for(int i = segment_number; i >= HEAD_SEG_NO; i--)
				{
					(*sam_fs) << r_block->get_revcom_seg_qual(i);
				}
			}	
		}
		(*sam_fs) << "\tNM:i:" << my_sam_info.mismatchs.size();
		if(my_sam_info.forward_strand_junc && !(my_sam_info.reverse_strand_junc))
			(*sam_fs) << "\t" << "XS:A:+";
		else if(!(my_sam_info.forward_strand_junc) && my_sam_info.reverse_strand_junc)
			(*sam_fs) << "\t" << "XS:A:-";
		else if(my_sam_info.forward_strand_junc && my_sam_info.reverse_strand_junc) 
			(*sam_fs) << "\t" << "XS:A:N";
		else if(my_sam_info.unknown_strand_junc)
			(*sam_fs) << "\t" << "XS:A:*";
		else if(my_sam_info.deletion)
			(*sam_fs) << "\t" << "XS:A:D";
		if(my_sam_info.junc_flank_seq.size() > 0)
		{
			(*sam_fs) << "\t" << "XF:Z:";
			for(size_t i = 0; i < my_sam_info.junc_flank_seq.size(); i++)
			{
					(*sam_fs) << my_sam_info.junc_flank_seq[i] << ",";
			}
		}
		(*sam_fs) << endl;
		my_sam_info.printed = true;
	}

	void print_unmapped_read_to_file(ofstream* unmapped_fs, Read_Block* r_block)
	{
		int segment_number = r_block->get_seg_num();
		if(fa)
		{
			(*unmapped_fs) << ">";
		}
		else
		{
			(*unmapped_fs) << "@";
		}
		(*unmapped_fs) << r_block->read_id << endl;
		for(int i = HEAD_SEG_NO; i <= segment_number; i++)
		{
			(*unmapped_fs) << r_block->get_seg_seq(i);
		}	
		(*unmapped_fs) << endl;
		if(!fa)
		{
			(*unmapped_fs) << "+" << endl;
			for(int i = HEAD_SEG_NO; i <= segment_number; i++)
			{
				(*unmapped_fs) << r_block->get_seg_qual(i);
			}	
			(*unmapped_fs) << endl;
		}	
	}

	void output_unspliced_read()
	{
		for(size_t i = 0; i < (*bwt_vector).size(); i++)
		{
			if((*bwt_vector)[i].splice_head.size() == 0 && (*bwt_vector)[i].splice_tail.size() == 0 && (*bwt_vector)[i].splice_internal.size() == 0)
			{
				int segment_number = (*bwt_vector)[i].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
				Read_Block* r_block = (*bwt_vector)[i].pair_no == 1 ? r_block_pe1 : r_block_pe2;
				Sam_Info new_sam_info((*bwt_vector)[i], i, r_block, segment_number);
				print_saminfo(new_sam_info);	
			}
		}
	}

	bool remap_to_juncdb(size_t bwt_index1, size_t bwt_index2, string chrom, string strand, int read_start, int read_end, int premap_len, int left_bound, int right_bound, vector<Splice_Info>& fixed_splice_vec, bool head_0M, bool tail_0M)
	{
		return  remap_to_juncdb(bwt_index1, bwt_index2, chrom, strand, read_start, read_end, premap_len, left_bound, right_bound, -1, -1, fixed_splice_vec, head_0M, tail_0M);
	}

	bool remap_to_juncdb(size_t bwt_index1, size_t bwt_index2, string chrom, string strand, int read_start, int read_end, int premap_len, int left_bound, int right_bound, int left_most, int right_most, vector<Splice_Info>& fixed_splice_vec, bool head_0M, bool tail_0M)
	{
		if(debug)
			cout << "remap start " << read_start << " " << read_end <<endl;
		sink->set_max_hits(max_remap_hits);		
		ReadBuf*	full_readbuf;
		if((!pair_end) || (*bwt_vector)[bwt_index1].pair_no == 1)
			full_readbuf = full_readbufa;
		else
			full_readbuf = full_readbufb;	
		remap_number++;
		if(strand == "+")
		{
			norc = true;
			nofw = false;
		}
		else
		{
			nofw = true;
			norc = false;
		}
		Ebwt< String<Dna, Alloc<> > >& ebwtFw = *juncdb_ebwt;
		Ebwt< String<Dna, Alloc<> > >& ebwtBw = *juncdb_ebwtBw; 
		GreedyDFSRangeSource& bt = *bt_ptr;
		EbwtSearchParams<String<Dna> >& params = *params_ptr;
		patsrc->set_readbufa_remap(*full_readbuf, read_start, read_end, seed, remap_number);
		SET_PARAM_PTR(patsrc);
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		bool only_run_once = true;
		#define DONEMASK_SET(p)
		while(only_run_once)
		{
			only_run_once = false;
			#include "search_1mm_phase1.c"
			#include "search_1mm_phase2.c"
		}
		#undef DONEMASK_SET
		//cout << "converting remapped" << endl;		
		bool remapped = sink->convert_remap(bwt_index1, bwt_index2, premap_len, fixed_splice_vec, chrom, left_bound, right_bound, left_most, right_most, ebwtFw.refnames(), head_0M, tail_0M, debug);
		sink->finishRead(*patsrc, false, false);	
		//cout << "remap finish" << endl;
		return remapped;	
	}

	bool remap_to_fusiondb(size_t index1, size_t index2, string chrom1, string chrom2, string strand1, string strand2, int read_start, int read_end, int premap_len, int buffer_len, int left_bound, int right_bound, vector<Fusion_Splice>& fixed_splice_vec)
	{
		if(debug)
			cout << "remap fusion start " << read_start << " " << read_end <<endl;
		sink->set_max_hits(max_remap_hits);		
		ReadBuf*	full_readbuf;
		if((!pair_end) || fusion_candidate[index1].pair_no == 1)
			full_readbuf = full_readbufa;
		else
			full_readbuf = full_readbufb;	
		remap_number++;
		norc = false;
		nofw = false;
		Ebwt< String<Dna, Alloc<> > >& ebwtFw = *fusiondb_ebwt;
		Ebwt< String<Dna, Alloc<> > >& ebwtBw = *fusiondb_ebwtBw; 
		GreedyDFSRangeSource& bt = *bt_ptr;
		EbwtSearchParams<String<Dna> >& params = *params_ptr;
		patsrc->set_readbufa_remap(*full_readbuf, read_start, read_end, seed, remap_number);
		SET_PARAM_PTR(patsrc);
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		bool only_run_once = true;
		#define DONEMASK_SET(p)
		while(only_run_once)
		{
			only_run_once = false;
			#include "search_1mm_phase1.c"
			#include "search_1mm_phase2.c"
		}
		#undef DONEMASK_SET
		bool remapped = sink->convert_fusion_remap(index1, index2, premap_len, buffer_len, fixed_splice_vec, chrom1, chrom2, strand1, strand2, left_bound, right_bound, ebwtFw.refnames(), *reference_sequence, debug);
		sink->finishRead(*patsrc, false, false);	
		return remapped;	
	}

	bool fix_double_anchor(size_t bwt_index1, size_t bwt_index2, int read_start, int read_end, string pending_seq, string chrom, string strand, int left_bound, int right_bound, bool insert_jumpcode)
	{
		Splice_Info fixed_splice;
		return fix_double_anchor(bwt_index1, bwt_index2, read_start, read_end, pending_seq, chrom, strand, left_bound, right_bound, fixed_splice, insert_jumpcode);
	}

	bool fix_double_anchor(size_t bwt_index1, size_t bwt_index2, int read_start, int read_end, string pending_seq, string chrom, string strand, int left_bound, int right_bound, Splice_Info& fixed_splice, bool insert_jumpcode)
	{
		bool return_fixed = false;
		bool remap_fixed = false;
		bool anchor_fixed = false;
		bool append_fixed = false;
		int segment_number = (*bwt_vector)[bwt_index1].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
		Read_Block* r_block = (*bwt_vector)[bwt_index1].pair_no == 1 ? r_block_pe1 : r_block_pe2;
		int total_read_len = r_block->get_seg_len(HEAD_SEG_NO, segment_number);
		int premap_len = (strand == "+" ? read_start : total_read_len - read_end);
		if(right_bound - left_bound == read_end - read_start)
		{  //appending double anchor
			size_t mismatch_bits = 0;
			string chrom_seq = (*reference_sequence)[chrom].substr(left_bound, right_bound - left_bound + 1);
			if(debug)
				cout << "appending double anchor for: " << (*bwt_vector)[bwt_index1].read_id << "\t" << chrom << "\t" << strand << "\t" << left_bound << "\t" << right_bound << "\t" << pending_seq << "\t" << chrom_seq << endl;
			append_fixed = score_string(pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);
			if(append_fixed)
			{
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				if(!insert_jumpcode)
				{
					if(new_splice.jump_code.size() == 0)
					{
						Jump_Code prefix_match(pending_seq.length(), "M");
						new_splice.jump_code.push_back(prefix_match);
					}
					else
					{
						new_splice.jump_code[new_splice.jump_code.size() - 1].len += pending_seq.length();
					}						
				}
				else
				{
					if(new_splice.jump_code.size() == 0)
					{
						Jump_Code suffix_match(pending_seq.length(), "M");
						new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
					}
					else
					{
						new_splice.jump_code[0].len += (pending_seq.length());
					}	
				}
				if(!check_duplicate_splice((*bwt_vector)[bwt_index1].splice_internal, new_splice))
					(*bwt_vector)[bwt_index1].splice_internal.push_back(new_splice);
			}
		}
		if(juncdb)
		{
			vector<Splice_Info> remap_vector;
			remap_fixed = remap_to_juncdb(bwt_index1, bwt_index2, chrom, strand, read_start, read_end, premap_len, left_bound, right_bound, remap_vector, true, true);
			if(debug)
				cout << "double anchor remap to juncdb complete, size is " <<remap_vector.size() << endl;
			for(size_t i = 0; i < remap_vector.size(); i++)
			{
				if(debug)
					cout << "generate splice for " << i << remap_vector[i].jump_code.size() << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				if(debug)
					cout << "inserting jumpcode for " << i << endl;
				if(!insert_jumpcode)
				{
					if(new_splice.jump_code.size() == 0)
					{
						new_splice.jump_code.push_back(remap_vector[i].jump_code[0]);
					}
					else
					{
						new_splice.jump_code[new_splice.jump_code.size() - 1].len += remap_vector[i].jump_code[0].len;
					}
					for(size_t j = 1; j < remap_vector[i].jump_code.size(); j++)
					{
						new_splice.jump_code.push_back(remap_vector[i].jump_code[j]);
					}
				}
				else
				{
					if(new_splice.jump_code.size() == 0)
					{
						new_splice.jump_code.insert(new_splice.jump_code.begin(), remap_vector[i].jump_code[remap_vector[i].jump_code.size() - 1]);
					}
					else
					{
						new_splice.jump_code[0].len += remap_vector[i].jump_code[remap_vector[i].jump_code.size() - 1].len;
					}
					for(int j = (int)(remap_vector[i].jump_code.size()) - 2; j >= 0; j--)
					{
						new_splice.jump_code.insert(new_splice.jump_code.begin(), remap_vector[i].jump_code[j]);
					}	
				}
				if(!check_duplicate_splice((*bwt_vector)[bwt_index1].splice_internal, new_splice))
					(*bwt_vector)[bwt_index1].splice_internal.push_back(new_splice);
			}	
			remap_vector.clear();		
		}
		if(debug)
			cout << "fixing double anchor" << endl;
		//else //if(!remap_fixed)
		{
			size_t prefix_length = 0;
			size_t mismatch_bits = 0;
			string left_chrom_seq = (*reference_sequence)[chrom].substr(left_bound, pending_seq.length() + 2);
			string right_chrom_seq = (*reference_sequence)[chrom].substr(right_bound - pending_seq.length() - 1, pending_seq.length() + 2);
			string flank_seq;
			bool small_deletion = ( (right_bound - left_bound) - (read_end - read_start) <= max_del );
			bool adjacent_segments = ((*bwt_vector)[bwt_index1].strand == "+" && ((*bwt_vector)[bwt_index2].start_seg_no - (*bwt_vector)[bwt_index1].end_seg_no == 1)) || ((*bwt_vector)[bwt_index1].strand == "-" && ((*bwt_vector)[bwt_index1].end_seg_no - (*bwt_vector)[bwt_index2].start_seg_no == 1));
			if(debug)
				cout << "fixing double anchor for " << (*bwt_vector)[bwt_index1].read_id << "\t" << (*bwt_vector)[bwt_index1].start_seg_no << "\t" << (*bwt_vector)[bwt_index1].end_seg_no << "\t" << (*bwt_vector)[bwt_index1].start << "\t" << (*bwt_vector)[bwt_index1].end << "\t" << (*bwt_vector)[bwt_index2].start_seg_no << "\t" << (*bwt_vector)[bwt_index2].end_seg_no << "\t" << (*bwt_vector)[bwt_index2].start << "\t" << (*bwt_vector)[bwt_index2].end << "\t" << chrom << "\t" << left_bound << "\t" << right_bound <<"\t" << (*bwt_vector)[bwt_index1].strand << "\t" << pending_seq << "\t" << left_chrom_seq << "\t" << right_chrom_seq << "\t" << (small_deletion ? "small_deleltion" : "not_small_deletion") << endl;
			if(small_deletion)	
				anchor_fixed = genome_scan.Double_anchored_score_least_mis(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, mismatch_bits, small_deletion);
			else
			{
				string flank_seq;
				anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, mismatch_bits, (!adjacent_segments) && double_anchor_noncanonical, flank_seq);
				if(anchor_fixed && (flank_seq == "ATAC" || flank_seq == "GTAT" || flank_seq == "CTGC" || flank_seq == "GCAG") && ((int)(right_bound - left_bound - pending_seq.length() + 1) > max_semicanonical_intron_length))
					anchor_fixed = false;
			}
			if(anchor_fixed)
			{
				if(debug)
					cout << "double anchor fixed, prefix length is " << prefix_length << endl;		
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				string chrom_seq = left_chrom_seq.substr(0, prefix_length) + right_chrom_seq.substr(right_chrom_seq.size() - pending_seq.length() + prefix_length);
				if(!insert_jumpcode)
				{
					if(right_bound - left_bound  + 1 - pending_seq.length()== 0)
					{
						if(new_splice.jump_code.size() == 0)
						{
							Jump_Code prefix_match(pending_seq.length(), "M");
							new_splice.jump_code.push_back(prefix_match);
						}
						else
						{
							new_splice.jump_code[new_splice.jump_code.size() - 1].len += pending_seq.length();
						}						
					}
					else
					{
						if(new_splice.jump_code.size() == 0)
						{
							Jump_Code prefix_match(prefix_length, "M");
							new_splice.jump_code.push_back(prefix_match);
						}
						else
						{
							new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
						}
						Jump_Code skipped_length(right_bound - left_bound - pending_seq.length() + 1, "N");
						Jump_Code suffix_match(pending_seq.length() -  prefix_length, "M");
						new_splice.jump_code.push_back(skipped_length);
						new_splice.jump_code.push_back(suffix_match);
					}
				}
				else
				{
					if(right_bound - left_bound + 1 - pending_seq.length() == 0)
					{
						if(new_splice.jump_code.size() == 0)
						{
							Jump_Code suffix_match(pending_seq.length(), "M");
							new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
						}
						else
						{
							new_splice.jump_code[0].len += (pending_seq.length());
						}	
					}
					else
					{
						if(new_splice.jump_code.size() == 0)
						{
							Jump_Code suffix_match(pending_seq.length() -  prefix_length, "M");
							new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
						}
						else
						{
							new_splice.jump_code[0].len += (pending_seq.length() -  prefix_length);
						}
						Jump_Code skipped_length(right_bound - left_bound - pending_seq.length() + 1, "N");
						Jump_Code prefix_match(prefix_length, "M");
						new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
						new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
					}
				}
				if(!check_duplicate_splice((*bwt_vector)[bwt_index1].splice_internal, new_splice))
					(*bwt_vector)[bwt_index1].splice_internal.push_back(new_splice);
			}
		}
		return_fixed = append_fixed || remap_fixed || anchor_fixed;
		/*if(!return_fixed && do_local_align && (right_bound - left_bound <= local_align_max_dist))
		{
				Splice_Info new_splice1(fixed_splice);
				Splice_Info new_splice2(fixed_splice);
				string chrom_seq = (*reference_sequence)[chrom].substr(left_bound, right_bound - left_bound + 1);
				int forward_mismatch = 0;
				int backward_mismatch = 0;
				if(debug)
				{
					cout << "local_align" << endl;
					cout << "chrom seq is " << chrom_seq << endl;
					cout << "read seq is " << pending_seq << endl;
				}
				int forward_score = local_align.align(chrom_seq, pending_seq, "GT", "AG", true, true, new_splice1, forward_mismatch, left_bound, premap_len, insert_jumpcode, debug);
				int backward_score = local_align.align(chrom_seq, pending_seq, "CT", "AC", true, true, new_splice2, backward_mismatch, left_bound, premap_len, insert_jumpcode, debug);
				if(debug)
					cout << "forward score is " << forward_score << ", backward score is " << backward_score << endl; 
 				Splice_Info* best_splice = (forward_score >= backward_score ? &new_splice1 : &new_splice2);
 				int best_mismatch = (forward_score >= backward_score ? forward_mismatch : backward_mismatch);
				if(best_mismatch <= local_align_max_mis && best_splice->Count_Gap() <= local_align_max_gap)
				{
					if(debug)
					{
						cout << "local_align fixed";
						for(size_t i = 0; i < best_splice->jump_code.size(); i++)
							cout << best_splice->jump_code[i].toString();
						cout << endl;
					}
					best_splice->start_contig = bwt_index1;
					best_splice->end_contig = bwt_index2;
					if(!check_duplicate_splice((*bwt_vector)[bwt_index1].splice_internal, *best_splice))
						(*bwt_vector)[bwt_index1].splice_internal.push_back(*best_splice);
					return_fixed = true;
				}
		}*/
		return return_fixed;
	}

	bool fix_double_anchor_insertion(size_t bwt_index1, size_t bwt_index2, int read_start, int read_end, string pending_seq, string chrom, string strand, int left_bound, int right_bound, bool insert_jumpcode)
	{
		Splice_Info fixed_splice;
		return fix_double_anchor_insertion(bwt_index1, bwt_index2, read_start, read_end, pending_seq, chrom, strand, left_bound, right_bound, fixed_splice, insert_jumpcode);
	}

	bool fix_double_anchor_insertion(size_t bwt_index1, size_t bwt_index2, int read_start, int read_end, string pending_seq, string chrom, string strand, int left_bound, int right_bound, Splice_Info& fixed_splice, bool insert_jumpcode)
	{
		bool insertion_fixed = false;
		size_t prefix_length = 0;
		size_t mismatch_bits = 0;
		string chrom_seq = (*reference_sequence)[chrom].substr(left_bound, right_bound - left_bound + 1);
		if(debug)
			cout << "fixing double anchor insertion for " << (*bwt_vector)[bwt_index1].read_id << "\t" << chrom << "\t" << left_bound << "\t" <<  right_bound << "\t" << (*bwt_vector)[bwt_index1].strand << "\t"	<< pending_seq << "\t" << chrom_seq << endl;;
		insertion_fixed = genome_scan.Double_anchored_score_ins(pending_seq, chrom_seq, max_ins_mismatch,  prefix_length, mismatch_bits);
		if(insertion_fixed)
		{
			if(debug)
				cout << "insertion fixed, prefix length is " << prefix_length << endl;		
			Splice_Info new_splice(fixed_splice);
			new_splice.start_contig = bwt_index1;
			new_splice.end_contig = bwt_index2;
			string mapped_read_seq;
			if(prefix_length > 0)
				mapped_read_seq = pending_seq.substr(0, prefix_length);
			if(prefix_length < chrom_seq.length())
				mapped_read_seq += pending_seq.substr(pending_seq.length() - chrom_seq.length() + prefix_length);		
			if(!insert_jumpcode)
			{
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code prefix_match(prefix_length, "M");
					new_splice.jump_code.push_back(prefix_match);
				}
				else
				{
					new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
				}
				Jump_Code skipped_length(pending_seq.length() - chrom_seq.length(), "I");
				Jump_Code suffix_match(chrom_seq.length() - prefix_length, "M");
				new_splice.jump_code.push_back(skipped_length);
				new_splice.jump_code.push_back(suffix_match);
			}
			else
			{
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code suffix_match(chrom_seq.length() - prefix_length, "M");
					new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
				}
				else
				{
					new_splice.jump_code[0].len += (chrom_seq.length() - prefix_length);
				}
				Jump_Code skipped_length(pending_seq.length() - chrom_seq.length(), "I");
				Jump_Code prefix_match(prefix_length, "M");
				new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
				new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
			}
			(*bwt_vector)[bwt_index1].splice_internal.push_back(new_splice);
		}
		return insertion_fixed;
	}

	void remap_repeats_segment(size_t repeat_index)
	{
		int mapping_pair_no = (sam_vec[repeat_index].pair_no == 1 ? 2 : 1);
		string mapping_strand = sam_vec[repeat_index].strand == "+" ? "-" : "+";
		Read_Block* r_block = mapping_pair_no == 1 ? r_block_pe1 : r_block_pe2;
		for(int i = 0; i < r_block->get_seg_num(); i++)
		{
			if(mapping_strand == "+")
			{
				int left_bound = max(0, sam_vec[repeat_index].start_pos - max_intron_length_double_anchor);
				int right_bound = max(0, sam_vec[repeat_index].end_pos);
				if(r_block->get_seg_len(i + 1) <= 32)
				{
					string segment_string = r_block->get_seg_seq(i + 1);
					if(debug)
						cout << "setting repeats segment" << (i+1) << endl;
					if(anchor_search_head.set(segment_string, (*reference_sequence)[sam_vec[repeat_index].chrom], mapping_strand == "+", mapping_strand == "-", max_repeat_hits))
					{
						if(debug)
							cout << "mapping repeats segment" << (i+1) << " " << segment_string << " " << left_bound << "\t" << right_bound <<endl;
						anchor_search_head.find_all_match((*reference_sequence)[sam_vec[repeat_index].chrom], left_bound, right_bound, *bwt_vector, r_block->get_read_id(), sam_vec[repeat_index].chrom, mapping_pair_no, i + 1);
					}
				}
			}
			else
			{
				int left_bound = min((int)((*chrom_map)[sam_vec[repeat_index].chrom] - 1), sam_vec[repeat_index].start_pos);
				int right_bound = min((int)((*chrom_map)[sam_vec[repeat_index].chrom] - 1), sam_vec[repeat_index].end_pos + max_intron_length_double_anchor);
				if(r_block->get_seg_len(i + 1) <= 32)
				{
					string segment_string = r_block->get_seg_seq(i + 1);
					if(debug)
						cout << "setting repeats segment tail" << (i+1) << endl;
					if(anchor_search_tail.set(segment_string, (*reference_sequence)[sam_vec[repeat_index].chrom], mapping_strand == "+", mapping_strand == "-", max_repeat_hits))
					{
						if(debug)
							cout << "mapping repeats segment tail" << (i+1) << " " << segment_string << " " << left_bound << "\t" << right_bound <<endl;
						anchor_search_tail.find_all_match((*reference_sequence)[sam_vec[repeat_index].chrom], left_bound, right_bound, *bwt_vector, r_block->get_read_id(), sam_vec[repeat_index].chrom, mapping_pair_no, i + 1);
					}
				}	
			}
		}
	}

	bool remap_repeats_full_read(int remap_pair_no)
	{
		int segment_number = remap_pair_no == 1 ? pe1_segment_number : pe2_segment_number;
		sink->set_max_hits(max_repeat_all_hits);
		ReadBuf*	full_readbuf;	
		if(remap_pair_no == 1)
			full_readbuf = full_readbufa;
		else
			full_readbuf = full_readbufb;	
		norc = false;
		nofw = false;
		GreedyDFSRangeSource& bt = *bt_ptr;
		EbwtSearchParams<String<Dna> >& params = *params_ptr;
		Ebwt<String<Dna> >& ebwtFw = *original_ebwt;
		Ebwt<String<Dna> >& ebwtBw = *original_ebwtBw;
		patsrc->set_readbufa(*full_readbuf);
		SET_PARAM_PTR(patsrc);
		uint32_t plen = length(patFw);
		uint32_t s = plen;
		uint32_t s3 = s >> 1; // length of 3' half of seed
		uint32_t s5 = (s >> 1) + (s & 1); // length of 5' half of seed
		bool continue_run = true;	
		#define DONEMASK_SET(p)
		while(continue_run)
		{
			continue_run = false;
		#include "search_1mm_phase1.c"
		#include "search_1mm_phase2.c"
		}
		#undef DONEMASK_SET
		bool full_read_remapped = sink->convert_full_read(*bwt_vector, pair_end, remap_pair_no, HEAD_SEG_NO, segment_number, debug);
		sink->finishRead(*patsrc, debug, false);
		return full_read_remapped;
	}

	bool fix_head(size_t bwt_index, size_t anchor_segment, Splice_Info& fixed_splice, string chrom, int pos, string strand, bool only_remap)
	{
		if(debug)
			cout << "fix head start" << endl;
		bool return_fixed = false;
		bool remap_fixed = false;
		bool append_fixed = false;
		bool anchor_fixed = false;
		int segment_number = (*bwt_vector)[bwt_index].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
		Read_Block* r_block = (*bwt_vector)[bwt_index].pair_no == 1 ? r_block_pe1 : r_block_pe2;
		vector<bool>* segment_no_repeats = (*bwt_vector)[bwt_index].pair_no == 1 ? pe1_segment_no_repeats : pe2_segment_no_repeats;
		int anchor_segment_length = r_block->get_seg_len(anchor_segment);
		size_t premap_len = (strand == "+" ? r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment - 1) : r_block->get_seg_len(segment_number, (int)anchor_segment + 1));
		int min_previous_len;
		if(strand == "+")
			min_previous_len = r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment);
		else
			min_previous_len = r_block->get_seg_len(segment_number, (int)anchor_segment);
		if(pos < min_previous_len)  // not enough previous length
			return return_fixed;
		vector<Splice_Info> candidate_splice_vec;
		string pending_seq;
		if(strand == "+")
			pending_seq = (*r_block).get_seg_seq(anchor_segment);
		else
			pending_seq = (*r_block).get_revcom_seg_seq(anchor_segment);
		string pending_seq_tmp = pending_seq;
		if(!only_remap)
		{     																						 // append
			if(debug)
				cout << chrom << "\t" << pos << "\t" << pending_seq.length() << endl;
			size_t mismatch_bits = 0;
			string chrom_seq = (*reference_sequence)[chrom].substr(pos - pending_seq.length(), pending_seq.length());
			if(debug)
				cout << "appending head for: " << (*bwt_vector)[bwt_index].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << "\t" << pending_seq << "\t" << chrom_seq << endl;
			append_fixed = score_string(pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);
			if(append_fixed)
			{
				if(debug)
					cout << "append head fixed" << endl;	
				Splice_Info new_splice(fixed_splice);
				new_splice.start_pos = pos - (int)pending_seq.length();
				new_splice.start_seg_no = anchor_segment;
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code prefix_match(anchor_segment_length, "M");
					new_splice.jump_code.push_back(prefix_match);
				}
				else
				{
					new_splice.jump_code[0].len += anchor_segment_length;
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
		}
		if(juncdb)
		{                                                  // remap
			if(debug)
				cout << "remapping head" << endl;	
			vector<Splice_Info> remap_vector;
			remap_fixed = remap_to_juncdb(bwt_index, bwt_index, chrom, strand, r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment -1), r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment) - 1, premap_len, -1, pos - 1, remap_vector, false, true);
			for(size_t i = 0; i < remap_vector.size(); i++)
			{
				if(debug)
					cout << "remap head fixed: " << i << endl;	
				Splice_Info new_splice(fixed_splice);
				new_splice.start_pos = remap_vector[i].start_pos;
				new_splice.start_seg_no = anchor_segment;
				if(new_splice.jump_code.size() == 0)
				{
					new_splice.jump_code.insert(new_splice.jump_code.begin(), remap_vector[i].jump_code[remap_vector[i].jump_code.size() - 1]);
				}
				else
				{
					new_splice.jump_code[0].len += remap_vector[i].jump_code[remap_vector[i].jump_code.size() - 1].len;
				}
				for(int j = (int)(remap_vector[i].jump_code.size()) - 2; j >= 0; j--)
				{
					new_splice.jump_code.insert(new_splice.jump_code.begin(), remap_vector[i].jump_code[j]);
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
			remap_vector.clear();			
		}
		if(!only_remap && try_hard_ins)
		{
			size_t prefix_length = 0;
			size_t insert_length = 0;
			size_t mismatch_bits = 0;
			int pending_chrom_seq_start = max(0, pos - anchor_segment_length);
			string pending_chrom_seq = (*reference_sequence)[chrom].substr(pending_chrom_seq_start, pos - pending_chrom_seq_start);
			if(try_hard_insertion(pending_seq, pending_chrom_seq, prefix_length, insert_length, mismatch_bits, true))
			{
				if(debug)
					cout << "try hard insertion fixed" << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_pos = pos - anchor_segment_length + insert_length;
				new_splice.start_seg_no = anchor_segment;
				string chrom_seq = pending_chrom_seq.substr(insert_length - (pending_seq.length() - pending_chrom_seq.length()) );
				string mapped_read_seq;
				if(prefix_length > 0)
					mapped_read_seq = pending_seq.substr(0, prefix_length);
				if(prefix_length < chrom_seq.length())
					mapped_read_seq += pending_seq.substr(insert_length + prefix_length);
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code suffix_match(chrom_seq.length() -  prefix_length, "M");
					new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
				}
				else
				{
					new_splice.jump_code[0].len += (chrom_seq.length() -  prefix_length);
				}
				Jump_Code skipped_length(insert_length, "I");					
				new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
				Jump_Code prefix_match(prefix_length, "M");
				new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
		}
		if(!only_remap)//else if(false/*!remap_fixed && !append_fixed*/)    						 //fix head
		{                               
			pending_seq = pending_seq.substr(anchor_length - extend_len);
			string anchor_string;
			if(strand == "+")
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(0, anchor_length);
			else
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(r_block->get_seg_len(anchor_segment) - anchor_length);
			if(debug)                              
				cout << "fixing head for: " << (*bwt_vector)[bwt_index].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << ""<<endl;
			int anchor_pos;
			bool mapped_fw;
			int left_bound = max(0, pos - max_intron_length_single_anchor - anchor_length);
      //int right_bound_forward = min( (int)((*chrom_map)[chrom]) - 1, pos - anchor_segment_length + max_ins + anchor_length - 1);
      int right_bound_forward = min( (int)((*chrom_map)[chrom]) - 1, pos - anchor_length + extend_len - 1);
			int right_bound_reverse = right_bound_forward;
			if((*segment_no_repeats)[anchor_segment - 1] && anchor_search_head.set(anchor_string, (*reference_sequence)[chrom], strand == "+", strand == "-", max_anchor_hits))
			{
				if(debug)
					cout << "searching anchor " << anchor_string << "\t" << left_bound << "\t" << right_bound_forward << "\t" << right_bound_reverse << endl;
				while(anchor_search_head.search_next_simple((*reference_sequence)[chrom], left_bound, right_bound_forward, right_bound_reverse, anchor_pos, mapped_fw))
				{
					if(debug)
						cout << "found head anchor at : " << anchor_pos << "\t" << anchor_string << "\t" << (*reference_sequence)[chrom].substr(anchor_pos, anchor_length) << endl; 
					size_t prefix_length = 0;
					size_t mismatch_bits = 0;		
					int dist = pos - anchor_pos;
					int gap = dist - anchor_segment_length;					
					if(gap == 0) 							// already did append, skip
						continue;
					else if(gap > max_del && gap < min_intron_length)
						continue;																	
					else if(gap < 0)  					// fix insertion
					{
						if(abs(gap) > max_ins)
							continue;
						bool insertion_fixed = false;
						int chorm_seq_start = anchor_pos + anchor_length - extend_len;
						string chrom_seq = (*reference_sequence)[chrom].substr(chorm_seq_start, pos - chorm_seq_start);
						if(debug)
							cout << "fixing insertion for " << (*bwt_vector)[bwt_index].read_id << "\t" << chrom << "\t" <<  (*bwt_vector)[bwt_index].strand << "\t" << pending_seq << "\t" << chrom_seq << endl;
						insertion_fixed = genome_scan.Double_anchored_score_ins(pending_seq, chrom_seq, max_ins_mismatch,  prefix_length, mismatch_bits);
						if(insertion_fixed)
						{
							if(debug)
								cout << "insertion fixed" << endl;
							Splice_Info new_splice(fixed_splice);
							new_splice.start_pos = anchor_pos;
							new_splice.start_seg_no = anchor_segment;
							string mapped_read_seq;
							if(prefix_length > 0)
								mapped_read_seq = pending_seq.substr(0, prefix_length);
							if(prefix_length < chrom_seq.length())
								mapped_read_seq += pending_seq.substr(pending_seq.length() - chrom_seq.length() + prefix_length);
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code suffix_match(chrom_seq.length() -  prefix_length, "M");
								new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
							}
							else
							{
								new_splice.jump_code[0].len += (chrom_seq.length() -  prefix_length);
							}
							Jump_Code skipped_length(pending_seq.length() - chrom_seq.length(), "I");					
							new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
							Jump_Code prefix_match(prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
						}
					}
					else                                            // fix splice
					{
						int donor_start = anchor_pos + anchor_length - extend_len;
						string left_chrom_seq = (*reference_sequence)[chrom].substr(donor_start, pending_seq.length() + 2);
						int acceptor_start;
						string right_chrom_seq;
						if(pos >= (int)(pending_seq.length()) + 2)
						{
							acceptor_start = pos - pending_seq.length() - 2;
							right_chrom_seq = (*reference_sequence)[chrom].substr(acceptor_start, pending_seq.length() + 2);
						}
						else
						{
							acceptor_start = pos - pending_seq.length();
							right_chrom_seq = "AA" + (*reference_sequence)[chrom].substr(acceptor_start, pending_seq.length());
						}
						bool small_deletion = (gap <= max_del);
						if(debug)
							cout << anchor_pos << "\t"<< pending_seq <<"\t" << left_chrom_seq << "\t" << right_chrom_seq << "\t" << (small_deletion ? "small_deletion" : "not_small_deletion") <<endl;
						if(small_deletion)
							anchor_fixed = genome_scan.Double_anchored_score_least_mis(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, small_deletion);
						else
						{
							string flank_seq;
							anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, single_anchor_noncanonical, flank_seq);
							if(anchor_fixed && (flank_seq == "ATAC" || flank_seq == "GTAT" || flank_seq == "CTGC" || flank_seq == "GCAG") && ((gap > max_semicanonical_intron_length) || ((int)(prefix_length + anchor_length - extend_len) < min_semicanonical_anchor_length)))
								anchor_fixed = false;
						}
						if(anchor_fixed)
						{
							if(debug)
								cout << "anchor fixed" << endl;
							Splice_Info new_splice(fixed_splice);
							new_splice.start_pos = anchor_pos;
							new_splice.start_seg_no = anchor_segment;
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code suffix_match(pending_seq.length() -  prefix_length, "M");
								new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
							}
							else
							{
								new_splice.jump_code[0].len += (pending_seq.length() -  prefix_length);
							}
							Jump_Code skipped_length(gap, "N");					
							Jump_Code prefix_match(prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
							new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
							break;
						}
					}
				}
			}
		}
		/*if(candidate_splice_vec.size() == 0 && do_local_align && (abs((*bwt_vector)[bwt_index].end_seg_no - (*bwt_vector)[bwt_index].start_seg_no) >= 1))
		{
				Splice_Info new_splice1(fixed_splice);
				Splice_Info new_splice2(fixed_splice);
				int left_bound = max(0, pos - local_align_max_dist);
				string chrom_seq = (*reference_sequence)[chrom].substr(left_bound, pos - left_bound);
				int forward_mismatch = 0;
				int backward_mismatch = 0;
				if(debug)
				{
					cout << "local_align" << endl;
					cout << "chrom seq is " << chrom_seq << endl;
					cout << "read seq is " << pending_seq << endl;
				}
				int forward_score = local_align.align(chrom_seq, pending_seq_tmp, "GT", "AG", false, true, new_splice1, forward_mismatch, left_bound, premap_len, true, debug);
				int backward_score = local_align.align(chrom_seq, pending_seq_tmp, "CT", "AC", false, true, new_splice2, backward_mismatch, left_bound, premap_len, true, debug);
				if(debug)
					cout << "forward score is " << forward_score << ", backward score is " << backward_score << endl; 
 				Splice_Info* best_splice = (forward_score >= backward_score ? &new_splice1 : &new_splice2);
 				int best_mismatch = (forward_score >= backward_score ? forward_mismatch : backward_mismatch);
				if(best_mismatch <= local_align_max_mis && best_splice->Count_Gap() <= local_align_max_gap)
				{
					best_splice->start_seg_no = anchor_segment;
					if(debug)
					{
						cout << "local_align fixed";
						for(size_t i = 0; i < best_splice->jump_code.size(); i++)
							cout << best_splice->jump_code[i].toString();
						cout << endl;
					}
					if(!check_duplicate_splice(candidate_splice_vec, *best_splice))
						candidate_splice_vec.push_back(*best_splice);
				}
		}*/
		for(size_t i = 0; i < candidate_splice_vec.size(); i++)
		{
			return_fixed = true;
			if((strand == "+" && anchor_segment == 1) || (strand == "-" && (int)anchor_segment == segment_number))
			{
				(*bwt_vector)[bwt_index].splice_head.push_back(candidate_splice_vec[i]);
				if(debug)
					cout << "head fixed for: " << (*bwt_vector)[bwt_index].read_id << endl;
			}
			else if(!fix_head(bwt_index, anchor_segment + (strand == "+"? -1 : 1), candidate_splice_vec[i], chrom, candidate_splice_vec[i].start_pos, strand, only_remap))
			{
				(*bwt_vector)[bwt_index].splice_head.push_back(candidate_splice_vec[i]);
				(*bwt_vector)[bwt_index].incomplete_head = true;
			}
		}
		candidate_splice_vec.clear();
		if(debug)
			cout << "fix head end" << endl;
		return return_fixed;
	}

	bool fix_tail(size_t bwt_index, size_t anchor_segment, Splice_Info& fixed_splice, string chrom, int pos, string strand, bool only_remap)
	{
		if(debug)
			cout << "fix tail start" << endl;
		bool return_fixed = false;
		bool remap_fixed = false;
		bool append_fixed = false;
		bool anchor_fixed = false;
		int segment_number = (*bwt_vector)[bwt_index].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
		Read_Block* r_block = (*bwt_vector)[bwt_index].pair_no == 1 ? r_block_pe1 : r_block_pe2;
		vector<bool>* segment_no_repeats = (*bwt_vector)[bwt_index].pair_no == 1 ? pe1_segment_no_repeats : pe2_segment_no_repeats;
		int anchor_segment_length = r_block->get_seg_len(anchor_segment);
		size_t premap_len = (strand == "+" ? r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment - 1) : r_block->get_seg_len(segment_number, (int)anchor_segment + 1));
		int min_next_len;
		if(strand == "+")
			min_next_len = r_block->get_seg_len(segment_number , (int)anchor_segment);
		else
			min_next_len = r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment);	
		if(pos + min_next_len >= (int)((*chrom_map)[chrom]))
			return return_fixed;
		vector<Splice_Info> candidate_splice_vec;
		string pending_seq;
		if(strand == "+")
			pending_seq = (*r_block).get_seg_seq(anchor_segment);
		else
			pending_seq = (*r_block).get_revcom_seg_seq(anchor_segment);
		string pending_seq_tmp = pending_seq;
		if(!only_remap)//if(!remap_fixed)
		{
			size_t mismatch_bits = 0;
			string chrom_seq = (*reference_sequence)[chrom].substr(pos + 1, anchor_segment_length);
			if(debug)
				cout << "appending tail for: " << (*bwt_vector)[bwt_index].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << "\t" << pending_seq << "\t" << chrom_seq << endl;
			append_fixed = score_string(pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);
			if(append_fixed)
			{
				if(debug)
					cout << "append tail fixed" << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_pos = pos + anchor_segment_length;
				new_splice.end_seg_no = anchor_segment;
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code prefix_match(anchor_segment_length, "M");
					new_splice.jump_code.push_back(prefix_match);
				}
				else
				{
					new_splice.jump_code[new_splice.jump_code.size() - 1].len += anchor_segment_length;
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
		}
		if(juncdb)
		{
			if(debug)
				cout << "remapping tail" << endl;
			vector<Splice_Info> remap_vector;
			remap_fixed = remap_to_juncdb(bwt_index, bwt_index, chrom, strand, r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment - 1), r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment) - 1, premap_len, pos + 1, -1, remap_vector, true, false);
			for(size_t i = 0; i < remap_vector.size(); i++)
			{
				if(debug)
					cout << "remap tail fixed" << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_pos = remap_vector[i].end_pos;
				new_splice.end_seg_no = anchor_segment;
				if(new_splice.jump_code.size() == 0)
				{
					new_splice.jump_code.push_back(remap_vector[i].jump_code[0]);
				}
				else
				{
					new_splice.jump_code[new_splice.jump_code.size() - 1].len += remap_vector[i].jump_code[0].len;
				}
				for(size_t j = 1; j < remap_vector[i].jump_code.size(); j++)
				{
					new_splice.jump_code.push_back(remap_vector[i].jump_code[j]);
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}	
			remap_vector.clear();		
		}		
		if(!only_remap && try_hard_ins)
		{
			size_t prefix_length = 0;
			size_t insert_length = 0;
			size_t mismatch_bits = 0;
			int pending_chrom_seq_end = min(pos + anchor_segment_length, (int)((*chrom_map)[chrom]) );
			string pending_chrom_seq = (*reference_sequence)[chrom].substr(pos + 1, pending_chrom_seq_end - pos);
			if(debug)
					cout << "try hard insertion" << endl;
			if(try_hard_insertion(pending_seq, pending_chrom_seq, prefix_length, insert_length, mismatch_bits, false))
			{
				if(debug)
					cout << "try hard insertion fixed" << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_pos = pos + anchor_segment_length - insert_length;
				new_splice.end_seg_no = anchor_segment;
				string chrom_seq = pending_chrom_seq.substr(0, anchor_segment_length - insert_length);
				string mapped_read_seq;
				if(prefix_length > 0)
					mapped_read_seq = pending_seq.substr(0, prefix_length);
				if(prefix_length < chrom_seq.length())
					mapped_read_seq += pending_seq.substr(pending_seq.length() - chrom_seq.length() + prefix_length);
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code prefix_match(prefix_length, "M");
					new_splice.jump_code.push_back(prefix_match);
				}
				else
				{
					new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
				}
				Jump_Code skipped_length(insert_length, "I");
				Jump_Code suffix_match(chrom_seq.length() - prefix_length, "M");
				new_splice.jump_code.push_back(skipped_length);
				new_splice.jump_code.push_back(suffix_match);
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
		}
		if(!only_remap)//else if(false/*!remap_fixed && !append_fixed*/)
		{
			pending_seq = pending_seq.substr(0, anchor_segment_length - anchor_length + extend_len);
			string anchor_string;
			if(strand == "+")
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(r_block->get_seg_len(anchor_segment) - anchor_length);
			else
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(0, anchor_length); 
			if(debug)
				cout << "searching tail anchor for: " << (*bwt_vector)[bwt_index].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << "\t" << anchor_string << endl;
			int anchor_pos;
			bool mapped_fw;
			//int left_bound_forward = max(0, pos + 1 + anchor_segment_length - anchor_length - max_ins);
			int left_bound_forward = max(0, pos - extend_len + 2);
			int left_bound_reverse = left_bound_forward;
			int right_bound = min( (int)((*chrom_map)[chrom]) - 1, pos + max_intron_length_single_anchor + anchor_segment_length - anchor_length);
			if((*segment_no_repeats)[anchor_segment - 1] && anchor_search_tail.set(anchor_string, (*reference_sequence)[chrom], strand == "+", strand == "-", max_anchor_hits))
			{
				if(debug)
					cout << "boundary is " << left_bound_forward << "\t" << right_bound << endl;
				while(anchor_search_tail.search_next_simple((*reference_sequence)[chrom], left_bound_forward, left_bound_reverse, right_bound, anchor_pos, mapped_fw))
				{
					if(debug)
						cout << "found tail anchor at : " << anchor_pos << "\t" << anchor_string << "\t" << (*reference_sequence)[chrom].substr(anchor_pos, anchor_length) << endl; 
					size_t prefix_length = 0;
					size_t mismatch_bits = 0;		
					int dist = anchor_pos + anchor_length - 1 - pos;
					int gap = dist - anchor_segment_length;		
					if(gap == 0)    
						continue;
					else if(gap > max_del && gap < min_intron_length)
						continue;		
					if(gap < 0)  						 // fix insertion
					{
						if(abs(gap) > max_ins)
							continue;
						bool insertion_fixed = false;
						string chrom_seq = (*reference_sequence)[chrom].substr(pos + 1, anchor_pos - pos - 1 + extend_len);
						if(debug)
							cout << "fixing insertion for " << (*bwt_vector)[bwt_index].read_id << "\t" << chrom << "\t" <<  (*bwt_vector)[bwt_index].strand << "\t" << pending_seq << "\t" << chrom_seq << endl;
						insertion_fixed = genome_scan.Double_anchored_score_ins(pending_seq, chrom_seq, max_ins_mismatch,  prefix_length, mismatch_bits);
						if(insertion_fixed)
						{
							if(debug)
								cout << "insertion fixed" << endl;
							Splice_Info new_splice(fixed_splice);
							new_splice.start_pos = anchor_pos + anchor_length - 1;
							new_splice.end_seg_no = anchor_segment;
							string mapped_read_seq;
							if(prefix_length > 0)
								mapped_read_seq = pending_seq.substr(0, prefix_length);
							if(prefix_length < chrom_seq.length())
								mapped_read_seq += pending_seq.substr(pending_seq.length() - chrom_seq.length() + prefix_length);
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code prefix_match(prefix_length, "M");
								new_splice.jump_code.push_back(prefix_match);
							}
							else
							{
								new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
							}
							Jump_Code skipped_length(pending_seq.length() - chrom_seq.length(), "I");
							Jump_Code suffix_match(chrom_seq.length() - prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.push_back(skipped_length);
							new_splice.jump_code.push_back(suffix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
						}
					}
					else                                           // fix splice
					{
						int donor_start = pos + 1;;
						int acceptor_start = anchor_pos + extend_len - pending_seq.length() - 2;
						string left_chrom_seq;
						string right_chrom_seq = (*reference_sequence)[chrom].substr(acceptor_start, pending_seq.length() + 2);
						if((int)(donor_start + pending_seq.length() + 1 < (*chrom_map)[chrom]))
							left_chrom_seq = (*reference_sequence)[chrom].substr(donor_start, pending_seq.length() + 2);
						else
							left_chrom_seq = (*reference_sequence)[chrom].substr(donor_start, pending_seq.length()) + "AA";
						bool small_deletion = (gap <= max_del);
						if(debug)
							cout << "fixing splice " << anchor_pos << "\t"<< pending_seq <<"\t" << left_chrom_seq << "\t" << right_chrom_seq << "\t" << (small_deletion ? "small_deletion" : "not_small_deletion") <<endl;
						if(small_deletion)
							anchor_fixed = genome_scan.Double_anchored_score_least_mis(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, small_deletion);
						else
						{
							string flank_seq;
							anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, single_anchor_noncanonical, flank_seq);
							if(anchor_fixed && (flank_seq == "ATAC" || flank_seq == "GTAT" || flank_seq == "CTGC" || flank_seq == "GCAG") && ((gap > max_semicanonical_intron_length) || ((int)(pending_seq.length() - prefix_length + anchor_length - extend_len) < min_semicanonical_anchor_length)))
								anchor_fixed = false;
						}
						if(anchor_fixed)
						{
							Splice_Info new_splice(fixed_splice);
							new_splice.start_pos = anchor_pos + anchor_length - 1;
							new_splice.end_seg_no = anchor_segment;
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code prefix_match(prefix_length, "M");
								new_splice.jump_code.push_back(prefix_match);
							}
							else
							{
								new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
							}
							Jump_Code skipped_length(gap, "N");
							Jump_Code suffix_match(pending_seq.length() - prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.push_back(skipped_length);
							new_splice.jump_code.push_back(suffix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
							break;
						}
					}
				}
			}
		}
		/*if(candidate_splice_vec.size() == 0 && do_local_align && (abs((*bwt_vector)[bwt_index].end_seg_no - (*bwt_vector)[bwt_index].start_seg_no) >= 1))
		{
				Splice_Info new_splice1(fixed_splice);
				Splice_Info new_splice2(fixed_splice);
				int right_bound = min((int)((*chrom_map)[chrom]) - 1, pos + local_align_max_dist);
				string chrom_seq = (*reference_sequence)[chrom].substr(pos + 1, right_bound - pos);
				int forward_mismatch = 0;
				int backward_mismatch = 0;
				if(debug)
				{
					cout << "local_align" << endl;
					cout << "chrom seq is " << chrom_seq << endl;
					cout << "read seq is " << pending_seq << endl;
				}
				int forward_score = local_align.align(chrom_seq, pending_seq_tmp, "GT", "AG", true, false, new_splice1, forward_mismatch, pos + 1, premap_len, false, debug);
				int backward_score = local_align.align(chrom_seq, pending_seq_tmp, "CT", "AC", true, false, new_splice2, backward_mismatch, pos + 1, premap_len, false, debug);
				if(debug)
					cout << "forward score is " << forward_score << ", backward score is " << backward_score << endl; 
 				Splice_Info* best_splice = (forward_score >= backward_score ? &new_splice1 : &new_splice2);
				int best_mismatch = (forward_score >= backward_score ? forward_mismatch : backward_mismatch);
				if(best_mismatch <= local_align_max_mis && best_splice->Count_Gap() <= local_align_max_gap)
				{
					best_splice->end_seg_no = anchor_segment;
					if(debug)
					{
						cout << "local_align fixed";
						for(size_t i = 0; i < best_splice->jump_code.size(); i++)
							cout << best_splice->jump_code[i].toString();
						cout << endl;
					}
					if(!check_duplicate_splice(candidate_splice_vec, *best_splice))
						candidate_splice_vec.push_back(*best_splice);
				}
		}*/		
		for(size_t i = 0; i < candidate_splice_vec.size(); i++)
		{
			return_fixed = true;
			if((strand == "+" && (int)anchor_segment == segment_number) || (strand == "-" && anchor_segment == 1))
			{
				(*bwt_vector)[bwt_index].splice_tail.push_back(candidate_splice_vec[i]);
				if(debug)
					cout << "tail fixed for: " << (*bwt_vector)[bwt_index].read_id << endl;
			}
			else if(!fix_tail(bwt_index, anchor_segment + (strand == "+" ? 1 : -1), candidate_splice_vec[i], chrom, candidate_splice_vec[i].start_pos, strand, only_remap))
			{
				(*bwt_vector)[bwt_index].splice_tail.push_back(candidate_splice_vec[i]);
				(*bwt_vector)[bwt_index].incomplete_tail = true;
			}
		}
		candidate_splice_vec.clear();
		if(debug)
			cout << "fix tail end" << endl;
		return return_fixed;
	}

	bool fix_hmer_fw(size_t bwt_index1, size_t bwt_index2, size_t anchor_segment, string chrom, int pos, string strand)
	{
		Splice_Info fixed_splice;
		return fix_hmer_fw(bwt_index1, bwt_index2, anchor_segment, chrom, pos, strand, fixed_splice);
	}

	bool fix_hmer_bw(size_t bwt_index1, size_t bwt_index2, size_t anchor_segment, string chrom, int pos, string strand)
	{
		Splice_Info fixed_splice;
		return fix_hmer_bw(bwt_index1, bwt_index2, anchor_segment, chrom, pos, strand, fixed_splice);
	}

	bool fix_hmer_fw(size_t bwt_index1, size_t bwt_index2, size_t anchor_segment, string chrom, int pos, string strand, Splice_Info& fixed_splice)
	{
		if(debug)
			cout << "fixing hmer fw start" << endl;
		bool return_fixed = false;
		bool remap_fixed = false;
		bool append_fixed = false;
		bool anchor_fixed = false;
		int segment_number = (*bwt_vector)[bwt_index1].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
		Read_Block* r_block = (*bwt_vector)[bwt_index1].pair_no == 1 ? r_block_pe1 : r_block_pe2;
		vector<bool>* segment_no_repeats = (*bwt_vector)[bwt_index1].pair_no == 1 ? pe1_segment_no_repeats : pe2_segment_no_repeats;
		int anchor_segment_length = r_block->get_seg_len(anchor_segment);
		size_t premap_len = (strand == "+" ? r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment - 1) : r_block->get_seg_len(segment_number, (int)anchor_segment + 1));
		int min_next_len;
		if(strand == "+")
			min_next_len = r_block->get_seg_len((int)anchor_segment, (*bwt_vector)[bwt_index2].start_seg_no - 1);
		else
			min_next_len = r_block->get_seg_len((int)anchor_segment, (*bwt_vector)[bwt_index2].start_seg_no + 1);	
		if(pos + min_next_len - max_ins >= (*bwt_vector)[bwt_index2].start)
		{
			if(debug)
				cout << "too short min next len" << pos << min_next_len << (*bwt_vector)[bwt_index2].start << endl;
			return return_fixed;
		}
		vector<Splice_Info> candidate_splice_vec;
		string pending_seq;
		if(strand == "+")
			pending_seq = (*r_block).get_seg_seq(anchor_segment);
		else
			pending_seq = (*r_block).get_revcom_seg_seq(anchor_segment);
		//if(!remap_fixed)
		{
			//cout << "appending tail for fw hmer: " << (*bwt_vector)[0].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << ""<<endl;
			//cout << pending_seq << "\t" << (*reference_sequence)[chrom].substr(pos + 1, segment_length) << endl;
			size_t mismatch_bits = 0;
			string chrom_seq = (*reference_sequence)[chrom].substr(pos + 1, anchor_segment_length);
			append_fixed = score_string(pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);
			if(append_fixed)
			{
				if(debug)
					cout << "fix hmer fw appended" << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				new_splice.start_pos = pos + anchor_segment_length;
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code prefix_match(anchor_segment_length, "M");
					new_splice.jump_code.push_back(prefix_match);
				}
				else
				{
					new_splice.jump_code[new_splice.jump_code.size() - 1].len += anchor_segment_length;
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
		}
		if(juncdb)
		{
			if(debug)
				cout << "remap hmer" << endl;
			vector<Splice_Info> remap_vector;
			remap_fixed = remap_to_juncdb(bwt_index1, bwt_index2, chrom, strand, r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment - 1), r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment) - 1, premap_len, pos + 1, -1, -1, (*bwt_vector)[bwt_index2].start - min_next_len + anchor_segment_length - 1 + max_ins, remap_vector, true, false);	
			for(size_t i = 0; i < remap_vector.size(); i++)
			{
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				new_splice.start_pos = remap_vector[i].end_pos;
				if(new_splice.jump_code.size() == 0)
				{
					new_splice.jump_code.push_back(remap_vector[i].jump_code[0]);
				}
				else
				{
					new_splice.jump_code[new_splice.jump_code.size() - 1].len += remap_vector[i].jump_code[0].len;
				}
				for(size_t j = 1; j < remap_vector[i].jump_code.size(); j++)
				{
					new_splice.jump_code.push_back(remap_vector[i].jump_code[j]);
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}	
			remap_vector.clear();		
		}	
		//if(!remap_fixed && !append_fixed)
		{
			if(debug)
				cout << "fixing tail for fw hmer: " << (*bwt_vector)[bwt_index1].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << ""<<endl;
			pending_seq = pending_seq.substr(0, anchor_segment_length - anchor_length + extend_len);
			string anchor_string;
			if(strand == "+")
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(r_block->get_seg_len(anchor_segment) - anchor_length);
			else
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(0, anchor_length); 
			int anchor_pos;
			bool mapped_fw;
			//int left_bound_forward = max(0, pos + 1 + anchor_segment_length - anchor_length - max_ins);
			int left_bound_forward = max(0, pos - extend_len + 2);
			int left_bound_reverse = left_bound_forward;
			int right_bound;
			if(strand == "+")
				right_bound = (*bwt_vector)[bwt_index2].start  + max_ins - r_block->get_seg_len((int)anchor_segment + 1, (*bwt_vector)[bwt_index1].start_seg_no - 1) - anchor_length;
			else
				right_bound = (*bwt_vector)[bwt_index2].start  + max_ins - r_block->get_seg_len((int)anchor_segment - 1, (*bwt_vector)[bwt_index1].start_seg_no + 1) - anchor_length;
			if((*segment_no_repeats)[anchor_segment - 1] && anchor_search_tail.set(anchor_string, (*reference_sequence)[chrom], strand == "+", strand == "-", max_anchor_hits))
			{
				if(debug)
					cout << "boundary is " << left_bound_forward << "\t" << right_bound << endl;
				while(anchor_search_tail.search_next_simple((*reference_sequence)[chrom], left_bound_forward, left_bound_reverse, right_bound, anchor_pos, mapped_fw))
				{
					int dist = anchor_pos + anchor_length - 1 - pos;
					int gap = dist - anchor_segment_length;		
					size_t prefix_length = 0;
					size_t mismatch_bits = 0;
					if(gap == 0)    
						continue;
					else if(gap > max_del && gap < min_intron_length)
						continue;		
					if(gap < 0)  						 // fix insertion
					{
						if(abs(gap) > max_ins)
							continue;
						bool insertion_fixed = false;
						string chrom_seq = (*reference_sequence)[chrom].substr(pos + 1, anchor_pos - 1 - pos + extend_len);
						if(debug)
							cout << "fixing insertion for " << (*bwt_vector)[bwt_index1].read_id << "\t" << chrom << "\t" <<  (*bwt_vector)[bwt_index1].strand << "\t" <<anchor_pos << "\t" << pending_seq << "\t" << chrom_seq << endl;
						insertion_fixed = genome_scan.Double_anchored_score_ins(pending_seq, chrom_seq, max_ins_mismatch,  prefix_length, mismatch_bits);
						if(insertion_fixed)
						{
							if(debug)
								cout << "insertion fixed"<< endl;
							Splice_Info new_splice(fixed_splice);
							new_splice.start_contig = bwt_index1;
							new_splice.end_contig = bwt_index2;
							new_splice.start_pos = anchor_pos + anchor_length - 1;
							string mapped_read_seq;
							if(prefix_length > 0)
								mapped_read_seq = pending_seq.substr(0, prefix_length);
							if(prefix_length < chrom_seq.length())
								mapped_read_seq += pending_seq.substr(pending_seq.length() - chrom_seq.length() + prefix_length);
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code prefix_match(prefix_length, "M");
								new_splice.jump_code.push_back(prefix_match);
							}
							else
							{
								new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
							}
							Jump_Code skipped_length(pending_seq.length() - chrom_seq.length(), "I");
							Jump_Code suffix_match(chrom_seq.length() -  prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.push_back(skipped_length);
							new_splice.jump_code.push_back(suffix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
						}
					}
					else
					{
						int donor_start = pos + 1;
						int acceptor_start = anchor_pos + extend_len - pending_seq.length() - 2;
						string left_chrom_seq = (*reference_sequence)[chrom].substr(donor_start, pending_seq.length() + 2);
						string right_chrom_seq = (*reference_sequence)[chrom].substr(acceptor_start, pending_seq.length() + 2);
						bool small_deletion = (gap <= max_del);
						if(debug)
							cout << "fixing normal splice " << anchor_pos << "\t"<< pending_seq <<"\t" << left_chrom_seq << "\t" << right_chrom_seq << "\t" << (small_deletion ? "small_deletion" : "not_small_deletion") <<endl;
						if(small_deletion)
							anchor_fixed = genome_scan.Double_anchored_score_least_mis(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, small_deletion);
						else
						{
							string flank_seq;
							anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, single_anchor_noncanonical, flank_seq);
							if(anchor_fixed && (flank_seq == "ATAC" || flank_seq == "GTAT" || flank_seq == "CTGC" || flank_seq == "GCAG") && (gap > max_semicanonical_intron_length))
								anchor_fixed = false;
						}
						if(anchor_fixed)
						{
							Splice_Info new_splice(fixed_splice);
							new_splice.start_contig = bwt_index1;
							new_splice.end_contig = bwt_index2;
							new_splice.start_pos = anchor_pos + anchor_length - 1;
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code prefix_match(prefix_length, "M");
								new_splice.jump_code.push_back(prefix_match);
							}
							else
							{
								new_splice.jump_code[new_splice.jump_code.size() - 1].len += prefix_length;
							}
							Jump_Code skipped_length(gap, "N");
							Jump_Code suffix_match(pending_seq.length() -  prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.push_back(skipped_length);
							new_splice.jump_code.push_back(suffix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
							break;
						}
					}
				}
			}
		}
		for(size_t i = 0; i < candidate_splice_vec.size(); i++)
		{
			if((strand == "+" && (int)anchor_segment + 2 == (*bwt_vector)[bwt_index2].start_seg_no) || (strand == "-" && (int)anchor_segment == (*bwt_vector)[bwt_index2].start_seg_no + 2 ) )
			{
				int double_anchor_segment = anchor_segment + (strand == "+" ? 1 : -1);
				string mid_seq;
				if(strand == "+")
					mid_seq = (*r_block).get_seg_seq(double_anchor_segment);
				else
					mid_seq = (*r_block).get_revcom_seg_seq(double_anchor_segment);
				int double_read_start = r_block->get_seg_len(HEAD_SEG_NO, double_anchor_segment - 1);
				int double_read_end = r_block->get_seg_len(HEAD_SEG_NO, double_anchor_segment) - 1;
				int double_left_bound = candidate_splice_vec[i].start_pos + 1;
				int double_right_bound = (*bwt_vector)[bwt_index2].start - 1;
				int gap = (double_right_bound - double_left_bound) - (double_read_end - double_read_start);
				if(gap < 0)
				{
					if(abs(gap) <= max_ins)
						return_fixed = fix_double_anchor_insertion(bwt_index1, bwt_index2, double_read_start, double_read_end, mid_seq, (*bwt_vector)[bwt_index1].chrom, (*bwt_vector)[bwt_index1].strand, double_left_bound, double_right_bound, candidate_splice_vec[i], false);
				}
				else
					return_fixed = fix_double_anchor(bwt_index1, bwt_index2, double_read_start, double_read_end, mid_seq, (*bwt_vector)[bwt_index1].chrom, (*bwt_vector)[bwt_index1].strand, double_left_bound, double_right_bound, candidate_splice_vec[i], false);
				if(debug && return_fixed )
						cout << "fw hmer fixed between index:" << candidate_splice_vec[i].start_contig << " " << candidate_splice_vec[i].end_contig << endl;	
			}
			else
				fix_hmer_fw(bwt_index1, bwt_index2, anchor_segment + (strand == "+"? 1 : -1), chrom, candidate_splice_vec[i].start_pos, strand, candidate_splice_vec[i]);	
		}
		candidate_splice_vec.clear();
		if(debug)
			cout << "fixing hmer fw end" << endl;
		return return_fixed;
	}

	bool fix_hmer_bw(size_t bwt_index1, size_t bwt_index2, size_t anchor_segment, string chrom, int pos, string strand, Splice_Info& fixed_splice)
	{
		if(debug)
			cout << "fixing hmer bw start " << bwt_index1 << "\t" << bwt_index2 << "\t" << anchor_segment << "\t" << (*bwt_vector)[bwt_index1].start_seg_no << "\t" << (*bwt_vector)[bwt_index1].end_seg_no << "\t" << (*bwt_vector)[bwt_index2].start_seg_no << "\t" << (*bwt_vector)[bwt_index2].end_seg_no <<endl;
		bool return_fixed = false;
		bool remap_fixed = false;
		bool append_fixed = false;
		bool anchor_fixed = false;
		int segment_number = (*bwt_vector)[bwt_index1].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
		Read_Block* r_block = (*bwt_vector)[bwt_index1].pair_no == 1 ? r_block_pe1 : r_block_pe2;
		vector<bool>* segment_no_repeats = (*bwt_vector)[bwt_index1].pair_no == 1 ? pe1_segment_no_repeats : pe2_segment_no_repeats;
		int anchor_segment_length = r_block->get_seg_len(anchor_segment);
		size_t premap_len = (strand == "+" ? r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment - 1) : r_block->get_seg_len(segment_number, (int)anchor_segment + 1));
		int min_previous_len;
		if(strand == "+")
			min_previous_len = r_block->get_seg_len((*bwt_vector)[bwt_index1].end_seg_no + 1, (int)anchor_segment);
		else
			min_previous_len = r_block->get_seg_len((*bwt_vector)[bwt_index1].end_seg_no - 1, (int)anchor_segment);
		if(pos - min_previous_len + max_ins < (*bwt_vector)[bwt_index1].end)  // not enough previous length
			return return_fixed;
		vector<Splice_Info> candidate_splice_vec;
		string pending_seq;
		if(strand == "+")
			pending_seq = (*r_block).get_seg_seq(anchor_segment);
		else
			pending_seq = (*r_block).get_revcom_seg_seq(anchor_segment);
		//if(!remap_fixed)
		{ // append
			//cout << "appending head for bw hmer: " << (*bwt_vector)[0].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << ""<<endl;
			//cout << pending_seq << "\t" << (*reference_sequence)[chrom].substr(pos - segment_length, segment_length) << endl;
			size_t mismatch_bits = 0;
			string chrom_seq = (*reference_sequence)[chrom].substr(pos - pending_seq.length(), pending_seq.length());
			append_fixed = score_string(pending_seq, chrom_seq, max_append_mismatch, mismatch_bits);		
			if(append_fixed)
			{
				if(debug)
					cout << "head appended for fix hmer bw" << endl;
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				new_splice.start_pos = pos - (int)pending_seq.length();
				if(new_splice.jump_code.size() == 0)
				{
					Jump_Code prefix_match(anchor_segment_length, "M");
					new_splice.jump_code.push_back(prefix_match);
				}
				else
				{
					new_splice.jump_code[0].len += anchor_segment_length;
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
		}
		if(juncdb)
		{
			if(debug)
				cout << "remap hmer" << endl;
			vector<Splice_Info> remap_vector;
			remap_fixed = remap_to_juncdb(bwt_index1, bwt_index2, chrom, strand, r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment -1), r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment) - 1, premap_len, -1, pos - 1, (*bwt_vector)[bwt_index1].end + min_previous_len - anchor_segment_length + 1 - max_ins, -1, remap_vector, false, true);
			if(debug && remap_fixed)
				cout << "remap fixed for read " << (*bwt_vector)[bwt_index1].read_id << endl;		
			for(size_t i = 0; i < remap_vector.size(); i++)
			{
				Splice_Info new_splice(fixed_splice);
				new_splice.start_contig = bwt_index1;
				new_splice.end_contig = bwt_index2;
				new_splice.start_pos = remap_vector[i].start_pos;
				if(new_splice.jump_code.size() == 0)
				{
					new_splice.jump_code.insert(new_splice.jump_code.begin(), remap_vector[i].jump_code[remap_vector[i].jump_code.size() - 1]);
				}
				else
				{
					new_splice.jump_code[0].len += remap_vector[i].jump_code[remap_vector[i].jump_code.size() - 1].len;
				}
				for(int j = (int)(remap_vector[i].jump_code.size()) - 2; j >= 0; j--)
				{
					new_splice.jump_code.insert(new_splice.jump_code.begin(), remap_vector[i].jump_code[j]);
				}
				if(!check_duplicate_splice(candidate_splice_vec, new_splice))
					candidate_splice_vec.push_back(new_splice);
			}
			remap_vector.clear();		
		}
		//if(!remap_fixed && !append_fixed)    //fix head
		{
			if(debug)
				cout << "fixing head for bw hmer: " << (*bwt_vector)[0].read_id << "\t" << anchor_segment << "\t" << chrom << "\t" << strand << "\t" << pos << ""<<endl;
			pending_seq = pending_seq.substr(anchor_length - extend_len);
			string anchor_string;
			if(strand == "+")
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(0, anchor_length);
			else
				anchor_string = r_block->get_seg_seq(anchor_segment).substr(r_block->get_seg_len(anchor_segment) - anchor_length);
			int anchor_pos;
			bool mapped_fw;
			int left_bound;
			if(strand == "+")
				left_bound = (*bwt_vector)[bwt_index1].end + r_block->get_seg_len((*bwt_vector)[bwt_index1].end_seg_no + 1, (int)anchor_segment - 1) + anchor_length - max_ins;
			else
				left_bound = (*bwt_vector)[bwt_index1].end + r_block->get_seg_len((*bwt_vector)[bwt_index1].end_seg_no - 1, (int)anchor_segment + 1) + anchor_length - max_ins;
			//int right_bound_forward = min( (int)((*chrom_map)[chrom]) - 1, pos - anchor_segment_length + max_ins + anchor_length - 1);
			int right_bound_forward = min( (int)((*chrom_map)[chrom]) - 1, pos - anchor_length + extend_len - 1);
			int right_bound_reverse = right_bound_forward;
			if((*segment_no_repeats)[anchor_segment - 1] && anchor_search_head.set(anchor_string, (*reference_sequence)[chrom], strand == "+", strand == "-", max_anchor_hits))
			{
				if(debug)
					cout << "searching anchor " << anchor_string << "\t" << left_bound << "\t" << right_bound_forward << "\t" << right_bound_reverse << endl;
				while(anchor_search_head.search_next_simple((*reference_sequence)[chrom], left_bound, right_bound_forward, right_bound_reverse, anchor_pos, mapped_fw))
				{
					size_t prefix_length = 0;
					size_t mismatch_bits = 0;	
					int dist = pos - anchor_pos;
					int gap = dist - anchor_segment_length;	
					if(gap == 0) 							// already did append, skip
						continue;
					else if(gap > max_del && gap < min_intron_length)
						continue;																	
					else if(gap < 0)  	 					// fix insertion
					{
						if(abs(gap) > max_ins)
							continue;
						bool insertion_fixed = false;
						int chorm_seq_start = anchor_pos + anchor_length - extend_len;
						string chrom_seq = (*reference_sequence)[chrom].substr(chorm_seq_start, pos - chorm_seq_start);
						if(debug)
							cout << "fixing insertion for " << (*bwt_vector)[bwt_index1].read_id << "\t" << chrom << "\t" <<  (*bwt_vector)[bwt_index1].strand << "\t" << pending_seq << "\t" << chrom_seq << endl;
						insertion_fixed = genome_scan.Double_anchored_score_ins(pending_seq, chrom_seq, max_ins_mismatch,  prefix_length, mismatch_bits);
						if(insertion_fixed)
						{
							if(debug)
								cout << "insertion fixed" << endl;
							Splice_Info new_splice(fixed_splice);
							new_splice.start_contig = bwt_index1;
							new_splice.end_contig = bwt_index2;
							new_splice.start_pos = anchor_pos;
							string mapped_read_seq;
							if(prefix_length > 0)
								mapped_read_seq = pending_seq.substr(0, prefix_length);
							if(prefix_length < chrom_seq.length())
								mapped_read_seq += pending_seq.substr(pending_seq.length() - chrom_seq.length() + prefix_length);
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code suffix_match(chrom_seq.length() -  prefix_length, "M");
								new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
							}
							else
							{
								new_splice.jump_code[0].len += (chrom_seq.length() -  prefix_length);
							}
							Jump_Code skipped_length(pending_seq.length() - chrom_seq.length(), "I");					
							Jump_Code prefix_match(prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
							new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
						}
					}
					else
					{
						int donor_start = anchor_pos + anchor_length - extend_len;
						int acceptor_start = pos - pending_seq.length() - 2;
						string left_chrom_seq = (*reference_sequence)[chrom].substr(donor_start, pending_seq.length() + 2);
						string right_chrom_seq = (*reference_sequence)[chrom].substr(acceptor_start, pending_seq.length() + 2);
						bool small_deletion = (gap <= max_del);
						if(small_deletion)
							anchor_fixed = genome_scan.Double_anchored_score_least_mis(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits,  small_deletion);
						else
						{
							string flank_seq;
							anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, single_anchor_noncanonical, flank_seq);
							if(anchor_fixed && (flank_seq == "ATAC" || flank_seq == "GTAT" || flank_seq == "CTGC" || flank_seq == "GCAG") && (gap > max_semicanonical_intron_length))
								anchor_fixed = false;
						}
						if(anchor_fixed)
						{
							Splice_Info new_splice(fixed_splice);
							new_splice.start_contig = bwt_index1;
							new_splice.end_contig = bwt_index2;
							new_splice.start_pos = anchor_pos;
							if(new_splice.jump_code.size() == 0)
							{
								Jump_Code suffix_match(pending_seq.length() -  prefix_length, "M");
								new_splice.jump_code.insert(new_splice.jump_code.begin(), suffix_match);	
							}
							else
							{
								new_splice.jump_code[0].len += (pending_seq.length() -  prefix_length);
							}
							Jump_Code skipped_length(gap, "N");					
							Jump_Code prefix_match(prefix_length + anchor_length - extend_len, "M");
							new_splice.jump_code.insert(new_splice.jump_code.begin(), skipped_length);
							new_splice.jump_code.insert(new_splice.jump_code.begin(), prefix_match);
							if(!check_duplicate_splice(candidate_splice_vec, new_splice))
								candidate_splice_vec.push_back(new_splice);
							break;
						}
					}
				}
			}
		}
		for(size_t i = 0; i < candidate_splice_vec.size(); i++)
		{
			if((strand == "+" && (int)anchor_segment == (*bwt_vector)[bwt_index1].end_seg_no + 2) || (strand == "-" && (int)anchor_segment + 2 == (*bwt_vector)[bwt_index1].end_seg_no) )
			{
				string mid_seq;
				int double_anchor_segment = anchor_segment + (strand == "+" ? -1 : 1);
				if(strand == "+")
					mid_seq = (*r_block).get_seg_seq(double_anchor_segment);
				else
					mid_seq = (*r_block).get_revcom_seg_seq(double_anchor_segment);
				int double_read_start = r_block->get_seg_len(HEAD_SEG_NO, (int)double_anchor_segment - 1);
				int double_read_end = r_block->get_seg_len(HEAD_SEG_NO, (int)double_anchor_segment) - 1;
				int double_left_bound = (*bwt_vector)[bwt_index1].end + 1;
				int double_right_bound = candidate_splice_vec[i].start_pos - 1;
				int gap = (double_right_bound - double_left_bound) - (double_read_end - double_read_start);
				if(gap < 0)
				{
					if(abs(gap) <= max_ins)
						return_fixed = fix_double_anchor_insertion(bwt_index1, bwt_index2, double_read_start, double_read_end, mid_seq, (*bwt_vector)[bwt_index1].chrom, (*bwt_vector)[bwt_index1].strand, double_left_bound, double_right_bound, candidate_splice_vec[i], true);
				}
				else
					return_fixed = fix_double_anchor(bwt_index1, bwt_index2, double_read_start, double_read_end, mid_seq, (*bwt_vector)[bwt_index1].chrom, (*bwt_vector)[bwt_index1].strand, double_left_bound, double_right_bound, candidate_splice_vec[i], true);
				if(debug && return_fixed)
						cout << "bw hmer fixed between index:" << candidate_splice_vec[i].start_contig << " " << candidate_splice_vec[i].end_contig << endl;			
			}
			else
				fix_hmer_bw(bwt_index1, bwt_index2, anchor_segment + (strand == "+"? -1 : 1), chrom, candidate_splice_vec[i].start_pos, strand, candidate_splice_vec[i]);
		}
		candidate_splice_vec.clear();
		if(debug)
			cout << "fixing hmer bw end" << endl;
		return return_fixed;
	}

	void search_double_anchor()
	{
		for (size_t i = 0; i < (*bwt_vector).size(); i++)
		{
			if((*bwt_vector)[i].splice_tail.size() != 0)
				continue;
			for (size_t j = i + 1; j < (*bwt_vector).size(); j++)
			{
				if((*bwt_vector)[j].splice_head.size() != 0)
					continue;
				Read_Block* r_block = ((*bwt_vector)[i].pair_no == 1) ? r_block_pe1 : r_block_pe2;
				int relation = check_relation((*bwt_vector)[i], (*bwt_vector)[j]);
				if(relation == FIX_NO_RELATIONSHIP || relation == FIX_TOO_FAR)            // not related or too far
					break;
				else if(relation == FIX_TOO_CLOSE)            // continue search
					continue;
				else                        // in range or small insertion
				{ 
					bool double_anchor_fixed = false;
					if ((*bwt_vector)[i].check_tail((*bwt_vector)[j].start_seg_no)&&(*bwt_vector)[i].strand =="+")   //strand + 
					{
						int seg_gap=(*bwt_vector)[j].start_seg_no - (*bwt_vector)[i].end_seg_no;
						if(relation == FIX_INSERTION)
						{
							if(seg_gap == 1)
							{
								int trunc_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - trunc_len;
								string mid_seq = r_block->get_seg_seq((*bwt_vector)[i].end_seg_no).substr(trunc_left_len) + r_block->get_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,trunc_len);
								double_anchor_fixed = fix_double_anchor_insertion(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no - 1) + trunc_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - trunc_len + 1, (*bwt_vector)[j].start + trunc_len - 1, false);	
							}
							else if(seg_gap == 2)
							{
								int extend_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - extend_len;
								string mid_seq = r_block->get_seg_seq((*bwt_vector)[i].end_seg_no).substr(extend_left_len) + r_block->get_seg_seq((*bwt_vector)[i].end_seg_no+1);
								mid_seq.append(r_block->get_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,extend_len));
								double_anchor_fixed = fix_double_anchor_insertion(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no - 1) + extend_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - extend_len + 1, (*bwt_vector)[j].start + extend_len - 1, false);
							}
						}
						else if(relation == FIX_DOUBLE_ANCHOR)
						{
							if(seg_gap==1)	  //continuous hole
							{
								int trunc_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - trunc_len;
								string mid_seq = r_block->get_seg_seq((*bwt_vector)[i].end_seg_no).substr(trunc_left_len) + r_block->get_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,trunc_len);
								double_anchor_fixed = fix_double_anchor(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no - 1) + trunc_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - trunc_len + 1, (*bwt_vector)[j].start+trunc_len - 1, false);	
							}
							else if(seg_gap== 2)  //gap hole
							{
								int extend_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - extend_len;
								string mid_seq = r_block->get_seg_seq((*bwt_vector)[i].end_seg_no).substr(extend_left_len) + r_block->get_seg_seq((*bwt_vector)[i].end_seg_no+1);
								mid_seq.append(r_block->get_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,extend_len));
								double_anchor_fixed = fix_double_anchor(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no - 1) + extend_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - extend_len + 1, (*bwt_vector)[j].start+extend_len - 1, false);
							}
							else if(seg_gap==0 && (*bwt_vector)[i].start_seg_no!=(*bwt_vector)[i].end_seg_no&&(*bwt_vector)[j].start_seg_no!=(*bwt_vector)[j].end_seg_no)
							{/////// special case
								int extend_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no - 1) - extend_len;
								string mid_seq = r_block->get_seg_seq((*bwt_vector)[i].end_seg_no-1).substr(extend_left_len) + r_block->get_seg_seq((*bwt_vector)[i].end_seg_no);
								mid_seq.append(r_block->get_seg_seq((*bwt_vector)[j].start_seg_no+1).substr(0,extend_len));
								double_anchor_fixed = fix_double_anchor(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 2) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no) + extend_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - extend_len + 1, (*bwt_vector)[j].start + r_block->get_seg_len((*bwt_vector)[j].start_seg_no) + extend_len - 1, false);
							}
							else if(seg_gap > 2 && seg_gap<=MAX_JUMP) //hmer
							{
								bool fw_fixed = false, bw_fixed = false;
								fw_fixed = fix_hmer_fw(i, j, (*bwt_vector)[i].end_seg_no + 1, (*bwt_vector)[i].chrom, (*bwt_vector)[i].end, (*bwt_vector)[i].strand);
								if(!fw_fixed)
									bw_fixed = fix_hmer_bw(i, j, (*bwt_vector)[j].start_seg_no - 1, (*bwt_vector)[j].chrom, (*bwt_vector)[j].start, (*bwt_vector)[j].strand);
								double_anchor_fixed = fw_fixed || bw_fixed;
							}
							else
								continue;
						}
						if(double_anchor_fixed)
						{
							if((*bwt_vector)[i].end_seg_no!=(*bwt_vector)[j].start_seg_no)  //not special case
							{
								(*bwt_vector)[i].tail=(*bwt_vector)[j].start_seg_no;
								(*bwt_vector)[j].head=(*bwt_vector)[i].end_seg_no;
							}
							else     //special case
							{
								(*bwt_vector)[i].tail=(*bwt_vector)[j].start_seg_no + 1;
								(*bwt_vector)[j].head=(*bwt_vector)[i].end_seg_no - 1;
							}
						}
					}
					else if((*bwt_vector)[i].check_tail((*bwt_vector)[j].start_seg_no)&&(*bwt_vector)[i].strand =="-")  //strand - hole
					{
						int seg_gap=(*bwt_vector)[i].end_seg_no - (*bwt_vector)[j].start_seg_no;
						if(relation == FIX_INSERTION)
						{
							if(seg_gap == 1)
							{
								int trunc_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - trunc_len;
								string mid_seq = r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no).substr(trunc_left_len) + r_block->get_revcom_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,trunc_len);
								double_anchor_fixed = fix_double_anchor_insertion(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no) - trunc_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + trunc_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - trunc_len + 1, (*bwt_vector)[j].start + trunc_len - 1, false);
							}
							else if(seg_gap == 2)
							{	
								int extend_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - extend_len;
								string mid_seq = r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no).substr(extend_left_len) + r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no-1);
								mid_seq.append(r_block->get_revcom_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,extend_len));
								double_anchor_fixed = fix_double_anchor_insertion(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no) - extend_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + extend_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - extend_len + 1, (*bwt_vector)[j].start+extend_len-1, false);
							}
						}
						else if(relation == FIX_DOUBLE_ANCHOR)
						{
							if(seg_gap== 1)  //continuous hole
							{
								int trunc_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - trunc_len;
								string mid_seq = r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no).substr(trunc_left_len) + r_block->get_revcom_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,trunc_len);
								double_anchor_fixed = fix_double_anchor(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no) - trunc_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + trunc_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - trunc_len + 1, (*bwt_vector)[j].start+trunc_len - 1, false);
							}
							else if(seg_gap== 2)  //gap hole
							{
								int extend_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - extend_len;
								string mid_seq = r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no).substr(extend_left_len) + r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no-1);
								mid_seq.append(r_block->get_revcom_seg_seq((*bwt_vector)[j].start_seg_no).substr(0,extend_len));
								double_anchor_fixed = fix_double_anchor(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no) - extend_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no - 1) + extend_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - extend_len + 1, (*bwt_vector)[j].start+extend_len - 1, false);
							}
							else if(seg_gap==0&&(*bwt_vector)[i].start_seg_no!=(*bwt_vector)[i].end_seg_no&&(*bwt_vector)[j].start_seg_no!=(*bwt_vector)[j].end_seg_no)
							{/////// special case
								int extend_left_len = r_block->get_seg_len((*bwt_vector)[i].end_seg_no + 1) - extend_len;
								string mid_seq = r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no+1).substr(extend_left_len) + r_block->get_revcom_seg_seq((*bwt_vector)[i].end_seg_no);
								mid_seq.append(r_block->get_revcom_seg_seq((*bwt_vector)[j].start_seg_no-1).substr(0,extend_len));
								double_anchor_fixed = fix_double_anchor(i, j, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[j].start_seg_no - 1) - extend_len, r_block->get_seg_len(HEAD_SEG_NO, (*bwt_vector)[i].end_seg_no) + extend_len - 1, mid_seq, (*bwt_vector)[i].chrom, (*bwt_vector)[i].strand, (*bwt_vector)[i].end - r_block->get_seg_len((*bwt_vector)[i].end_seg_no) - extend_len + 1, (*bwt_vector)[j].start + r_block->get_seg_len((*bwt_vector)[j].start_seg_no) + extend_len - 1, false);
							}
							else if(seg_gap>2&&seg_gap<=MAX_JUMP) //hmer
							{
								bool fw_fixed, bw_fixed;
								fw_fixed = fix_hmer_fw(i, j, (*bwt_vector)[i].end_seg_no - 1, (*bwt_vector)[i].chrom, (*bwt_vector)[i].end, (*bwt_vector)[i].strand);
								if(!fw_fixed)
									bw_fixed = fix_hmer_bw(i, j, (*bwt_vector)[j].start_seg_no + 1, (*bwt_vector)[j].chrom, (*bwt_vector)[j].start, (*bwt_vector)[j].strand);
								double_anchor_fixed = fw_fixed || bw_fixed;
							}
							else 
								continue;
						}
						if(double_anchor_fixed)
						{
							if((*bwt_vector)[i].end_seg_no!=(*bwt_vector)[j].start_seg_no)   //not special case
							{
								(*bwt_vector)[i].tail = (*bwt_vector)[j].start_seg_no;
								(*bwt_vector)[j].head = (*bwt_vector)[i].end_seg_no;
							}
							else     //special case
							{
								(*bwt_vector)[i].tail = (*bwt_vector)[j].start_seg_no - 1;
								(*bwt_vector)[j].head = (*bwt_vector)[i].end_seg_no + 1;
							}
						}
					}	
				}
			}
		}	
	}

	void search_single_anchor()
	{
		for (unsigned int i = 0; i < (*bwt_vector).size(); i++)   //single anchored
		{
			int segment_number = (*bwt_vector)[i].pair_no == 1 ? pe1_segment_number : pe2_segment_number;
			if((*bwt_vector)[i].strand=="+")
			{
				if ((*bwt_vector)[i].head == -1 && (*bwt_vector)[i].start_seg_no != 1 && ((*bwt_vector)[i].start_seg_no - 1 <= MAX_JUMP))
				{// fix head	
					bool only_remap = ((*bwt_vector)[i].splice_head.size() != 0);
					Splice_Info fixed_splice;
					fix_head(i, (*bwt_vector)[i].start_seg_no - 1, fixed_splice, (*bwt_vector)[i].chrom, (*bwt_vector)[i].start, (*bwt_vector)[i].strand, only_remap);
				}
				if ((*bwt_vector)[i].tail== -1 && (*bwt_vector)[i].end_seg_no!= segment_number&& (segment_number - (*bwt_vector)[i].end_seg_no <= MAX_JUMP))
				{// fix tail
					bool only_remap = ((*bwt_vector)[i].splice_tail.size() != 0);
					Splice_Info fixed_splice;
					fix_tail(i, (*bwt_vector)[i].end_seg_no + 1, fixed_splice, (*bwt_vector)[i].chrom, (*bwt_vector)[i].end, (*bwt_vector)[i].strand, only_remap);
				}
			}
			else if((*bwt_vector)[i].strand=="-")
			{
				if ((*bwt_vector)[i].head== -1 && (*bwt_vector)[i].start_seg_no!=segment_number&& (segment_number - (*bwt_vector)[i].start_seg_no <= MAX_JUMP))
				{// fix head
					bool only_remap = ((*bwt_vector)[i].splice_head.size() != 0);
					Splice_Info fixed_splice;
					fix_head(i, (*bwt_vector)[i].start_seg_no + 1, fixed_splice, (*bwt_vector)[i].chrom, (*bwt_vector)[i].start, (*bwt_vector)[i].strand, only_remap);
				}
				if ((*bwt_vector)[i].tail== -1 && (*bwt_vector)[i].end_seg_no!=1 && ((*bwt_vector)[i].end_seg_no - 1 <= MAX_JUMP)) 
				{// fix tail
					bool only_remap = ((*bwt_vector)[i].splice_tail.size() != 0);
					Splice_Info fixed_splice;
					fix_tail(i, (*bwt_vector)[i].end_seg_no - 1, fixed_splice, (*bwt_vector)[i].chrom, (*bwt_vector)[i].end, (*bwt_vector)[i].strand, only_remap);
				}
			}
		}
	}
	
	bool check_valid_overlapping(vector<Sam_Info*>::iterator it1, vector<Sam_Info*>::iterator it2)
	{
		//cout << "checking valid overlapping" << endl;
		Sam_Info *sam1, *sam2;
		if((*it1)->start_pos < (*it2)->start_pos)
		{
			sam1 = (*it1);
			sam2 = (*it2);
		}
		else
		{
			sam1 = (*it2);
			sam2 = (*it1);
		}
		vector<Jump_Code> sam1_jumpcode;
		int sam1_cur_pos = sam1->start_pos - 1;		
		for(size_t i = 0; i < sam1->jump_code.size(); i++)
		{
			if(sam1->jump_code[i].type != "I")
				sam1_cur_pos += sam1->jump_code[i].len;
			if(sam1_cur_pos >= sam2->start_pos)
			{
				if(sam1->jump_code[i].type == "M")
				{
					for(size_t j = i; j < sam1->jump_code.size(); j++)
						sam1_jumpcode.push_back(sam1->jump_code[j]);	
					sam1_jumpcode[0].len = (sam1_cur_pos - sam2->start_pos + 1);
				}
				else if(sam1->jump_code[i].type == "N")
				{
					return false;
				}
			}
		}
		for(size_t i = 0; i < sam1_jumpcode.size() && i < sam2->jump_code.size(); i++)
		{
			if(sam1_jumpcode[i].type != sam2->jump_code[i].type)
				return false;
			if((i != sam1_jumpcode.size() - 1) && (i != sam2->jump_code.size() - 1))
			{
				if(sam1_jumpcode[i].len != sam2->jump_code[i].len)
					return false;
			}
			else if((i == sam1_jumpcode.size() - 1) && (i != sam2->jump_code.size() - 1))
			{
				if(sam1_jumpcode[i].len > sam2->jump_code[i].len)
					return false;	
			}
			else if((i != sam1_jumpcode.size() - 1) && (i == sam2->jump_code.size() - 1))
			{
				if(sam1_jumpcode[i].len < sam2->jump_code[i].len)
					return false;
			}	
		}
		return true;
	}
	
	bool check_valid_pairing(vector<Sam_Info*>::iterator it1, vector<Sam_Info*>::iterator it2)
	{
			if(((*it1)->strand == "+" && (*it2)->start_pos < (*it1)->start_pos) || ((*it1)->strand == "-" && (*it1)->start_pos < (*it2)->start_pos))
				return false;
			int mate_dist = ((*it1)->strand == "+") ? ((*it2)->start_pos - (*it1)->end_pos) : ((*it1)->start_pos - (*it2)->end_pos);
			if((*it1)->chrom == (*it2)->chrom && (*it1)->strand != (*it2)->strand && mate_dist < max_pair_dist)
			{
				if(mate_dist > 0 || check_valid_overlapping(it1, it2))
					return true;	
			}
			return false;
	}
	
	void pair_alignment()
	{
		for(size_t i =0; i < sam_vec.size(); i++)
		{
			if(sam_vec[i].pair_no == 1)
				pe1_sam.push_back(&sam_vec[i]);
			else
				pe2_sam.push_back(&sam_vec[i]);	
		}
		if(debug)
			cout << "end 1 alignment size " << pe1_sam.size() << " end 2 alignment size " << pe2_sam.size() << endl;
		//if(pe1_sam.size() < 1 || pe2_sam.size() < 1)
		//	return;
		for(vector<Sam_Info*>::iterator it1 = pe1_sam.begin(); it1 != pe1_sam.end(); it1 ++)
		{
			for(vector<Sam_Info*>::iterator it2 = pe2_sam.begin(); it2 != pe2_sam.end(); it2 ++)
			{
				if(check_valid_pairing(it1, it2))
				{
					Paired_Sam_Info new_pair_info(*it1, *it2, adaptive_mate_pair_dist);
					paired_sam_vec.push_back(new_pair_info);
					read_paired = true;
				}
			}	
		}
		if(debug)
			cout << "read pairing complete" <<endl;
	}

	bool check_pairing()
	{
		vector<Sam_Info*> pe1_sam;
		vector<Sam_Info*> pe2_sam;
		for(size_t i =0; i < sam_vec.size(); i++)
		{
			if(sam_vec[i].pair_no == 1)
				pe1_sam.push_back(&sam_vec[i]);
			else
				pe2_sam.push_back(&sam_vec[i]);	
		}
		for(vector<Sam_Info*>::iterator it1 = pe1_sam.begin(); it1 != pe1_sam.end(); it1 ++)
		{
			for(vector<Sam_Info*>::iterator it2 = pe2_sam.begin(); it2 != pe2_sam.end(); it2 ++)
			{
				if(check_valid_pairing(it1, it2))
				{
					return true;
				}
			}	
		}
		return false;
	}

	void train_adaptive_mate_dist(int num_output)
	{
		//cout << r_block_pe1-> read_id << " " << read_paired << " " << mate_distance_training_count << " " << num_output << " " << endl;
		if(mate_distance_training_count >= mate_distance_training_set)
			return;
		if(read_paired && mate_distance_training_count < mate_distance_training_set && num_output == 1 && paired_sam_vec[0].mate_dist > 0)
		{
			//cout << mate_distance_training_count << " th read " << r_block_pe1-> read_id << " trained" << ", mate distance is " << paired_sam_vec[0].mate_dist << endl;
			mate_distance_training_count ++;
			mate_distance_training_array.push_back(paired_sam_vec[0].mate_dist);
		}
		if(mate_distance_training_count == mate_distance_training_set)
		{
			sort(mate_distance_training_array.begin(), mate_distance_training_array.end());
			adaptive_mate_pair_dist = mate_distance_training_array[mate_distance_training_set / 2];
			mate_distance_training_array.clear();
			//if(debug)
				//cout << "adaptive mate pair distance is " << adaptive_mate_pair_dist << endl;
		}	
	}

	void mapsplice_search()
	{
		search_double_anchor();
		search_single_anchor();
	}

	void mapsplice_report()
	{
		output_unspliced_read();
		if(debug)
			cout << "preparing unspliced read done" << endl;
		for(size_t i = 0; i < (*bwt_vector).size(); i++)
		{
			if((*bwt_vector)[i].head == -1 && ((*bwt_vector)[i].splice_internal.size() != 0 || (*bwt_vector)[i].splice_head.size() != 0 || (*bwt_vector)[i].splice_tail.size() != 0))
			{
				if(debug)
					cout << "generating sam for " << i << " of " << (*bwt_vector).size() << endl;
				size_t num_current_alignments = ((*bwt_vector)[i].pair_no == 1 ? num_pe1_alignment : num_pe2_alignment);
				if(num_current_alignments < maximum_alignments)
					generate_sam_head(i);
			}
		}
	}

	int output_sam_single(vector<Sam_Info*>& single_sam_vec)
	{
		int num_output = 0;
		vector<Sam_Info*>::iterator it, it_end;
		sort_single_alignment(single_sam_vec.begin(), single_sam_vec.end(), it_end, pairing_intron_diff);
		for(it = single_sam_vec.begin(); it != it_end; it++)
		{
			if(it - single_sam_vec.begin() < maximum_alignments_output)
			{
				print_sam_info_to_file((*(*it)));
				num_output ++;
			}
		}
		return num_output;
	}

	int output_sam_paired()
	{
		int num_output = 0;
		if(debug)
			cout << "sorting paired alignments, size is " << paired_sam_vec.size() << endl;
		vector<Paired_Sam_Info>::iterator it, it_end;
		sort_pair_alignment(paired_sam_vec.begin(), paired_sam_vec.end(), it_end, pairing_mate_dist_diff, pairing_intron_diff);	
		map<Sam_Info*, int> printed_sam;
		for(it = paired_sam_vec.begin(); it != it_end; it++)
		{
				if((printed_sam.find(it->paired_sam.first) == printed_sam.end()) && (printed_sam.find(it->paired_sam.second) == printed_sam.end()) && (it - paired_sam_vec.begin() < maximum_alignments_output))
				{
					print_sam_info_to_file(*(it->paired_sam.first));
					print_sam_info_to_file(*(it->paired_sam.second));
					printed_sam.insert(make_pair(it->paired_sam.first, 1));
					printed_sam.insert(make_pair(it->paired_sam.second, 1));
					num_output ++;
				}
		}
		printed_sam.clear();
		return num_output;
	}

	int output_sam()
	{
		int num_output = 0;
		if(read_paired)
			num_output = output_sam_paired();
		else
		{
			if(pe1_sam.size() > 0)
				num_output += output_sam_single(pe1_sam);
			if(pe2_sam.size() > 0)
				num_output += output_sam_single(pe2_sam);	
		}
		return num_output;
	}

	void output_unmapped_read()
	{
		if(!pe1_aligned || output_unmapped_pe)
			print_unmapped_read_to_file(unmapped_fs1, r_block_pe1);
		if(pair_end && (!pe2_aligned || output_unmapped_pe))
			print_unmapped_read_to_file(unmapped_fs2, r_block_pe2);
	}

	void optimize_repeats()
	{
		bool temp_double_anchor_noncanonical = double_anchor_noncanonical;
		bool temp_single_anchor_noncanonical = single_anchor_noncanonical;
		bool temp_fusion = fusion;
		bool temp_try_hard_ins = try_hard_ins;
		bool temp_do_local_align = do_local_align;
		bool temp_juncdb = juncdb;
		double_anchor_noncanonical = false;
		single_anchor_noncanonical = false;
		try_hard_ins = false;
		fusion = false;
		do_local_align = false;
		juncdb = false;
		if(!pair_end && !pe1_aligned && !(*pe1_no_repeats))
		{
			bwt_vector->clear();
			remap_repeats_full_read(1);
			output_unspliced_read();
		}
		if(pair_end)
		{
 			if(debug)
				cout << "full read repeat fix" << endl;
 			bwt_vector->clear();
 			if(!pe1_aligned && !pe2_aligned)
 			{
 				if(!(*pe1_no_repeats))
					remap_repeats_full_read(1);
				if(!(*pe2_no_repeats))
					remap_repeats_full_read(2);	
 			}
			else if(pe2_aligned && !pe1_aligned && !(*pe1_no_repeats))
				remap_repeats_full_read(1);
			else if(pe1_aligned && !pe2_aligned && !(*pe2_no_repeats))
				remap_repeats_full_read(2);	
			else if(pe1_aligned && pe2_aligned && !check_pairing())
			{
				if(!(*pe1_no_repeats))
					remap_repeats_full_read(1);
				if(!(*pe2_no_repeats))
					remap_repeats_full_read(2);
			}
			output_unspliced_read();	
			if( ((pe1_aligned && !pe2_aligned && !(*pe2_no_repeats)) || (pe2_aligned && !pe1_aligned && !(*pe1_no_repeats)) || (pe1_aligned && pe2_aligned && !check_pairing()) ) )
			{
				if(debug)
					cout << "segment repeat fix" << endl;
				size_t repeat_index_bound = sam_vec.size();
				for(size_t i = 0; i < repeat_index_bound; i++)
				{
					if((sam_vec[i].pair_no == 1 && (*pe2_no_repeats)) || (sam_vec[i].pair_no == 2 && (*pe1_no_repeats)))
						continue;
					if(debug)
						cout << "fixing repeat for alignment " << sam_vec[i].chrom << "\t" << sam_vec[i].start_pos << endl;
					bwt_vector->clear();
					remap_repeats_segment(i);
					if(debug)
					{
						for(size_t j = 0; j < bwt_vector->size(); j++)
							cout << (*bwt_vector)[j].start_seg_no << "\t" << (*bwt_vector)[j].chrom << "\t" << (*bwt_vector)[j].strand << "\t" << (*bwt_vector)[j].start << "\t" << (*bwt_vector)[j].end << endl;
					}
					if(bwt_vector->size() > 0)
					{
						sort((*bwt_vector).begin(), (*bwt_vector).end(), bwtmap_compare);
						combine_contig();
						mapsplice_search();
						mapsplice_report();
					}
				}
			}
		}
		double_anchor_noncanonical = temp_double_anchor_noncanonical;
		single_anchor_noncanonical = temp_single_anchor_noncanonical;
		try_hard_ins = temp_try_hard_ins;
		do_local_align = temp_do_local_align;
		fusion = temp_fusion;
		juncdb = temp_juncdb;
		bwt_vector->clear();
	}

	void full_read_local_align()
	{
		if(debug)
			cout << "doing full read local alignment" << endl;
		if(sam_vec.size() != 1)
			return;
		int left_bound, right_bound;
		string r_strand;
		if(sam_vec[0].strand == "+")
		{
			left_bound = sam_vec[0].start_pos + 1; ;
			right_bound = min((int)((*chrom_map)[sam_vec[0].chrom]) - 1, sam_vec[0].end_pos + local_align_max_dist);
			r_strand = "-";
		}
		else
		{
			left_bound = max(0, sam_vec[0].start_pos - local_align_max_dist);
			right_bound = sam_vec[0].start_pos - 1;
			r_strand = "+";
		}
		Read_Block* r_block = sam_vec[0].pair_no == 1 ? r_block_pe2 : r_block_pe1;
		int start_seg_no = HEAD_SEG_NO;
		int end_seg_no = (sam_vec[0].pair_no == 1 ? pe2_segment_number : pe1_segment_number);
		string pending_seq;
		if(r_strand == "+")
		{
			for(int i = start_seg_no; i <= end_seg_no; i++)
			{
				pending_seq += r_block->get_seg_seq(i);
			}
		}
		else                   //out put - sequence
		{
			for(int i = end_seg_no; i >= start_seg_no; i--)
			{
				pending_seq += r_block->get_revcom_seg_seq(i);
			}
		}
		string chrom_seq = (*reference_sequence)[sam_vec[0].chrom].substr(left_bound, right_bound - left_bound + 1);
		Sam_Info forward_strand_sam_info, reverse_strand_sam_info;
		int forward_mismatch = 0, backward_mismatch = 0;
		int forward_score = local_align.align(chrom_seq, pending_seq, "GT", "AG", false, false, forward_strand_sam_info, forward_mismatch, left_bound, debug);
		int backward_score = local_align.align(chrom_seq, pending_seq, "CT", "AC", false, false, reverse_strand_sam_info, backward_mismatch, left_bound, debug);	
		Sam_Info* best_sam = (forward_score >= backward_score ? &forward_strand_sam_info : &reverse_strand_sam_info);
 		int best_mismatch = (forward_score >= backward_score ? forward_mismatch : backward_mismatch);
		if(best_mismatch <= local_align_max_mis && best_sam->count_gap() <= local_align_max_gap)
		{
			best_sam->chrom = sam_vec[0].chrom;
			best_sam->strand = (sam_vec[0].strand == "+" ? "-" : "+");
			best_sam->pair_no = (sam_vec[0].pair_no == 1 ? 2 : 1);
			best_sam->start_seg_no = (best_sam->strand == "+" ? start_seg_no : end_seg_no);
			best_sam->end_seg_no = (best_sam->strand == "+" ? end_seg_no : start_seg_no);
			print_saminfo((*best_sam));
		}
	}
	
	bool fix_double_anchor_fusion(size_t bwt_index1, size_t bwt_index2, int read_start, int read_end, string pending_seq, int left_bound, int right_bound, int buffer_len)
	{
		if(debug)
			cout << "fixing fusion double anchor" <<endl;
		bool return_fixed = false;
		int premap_len = 0; //(strand == "+" ? read_start : total_read_len - read_end);
		{
			size_t prefix_length = 0;
			size_t mismatch_bits = 0;
			string left_chrom_seq;
			string right_chrom_seq;
			string flank_seq;
			int left_anchor_len = 0;
			int right_anchor_len = 0;
			if(fusion_candidate[bwt_index1].strand == "+")
			{
				if(left_bound + (int)pending_seq.length() + 2 > (int)((*chrom_map)[fusion_candidate[bwt_index1].chrom]))
					return return_fixed;
				if(debug)
					cout << "+ left chrom seq start from " << fusion_candidate[bwt_index1].chrom << " " << left_bound << ", length is " << pending_seq.length() + 2 << endl;
				left_chrom_seq = (*reference_sequence)[fusion_candidate[bwt_index1].chrom].substr(left_bound, pending_seq.length() + 2);
				left_anchor_len = fusion_candidate[bwt_index1].jump_code[fusion_candidate[bwt_index1].jump_code.size() - 1].len;
			}
			else
			{
				if(left_bound < (int)pending_seq.length() + 1)
					return return_fixed;
				if(debug)
					cout << "- left chrom seq start from " << fusion_candidate[bwt_index1].chrom << " " << left_bound - pending_seq.length() - 1 << ", length is " << pending_seq.length() + 2 << endl;
				left_chrom_seq = revcomp((*reference_sequence)[fusion_candidate[bwt_index1].chrom].substr(left_bound - pending_seq.length() - 1, pending_seq.length() + 2));
				left_anchor_len = fusion_candidate[bwt_index1].jump_code[0].len;
			}
			if(fusion_candidate[bwt_index2].strand == "+")
			{
				if(right_bound < (int)pending_seq.length() + 1)
					return return_fixed;
				if(debug)
					cout << "+ right chrom seq start from " << fusion_candidate[bwt_index2].chrom << " " << right_bound - pending_seq.length() - 1 << ", length is " << pending_seq.length() + 2 << endl;
				right_chrom_seq = (*reference_sequence)[fusion_candidate[bwt_index2].chrom].substr(right_bound - pending_seq.length() - 1, pending_seq.length() + 2);	
				right_anchor_len = fusion_candidate[bwt_index2].jump_code[0].len;
			}
			else
			{
				if(right_bound + (int)pending_seq.length() + 2 > (int)((*chrom_map)[fusion_candidate[bwt_index2].chrom]))
					return return_fixed;
				if(debug)
					cout << "- right chrom seq start from" << fusion_candidate[bwt_index2].chrom << " " << right_bound << ", length is " << pending_seq.length() + 2 << endl;
				right_chrom_seq = revcomp((*reference_sequence)[fusion_candidate[bwt_index2].chrom].substr(right_bound, pending_seq.length() + 2));
				right_anchor_len = fusion_candidate[bwt_index2].jump_code[fusion_candidate[bwt_index2].jump_code.size() - 1].len;
			}
			if(left_anchor_len <= buffer_len || right_anchor_len <= buffer_len)
				return return_fixed;
			if(debug)
				cout << left_chrom_seq << "\t" << right_chrom_seq << "\t" << pending_seq << endl;
			return_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_double_splice_mismatch, mismatch_bits, fusion_double_anchor_noncanonical, flank_seq);
			size_t suffix_length = pending_seq.length() -  prefix_length;
			/*bool no_cross_junc = true;
			if((int)prefix_length <= buffer_len - left_anchor_len)
				no_cross_junc = false;
			if((int)suffix_length <= buffer_len - right_anchor_len)
				no_cross_junc = false;*/
			if(return_fixed)
			{
				if(debug)
					cout << "double anchor fusion fixed, prefix length is " << prefix_length << endl;		
				Fusion_Splice new_fusion_splice;
				new_fusion_splice.flank_seq = flank_seq;
				new_fusion_splice.first_splice.start_contig = bwt_index1;
				new_fusion_splice.first_splice.end_contig = bwt_index2;
				new_fusion_splice.second_splice.start_contig = bwt_index1;
				new_fusion_splice.second_splice.end_contig = bwt_index2;
				new_fusion_splice.first_splice.chrom = fusion_candidate[bwt_index1].chrom;
				new_fusion_splice.second_splice.chrom = fusion_candidate[bwt_index2].chrom;
				new_fusion_splice.first_splice.strand = fusion_candidate[bwt_index1].strand;
				new_fusion_splice.second_splice.strand = fusion_candidate[bwt_index2].strand;
				new_fusion_splice.first_splice.buffer_len = buffer_len;
				new_fusion_splice.second_splice.buffer_len = buffer_len;
				if(fusion_candidate[bwt_index1].strand == "+")
					new_fusion_splice.first_splice.start_pos = fusion_candidate[bwt_index1].end_pos - buffer_len + 1;
				else
					new_fusion_splice.first_splice.start_pos = fusion_candidate[bwt_index1].start_pos + buffer_len - prefix_length;
				if(fusion_candidate[bwt_index2].strand == "+")
					new_fusion_splice.second_splice.start_pos = fusion_candidate[bwt_index2].start_pos + buffer_len - suffix_length;
				else
					new_fusion_splice.second_splice.start_pos = fusion_candidate[bwt_index2].end_pos - buffer_len + 1;
				string chrom_seq = left_chrom_seq.substr(0, prefix_length) + right_chrom_seq.substr(right_chrom_seq.size() - pending_seq.length() + prefix_length);
				Jump_Code left_match(prefix_length, "M");
				new_fusion_splice.first_splice.jump_code.push_back(left_match);
				Jump_Code right_match(suffix_length, "M");
				new_fusion_splice.second_splice.jump_code.push_back(right_match);
				fusion_candidate[bwt_index1].fusion_internal.push_back(new_fusion_splice);
			}
		}
		if(fusiondb)
		{
			int remap_left_bound, remap_right_bound;
			if(fusion_candidate[bwt_index1].strand == "+")
				remap_left_bound = max(0, fusion_candidate[bwt_index1].end_pos + 1 - buffer_len);
			else
				remap_left_bound = max(0, fusion_candidate[bwt_index1].start_pos - 1 + buffer_len);
			if(fusion_candidate[bwt_index2].strand == "+")
				remap_right_bound = max(0, fusion_candidate[bwt_index2].start_pos - 1 + buffer_len);
			else
				remap_right_bound = max(0, fusion_candidate[bwt_index2].end_pos + 1 - buffer_len);
			if(remap_to_fusiondb(bwt_index1, bwt_index2, fusion_candidate[bwt_index1].chrom, fusion_candidate[bwt_index2].chrom, fusion_candidate[bwt_index1].strand, fusion_candidate[bwt_index2].strand, read_start, read_end, premap_len, buffer_len, remap_left_bound, remap_right_bound, fusion_candidate[bwt_index1].fusion_internal))
			{
				if(debug)
					cout << "double anchor fusion remap fixed" << endl;
				return_fixed = true;
			}
		}
		if(debug)
			cout << "fixing fusion double anchor complete" <<endl;
		return return_fixed;
	}

	bool remap_head_fusion(size_t bwt_index, size_t anchor_segment)
	{
		if(debug)
			cout << "remap fusion head " << fusion_candidate[bwt_index].start_seg_no << "\t" << fusion_candidate[bwt_index].end_seg_no << "\t" << fusion_candidate[bwt_index].chrom << "\t" << fusion_candidate[bwt_index].strand << "\t" << fusion_candidate[bwt_index].start_pos << "\t" << fusion_candidate[bwt_index].end_pos << endl;
		bool head_fixed = false;
		Read_Block* r_block = (fusion_candidate[bwt_index].pair_no == 1) ? r_block_pe1 : r_block_pe2;
		size_t premap_len = 0;
		if(fusiondb)
		{
			int remap_left_bound = -1, remap_right_bound;
			if(fusion_candidate[bwt_index].strand == "+")
				remap_right_bound = fusion_candidate[bwt_index].start_pos - 1;
			else
				remap_right_bound = fusion_candidate[bwt_index].end_pos + 1;
			if(remap_to_fusiondb(bwt_index, bwt_index, "*", fusion_candidate[bwt_index].chrom, "*", fusion_candidate[bwt_index].strand, r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment -1), r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment) - 1, premap_len, 0, remap_left_bound, remap_right_bound, fusion_candidate[bwt_index].fusion_head))
			{
				if(debug)
					cout << "head fusion remap fixed" << endl;
				head_fixed = true;
			}
		}
		return head_fixed;
	}


	bool fix_head_fusion(size_t bwt_index, size_t anchor_segment, string chrom, int left_bound, int right_bound, bool same_strand)
	{
		if(debug)
			cout << "fixing fusion head" <<endl;	
		bool head_fixed = false;	
		Read_Block* r_block = (fusion_candidate[bwt_index].pair_no == 1) ? r_block_pe1 : r_block_pe2;
		string pending_seq = r_block->get_seg_seq(anchor_segment).substr(anchor_length - extend_len);
		if(fusion_candidate[bwt_index].strand == "+")
		{
			if(fusion_candidate[bwt_index].start_pos < (int)pending_seq.length() + 2)
				return head_fixed;
		}
		else
		{
			if(fusion_candidate[bwt_index].end_pos + 1 + (int)pending_seq.length() + 2 > (int)((*chrom_map)[fusion_candidate[bwt_index].chrom]))
				return head_fixed;
		}
		{
			string anchor_strand = (fusion_candidate[bwt_index].strand == "+") == same_strand ? "+" : "-";
			string anchor_string = r_block->get_seg_seq(anchor_segment).substr(0, anchor_length);
			int anchor_pos;
			bool mapped_fw;
			int left_bound_forward = left_bound;
			int right_bound_forward = right_bound;
			if(fusion_candidate[bwt_index].chrom == chrom)
			{
				if(left_bound_forward > fusion_candidate[bwt_index].start_pos)
				{
					if(left_bound_forward - fusion_candidate[bwt_index].end_pos < min_fusion_distance)
						left_bound_forward = fusion_candidate[bwt_index].end_pos + min_fusion_distance; 
				}
				else if(right_bound_forward < fusion_candidate[bwt_index].end_pos)
				{
					if(fusion_candidate[bwt_index].start_pos - right_bound_forward < min_fusion_distance)
						right_bound_forward = fusion_candidate[bwt_index].start_pos - min_fusion_distance;
				}
				if(left_bound_forward >= right_bound_forward)
					return head_fixed;
			}
			int right_bound_reverse = right_bound;
			if(debug)
				cout << "searching anchor " << anchor_string << " " << anchor_strand << " " << chrom << left_bound << right_bound << endl;
			if(anchor_search_head.set(anchor_string, (*reference_sequence)[chrom], anchor_strand == "+", anchor_strand == "-", max_anchor_hits))
			{
				while(anchor_search_head.search_next_simple((*reference_sequence)[chrom], left_bound_forward, right_bound_forward, right_bound_reverse, anchor_pos, mapped_fw))
				{
					if(debug)
						cout << "found head anchor at : " << anchor_pos << "\t" << anchor_string << "\t" << (*reference_sequence)[chrom].substr(anchor_pos, anchor_length) << endl; 
					size_t prefix_length = 0;
					size_t mismatch_bits = 0;
					string left_chrom_seq;
					string right_chrom_seq;
					string flank_seq;
					if(anchor_strand == "+")
					{
						if(anchor_pos + anchor_length - extend_len +(int)pending_seq.length() + 2 > (int)((*chrom_map)[chrom]))
							continue;
						left_chrom_seq = (*reference_sequence)[chrom].substr(anchor_pos + anchor_length - extend_len, pending_seq.length() + 2);	
					}
					else
					{
						if(anchor_pos - 1 + extend_len < (int)pending_seq.length() + 1)
							continue;
						left_chrom_seq = revcomp((*reference_sequence)[chrom].substr(anchor_pos + extend_len - pending_seq.length() - 2, pending_seq.length() + 2));	
					}
					if(fusion_candidate[bwt_index].strand == "+")
					{
						right_chrom_seq = (*reference_sequence)[fusion_candidate[bwt_index].chrom].substr(fusion_candidate[bwt_index].start_pos - pending_seq.length() - 2, pending_seq.length() + 2);
					}
					else
					{
						right_chrom_seq = revcomp((*reference_sequence)[fusion_candidate[bwt_index].chrom].substr(fusion_candidate[bwt_index].end_pos + 1, pending_seq.length() + 2));
					}
					if(debug)
						cout << left_chrom_seq << "\t" << right_chrom_seq << "\t" << pending_seq <<endl;
					bool anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, fusion_single_anchor_noncanonical, flank_seq);
					if(anchor_fixed && prefix_length != 0)
					{
						if(debug)
							cout << "head fusion fixed, prefix length is " << prefix_length << endl;		
						size_t suffix_length = pending_seq.length() -  prefix_length;
						Fusion_Splice new_fusion_splice;
						new_fusion_splice.flank_seq = flank_seq;
						new_fusion_splice.first_splice.chrom = chrom;
						new_fusion_splice.second_splice.chrom = fusion_candidate[bwt_index].chrom;
						new_fusion_splice.first_splice.strand = anchor_strand;
						new_fusion_splice.second_splice.strand = fusion_candidate[bwt_index].strand;
						if(anchor_strand == "+")
							new_fusion_splice.first_splice.start_pos = anchor_pos;
						else
							new_fusion_splice.first_splice.start_pos = anchor_pos + extend_len - prefix_length;
						if(fusion_candidate[bwt_index].strand == "+")
							new_fusion_splice.second_splice.start_pos = fusion_candidate[bwt_index].start_pos - suffix_length;
						else
							new_fusion_splice.second_splice.start_pos = fusion_candidate[bwt_index].end_pos + 1;
						Jump_Code left_match(prefix_length + anchor_length - extend_len, "M");
						new_fusion_splice.first_splice.jump_code.push_back(left_match);
						Jump_Code right_match(suffix_length, "M");
						new_fusion_splice.second_splice.jump_code.push_back(right_match);
						fusion_candidate[bwt_index].fusion_head.push_back(new_fusion_splice);	
						head_fixed = true;
					}
				}
			}
		}
		if(debug)
			cout << "fixing fusion head fixed: " << head_fixed << endl;
		return head_fixed;
	}

	bool remap_tail_fusion(size_t bwt_index, size_t anchor_segment)
	{
		if(debug)
			cout << "remap fusion tail " << fusion_candidate[bwt_index].start_seg_no << "\t" << fusion_candidate[bwt_index].end_seg_no << "\t" << fusion_candidate[bwt_index].chrom << "\t" << fusion_candidate[bwt_index].strand << "\t" << fusion_candidate[bwt_index].start_pos << "\t" << fusion_candidate[bwt_index].end_pos << endl;
		bool tail_fixed = false;
		Read_Block* r_block = (fusion_candidate[bwt_index].pair_no == 1) ? r_block_pe1 : r_block_pe2;
		int premap_len = 0;
		if(fusiondb)
		{
			int remap_left_bound, remap_right_bound = -1;
			if(fusion_candidate[bwt_index].strand == "+")
				remap_left_bound = fusion_candidate[bwt_index].end_pos + 1;
			else
				remap_left_bound = fusion_candidate[bwt_index].start_pos - 1;
			if(remap_to_fusiondb(bwt_index, bwt_index, fusion_candidate[bwt_index].chrom, "*", fusion_candidate[bwt_index].strand, "*", r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment -1), r_block->get_seg_len(HEAD_SEG_NO, (int)anchor_segment) - 1, premap_len, 0, remap_left_bound, remap_right_bound, fusion_candidate[bwt_index].fusion_tail))
			{
				if(debug)
					cout << "tail fusion remap fixed" << endl;
				tail_fixed = true;
			}
		}
		return tail_fixed;
	}
	
	bool fix_tail_fusion(size_t bwt_index, size_t anchor_segment, string chrom, int left_bound, int right_bound, bool same_strand)
	{
		if(debug)
			cout << "fixing fusion tail" <<endl;
		bool tail_fixed = false;
		if(debug)
			cout << "bwt vector strand is " << fusion_candidate[bwt_index].strand << ", same strand is " << same_strand <<endl;
		string anchor_strand = (fusion_candidate[bwt_index].strand == "+") == same_strand ? "+" : "-";
		Read_Block* r_block = (fusion_candidate[bwt_index].pair_no == 1) ? r_block_pe1 : r_block_pe2;
		int anchor_segment_length = r_block->get_seg_len(anchor_segment);
		string pending_seq = (*r_block).get_seg_seq(anchor_segment).substr(0, anchor_segment_length - anchor_length + extend_len);
		if(fusion_candidate[bwt_index].strand == "+")
		{
			if(fusion_candidate[bwt_index].end_pos + 1 + (int)pending_seq.length() + 2 > (int)((*chrom_map)[fusion_candidate[bwt_index].chrom]))
				return tail_fixed;
		}
		else
		{
			if(fusion_candidate[bwt_index].start_pos < (int)pending_seq.length() + 2)
				return tail_fixed;
		}
		{
			string anchor_string = r_block->get_seg_seq(anchor_segment).substr(r_block->get_seg_len(anchor_segment) - anchor_length);
			int anchor_pos;
			bool mapped_fw;
			int left_bound_forward = left_bound;
			int right_bound_forward = right_bound;
			if(fusion_candidate[bwt_index].chrom == chrom)
			{
				if(left_bound_forward > fusion_candidate[bwt_index].start_pos)
				{
					if(left_bound_forward - fusion_candidate[bwt_index].end_pos < min_fusion_distance)
						left_bound_forward = fusion_candidate[bwt_index].end_pos + min_fusion_distance; 
				}
				else if(right_bound_forward < fusion_candidate[bwt_index].end_pos)
				{
					if(fusion_candidate[bwt_index].start_pos - right_bound_forward < min_fusion_distance)
						right_bound_forward = fusion_candidate[bwt_index].start_pos - min_fusion_distance;
				}
				if(left_bound_forward >= right_bound_forward)
					return tail_fixed;
			}
			int left_bound_reverse = left_bound_forward;
			if(debug)
				cout << "searching anchor " << anchor_string << " " << anchor_strand << " " << chrom << left_bound << right_bound << endl;
			if(anchor_search_tail.set(anchor_string, (*reference_sequence)[chrom], anchor_strand == "+", anchor_strand == "-", max_anchor_hits))
			{
				while(anchor_search_tail.search_next_simple((*reference_sequence)[chrom], left_bound_forward, left_bound_reverse, right_bound_forward, anchor_pos, mapped_fw))
				{
					if(debug)
						cout << "found tail anchor at : " << anchor_pos << "\t" << anchor_string << "\t" << (*reference_sequence)[chrom].substr(anchor_pos, anchor_length) << endl; 
					size_t prefix_length = 0;
					size_t mismatch_bits = 0;
					string left_chrom_seq;
					string right_chrom_seq;
					string flank_seq;
					if(fusion_candidate[bwt_index].strand == "+")
					{
						left_chrom_seq = (*reference_sequence)[fusion_candidate[bwt_index].chrom].substr(fusion_candidate[bwt_index].end_pos + 1, pending_seq.length() + 2);
					}
					else
					{
						left_chrom_seq = revcomp((*reference_sequence)[fusion_candidate[bwt_index].chrom].substr(fusion_candidate[bwt_index].start_pos - pending_seq.length() - 2, pending_seq.length() + 2));
					}
					if(anchor_strand == "+")
					{
						if(anchor_pos + extend_len < (int)pending_seq.length() + 2)
							continue;
						right_chrom_seq = (*reference_sequence)[chrom].substr(anchor_pos + extend_len - pending_seq.length() - 2, pending_seq.length() + 2);
					}
					else
					{
						if(anchor_pos + anchor_length - extend_len + (int)pending_seq.length() + 2 > (int)((*chrom_map)[chrom]))
							continue;
						right_chrom_seq = revcomp((*reference_sequence)[chrom].substr(anchor_pos + anchor_length - extend_len, pending_seq.length() + 2));	
					}
					if(debug)
						cout << left_chrom_seq << "\t" << right_chrom_seq << "\t" << pending_seq <<endl;
					bool anchor_fixed = genome_scan.Double_anchored_score(pending_seq, left_chrom_seq, right_chrom_seq, prefix_length, max_single_splice_mismatch, mismatch_bits, fusion_single_anchor_noncanonical, flank_seq);
					if(anchor_fixed && prefix_length != pending_seq.length())
					{
						if(debug)
							cout << "tail fusion fixed, prefix length is " << prefix_length << endl;		
						size_t suffix_length = pending_seq.length() -  prefix_length;
						Fusion_Splice new_fusion_splice;
						new_fusion_splice.flank_seq = flank_seq;
						new_fusion_splice.first_splice.chrom = fusion_candidate[bwt_index].chrom;;
						new_fusion_splice.second_splice.chrom = chrom;
						new_fusion_splice.first_splice.strand = fusion_candidate[bwt_index].strand;
						new_fusion_splice.second_splice.strand = anchor_strand;
						if(fusion_candidate[bwt_index].strand == "+")
							new_fusion_splice.first_splice.start_pos = fusion_candidate[bwt_index].end_pos + 1;
						else
							new_fusion_splice.first_splice.start_pos = fusion_candidate[bwt_index].start_pos - prefix_length;
						if(anchor_strand == "+")
							new_fusion_splice.second_splice.start_pos = anchor_pos + extend_len - suffix_length;
						else
						new_fusion_splice.second_splice.start_pos = anchor_pos;
						Jump_Code left_match(prefix_length, "M");
						new_fusion_splice.first_splice.jump_code.push_back(left_match);
						Jump_Code right_match(suffix_length + anchor_length - extend_len, "M");
						new_fusion_splice.second_splice.jump_code.push_back(right_match);
						fusion_candidate[bwt_index].fusion_tail.push_back(new_fusion_splice);	
						tail_fixed = true;
					}
				}
			}
		}
		if(debug)
			cout << "fixing fusion tail fixed :" << tail_fixed <<endl;
		return tail_fixed;
	}

	void get_fusion_region(string chrom, int start, int end, vector< pair<Cluster*, bool> >& fusion_region)
	{
		for(size_t i = 0; i< (*fusion_cluster)[chrom].size(); i++)
		{
			//if(debug)
			//	cout << (*fusion_cluster)[chrom][i].start1 << "\t" << (*fusion_cluster)[chrom][i].end1 << endl;
			if((*fusion_cluster)[chrom][i].end1 < start)
				continue;
			else if((*fusion_cluster)[chrom][i].start1 > end)
				break;
			else if((*fusion_cluster)[chrom][i].start1 < start && (*fusion_cluster)[chrom][i].end1 > end)
			{
				if((*fusion_cluster)[chrom][i].chrom1 == (*fusion_cluster)[chrom][i].chrom2)
				{	
					if(max(abs((*fusion_cluster)[chrom][i].start1 - (*fusion_cluster)[chrom][i].end2), abs((*fusion_cluster)[chrom][i].start2 - (*fusion_cluster)[chrom][i].end1)) < min_fusion_distance)
						continue;
					if((*fusion_cluster)[chrom][i].start2 < start && (*fusion_cluster)[chrom][i].end2 > end)
						continue;
					/*if((*fusion_cluster)[chrom][i].start1 <= (*fusion_cluster)[chrom][i].end2 && (*fusion_cluster)[chrom][i].start2 <= (*fusion_cluster)[chrom][i].end1)
						continue;
					if(min(abs((*fusion_cluster)[chrom][i].start2 - (*fusion_cluster)[chrom][i].end1), abs((*fusion_cluster)[chrom][i].start1 - (*fusion_cluster)[chrom][i].end2)) < min_fusion_distance)
						continue;*/
				}
				fusion_region.push_back(make_pair(&((*fusion_cluster)[chrom][i]), false));
				if(debug)
					cout << (*fusion_cluster)[chrom][i].chrom1 << "\t" << (*fusion_cluster)[chrom][i].start1 << "\t" << (*fusion_cluster)[chrom][i].end1 << "\t" << (*fusion_cluster)[chrom][i].chrom2 << "\t" << (*fusion_cluster)[chrom][i].start2 << "\t" << (*fusion_cluster)[chrom][i].end2 << endl;
			}
		}
	}

	bool check_overlap(int start1, int end1, int start2, int end2)
	{
		if(start1 > end2 || start2 > end1)
			return false;
		else
			return true;
	}

	bool in_fusion_cluster(Bwtmap_Info& bwt1, Bwtmap_Info& bwt2)
	{
		bool segment_same_strand = (bwt1.strand == bwt2.strand);
		for(size_t i = 0; i < (*fusion_cluster)[bwt1.chrom].size(); i++)
		{
			if(segment_same_strand == (*fusion_cluster)[bwt1.chrom][i].same_strand \
				&& in_range(bwt1.chrom, bwt1.start, bwt1.end, (*fusion_cluster)[bwt1.chrom][i].chrom1, (*fusion_cluster)[bwt1.chrom][i].start1, (*fusion_cluster)[bwt1.chrom][i].end1) \
				&& in_range(bwt2.chrom, bwt2.start, bwt2.end, (*fusion_cluster)[bwt1.chrom][i].chrom2, (*fusion_cluster)[bwt1.chrom][i].start2, (*fusion_cluster)[bwt1.chrom][i].end2))
			{
				return true;	
			}	
		}
		return false;
	}

	bool near_normal_alignment(Sam_Info& fusion_can)
	{
		for(size_t i = 0; i < sam_vec.size(); i++)
		{
			if(debug)
				cout << "checking near alignment: " << sam_vec[i].chrom << "\t" << sam_vec[i].start_pos << "\t" << sam_vec[i].end_pos << endl;
			if(fusion_can.chrom == sam_vec[i].chrom && (abs(fusion_can.start_pos - sam_vec[i].end_pos) < max_pair_dist || abs(fusion_can.end_pos - sam_vec[i].start_pos) < max_pair_dist))
				return true;
		}
		return false;
	}

	bool near_normal_alignment(vector< pair<Cluster*, bool> >& fusion_region)
	{
		bool near_normal_alignment = false;
		for(size_t i = 0; i < sam_vec.size(); i++)
		{
			for(size_t j = 0; j < fusion_region.size(); j++)
			{
				if(debug)
					cout << "checking near alignment: " << sam_vec[i].chrom << "\t" << sam_vec[i].start_pos << "\t" << sam_vec[i].end_pos << "\t" << fusion_region[j].first->chrom2 << "\t" << fusion_region[j].first->start2 << "\t" << fusion_region[j].first->end2 << endl;
				if(in_range(sam_vec[i].chrom, sam_vec[i].start_pos, sam_vec[i].end_pos, fusion_region[j].first->chrom2, fusion_region[j].first->start2, fusion_region[j].first->end2))
				{
					near_normal_alignment = true;
					fusion_region[j].second = true;
				}
			}
		}
		return near_normal_alignment;
	}

	void search_fusion_double_anchor(bool search_pe1, bool search_pe2)
	{
		if(debug)
			cout << "searching fusion double anchor" << endl;
		for (size_t i = 0; i < fusion_candidate.size(); i++)
		{
			if((!search_pe1 && fusion_candidate[i].pair_no == 1) || (!search_pe2 && fusion_candidate[i].pair_no == 2))
				continue;
			for (size_t j = i + 1; j < fusion_candidate.size(); j++)
			{
				if(fusion_candidate[i].pair_no != fusion_candidate[j].pair_no)
					continue;
				Read_Block* r_block = (fusion_candidate[i].pair_no == 1) ? r_block_pe1 : r_block_pe2;
				int segment_number = (fusion_candidate[i].pair_no == 1 ? pe1_segment_number : pe2_segment_number);
				if(debug)
				{
					cout << fusion_candidate[i].start_seg_no << "\t" << fusion_candidate[i].end_seg_no << "\t" << fusion_candidate[i].chrom << "\t" << fusion_candidate[i].strand << "\t" << fusion_candidate[i].start_pos << "\t" << fusion_candidate[i].end_pos \
					<< "\t" << fusion_candidate[j].start_seg_no << "\t" << fusion_candidate[j].end_seg_no << "\t" << fusion_candidate[j].chrom << "\t" << fusion_candidate[j].strand << "\t" << fusion_candidate[j].start_pos << "\t" << fusion_candidate[j].end_pos << endl;
					for(size_t x = 0; x < fusion_candidate[i].jump_code.size(); x++)
						cout << fusion_candidate[i].jump_code[x].toString();
					cout << "\t";
					for(size_t x = 0; x < fusion_candidate[j].jump_code.size(); x++)
						cout << fusion_candidate[j].jump_code[x].toString();
					cout << endl;
				}
				if(fusion_candidate[i].strand == "+" && fusion_candidate[j].strand == "+")
				{
					size_t index1;
					size_t index2;
					bool candidate = false;
					if(fusion_candidate[i].end_seg_no < fusion_candidate[j].start_seg_no && fusion_candidate[i].start_seg_no == HEAD_SEG_NO && fusion_candidate[j].end_seg_no == segment_number)
					{
						int max_gap = (fusion_candidate[j].start_seg_no - fusion_candidate[i].end_seg_no - 1) * segment_length + max_intron_length_double_anchor;
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > max_gap) // || fusion_candidate[j].start_pos <= fusion_candidate[i].end_pos
						{
							candidate = true;
							index1 = i;
							index2 = j;
						}
					}
					else if(fusion_candidate[j].end_seg_no < fusion_candidate[i].start_seg_no  && fusion_candidate[j].start_seg_no == HEAD_SEG_NO && fusion_candidate[i].end_seg_no == segment_number)
					{
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > min_fusion_distance)
						{
							candidate = true;
							index1 = j;
							index2 = i;
						}
					}
					if(candidate)
					{
						int seg_gap = fusion_candidate[index2].start_seg_no - fusion_candidate[index1].end_seg_no;
						if(seg_gap == 1)
						{
							int trunc_left_len = r_block->get_seg_len(fusion_candidate[index1].end_seg_no) - trunc_len;
							string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].end_seg_no).substr(trunc_left_len) + r_block->get_seg_seq(fusion_candidate[index2].start_seg_no).substr(0,trunc_len);
							fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].end_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].start_seg_no - 1) + trunc_len - 1, mid_seq, fusion_candidate[index1].end_pos - trunc_len + 1, fusion_candidate[index2].start_pos + trunc_len - 1, trunc_len);		
						}
						else if(seg_gap == 2)
						{
							int extend_left_len = r_block->get_seg_len(fusion_candidate[index1].end_seg_no) - extend_len;
							string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].end_seg_no).substr(extend_left_len) + r_block->get_seg_seq(fusion_candidate[index1].end_seg_no + 1);
							mid_seq.append(r_block->get_seg_seq(fusion_candidate[index2].start_seg_no).substr(0,extend_len));
							fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].end_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].start_seg_no - 1) + extend_len - 1, mid_seq, fusion_candidate[index1].end_pos - extend_len + 1, fusion_candidate[index2].start_pos + extend_len - 1, extend_len);
						}
					}
				}
				else if(fusion_candidate[i].strand == "-" && fusion_candidate[j].strand == "-")
				{
					size_t index1;
					size_t index2;
					bool candidate = false;
					if(fusion_candidate[j].end_seg_no > fusion_candidate[i].start_seg_no && fusion_candidate[i].end_seg_no == HEAD_SEG_NO && fusion_candidate[j].start_seg_no == segment_number)
					{
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > min_fusion_distance)
						{
							candidate = true;
							index1 = i;
							index2 = j;
						}
					}
					else if(fusion_candidate[i].end_seg_no > fusion_candidate[j].start_seg_no && fusion_candidate[j].end_seg_no == HEAD_SEG_NO && fusion_candidate[i].start_seg_no == segment_number)
					{
						int max_gap = (fusion_candidate[i].end_seg_no - fusion_candidate[j].start_seg_no - 1) * segment_length + max_intron_length_double_anchor;
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > max_gap)  //|| fusion_candidate[j].start_pos <= fusion_candidate[i].end_pos
						{
							candidate = true;
							index1 = j;
							index2 = i;
						}
					}
					if(candidate)
					{
						int seg_gap = fusion_candidate[index2].end_seg_no - fusion_candidate[index1].start_seg_no;
						if(seg_gap == 1)
						{
							int trunc_left_len = r_block->get_seg_len(fusion_candidate[index1].start_seg_no) - trunc_len;
							string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].start_seg_no).substr(trunc_left_len) + r_block->get_seg_seq(fusion_candidate[index2].end_seg_no).substr(0,trunc_len);
							fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].start_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].end_seg_no - 1) + trunc_len - 1, mid_seq, fusion_candidate[index1].start_pos + trunc_len - 1, fusion_candidate[index2].end_pos - trunc_len + 1, trunc_len);		
						}
						else if(seg_gap == 2)
						{
							int extend_left_len = r_block->get_seg_len(fusion_candidate[index1].start_seg_no) - extend_len;
							string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].start_seg_no).substr(extend_left_len) + r_block->get_seg_seq(fusion_candidate[index1].start_seg_no + 1);
							mid_seq.append(r_block->get_seg_seq(fusion_candidate[index2].end_seg_no).substr(0,extend_len));
							fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].start_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].end_seg_no - 1) + extend_len - 1, mid_seq, fusion_candidate[index1].start_pos + extend_len - 1, fusion_candidate[index2].end_pos - extend_len + 1, extend_len);
						}
					}	
				}
				else if(fusion_candidate[i].strand == "+" && fusion_candidate[j].strand == "-")
				{
					size_t index1;
					size_t index2;
					int seg_gap;
					if(fusion_candidate[i].end_seg_no < fusion_candidate[j].end_seg_no && fusion_candidate[i].start_seg_no == HEAD_SEG_NO && fusion_candidate[j].start_seg_no == segment_number)
					{
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > min_fusion_distance)
						{
							index1 = i;
							index2 = j;
							seg_gap = fusion_candidate[j].end_seg_no - fusion_candidate[i].end_seg_no;
							if(seg_gap == 1)
							{
								int trunc_left_len = r_block->get_seg_len(fusion_candidate[index1].end_seg_no) - trunc_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].end_seg_no).substr(trunc_left_len) + r_block->get_seg_seq(fusion_candidate[index2].end_seg_no).substr(0, trunc_len);
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].end_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].end_seg_no - 1) + trunc_len - 1, mid_seq, fusion_candidate[index1].end_pos - trunc_len + 1, fusion_candidate[index2].end_pos - trunc_len + 1, trunc_len);		
							}
							else if(seg_gap == 2)
							{
								int extend_left_len = r_block->get_seg_len(fusion_candidate[index1].end_seg_no) - extend_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].end_seg_no).substr(extend_left_len) + r_block->get_seg_seq(fusion_candidate[index1].end_seg_no + 1);
								mid_seq.append(r_block->get_seg_seq(fusion_candidate[index2].end_seg_no).substr(0,extend_len));
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].end_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].end_seg_no - 1) + extend_len - 1, mid_seq, fusion_candidate[index1].end_pos - extend_len + 1, fusion_candidate[index2].end_pos - extend_len + 1, extend_len);
							}
						}
					}
					else if(fusion_candidate[j].start_seg_no < fusion_candidate[i].start_seg_no && fusion_candidate[j].end_seg_no == HEAD_SEG_NO && fusion_candidate[i].end_seg_no == segment_number)
					{
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > min_fusion_distance)
						{
							index1 = j;
							index2 = i;
							seg_gap = fusion_candidate[i].start_seg_no - fusion_candidate[j].start_seg_no;
							if(seg_gap == 1)
							{
								int trunc_left_len = r_block->get_seg_len(fusion_candidate[index1].start_seg_no) - trunc_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].start_seg_no).substr(trunc_left_len) + r_block->get_seg_seq(fusion_candidate[index2].start_seg_no).substr(0, trunc_len);
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].start_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].start_seg_no - 1) + trunc_len - 1, mid_seq, fusion_candidate[index1].start_pos + trunc_len - 1, fusion_candidate[index2].start_pos + trunc_len - 1, trunc_len);		
							}
							else if(seg_gap == 2)
							{
								int extend_left_len = r_block->get_seg_len(fusion_candidate[index1].start_seg_no) - extend_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].start_seg_no).substr(extend_left_len) + r_block->get_seg_seq(fusion_candidate[index1].start_seg_no+1);
								mid_seq.append(r_block->get_seg_seq(fusion_candidate[index2].start_seg_no).substr(0,extend_len));
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].start_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].start_seg_no - 1) + extend_len - 1, mid_seq, fusion_candidate[index1].start_pos + extend_len - 1, fusion_candidate[index2].start_pos + extend_len - 1, extend_len);
							}
						}
					}
				}
				else
				{
					size_t index1;
					size_t index2;
					int seg_gap;
					if(fusion_candidate[i].start_seg_no < fusion_candidate[j].start_seg_no && fusion_candidate[i].end_seg_no == HEAD_SEG_NO && fusion_candidate[j].end_seg_no == segment_number)
					{
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > min_fusion_distance)
						{
							index1 = i;
							index2 = j;
							seg_gap = fusion_candidate[j].start_seg_no - fusion_candidate[i].start_seg_no;
							if(seg_gap == 1)
							{
								int trunc_left_len = r_block->get_seg_len(fusion_candidate[index1].start_seg_no) - trunc_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].start_seg_no).substr(trunc_left_len) + r_block->get_seg_seq(fusion_candidate[index2].start_seg_no).substr(0, trunc_len);
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].start_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].start_seg_no - 1) + trunc_len - 1, mid_seq, fusion_candidate[index1].start_pos + trunc_len - 1, fusion_candidate[index2].start_pos + trunc_len - 1, trunc_len);		
							}
							else if(seg_gap == 2)
							{
								int extend_left_len = r_block->get_seg_len(fusion_candidate[index1].start_seg_no) - extend_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].start_seg_no).substr(extend_left_len) + r_block->get_seg_seq(fusion_candidate[index1].start_seg_no + 1);
								mid_seq.append(r_block->get_seg_seq(fusion_candidate[index2].start_seg_no).substr(0,extend_len));
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].start_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].start_seg_no - 1) + extend_len - 1, mid_seq, fusion_candidate[index1].start_pos + extend_len - 1, fusion_candidate[index2].start_pos + extend_len - 1, extend_len);
							}
						}
					}
					else if(fusion_candidate[j].end_seg_no < fusion_candidate[i].end_seg_no && fusion_candidate[j].start_seg_no == HEAD_SEG_NO && fusion_candidate[i].start_seg_no == segment_number)
					{
						if(fusion_candidate[i].chrom != fusion_candidate[j].chrom || fusion_candidate[j].start_pos - fusion_candidate[i].end_pos - 1 > min_fusion_distance)
						{
							index1 = j;
							index2 = i;
							seg_gap = fusion_candidate[i].end_seg_no - fusion_candidate[j].end_seg_no;
							if(seg_gap == 1)
							{
								int trunc_left_len = r_block->get_seg_len(fusion_candidate[index1].end_seg_no) - trunc_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].end_seg_no).substr(trunc_left_len) + r_block->get_seg_seq(fusion_candidate[index2].end_seg_no).substr(0, trunc_len);
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].end_seg_no - 1) + trunc_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].end_seg_no - 1) + trunc_len - 1, mid_seq, fusion_candidate[index1].end_pos - trunc_len + 1, fusion_candidate[index2].end_pos - trunc_len + 1, trunc_len);		
							}
							else if(seg_gap == 2)
							{
								int extend_left_len = r_block->get_seg_len(fusion_candidate[index1].end_seg_no) - extend_len;
								string mid_seq = r_block->get_seg_seq(fusion_candidate[index1].end_seg_no).substr(extend_left_len) + r_block->get_seg_seq(fusion_candidate[index1].end_seg_no + 1);
								mid_seq.append(r_block->get_seg_seq(fusion_candidate[index2].end_seg_no).substr(0,extend_len));
								fix_double_anchor_fusion(index1, index2, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index1].end_seg_no - 1) + extend_left_len, r_block->get_seg_len(HEAD_SEG_NO, fusion_candidate[index2].end_seg_no - 1) + extend_len - 1, mid_seq, fusion_candidate[index1].end_pos - extend_len + 1, fusion_candidate[index2].end_pos - extend_len + 1, extend_len);
							}
						}
					}
				}
			}
		}
	}

	void search_fusion_single_anchor(bool remap_pe1, bool remap_pe2, bool search_pe1, bool search_pe2)
	{
		for (size_t i = 0; i < fusion_candidate.size(); i++)
		{
			if(debug)
				cout << "searching single anchor fusion for: " << fusion_candidate[i].chrom << "\t" << fusion_candidate[i].strand << "\t" << fusion_candidate[i].start_pos << "\t" << fusion_candidate[i].end_pos << "\t" << fusion_candidate[i].start_seg_no << "\t" << fusion_candidate[i].end_seg_no << endl;
			int segment_number = (fusion_candidate[i].pair_no == 1 ? pe1_segment_number : pe2_segment_number);		
			if((fusion_candidate[i].pair_no == 1 && remap_pe1) || (fusion_candidate[i].pair_no == 2 && remap_pe2))
			{
				if(fusion_candidate[i].strand == "+")
				{
					if(fusion_candidate[i].start_seg_no == HEAD_SEG_NO && fusion_candidate[i].end_seg_no == segment_number - 1) // remap tail +
						remap_tail_fusion(i, fusion_candidate[i].end_seg_no + 1);
					else if(fusion_candidate[i].start_seg_no == HEAD_SEG_NO + 1 && fusion_candidate[i].end_seg_no == segment_number) // remap head +
						remap_head_fusion(i, fusion_candidate[i].start_seg_no - 1);
				}
				else
				{
					if(fusion_candidate[i].start_seg_no == segment_number && fusion_candidate[i].end_seg_no == HEAD_SEG_NO + 1) // remap head -
						remap_head_fusion(i, fusion_candidate[i].end_seg_no - 1);
					else if(fusion_candidate[i].start_seg_no == segment_number - 1 && fusion_candidate[i].end_seg_no == HEAD_SEG_NO) // remap tail -
						remap_tail_fusion(i, fusion_candidate[i].start_seg_no + 1);	
				}
			}
			if((fusion_candidate[i].pair_no == 1 && search_pe1) || (fusion_candidate[i].pair_no == 2 && search_pe2))
			{
				vector< pair<Cluster*, bool> > fusion_region;
				get_fusion_region(fusion_candidate[i].chrom, fusion_candidate[i].start_pos, fusion_candidate[i].end_pos, fusion_region);
				if(debug)
					cout << "fusion cluster size is " << fusion_region.size() << endl;
				if(fusion_region.size() == 0)
					return;
				bool segment_near_normal_alignment = near_normal_alignment(fusion_candidate[i]);
				bool fusion_region_near_normal_alignent = false;
				if(!segment_near_normal_alignment)
				 fusion_region_near_normal_alignent = near_normal_alignment(fusion_region);
				if(!segment_near_normal_alignment && !fusion_region_near_normal_alignent)
					continue;
				for(size_t j = 0; j < fusion_region.size(); j++)
				{
					if(segment_near_normal_alignment || fusion_region[j].second)
					{
						if(debug)
							cout << "fixing single anchor for fusion region " << fusion_region[j].first->chrom1 << "\t" << fusion_region[j].first->start1 << "\t" << fusion_region[j].first->end1 << "\t" << fusion_region[j].first->chrom2 << "\t" << fusion_region[j].first->start2 << "\t" << fusion_region[j].first->end2<< endl;
						bool same_strand = (fusion_region[j].first->strand1 == fusion_region[j].first->strand2);
						if(fusion_candidate[i].strand == "+")
						{
							if(fusion_candidate[i].start_seg_no == HEAD_SEG_NO && fusion_candidate[i].end_seg_no == segment_number - 1) // fix tail +
								fix_tail_fusion(i, fusion_candidate[i].end_seg_no + 1, fusion_region[j].first->chrom2, fusion_region[j].first->start2, fusion_region[j].first->end2, same_strand);
							else if(fusion_candidate[i].start_seg_no == HEAD_SEG_NO + 1 && fusion_candidate[i].end_seg_no == segment_number) // fix head +
								fix_head_fusion(i, fusion_candidate[i].start_seg_no - 1, fusion_region[j].first->chrom2, fusion_region[j].first->start2, fusion_region[j].first->end2, same_strand);
						}	
						else
						{
							if(fusion_candidate[i].start_seg_no == segment_number && fusion_candidate[i].end_seg_no == HEAD_SEG_NO + 1) // fix head -
								fix_head_fusion(i, fusion_candidate[i].end_seg_no - 1, fusion_region[j].first->chrom2, fusion_region[j].first->start2, fusion_region[j].first->end2, same_strand);
							else if(fusion_candidate[i].start_seg_no == segment_number - 1 && fusion_candidate[i].end_seg_no == HEAD_SEG_NO) // fix tail -
								fix_tail_fusion(i, fusion_candidate[i].start_seg_no + 1, fusion_region[j].first->chrom2, fusion_region[j].first->start2, fusion_region[j].first->end2, same_strand);	
						}
					}
				}
				fusion_region.clear();
			}
		}
	}

	void combine_jumpcode(vector<Jump_Code>& jump_code1, vector<Jump_Code>& jump_code2, vector<Jump_Code>& jump_code) 
	{
		for(size_t i = 0; i < jump_code1.size(); i++)
		{
			if(i != jump_code1.size() - 1)	
				jump_code.push_back(jump_code1[i]);
			else
			{
				Jump_Code new_jc(jump_code1[i].len + jump_code2[0].len, jump_code1[i].type);
				jump_code.push_back(new_jc);
			}
		}
		for(size_t i = 1; i < jump_code2.size(); i++)
				jump_code.push_back(jump_code2[i]);
	}

	void combine_jumpcode(vector<Jump_Code>& jump_code1, vector<Jump_Code>& jump_code2, int buffer_len, vector<Jump_Code>& jump_code) 
	{
		for(size_t i = 0; i < jump_code1.size(); i++)
		{
			if(i != jump_code1.size() - 1)	
				jump_code.push_back(jump_code1[i]);
			else
			{
				Jump_Code new_jc(jump_code1[i].len + jump_code2[0].len - buffer_len, jump_code1[i].type);
				jump_code.push_back(new_jc);
			}
		}
		for(size_t i = 1; i < jump_code2.size(); i++)
				jump_code.push_back(jump_code2[i]);
	}
	
	string jump_code_to_str(vector<Jump_Code>& jump_code)
	{
		string jc_str;
		for(size_t i = 0; i < jump_code.size(); i++)
			jc_str += jump_code[i].toStringConvertDeletion(max_del);
		return jc_str;
	}

	string get_chrom_sequence(string& chrom_seq, string& read_seq, int start_pos, vector<Jump_Code>& jump_code)
	{
		int end_pos = start_pos - 1;
		int mapped_len = 0;
		string map_seq;
		for(size_t i = 0; i < jump_code.size(); i++)
		{
			if(jump_code[i].type == "M")
			{
				map_seq.append(chrom_seq.substr(end_pos + 1, jump_code[i].len));
				end_pos += jump_code[i].len;
				mapped_len += jump_code[i].len;
			}
			else if(jump_code[i].type == "N")
			{
				end_pos += jump_code[i].len;
			}
			else   // "I"
			{
				map_seq.append(read_seq.substr(mapped_len, jump_code[i].len));
				mapped_len += jump_code[i].len;
			}
		}
		return map_seq;
	}
					
	void report_fusion()
	{
		if(debug)
			cout << "report fusion" << endl;
		for(size_t i = 0; i < fusion_candidate.size(); i++)
		{
			if(fusion_candidate[i].fusion_internal.size() == 0 && fusion_candidate[i].fusion_head.size() == 0 && fusion_candidate[i].fusion_tail.size() == 0)
				continue;
			if(fusion_candidate[i].pair_no == 1)
			{
				pe1_aligned = true;
				pe1_fusion = true;
			}
			else
			{
				pe2_aligned = true;
				pe2_fusion = true;
			}
			Read_Block* r_block = (fusion_candidate[i].pair_no == 1) ? r_block_pe1 : r_block_pe2;
			int read_length = fusion_candidate[i].pair_no == 1 ? pe1_read_length : pe2_read_length;
			int segment_number = (fusion_candidate[i].pair_no == 1 ? pe1_segment_number : pe2_segment_number);
			string sequence;
			for(size_t j = HEAD_SEG_NO; j <= (size_t)segment_number; j++)
			{
				sequence.append(r_block->get_seg_seq((int)j));
			}
			string quality_string;
			if(!fa)
			{
				for(int j = HEAD_SEG_NO; j <= segment_number; j++)
				{
					quality_string.append(r_block->get_seg_qual(j));
				}
			}
			else
				quality_string = global_qual_string.substr(0, read_length);
			for(size_t j = 0; j < fusion_candidate[i].fusion_internal.size(); j++)
			{
				if(debug)
					cout << "print fusion for fusion candidate " << i << ", fusion internal " << j << endl;
				int index1 = i;
				int index2 = fusion_candidate[i].fusion_internal[j].first_splice.end_contig;
				int start1 = (fusion_candidate[index1].strand == "+" ? fusion_candidate[index1].start_pos : fusion_candidate[index1].fusion_internal[j].first_splice.start_pos);
				int match1 = fusion_candidate[index1].mapped_len + fusion_candidate[index1].fusion_internal[j].first_splice.GetMapLen() - fusion_candidate[index1].fusion_internal[j].first_splice.buffer_len;
				int start2 = (fusion_candidate[index2].strand == "+" ? fusion_candidate[index1].fusion_internal[j].second_splice.start_pos : fusion_candidate[index2].start_pos);
				size_t total_mismatch = 0; 
				string map_sequence1;
				string map_sequence2;
				string map_sequence;
				string map_qual1;
				string map_qual2;
				string map_qual;
				string chrom_sequence;
				vector<Jump_Code> jump_code1;
				vector<Jump_Code> jump_code2;
				string jump_code1_str;
				string jump_code2_str;
				if(fusion_candidate[index1].strand == "+" && fusion_candidate[index2].strand == "+")
				{
					map_sequence1 = sequence.substr(0, match1);
					map_sequence2 = sequence.substr(match1);
					if(!fa)
					{
						map_qual1 = quality_string.substr(0, match1);
						map_qual2 = quality_string.substr(match1);
					}
					combine_jumpcode(fusion_candidate[index1].jump_code, fusion_candidate[index1].fusion_internal[j].first_splice.jump_code, fusion_candidate[index1].fusion_internal[j].first_splice.buffer_len, jump_code1);
					combine_jumpcode(fusion_candidate[index1].fusion_internal[j].second_splice.jump_code, fusion_candidate[index2].jump_code, fusion_candidate[index1].fusion_internal[j].second_splice.buffer_len, jump_code2);	
				}
				else if(fusion_candidate[index1].strand == "-" && fusion_candidate[index2].strand == "-")
				{
					map_sequence1 = revcomp(sequence.substr(0, match1));
					map_sequence2 = revcomp(sequence.substr(match1));
					if(!fa)
					{
						map_qual1 = revqual(quality_string.substr(0, match1));
						map_qual2 = revqual(quality_string.substr(match1));
					}
					combine_jumpcode(fusion_candidate[index1].fusion_internal[j].first_splice.jump_code, fusion_candidate[index1].jump_code, fusion_candidate[index1].fusion_internal[j].first_splice.buffer_len, jump_code1);
					combine_jumpcode(fusion_candidate[index2].jump_code, fusion_candidate[index1].fusion_internal[j].second_splice.jump_code, fusion_candidate[index1].fusion_internal[j].second_splice.buffer_len, jump_code2);	
				}
				else if(fusion_candidate[index1].strand == "+" && fusion_candidate[index2].strand == "-")
				{
					map_sequence1 = sequence.substr(0, match1);
					map_sequence2 = revcomp(sequence.substr(match1));
					if(!fa)
					{
						map_qual1 = quality_string.substr(0, match1);
						map_qual2 = revqual(quality_string.substr(match1));
					}
					combine_jumpcode(fusion_candidate[index1].jump_code, fusion_candidate[index1].fusion_internal[j].first_splice.jump_code, fusion_candidate[index1].fusion_internal[j].first_splice.buffer_len, jump_code1);
					combine_jumpcode(fusion_candidate[index2].jump_code, fusion_candidate[index1].fusion_internal[j].second_splice.jump_code, fusion_candidate[index1].fusion_internal[j].second_splice.buffer_len, jump_code2);	
				}
				else
				{
					map_sequence1 = revcomp(sequence.substr(0, match1));
					map_sequence2 = sequence.substr(match1);
					if(!fa)
					{
						map_qual1 = revqual(quality_string.substr(0, match1));
						map_qual2 = quality_string.substr(match1);
					}
					combine_jumpcode(fusion_candidate[index1].fusion_internal[j].first_splice.jump_code, fusion_candidate[index1].jump_code, fusion_candidate[index1].fusion_internal[j].first_splice.buffer_len, jump_code1);
					combine_jumpcode(fusion_candidate[index1].fusion_internal[j].second_splice.jump_code, fusion_candidate[index2].jump_code, fusion_candidate[index1].fusion_internal[j].second_splice.buffer_len, jump_code2);	
				}
				map_sequence = map_sequence1 + map_sequence2;
				if(!fa)
					map_qual = map_qual1 + map_qual2;
				else
					map_qual = quality_string;
				jump_code1_str = jump_code_to_str(jump_code1);
				jump_code2_str = jump_code_to_str(jump_code2);
				chrom_sequence = get_chrom_sequence((*reference_sequence)[fusion_candidate[i].fusion_internal[j].first_splice.chrom], map_sequence1, start1, jump_code1) + get_chrom_sequence((*reference_sequence)[fusion_candidate[i].fusion_internal[j].second_splice.chrom], map_sequence2, start2, jump_code2);
				for(size_t k = 0; k < map_sequence.length(); k++)
				{
					if(map_sequence[k] != chrom_sequence[k])
						total_mismatch ++ ;
				}
				string flank_seq;
				if(fusion_candidate[index1].strand == "+")
				{
					for(size_t x = 0; x < fusion_candidate[index1].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index1].junc_flank_seq[x];
						flank_seq += ",";
					}
					for(size_t x = 0; x < fusion_candidate[index1].fusion_internal[j].first_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index1].fusion_internal[j].first_splice.junc_flank_seq[x];
						flank_seq += ",";
					}
				}
				else
				{
					for(size_t x = 0; x < fusion_candidate[index1].fusion_internal[j].first_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index1].fusion_internal[j].first_splice.junc_flank_seq[x];
						flank_seq += ",";
					}	
					for(size_t x = 0; x < fusion_candidate[index1].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index1].junc_flank_seq[x];
						flank_seq += ",";
					}	
				}
				flank_seq += fusion_candidate[index1].fusion_internal[j].flank_seq;
				flank_seq += ",";
				if(fusion_candidate[index2].strand == "+")
				{
					for(size_t x = 0; x < fusion_candidate[index1].fusion_internal[j].second_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index1].fusion_internal[j].second_splice.junc_flank_seq[x];
						flank_seq += ",";
					}	
					for(size_t x = 0; x < fusion_candidate[index2].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index2].junc_flank_seq[x];
						flank_seq += ",";
					}
				}
				else
				{
					for(size_t x = 0; x < fusion_candidate[index2].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index2].junc_flank_seq[x];
						flank_seq += ",";
					}
					for(size_t x = 0; x < fusion_candidate[index1].fusion_internal[j].second_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[index1].fusion_internal[j].second_splice.junc_flank_seq[x];
						flank_seq += ",";
					}	
				}
				char output_buf[5000];
				sprintf(output_buf, "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\tNM:i:%zu\tXF:Z:%s",
					r_block->read_id.c_str(), fusion_candidate[i].fusion_internal[j].first_splice.chrom.c_str(), fusion_candidate[i].fusion_internal[j].first_splice.strand.c_str(), start1 + 1, jump_code1_str.c_str(), \
					fusion_candidate[i].fusion_internal[j].second_splice.chrom.c_str(), fusion_candidate[i].fusion_internal[j].second_splice.strand.c_str(), start2 + 1, jump_code2_str.c_str(), \
					map_sequence.c_str(), map_qual.c_str(), total_mismatch, flank_seq.c_str());
				pthread_mutex_lock (fusion_fs_lock);
				(*fusion_fs) << output_buf << endl;
				pthread_mutex_unlock (fusion_fs_lock);
			}
			for(size_t j = 0; j <fusion_candidate[i].fusion_head.size(); j++)
			{		
				if(debug)
					cout << "print fusion for fusion candidate " << i << ", fusion head " << j << endl;
				int start1 = fusion_candidate[i].fusion_head[j].first_splice.start_pos;
				int match1 = fusion_candidate[i].fusion_head[j].first_splice.GetMapLen();
				int start2 = (fusion_candidate[i].strand == "+" ? fusion_candidate[i].fusion_head[j].second_splice.start_pos : fusion_candidate[i].start_pos);
				size_t total_mismatch = 0;
				string map_sequence1;
				string map_sequence2;
				string map_sequence;
				string map_qual1;
				string map_qual2;
				string map_qual;
				string chrom_sequence;
				vector<Jump_Code> jump_code1;
				vector<Jump_Code> jump_code2;
				string jump_code1_str;
				string jump_code2_str;
				for(size_t x = 0; x < fusion_candidate[i].fusion_head[j].first_splice.jump_code.size(); x++)
						jump_code1.push_back(fusion_candidate[i].fusion_head[j].first_splice.jump_code[x]);
				if(fusion_candidate[i].fusion_head[j].first_splice.strand == "+" && fusion_candidate[i].fusion_head[j].second_splice.strand == "+")
				{
					map_sequence1 = sequence.substr(0, match1);
					map_sequence2 = sequence.substr(match1);
					if(!fa)
					{
						map_qual1 = quality_string.substr(0, match1);
						map_qual2 = quality_string.substr(match1);
					}
					combine_jumpcode(fusion_candidate[i].fusion_head[j].second_splice.jump_code, fusion_candidate[i].jump_code, jump_code2);
				}
				else if(fusion_candidate[i].fusion_head[j].first_splice.strand == "-" && fusion_candidate[i].fusion_head[j].second_splice.strand == "-")
				{
					map_sequence1 = revcomp(sequence.substr(0, match1));
					map_sequence2 = revcomp(sequence.substr(match1));
					if(!fa)
					{
						map_qual1 = revqual(quality_string.substr(0, match1));
						map_qual2 = revqual(quality_string.substr(match1));
					}
					combine_jumpcode(fusion_candidate[i].jump_code,  fusion_candidate[i].fusion_head[j].second_splice.jump_code, jump_code2);
				}
				else if(fusion_candidate[i].fusion_head[j].first_splice.strand == "+" && fusion_candidate[i].fusion_head[j].second_splice.strand == "-")
				{
					map_sequence1 = sequence.substr(0, match1);
					map_sequence2 = revcomp(sequence.substr(match1));
					if(!fa)
					{
						map_qual1 = quality_string.substr(0, match1);
						map_qual2 = revqual(quality_string.substr(match1));
					}
					combine_jumpcode(fusion_candidate[i].jump_code, fusion_candidate[i].fusion_head[j].second_splice.jump_code, jump_code2);	
				}
				else
				{
					map_sequence1 = revcomp(sequence.substr(0, match1));
					map_sequence2 = sequence.substr(match1);
					if(!fa)
					{
						map_qual1 = revqual(quality_string.substr(0, match1));
						map_qual2 = quality_string.substr(match1);
					}
					combine_jumpcode(fusion_candidate[i].fusion_head[j].second_splice.jump_code, fusion_candidate[i].jump_code, jump_code2);
				}
				map_sequence = map_sequence1 + map_sequence2;
				if(!fa)
					map_qual = map_qual1 + map_qual2;
				else
					map_qual = quality_string;
				jump_code1_str = jump_code_to_str(jump_code1);
				jump_code2_str = jump_code_to_str(jump_code2);
				chrom_sequence = get_chrom_sequence((*reference_sequence)[fusion_candidate[i].fusion_head[j].first_splice.chrom], map_sequence1, start1, jump_code1) + get_chrom_sequence((*reference_sequence)[fusion_candidate[i].fusion_head[j].second_splice.chrom], map_sequence2, start2, jump_code2);
				for(size_t k = 0; k < map_sequence.length(); k++)
				{
					if(map_sequence[k] != chrom_sequence[k])
						total_mismatch ++ ;
				}
				string flank_seq;
				for(size_t x = 0; x < fusion_candidate[i].fusion_head[j].first_splice.junc_flank_seq.size(); x++)
				{
					flank_seq += fusion_candidate[i].fusion_head[j].first_splice.junc_flank_seq[x];
					flank_seq += ",";
				}
				flank_seq += fusion_candidate[i].fusion_head[j].flank_seq + ",";
				if( fusion_candidate[i].strand == "+")
				{
					for(size_t x = 0; x < fusion_candidate[i].fusion_head[j].second_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].fusion_head[j].second_splice.junc_flank_seq[x];
						flank_seq += ",";
					}	
					for(size_t x = 0; x < fusion_candidate[i].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].junc_flank_seq[x];
						flank_seq += ",";
					}
				}
				else
				{
					for(size_t x = 0; x < fusion_candidate[i].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].junc_flank_seq[x];
						flank_seq += ",";
					}	
					for(size_t x = 0; x < fusion_candidate[i].fusion_head[j].second_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].fusion_head[j].second_splice.junc_flank_seq[x];
						flank_seq += ",";
					}	
				}
				char output_buf[5000];
				sprintf(output_buf, "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\tNM:i:%zu\tXF:Z:%s",
					r_block->read_id.c_str(), fusion_candidate[i].fusion_head[j].first_splice.chrom.c_str(), fusion_candidate[i].fusion_head[j].first_splice.strand.c_str(), start1 + 1, jump_code1_str.c_str(), \
					fusion_candidate[i].fusion_head[j].second_splice.chrom.c_str(), fusion_candidate[i].fusion_head[j].second_splice.strand.c_str(), start2 + 1, jump_code2_str.c_str(), \
					map_sequence.c_str(), map_qual.c_str(), total_mismatch, flank_seq.c_str());
				pthread_mutex_lock (fusion_fs_lock);
				(*fusion_fs) << output_buf << endl;
				pthread_mutex_unlock (fusion_fs_lock);
			}
			for(size_t j = 0; j < fusion_candidate[i].fusion_tail.size(); j++)
			{
				if(debug)
					cout << "print fusion for bwt vector " << i << ", fusion tail " << j << endl;
				int start1 = (fusion_candidate[i].strand == "+" ? fusion_candidate[i].start_pos : fusion_candidate[i].fusion_tail[j].first_splice.start_pos);
				int match1 = fusion_candidate[i].mapped_len + fusion_candidate[i].fusion_tail[j].first_splice.GetMapLen();
				int start2 = fusion_candidate[i].fusion_tail[j].second_splice.start_pos;
				size_t total_mismatch = 0;
				string map_sequence1;
				string map_sequence2;
				string map_sequence;
				string map_qual1;
				string map_qual2;
				string map_qual;
				string chrom_sequence;
				vector<Jump_Code> jump_code1;
				vector<Jump_Code> jump_code2;
				string jump_code1_str;
				string jump_code2_str;
				if(fusion_candidate[i].fusion_tail[j].first_splice.strand == "+" && fusion_candidate[i].fusion_tail[j].second_splice.strand == "+")
				{
					map_sequence1 = sequence.substr(0, match1);
					map_sequence2 = sequence.substr(match1);
					if(!fa)
					{
						map_qual1 = quality_string.substr(0, match1);
						map_qual2 = quality_string.substr(match1);
					}
					combine_jumpcode(fusion_candidate[i].jump_code, fusion_candidate[i].fusion_tail[j].first_splice.jump_code, jump_code1);
				}
				else if(fusion_candidate[i].fusion_tail[j].first_splice.strand == "-" && fusion_candidate[i].fusion_tail[j].second_splice.strand == "-")
				{
					map_sequence1 = revcomp(sequence.substr(0, match1));
					map_sequence2 = revcomp(sequence.substr(match1));
					if(!fa)
					{
						map_qual1 = revqual(quality_string.substr(0, match1));
						map_qual2 = revqual(quality_string.substr(match1));
					}
					combine_jumpcode(fusion_candidate[i].fusion_tail[j].first_splice.jump_code, fusion_candidate[i].jump_code, jump_code1);
				}
				else if(fusion_candidate[i].fusion_tail[j].first_splice.strand == "+" && fusion_candidate[i].fusion_tail[j].second_splice.strand == "-")
				{
					map_sequence1 = sequence.substr(0, match1);
					map_sequence2 = revcomp(sequence.substr(match1));
					if(!fa)
					{
						map_qual1 = quality_string.substr(0, match1);
						map_qual2 = revqual(quality_string.substr(match1));
					}
					combine_jumpcode(fusion_candidate[i].jump_code, fusion_candidate[i].fusion_tail[j].first_splice.jump_code, jump_code1);
				}
				else
				{
					map_sequence1 = revcomp(sequence.substr(0, match1));
					map_sequence2 = sequence.substr(match1);
					if(!fa)
					{
						map_qual1 = revqual(quality_string.substr(0, match1));
						map_qual2 = quality_string.substr(match1);					
					}
					combine_jumpcode(fusion_candidate[i].fusion_tail[j].first_splice.jump_code, fusion_candidate[i].jump_code, jump_code1);
				}
				for(size_t x = 0; x < fusion_candidate[i].fusion_tail[j].second_splice.jump_code.size(); x++)
					jump_code2.push_back(fusion_candidate[i].fusion_tail[j].second_splice.jump_code[x]);
				map_sequence = map_sequence1 + map_sequence2;
				if(!fa)
					map_qual = map_qual1 + map_qual2;
				else
					map_qual = quality_string;
				jump_code1_str = jump_code_to_str(jump_code1);
				jump_code2_str = jump_code_to_str(jump_code2);
				chrom_sequence = get_chrom_sequence((*reference_sequence)[fusion_candidate[i].fusion_tail[j].first_splice.chrom], map_sequence1, start1, jump_code1) + get_chrom_sequence((*reference_sequence)[fusion_candidate[i].fusion_tail[j].second_splice.chrom], map_sequence2, start2, jump_code2);
				for(size_t k = 0; k < map_sequence.length(); k++)
				{
					if(map_sequence[k] != chrom_sequence[k])
						total_mismatch ++ ;
				}	
				string flank_seq;
				if(fusion_candidate[i].strand == "+")
				{
					for(size_t x = 0; x < fusion_candidate[i].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].junc_flank_seq[x];
						flank_seq += ",";
					}
					for(size_t x = 0; x < fusion_candidate[i].fusion_tail[j].first_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].fusion_tail[j].first_splice.junc_flank_seq[x];
						flank_seq += ",";
					}
				}
				else
				{
					for(size_t x = 0; x < fusion_candidate[i].fusion_tail[j].first_splice.junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].fusion_tail[j].first_splice.junc_flank_seq[x];
						flank_seq += ",";
					}
					for(size_t x = 0; x < fusion_candidate[i].junc_flank_seq.size(); x++)
					{
						flank_seq += fusion_candidate[i].junc_flank_seq[x];
						flank_seq += ",";
					}
				}
				flank_seq += fusion_candidate[i].fusion_tail[j].flank_seq;
				flank_seq += ",";
				for(size_t x = 0; x < fusion_candidate[i].fusion_tail[j].second_splice.junc_flank_seq.size(); x++)
				{
					flank_seq += fusion_candidate[i].fusion_tail[j].second_splice.junc_flank_seq[x];
					flank_seq += ",";
				}
				char output_buf[5000];
				sprintf(output_buf, "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\tNM:i:%zu\tXF:Z:%s",
					r_block->read_id.c_str(), fusion_candidate[i].fusion_tail[j].first_splice.chrom.c_str(), fusion_candidate[i].fusion_tail[j].first_splice.strand.c_str(), start1 + 1, jump_code1_str.c_str(), \
					fusion_candidate[i].fusion_tail[j].second_splice.chrom.c_str(), fusion_candidate[i].fusion_tail[j].second_splice.strand.c_str(), start2 + 1, jump_code2_str.c_str(), \
					map_sequence.c_str(), map_qual.c_str(), total_mismatch, flank_seq.c_str());
				pthread_mutex_lock (fusion_fs_lock);
				(*fusion_fs) << output_buf << endl;
				pthread_mutex_unlock (fusion_fs_lock);
			}
		}
		if(debug)
			cout << "report fusion complete" << endl;
	}

	void mapsplice_fusion()
	{
		if(debug)
			cout << "pe1 aligned " << pe1_aligned << " pe2 aligned " << pe2_aligned << endl;
		if(pair_end)
		{
			if(!pe1_aligned || !pe2_aligned)
			{
				sort(fusion_candidate.begin(), fusion_candidate.end(), sort_sam_by_pos);
				search_fusion_double_anchor(!pe1_aligned, !pe2_aligned);
				bool remap_pe1 = !pe1_aligned && fusiondb;
				bool remap_pe2 = !pe2_aligned && fusiondb;
				bool search_pe1 = !pe1_aligned && pe2_aligned;
				bool search_pe2 = !pe2_aligned && pe1_aligned;
				search_fusion_single_anchor(remap_pe1, remap_pe2, search_pe1, search_pe2);
				report_fusion();
			}
		}
		else
		{
			if(!pe1_aligned)
			{
				sort(fusion_candidate.begin(), fusion_candidate.end(), sort_sam_by_pos);
				search_fusion_double_anchor(!pe1_aligned, !pe2_aligned);
				if(fusiondb)
					search_fusion_single_anchor(true, false, false, false);	
				report_fusion();
			}
		}
	}

	void update_stats()
	{
		if(pair_end)
			mapping_stats->read_processed += 2;
		else
			mapping_stats->read_processed++;
		if(pe1_aligned)
			mapping_stats->read_aligned ++;
		if(pe2_aligned)
			mapping_stats->read_aligned ++;
		if(pe1_spliced)
			mapping_stats->read_spliced ++;
		if(pe2_spliced)
			mapping_stats->read_spliced ++;
		if(pe1_fusion)
			mapping_stats->read_fusion ++;
		if(pe2_fusion)
			mapping_stats->read_fusion ++;			
	}

	void align()
	{
		nofw = false;
		norc = false;
		pe1_segment_number = r_block_pe1->get_seg_num();
		pe2_segment_number = r_block_pe2->get_seg_num();
		pe1_read_length = r_block_pe1->get_seg_len(HEAD_SEG_NO, pe1_segment_number);
		pe2_read_length = r_block_pe2->get_seg_len(HEAD_SEG_NO, pe2_segment_number);
		anchor_number = 0;
		remap_number = 0;
		pe1_spliced = false;
		pe1_unspliced = false;
		pe2_spliced = false;
		pe2_unspliced = false;
		pe1_aligned = false;
		pe2_aligned = false;
		pe1_fusion = false;
		pe2_fusion = false;
		read_paired = false;
		num_pe1_alignment = 0;
		num_pe2_alignment = 0;
		if(bwt_vector->size() > 0)
		{
			sort((*bwt_vector).begin(), (*bwt_vector).end(), bwtmap_compare);
			combine_contig();
			mapsplice_search();
			if(debug)
				cout << "mapsplice search done" << endl;
			mapsplice_report();
			if(debug)
				cout << "mapsplice report done" << endl;
		}
		if(optimize_for_repeats && (!(*pe1_no_repeats) || !(*pe2_no_repeats)) )
		{

			if(debug)
				cout << "optimizing repeats" << endl;
			set_bwt_vector_repeats();
			optimize_repeats();
			set_bwt_vector_normal();
			if(debug)
				cout << "optimizing repeats done" << endl;
		}
		if(pair_end && do_local_align && ((pe1_aligned && !pe2_aligned) || (!pe1_aligned && pe2_aligned)))
		{
			full_read_local_align();	
		}
		if(fusion)
		{
			if(debug)
				cout << "fixing fusion" << endl;
			mapsplice_fusion();
			if(debug)
			{
				output_sam();
				cout << "fixing fusion done" << endl;
			}
		}
		else
		{
			//if(do_pairing && pe1_aligned && pe2_aligned)
			pair_alignment();
			int num_output = output_sam();
			train_adaptive_mate_dist(num_output);
		}
		if(output_unmapped && ((pair_end && (!pe1_aligned || !pe2_aligned)) || (!pair_end && !pe1_aligned)))
			output_unmapped_read();
		{
			vector<Paired_Sam_Info>().swap(paired_sam_vec);
		}
		pe1_sam.clear();
		pe2_sam.clear();
		sam_vec.clear();
		fusion_candidate.clear();
		update_stats();
		if(debug)
			cout << "end 1 aligned " << pe1_aligned << " end 2 aligned " << pe2_aligned << endl;
	}
};

#endif

