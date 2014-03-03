#include "AlignmentHandler.h"

bool comp_sam_rec_maplen(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.mappedlen > rhs.mappedlen;
}

bool comp_min_anchor(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.min_anchor > rhs.min_anchor;
}

bool comp_mate_dist(size_t lhs, size_t rhs)
{
	if (lhs >= rhs && lhs - rhs <= mate_dist_sd)
		return false;
	if (rhs > lhs  && rhs - lhs <= mate_dist_sd)
		return false;
	return lhs < rhs;
}

bool comp_filter_score_range(double lhs, double rhs)
{
	if (lhs >= rhs && lhs - rhs <= 0.00001)
		return false;
	if (rhs > lhs  && rhs - lhs <= 0.00001)
		return false;
	return lhs < rhs;
}

bool better_than(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (comp_mate_dist(lhs.mate_dist, rhs.mate_dist) && comp_mate_dist(rhs.mate_dist, lhs.mate_dist))
	{
		if (lhs.total_mismatch == rhs.total_mismatch)
		{
			if (lhs.mappedlen == rhs.mappedlen)
			{
				if (lhs.intron_size == rhs.intron_size)
					return false;
				else 
					return lhs.intron_size < rhs.intron_size;
			}
			else
				return lhs.mappedlen > rhs.mappedlen;
		}			
		else
		{
			return lhs.total_mismatch < rhs.total_mismatch;
		}
	}
	else
		return lhs.mate_dist < rhs.mate_dist;
}
bool comp_dist(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.mappedlen == rhs.mappedlen)
	{
		if (lhs.total_mismatch == rhs.total_mismatch)
		{
			if (!comp_mate_dist(lhs.mate_dist, rhs.mate_dist) && !comp_mate_dist(rhs.mate_dist, lhs.mate_dist))
			{
				//if (lhs.total_ave_mismatch == rhs.total_ave_mismatch)
				//{
					if (!comp_intron_dist(lhs.intron_size, rhs.intron_size) && !comp_intron_dist(rhs.intron_size, lhs.intron_size))
					{
						if (lhs.total_hits == rhs.total_hits)
						{
							if (lhs.total_filter_score == rhs.total_filter_score)
							{
								if (lhs.contiglen == rhs.contiglen)
									return lhs.min_anchor > rhs.min_anchor;
								else;
								return lhs.contiglen > rhs.contiglen;
							}
							else
								return lhs.total_filter_score > rhs.total_filter_score;
						}
						else
							return lhs.total_hits > rhs.total_hits;
					}
					else
						return lhs.intron_size < rhs.intron_size;
				//}
				//else
				//	return lhs.total_ave_mismatch < rhs.total_ave_mismatch;

			}
			else
				return lhs.mate_dist < rhs.mate_dist;
		}			
		else
		{
			return lhs.total_mismatch < rhs.total_mismatch;
		}
	}
	else
		return lhs.mappedlen > rhs.mappedlen;
}

bool comp_dist_0(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.mappedlen == rhs.mappedlen)
	{
		if (lhs.total_mismatch == rhs.total_mismatch)
		{
			if (lhs.total_canon_rate == rhs.total_canon_rate)
			{
				if (lhs.mate_dist == rhs.mate_dist)
				{
					//if (lhs.total_ave_mismatch == rhs.total_ave_mismatch)
					//{
					if (lhs.contiglen == rhs.contiglen)
					{
						if (lhs.intron_size == rhs.intron_size)
						{
							if (lhs.total_hits == rhs.total_hits)
							{
								if (lhs.total_filter_score == rhs.total_filter_score)
								{

									if (lhs.min_anchor == rhs.min_anchor)
										return lhs.total_confid > rhs.total_confid;
									else
										return lhs.min_anchor > rhs.min_anchor;

								}
								else
									return lhs.total_filter_score > rhs.total_filter_score;
							}
							else
								return lhs.total_hits > rhs.total_hits;

							//if (lhs.total_filter_score == rhs.total_filter_score)
							//{
							//	if (lhs.contiglen == rhs.contiglen)
							//		return lhs.min_anchor > rhs.min_anchor;
							//	else;
							//	return lhs.contiglen > rhs.contiglen;
							//}
							//else
							//	return lhs.total_filter_score > rhs.total_filter_score;
						}
						else
							return lhs.intron_size < rhs.intron_size;
						//}
						//else
						//	return lhs.total_ave_mismatch < rhs.total_ave_mismatch;

					}
					else
						return lhs.contiglen > rhs.contiglen;
				}
				else
					return lhs.mate_dist < rhs.mate_dist;
			}
			else
				return lhs.total_canon_rate > rhs.total_canon_rate;
		}			
		else
		{
			return lhs.total_mismatch < rhs.total_mismatch;
		}
	}
	else
		return lhs.mappedlen > rhs.mappedlen;
}

bool comp_dist_1(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.intron_size == rhs.intron_size)
	{
		if (lhs.total_hits == rhs.total_hits)
		{
			if (lhs.total_filter_score == rhs.total_filter_score)
			{
				//if (lhs.contiglen == rhs.contiglen)
				//{
					if (lhs.min_anchor == rhs.min_anchor)
						return lhs.total_confid > rhs.total_confid;
					else
						return lhs.min_anchor > rhs.min_anchor;
				//}
				//else
				//	return lhs.contiglen > rhs.contiglen;
			}
			else
				return lhs.total_filter_score > rhs.total_filter_score;
		}
		else
			return lhs.total_hits > rhs.total_hits;
		//if (lhs.total_filter_score == rhs.total_filter_score)
		//{
		//	if (lhs.contiglen == rhs.contiglen)
		//		return lhs.min_anchor > rhs.min_anchor;
		//	else;
		//	return lhs.contiglen > rhs.contiglen;
		//}
		//else
		//	return lhs.total_filter_score > rhs.total_filter_score;
	}
	else
		return lhs.intron_size < rhs.intron_size;
}

bool comp_dist_2(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.total_hits == rhs.total_hits)
	{
		if (lhs.total_filter_score == rhs.total_filter_score)
		{
			//if (lhs.contiglen == rhs.contiglen)
			//{
				if (lhs.min_anchor == rhs.min_anchor)
					return lhs.total_confid > rhs.total_confid;
				else
					return lhs.min_anchor > rhs.min_anchor;
			//}
			//else
			//	return lhs.contiglen > rhs.contiglen;
		}
		else
			return lhs.total_filter_score > rhs.total_filter_score;
	}
	else
		return lhs.total_hits > rhs.total_hits;
	/*if (lhs.total_filter_score == rhs.total_filter_score)
	{
		if (lhs.contiglen == rhs.contiglen)
			return lhs.min_anchor > rhs.min_anchor;
		else;
		return lhs.contiglen > rhs.contiglen;
	}
	else
		return lhs.total_filter_score > rhs.total_filter_score;*/
}

bool comp_intron_dist(size_t lhs, size_t rhs)
{
	if (lhs >= rhs && lhs - rhs <= intron_dist_sd)
		return false;
	if (rhs > lhs  && rhs - lhs <= intron_dist_sd)
		return false;
	return lhs < rhs;
}

bool comp_mis(const SamRec* lhs, const SamRec* rhs)
{
	return lhs->mis_match < rhs->mis_match;
}

bool comp_dup(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->chrom_name == rhs->chrom_name)
	{
		if (lhs->start == rhs->start)
		{
			if (lhs->splice_way == rhs->splice_way)
			{
				if (lhs->strand_t == rhs->strand_t)
				{
					if (lhs->chrom_name2 == rhs->chrom_name2)
					{
						if (lhs->start2 == rhs->start2)
						{
							if (lhs->splice_way2 == rhs->splice_way2)
								return lhs->strand_t2  < rhs->strand_t2;
							else
								return lhs->splice_way2 < rhs->splice_way2;
						}
						else
							return lhs->start2 < rhs->start2;
					}
					else
						return lhs->chrom_name2 < rhs->chrom_name2;
				}
				else
					return lhs->strand_t  < rhs->strand_t;
			}
			else
				return lhs->splice_way < rhs->splice_way;
		}
		else
			return lhs->start < rhs->start;
	}
	else
		return lhs->chrom_name < rhs->chrom_name;
}

bool
comp_filterscore(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->isexonic && !rhs->isexonic)
		return true;
	else if (!lhs->isexonic && rhs->isexonic)
		return false;

	//if (lhs->mappedlen == rhs->mappedlen)
	//{

	//return lhs->ave_intron_len < rhs->ave_intron_len;

	;
	if (lhs->mappedlen == rhs->mappedlen)
	{
		if (lhs->filter_score == rhs->filter_score)
		{
			if (lhs->ave_intron_len == rhs->ave_intron_len)
			{
				if (lhs->ave_junc_mis == rhs->ave_junc_mis)
					return lhs->m_contig_len > rhs->m_contig_len;
				else
					return lhs->ave_junc_mis < rhs->ave_junc_mis;

			}
			else
				return lhs->ave_intron_len < rhs->ave_intron_len;
		}
		else
			return lhs->filter_score > rhs->filter_score;
	}
	else
		return lhs->mappedlen > rhs->mappedlen;

	;
	//}
	//else
	//	return lhs->mappedlen > rhs->mappedlen;
}

bool
comp_filterscore1(const SamRec* lhs, const SamRec* rhs)
{
	//if (lhs->isexonic && !rhs->isexonic)
	//	return true;
	//else if (!lhs->isexonic && rhs->isexonic)
	//	return false;

	//if (lhs->mappedlen == rhs->mappedlen)
	//{

	//return lhs->ave_intron_len < rhs->ave_intron_len;

	;

	if (lhs->ave_intron_len == rhs->ave_intron_len)
	{
		if (lhs->ave_junc_mis == rhs->ave_junc_mis)
			return lhs->m_contig_len > rhs->m_contig_len;
		else
			return lhs->ave_junc_mis < rhs->ave_junc_mis;

	}
	else
		return lhs->ave_intron_len < rhs->ave_intron_len;

	;
	//}
	//else
	//	return lhs->mappedlen > rhs->mappedlen;
}

bool
comp_filterscore2(const SamRec* lhs, const SamRec* rhs)
{
	//if (lhs->isexonic && !rhs->isexonic)
	//	return true;
	//else if (!lhs->isexonic && rhs->isexonic)
	//	return false;

	//if (lhs->mappedlen == rhs->mappedlen)
	//{

	//return lhs->ave_intron_len < rhs->ave_intron_len;

	;

	if (lhs->ave_junc_mis == rhs->ave_junc_mis)
		return lhs->m_contig_len > rhs->m_contig_len;
	else
		return lhs->ave_junc_mis < rhs->ave_junc_mis;

	;
	//}
	//else
	//	return lhs->mappedlen > rhs->mappedlen;
}

bool comp_fusion_filterscore(const SamRec* lhs, const SamRec* rhs)
{
	if (lhs->isexonic && !lhs->is_fusion && !rhs->isexonic)
		return true;
	else if (!lhs->isexonic && rhs->isexonic && !rhs->is_fusion)
		return false;

	if (lhs->mappedlen == rhs->mappedlen)
	{
		if (lhs->filter_score == rhs->filter_score)
		{
			if (lhs->ave_junc_mis == rhs->ave_junc_mis)
				return lhs->m_contig_len > rhs->m_contig_len;
			else
				return lhs->ave_junc_mis < rhs->ave_junc_mis;
		}
		else
			return lhs->filter_score > rhs->filter_score;
	}
	else
		return lhs->mappedlen > rhs->mappedlen;
	//}
	//else
	//	return lhs->mappedlen > rhs->mappedlen;
}

bool check_overlap(vector<pair<size_t, int> >& spliceway_vec_1, vector<pair<size_t, int> >& spliceway_vec_2)
{
	vector<pair<size_t, int> >::iterator sp_vec_iter1, sp_vec_iter2;

	for (sp_vec_iter1 = spliceway_vec_1.begin(); sp_vec_iter1 != spliceway_vec_1.end()/* - 1*/; ++sp_vec_iter1)
	{
		//if (sp_vec_iter1->second > 0/* && (sp_vec_iter1 + 1)->second > 0*/)
		{
			for (sp_vec_iter2 = spliceway_vec_2.begin(); sp_vec_iter2 != spliceway_vec_2.end() - 1; ++sp_vec_iter2)
			{
				if (sp_vec_iter2->second > 0 && (sp_vec_iter2 + 1)->second > 0)
				{
					if (sp_vec_iter1->first > sp_vec_iter2->first + sp_vec_iter2->second && 
						sp_vec_iter1->first < (sp_vec_iter2 + 1)->first)
						return true;
					if (sp_vec_iter1->second > 0 &&
						sp_vec_iter1->first + sp_vec_iter1->second > sp_vec_iter2->first + sp_vec_iter2->second && 
						sp_vec_iter1->first + sp_vec_iter1->second - 1 < (sp_vec_iter2 + 1)->first)
						return true;
				}
			}
		}
	}

	return false;
}

struct EstPairingStruct {
	pair<vector<SamRec*>, vector<SamRec*> >* sam_rec_pe_ptr;

	vector<PairedSamRec>* fusion_paired_reads_ptr;

	int do_filter;

	EstPairingStruct(pair<vector<SamRec*>, vector<SamRec*> >* _sam_rec_pe_ptr, vector<PairedSamRec>* _fusion_paired_reads_ptr, int _do_filter)
		: sam_rec_pe_ptr(_sam_rec_pe_ptr), fusion_paired_reads_ptr(_fusion_paired_reads_ptr), do_filter(_do_filter) {}

};

struct FilterPairingStruct {

	vector<PairedSamRec>* pairedsamrec_vec;

	pair<vector<SamRec*>, vector<SamRec*> >* sam_rec_pe_ptr;

	int do_filter;

	FilterPairingStruct(vector<PairedSamRec>* _pairedsamrec_vec, pair<vector<SamRec*>, vector<SamRec*> >* _sam_rec_pe_ptr, int _do_filter)
		: pairedsamrec_vec(_pairedsamrec_vec), sam_rec_pe_ptr(_sam_rec_pe_ptr), do_filter(_do_filter) {}
};

struct FindFusionJuncRegionStruct {
	FusionJuncRegion* cur_fusion_region;

	PairedSamRec* cur_paired_sam_rec;

	vector<FusionJuncRegion>* sorted_fusion_regions;

	FindFusionJuncRegionStruct(FusionJuncRegion* _cur_fusion_region, PairedSamRec* _cur_paired_sam_rec, vector<FusionJuncRegion>* _sorted_fusion_regions)
		: cur_fusion_region(_cur_fusion_region), cur_paired_sam_rec(_cur_paired_sam_rec), sorted_fusion_regions(_sorted_fusion_regions)
	{
	}
};

string AlignmentHandler::m_chrom_dir;

size_t AlignmentHandler::m_cur_read_alignments_index = 0;

size_t AlignmentHandler::m_cur_line_index = 0;

vector<string> AlignmentHandler::m_stored_lines(buf_size);

vector<SamRec> AlignmentHandler::m_stored_sam_recs(buf_size);

vector<pair<vector<SamRec>::iterator, vector<SamRec>::iterator> > AlignmentHandler::m_same_reads_alignments;

vector<pair<vector<SamRec*>, vector<SamRec*> > > AlignmentHandler::m_same_reads_alignments_paired;

vector<bool> AlignmentHandler::m_same_reads_alignments_paired_is_unspliced;

int AlignmentHandler::m_min_ins_static;

int AlignmentHandler::m_do_filter;

size_t AlignmentHandler::m_max_pair_dist;

JunctionHandler AlignmentHandler::m_junction_handler;

JunctionHandler AlignmentHandler::m_junction_handler_filtered;

JunctionHandler AlignmentHandler::m_junction_handler_intermediate;

JunctionHandler AlignmentHandler::m_junction_handler_annotated;

JunctionHandler* AlignmentHandler::m_junction_handler_prev_ptr;

double AlignmentHandler::m_entrpy_weight;

double AlignmentHandler::m_pqlen_weight;

double AlignmentHandler::m_ave_mis_weight;

bool AlignmentHandler::m_add_S;

size_t AlignmentHandler::m_max_hits;

hash_map<string, vector<pair<int, int> > > AlignmentHandler::m_normal_paired_locs;

hash_map<string, vector<pair<int, int> > > AlignmentHandler::m_original_junc_locs;

size_t AlignmentHandler::m_spliced;
size_t AlignmentHandler::m_unspliced;
size_t AlignmentHandler::m_insertion;
size_t AlignmentHandler::m_deletion;
size_t AlignmentHandler::m_unique;
size_t AlignmentHandler::m_multiple;
size_t AlignmentHandler::m_canoical;
size_t AlignmentHandler::m_semi_canonical;
size_t AlignmentHandler::m_non_canonical;
size_t AlignmentHandler::m_paired;
size_t AlignmentHandler::m_single;
size_t AlignmentHandler::m_fusion_paired;
size_t AlignmentHandler::m_unmapped;
size_t AlignmentHandler::m_clipped;
size_t AlignmentHandler::m_fusion;

size_t AlignmentHandler::m_paired_multiple;
size_t AlignmentHandler::m_paired_unique;
size_t AlignmentHandler::m_single_multiple;
size_t AlignmentHandler::m_single_unique;
size_t AlignmentHandler::m_fusion_paired_multiple;
size_t AlignmentHandler::m_fusion_paired_unique;

size_t AlignmentHandler::m_both_unspliced_multiple;
size_t AlignmentHandler::m_both_unspliced_unique;
size_t AlignmentHandler::m_both_unspliced_paired;
size_t AlignmentHandler::m_both_unspliced_single;
size_t AlignmentHandler::m_both_unspliced_fusion_paired;

size_t AlignmentHandler::m_both_unspliced_paired_multiple;
size_t AlignmentHandler::m_both_unspliced_paired_unique;
size_t AlignmentHandler::m_both_unspliced_single_multiple;
size_t AlignmentHandler::m_both_unspliced_single_unique;
size_t AlignmentHandler::m_both_unspliced_fusion_paired_multiple;
size_t AlignmentHandler::m_both_unspliced_fusion_paired_unique;



size_t AlignmentHandler::m_spliced_paired;
size_t AlignmentHandler::m_spliced_single;
size_t AlignmentHandler::m_spliced_fusion_paired;
	
size_t AlignmentHandler::m_unspliced_paired;
size_t AlignmentHandler::m_unspliced_single;
size_t AlignmentHandler::m_unspliced_fusion_paired;

size_t AlignmentHandler::m_insertion_paired;
size_t AlignmentHandler::m_insertion_single;
size_t AlignmentHandler::m_insertion_fusion_paired;

size_t AlignmentHandler::m_deletion_paired;
size_t AlignmentHandler::m_deletion_single;
size_t AlignmentHandler::m_deletion_fusion_paired;

size_t AlignmentHandler::m_canoical_paired;
size_t AlignmentHandler::m_canoical_single;
size_t AlignmentHandler::m_canoical_fusion_paired;

size_t AlignmentHandler::m_semi_canonical_paired;
size_t AlignmentHandler::m_semi_canonical_single;
size_t AlignmentHandler::m_semi_canonical_fusion_paired;

size_t AlignmentHandler::m_non_canonical_paired;
size_t AlignmentHandler::m_non_canonical_single;
size_t AlignmentHandler::m_non_canonical_fusion_paired;

size_t AlignmentHandler::m_clipped_paired;
size_t AlignmentHandler::m_clipped_single;
size_t AlignmentHandler::m_clipped_fusion_paired;

size_t AlignmentHandler::m_fusion_canonical;
size_t AlignmentHandler::m_fusion_semi_canonical;
size_t AlignmentHandler::m_fusion_non_canonical;

size_t AlignmentHandler::m_fusion_multiple;
size_t AlignmentHandler::m_fusion_unique;

string AlignmentHandler::m_prev_tag;

bool AlignmentHandler::m_is_paired;


hash_map<string, vector<bool> > AlignmentHandler::m_mapped_pos;

hash_map<string, vector<pair<size_t, size_t> > > AlignmentHandler::m_expressed_regions;

hash_map<string, DisjointSet> AlignmentHandler::m_disjointset_map;

hash_map<string, hash_map<size_t, UnionExpressedRegions > > AlignmentHandler::m_unioned_expressed_regions_map;

AlignmentHandler::AlignmentHandler(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, string filtered_alignment_file,
	size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor, char* chr_sz_file,
	size_t min_fusion_coverage, double entrpy_weight, double pqlen_weight, double ave_mis_weight) : /* m_sam_file_handler(alignment_file.c_str(), min_ins), */
	m_junction_file(junction_file), m_sam_files(alignment_files), 
	 m_filtered_alignment_file(filtered_alignment_file),
	/*m_ofs_filtered_alignment(m_filtered_alignment_file.c_str()),*/ m_max_read_width(max_read_width),  
	m_min_ins(min_ins), m_max_del(max_del),
	m_min_anchor(min_anchor), m_min_mismatch(min_mismatch), 
	m_min_junc_anchor(min_junc_anchor), m_chrom_size_file(chr_sz_file), m_min_fusion_coverage(min_fusion_coverage)
	//m_cur_read_alignments_index(0), m_cur_line_index(0)//, m_stored_lines(buf_size), m_stored_sam_recs(buf_size)
{
	string to_mapper = m_filtered_alignment_file; to_mapper.append(".tomapper");

	m_do_filter = do_filter;

	m_min_ins_static = m_min_ins;

	m_max_pair_dist = max_pair_dist;

	m_entrpy_weight = entrpy_weight;
	
	m_pqlen_weight = pqlen_weight;
	
	m_ave_mis_weight = ave_mis_weight;

	m_chrom_dir = chrom_dir;

	m_add_S = add_S;

	m_max_hits = max_hits;

	m_is_paired = is_paired; 

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0;
	
	m_multiple = 0;
	
	m_canoical = 0;
	
	m_semi_canonical = 0;
	
	m_non_canonical = 0;
	
	m_paired = 0;
	
	m_single = 0;

	m_fusion_paired = 0;

	m_paired_multiple = 0;

	m_paired_unique = 0;

	m_single_multiple = 0;

	m_single_unique = 0;

	m_fusion_paired_multiple = 0;

	m_fusion_paired_unique = 0;

	m_unmapped = 0;
	
	m_fusion = 0;

	m_clipped = 0;

	m_both_unspliced_multiple = 0;
	
	m_both_unspliced_unique = 0;
	
	m_both_unspliced_paired = 0;
	
	m_both_unspliced_single = 0;

	m_both_unspliced_fusion_paired = 0;

	m_both_unspliced_paired_multiple = 0;

	m_both_unspliced_paired_unique = 0;

	m_both_unspliced_single_multiple = 0;

	m_both_unspliced_single_unique = 0;

	m_both_unspliced_fusion_paired_multiple = 0;

	m_both_unspliced_fusion_paired_unique = 0;

	/////////////////

	m_spliced_paired = 0;

	m_spliced_single = 0;

	m_spliced_fusion_paired = 0;

	m_unspliced_paired = 0;

	m_unspliced_single = 0;

	m_unspliced_fusion_paired = 0;

	m_insertion_paired = 0;

	m_insertion_single = 0;

	m_insertion_fusion_paired = 0;

	m_deletion_paired = 0;

	m_deletion_single = 0;

	m_deletion_fusion_paired = 0;

	m_canoical_paired = 0;

	m_canoical_single = 0;

	m_canoical_fusion_paired = 0;

	m_semi_canonical_paired = 0;

	m_semi_canonical_single = 0;

	m_semi_canonical_fusion_paired = 0;

	m_non_canonical_paired = 0;

	m_non_canonical_single = 0;

	m_non_canonical_fusion_paired = 0;

	m_clipped_paired = 0;

	m_clipped_single = 0;

	m_clipped_fusion_paired = 0;

	m_fusion_canonical = 0;

	m_fusion_semi_canonical = 0;

	m_fusion_non_canonical = 0;

	m_fusion_multiple = 0;

	m_fusion_unique = 0;
	
	//m_ofs_to_mapper.open(to_mapper.c_str());
}

bool AlignmentHandler::Init(vector<string> alignment_files, string junction_file, bool is_paired, bool add_S, size_t max_pair_dist, size_t max_hits, 
	string filtered_alignment_file, size_t max_read_width, string chrom_dir, int do_filter, int min_ins, int max_del, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor,
	double entrpy_weight, double pqlen_weight, double ave_mis_weight) 
{
	Clear();

	m_filtered_alignment_file = filtered_alignment_file;

	string to_mapper = m_filtered_alignment_file; to_mapper.append(".tomapper");

	m_entrpy_weight = entrpy_weight;

	m_pqlen_weight = pqlen_weight;

	m_ave_mis_weight = ave_mis_weight;

	m_is_paired = is_paired;

	m_add_S = add_S;

	m_max_pair_dist = max_pair_dist;

	m_max_hits = max_hits;

	m_sam_files = alignment_files;

	m_max_read_width = max_read_width;

	m_chrom_dir = chrom_dir;

	m_do_filter = do_filter;

	m_min_ins = min_ins;

	m_max_del = max_del;

	m_min_anchor = min_anchor;

	m_min_junc_anchor = min_junc_anchor;

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0;
	
	m_multiple = 0;
	
	m_canoical = 0;

	m_semi_canonical = 0;
	
	m_non_canonical = 0;
	
	m_paired = 0;
	
	m_single = 0;

	m_fusion_paired = 0;

	m_paired_multiple = 0;

	m_paired_unique = 0;

	m_single_multiple = 0;

	m_single_unique = 0;

	m_fusion_paired_multiple = 0;

	m_fusion_paired_unique = 0;

	m_unmapped = 0;

	m_fusion = 0;

	m_clipped = 0;

	m_both_unspliced_multiple = 0;
	
	m_both_unspliced_unique = 0;
	
	m_both_unspliced_paired = 0;
	
	m_both_unspliced_single = 0;

	m_both_unspliced_fusion_paired = 0;

	m_both_unspliced_paired_multiple = 0;

	m_both_unspliced_paired_unique = 0;

	m_both_unspliced_single_multiple = 0;

	m_both_unspliced_single_unique = 0;

	m_both_unspliced_fusion_paired_multiple = 0;

	m_both_unspliced_fusion_paired_unique = 0;

	///////////////////////////

	m_spliced_paired = 0;

	m_spliced_single = 0;

	m_spliced_fusion_paired = 0;

	m_unspliced_paired = 0;

	m_unspliced_single = 0;

	m_unspliced_fusion_paired = 0;

	m_insertion_paired = 0;

	m_insertion_single = 0;

	m_insertion_fusion_paired = 0;

	m_deletion_paired = 0;

	m_deletion_single = 0;

	m_deletion_fusion_paired = 0;

	m_canoical_paired = 0;

	m_canoical_single = 0;

	m_canoical_fusion_paired = 0;

	m_semi_canonical_paired = 0;

	m_semi_canonical_single = 0;

	m_semi_canonical_fusion_paired = 0;

	m_non_canonical_paired = 0;

	m_non_canonical_single = 0;

	m_non_canonical_fusion_paired = 0;

	m_clipped_paired = 0;

	m_clipped_single = 0;

	m_clipped_fusion_paired = 0;

	m_fusion_canonical = 0;

	m_fusion_semi_canonical = 0;

	m_fusion_non_canonical = 0;

	m_fusion_multiple = 0;

	m_fusion_unique = 0;

	return true;
}

bool AlignmentHandler::Clear()
{
	m_is_paired = false;

	m_add_S = false;

	m_max_pair_dist = 0;

	m_max_pair_dist = 0;

	m_entrpy_weight = 0;

	m_pqlen_weight = 0;

	m_ave_mis_weight = 0;

	m_sam_files.clear();

	m_junction_handler.Clear();

	m_junc_hash_ptr = 0;

	m_sam_rec_pe.first.clear();

	m_sam_rec_pe.second.clear();

	m_sam_rec_pe_ptr.first.clear();

	m_sam_rec_pe_ptr.second.clear();

	m_filtered_alignment_file.clear();

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_bothunspliced.close();

	m_ofs_onespliced.close();

	m_chrom_dir.clear();

	m_do_filter = 0;

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0; 
	
	m_multiple = 0;
	
	m_canoical = 0;

	m_semi_canonical = 0;
	
	m_non_canonical = 0;

	return true;
}

bool AlignmentHandler::ClearStats()
{
	m_unmapped = 0;

	m_fusion = 0;

	m_clipped = 0;

	m_spliced = 0;
	
	m_unspliced = 0;
	
	m_insertion = 0;
	
	m_deletion = 0;
	
	m_unique = 0; 
	
	m_multiple = 0;
	
	m_canoical = 0;

	m_semi_canonical = 0;
	
	m_non_canonical = 0;

	m_paired = 0;

	m_single = 0;

	///////////////////

	m_fusion_paired = 0;

	m_paired_multiple = 0;

	m_paired_unique = 0;

	m_single_multiple = 0;

	m_single_unique = 0;

	m_fusion_paired_multiple = 0;

	m_fusion_paired_unique = 0;

	m_spliced_paired = 0;

	m_spliced_single = 0;

	m_spliced_fusion_paired = 0;

	m_unspliced_paired = 0;

	m_unspliced_single = 0;

	m_unspliced_fusion_paired = 0;

	m_insertion_paired = 0;

	m_insertion_single = 0;

	m_insertion_fusion_paired = 0;

	m_deletion_paired = 0;

	m_deletion_single = 0;

	m_deletion_fusion_paired = 0;

	m_canoical_paired = 0;

	m_canoical_single = 0;

	m_canoical_fusion_paired = 0;

	m_semi_canonical_paired = 0;

	m_semi_canonical_single = 0;

	m_semi_canonical_fusion_paired = 0;

	m_non_canonical_paired = 0;

	m_non_canonical_single = 0;

	m_non_canonical_fusion_paired = 0;

	m_clipped_paired = 0;

	m_clipped_single = 0;

	m_clipped_fusion_paired = 0;

	m_fusion_canonical = 0;

	m_fusion_semi_canonical = 0;

	m_fusion_non_canonical = 0;

	m_fusion_multiple = 0;

	m_fusion_unique = 0;

	return true;
}

void 
AlignmentHandler::FilterByMinAnchor(vector<SamRec>& read_sam)
{
	sort(read_sam.begin(), read_sam.end(), comp_min_anchor);

	vector<SamRec >::iterator sam_rec_iter1;

	for (sam_rec_iter1 = read_sam.begin(); sam_rec_iter1 != read_sam.end(); ++sam_rec_iter1)
	{
		if (sam_rec_iter1->min_anchor < m_min_anchor)
		{
			read_sam.resize(sam_rec_iter1 - read_sam.begin());

			break;
		}
	}
}

void 
AlignmentHandler::FilterByMapLen(vector<SamRec>& read_sam)
{
	sort(read_sam.begin(), read_sam.end(), comp_sam_rec_maplen);

	vector<SamRec >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin() + 1; sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		if (sam_rec_iter->mappedlen != (sam_rec_iter - 1)->mappedlen)
		{
			read_sam.resize(sam_rec_iter - read_sam.begin());

			break;
		}
	}
}

void
AlignmentHandler::MarkCanonNoncanonByReads(vector<SamRec>& read_sam, JunctionHandler* junc_hash_ptr)
{
	if (read_sam.empty())
		return;

	vector<SamRec >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin(); sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		if ((sam_rec_iter)->spliceway_vec.size() <= 1)
			continue;

		vector<pair<size_t, int> >::iterator oft_mpl_iter;

		double sumscore = 0;

		double sum_intron_len = 0;

		double junc_anchor_len = 0;

		double sum_ave_mis = 0;

		double pair_rate = 0;

		for (oft_mpl_iter = (sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
		{
			if ((oft_mpl_iter + 1)->second < 0)
			{
				continue;
			}

			size_t comb_offset;

			if (oft_mpl_iter->second < 0)
			{
				comb_offset = ((oft_mpl_iter->first - 1) << THIRTY_TWO) + oft_mpl_iter->first - 1;
			}
			else
				comb_offset = ((oft_mpl_iter->first + oft_mpl_iter->second - 1) << THIRTY_TWO) + (oft_mpl_iter + 1)->first;

			CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = junc_hash_ptr->m_junc_hash.find(sam_rec_iter->chrom_name);

			if (chrom_junc_hash_iter == junc_hash_ptr->m_junc_hash.end())
			{
				sam_rec_iter->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				continue;
			}

			JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

			JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(comb_offset);

			if (junc_hash_comb_iter == junc_hash_comb.end())
			{
				sam_rec_iter->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				continue;
			}

			sam_rec_iter->filter_type= (junc_hash_comb_iter->second.m_filtered_type);

			sumscore += (junc_hash_comb_iter->second.m_entropy * m_entrpy_weight) + (junc_hash_comb_iter->second.m_lpq * m_pqlen_weight) 
				+ (junc_hash_comb_iter->second.m_ave_mismatch * m_ave_mis_weight);

			sum_intron_len += junc_hash_comb_iter->second.m_intronlen;

			junc_anchor_len += junc_hash_comb_iter->second.m_max_prefix_len + junc_hash_comb_iter->second.m_max_suffix_len;

			sum_ave_mis += junc_hash_comb_iter->second.m_ave_mismatch;

			pair_rate += (double)junc_hash_comb_iter->second.m_paired_count;
			
			if (junc_hash_comb_iter->second.m_flankcase >= 5)
			{
				(sam_rec_iter)->canon_count++;
				(sam_rec_iter)->iscanonical = true;
			}
			else if (junc_hash_comb_iter->second.m_flankcase >= 1)
			{
				(sam_rec_iter)->canon_count++;

				(sam_rec_iter)->noncanon_count++;

				(sam_rec_iter)->issemicanonical = true;
			}
			else
			{
				(sam_rec_iter)->noncanon_count++;

				(sam_rec_iter)->isnoncanoical = true;
			}
		}

		sumscore = sumscore / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		sum_intron_len =  sum_intron_len / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		junc_anchor_len = junc_anchor_len / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		sum_ave_mis = sum_ave_mis / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		pair_rate = pair_rate / double ((sam_rec_iter)->spliceway_vec.size() - 1);

		(sam_rec_iter)->filter_score = sumscore;

		(sam_rec_iter)->junc_anchor_len = junc_anchor_len;

		(sam_rec_iter)->ave_intron_len = sum_intron_len;

		(sam_rec_iter)->ave_junc_mis = sum_ave_mis;

		(sam_rec_iter)->pair_rate = pair_rate;

		(sam_rec_iter)->canon_rate = (double)(sam_rec_iter)->canon_count / (double)((sam_rec_iter)->canon_count + (sam_rec_iter)->noncanon_count);
	}
}


void*
AlignmentHandler::MarkCanonNoncanonByReadsStatic(void * str)
{
	vector<SamRec*>& read_sam = *((vector<SamRec*>*)str);

	JunctionHandler* junc_hash_ptr = m_junction_handler_prev_ptr;

	if (read_sam.empty())
		return 0;

	vector<SamRec* >::iterator sam_rec_iter;

	for (sam_rec_iter = read_sam.begin(); sam_rec_iter != read_sam.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
			continue;

		if ((*sam_rec_iter)->is_fusion)
		{
			double sumscore = 0;

			double sum_intron_len = 0;

			double junc_anchor_len = 0;

			double junc_hits = 0;

			double sum_ave_mis = 0;

			double pair_rate = 0;

			size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

			string ins_str;

			string chrom_name;
			//not reverse alignment 2

			if ((*sam_rec_iter)->need_swap == false)
			{
				prefixst = (*sam_rec_iter)->fusion_prefix_st;

				prefixend = (*sam_rec_iter)->fusion_prefix_end;

				prefixlen = (*sam_rec_iter)->fusion_prefix_len;

				suffixst = (*sam_rec_iter)->fusion_suffix_st;

				suffixlen = (*sam_rec_iter)->fusion_suffix_len;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;

				chrom_name = (*sam_rec_iter)->chrom_name + "~" + (*sam_rec_iter)->chrom_name2;
			}
			else 
			{
				prefixst = (*sam_rec_iter)->fusion_prefix_st;

				suffixst = (*sam_rec_iter)->fusion_prefix_end;

				suffixlen = (*sam_rec_iter)->fusion_prefix_len;

				prefixend = (*sam_rec_iter)->fusion_suffix_st;

				prefixlen = (*sam_rec_iter)->fusion_suffix_len;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;

				chrom_name = (*sam_rec_iter)->chrom_name2 + "~" + (*sam_rec_iter)->chrom_name;
			}

			CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = junc_hash_ptr->m_junc_hash.find(chrom_name);

			if (chrom_junc_hash_iter == junc_hash_ptr->m_junc_hash.end())
			{
				(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

#ifdef DEBUG
				cout << "FILTERED_BY_SMALL_ANCHOR"<<endl;

				cout << (*sam_rec_iter)->tostring(0, 0) << endl;
#endif

				continue;
			}

			JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

			JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

			if (junc_hash_comb_iter == junc_hash_comb.end())
			{
				(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

				#ifdef DEBUG

				cout << "FILTERED_BY_SMALL_ANCHOR"<<endl;

				cout << (*sam_rec_iter)->tostring(0, 0) << endl;
				#endif

				continue;
			}

			if ((*sam_rec_iter)->filter_type == NOT_FILTERED)
				(*sam_rec_iter)->filter_type = (junc_hash_comb_iter->second.m_filtered_type);

			sumscore += (junc_hash_comb_iter->second.m_entropy * m_entrpy_weight) + (junc_hash_comb_iter->second.m_lpq * m_pqlen_weight) 
				+ (junc_hash_comb_iter->second.m_ave_mismatch * m_ave_mis_weight);

			sum_intron_len += junc_hash_comb_iter->second.m_intronlen;

			junc_anchor_len += junc_hash_comb_iter->second.m_max_prefix_len + junc_hash_comb_iter->second.m_max_suffix_len;

			junc_hits += junc_hash_comb_iter->second.m_positive_count + junc_hash_comb_iter->second.m_negative_count;

			sum_ave_mis += junc_hash_comb_iter->second.m_ave_mismatch;

			pair_rate += (double)junc_hash_comb_iter->second.m_paired_count;

			if (junc_hash_comb_iter->second.m_flankcase >= 5)
			{
				(*sam_rec_iter)->canon_count++;
				(*sam_rec_iter)->iscanonical = true;
			}
			else if (junc_hash_comb_iter->second.m_flankcase >= 1)
			{
				(*sam_rec_iter)->canon_count++;

				(*sam_rec_iter)->noncanon_count++;

				(*sam_rec_iter)->issemicanonical = true;
			}
			else
			{
				(*sam_rec_iter)->noncanon_count++;

				(*sam_rec_iter)->isnoncanoical = true;
			}

			(*sam_rec_iter)->filter_score = sumscore;

			(*sam_rec_iter)->junc_anchor_len = junc_anchor_len;

			(*sam_rec_iter)->ave_intron_len = 0;

			(*sam_rec_iter)->ave_junc_mis = sum_ave_mis;

			(*sam_rec_iter)->junc_hits = junc_hits;

			(*sam_rec_iter)->pair_rate = pair_rate;

			(*sam_rec_iter)->canon_rate = (double)(*sam_rec_iter)->canon_count / (double)((*sam_rec_iter)->canon_count + (*sam_rec_iter)->noncanon_count);

			if (!m_chrom_dir.empty() && (m_do_filter & 16))
			{
				(*sam_rec_iter)->flankstrings.push_back(junc_hash_comb_iter->second.m_flankstring);

				char cur_xs_tag;
				switch (junc_hash_comb_iter->second.m_flankcase) 
				{
				case 0:
					cur_xs_tag = '*';
					break;
				case 1:
					cur_xs_tag = '+';
					break;
				case 2:
					cur_xs_tag = '-';
					break;
				case 3:
					cur_xs_tag = '-';
					break;
				case 4:
					cur_xs_tag = '+';
					break;
				case 5:
					cur_xs_tag = '+';
					break;
				case 6:
					cur_xs_tag = '-';
					break;
				default:
					cout << "none specified case: "<<junc_hash_comb_iter->second.m_flankcase << endl;
					break;
				}

				if ((*sam_rec_iter)->xs_tag == cur_xs_tag)
					;
				else if ((*sam_rec_iter)->xs_tag != '-' && (*sam_rec_iter)->xs_tag != '+')
				{
					(*sam_rec_iter)->xs_tag = cur_xs_tag;
				}
				else if ((*sam_rec_iter)->xs_tag == '-' && cur_xs_tag == '+')
				{
					(*sam_rec_iter)->xs_tag = 'C';
				}
				else if ((*sam_rec_iter)->xs_tag == '+' && cur_xs_tag == '-')
				{
					(*sam_rec_iter)->xs_tag = 'C';
				}
				else if ((*sam_rec_iter)->xs_tag == '-' && cur_xs_tag == '*')
				{
				}
				else if ((*sam_rec_iter)->xs_tag == '+' && cur_xs_tag == '*')
				{
				}
				else
					(*sam_rec_iter)->xs_tag = cur_xs_tag;

			}
		}
		else
		{
			if ((*sam_rec_iter)->spliceway_vec.size() <= 1 || (*sam_rec_iter)->is_fusion)
				continue;

			vector<pair<size_t, int> >::iterator oft_mpl_iter;

			double sumscore = 0;

			double sum_intron_len = 0;

			double junc_anchor_len = 0;

			double junc_hits = 0;

			double sum_ave_mis = 0;

			double pair_rate = 0;

			for (oft_mpl_iter = (*sam_rec_iter)->spliceway_vec.begin(); oft_mpl_iter != (*sam_rec_iter)->spliceway_vec.end() - 1; ++oft_mpl_iter)
			{
				if ((oft_mpl_iter + 1)->second < 0)
				{
					continue;
				}

				size_t comb_offset;

				if (oft_mpl_iter->second < 0)
				{
					comb_offset = ((oft_mpl_iter->first - 1) << THIRTY_TWO) + oft_mpl_iter->first - 1;
				}
				else
					comb_offset = ((oft_mpl_iter->first + oft_mpl_iter->second - 1) << THIRTY_TWO) + (oft_mpl_iter + 1)->first;

				CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = junc_hash_ptr->m_junc_hash.find((*sam_rec_iter)->chrom_name);

				if (chrom_junc_hash_iter == junc_hash_ptr->m_junc_hash.end())
				{
					(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

					continue;
				}

				JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

				JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(comb_offset);

				if (junc_hash_comb_iter == junc_hash_comb.end())
				{

					(*sam_rec_iter)->filter_type = (FILTERED_BY_SMALL_ANCHOR);

					continue;
				}

				if ((*sam_rec_iter)->filter_type == NOT_FILTERED)
					(*sam_rec_iter)->filter_type = (junc_hash_comb_iter->second.m_filtered_type);

				sumscore += (junc_hash_comb_iter->second.m_entropy * m_entrpy_weight) + (junc_hash_comb_iter->second.m_lpq * m_pqlen_weight) 
					+ (junc_hash_comb_iter->second.m_ave_mismatch * m_ave_mis_weight);

				sum_intron_len += junc_hash_comb_iter->second.m_intronlen;

				junc_anchor_len += junc_hash_comb_iter->second.m_max_prefix_len + junc_hash_comb_iter->second.m_max_suffix_len;

				junc_hits += junc_hash_comb_iter->second.m_positive_count + junc_hash_comb_iter->second.m_negative_count;

				sum_ave_mis += junc_hash_comb_iter->second.m_ave_mismatch;

				pair_rate += (double)junc_hash_comb_iter->second.m_paired_count;

				if (junc_hash_comb_iter->second.m_flankcase >= 5)
				{
					(*sam_rec_iter)->canon_count++;
					(*sam_rec_iter)->iscanonical = true;
				}
				else if (junc_hash_comb_iter->second.m_flankcase >= 1)
				{
					(*sam_rec_iter)->canon_count++;

					(*sam_rec_iter)->noncanon_count++;

					(*sam_rec_iter)->issemicanonical = true;
				}
				else
				{
					(*sam_rec_iter)->noncanon_count++;

					(*sam_rec_iter)->isnoncanoical = true;
				}

				if (!m_chrom_dir.empty() && (m_do_filter & 16))
				{
					(*sam_rec_iter)->flankstrings.push_back(junc_hash_comb_iter->second.m_flankstring);

					char cur_xs_tag;
					switch (junc_hash_comb_iter->second.m_flankcase) 
					{
					case 0:
						cur_xs_tag = '*';
						break;
					case 1:
						cur_xs_tag = '+';
						break;
					case 2:
						cur_xs_tag = '-';
						break;
					case 3:
						cur_xs_tag = '-';
						break;
					case 4:
						cur_xs_tag = '+';
						break;
					case 5:
						cur_xs_tag = '+';
						break;
					case 6:
						cur_xs_tag = '-';
						break;
					default:
						cout << "none specified case: "<<junc_hash_comb_iter->second.m_flankcase << endl;
						break;
					}

					if ((*sam_rec_iter)->xs_tag == cur_xs_tag)
						;
					else if ((*sam_rec_iter)->xs_tag != '-' && (*sam_rec_iter)->xs_tag != '+')
					{
						(*sam_rec_iter)->xs_tag = cur_xs_tag;
					}
					else if ((*sam_rec_iter)->xs_tag == '-' && cur_xs_tag == '+')
					{
						(*sam_rec_iter)->xs_tag = 'C';
					}
					else if ((*sam_rec_iter)->xs_tag == '+' && cur_xs_tag == '-')
					{
						(*sam_rec_iter)->xs_tag = 'C';
					}
					else if ((*sam_rec_iter)->xs_tag == '-' && cur_xs_tag == '*')
					{
					}
					else if ((*sam_rec_iter)->xs_tag == '+' && cur_xs_tag == '*')
					{
					}
					else
						(*sam_rec_iter)->xs_tag = cur_xs_tag;
				}
			}

			sumscore = sumscore / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			sum_intron_len =  sum_intron_len / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			junc_anchor_len = junc_anchor_len / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			sum_ave_mis = sum_ave_mis / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			junc_hits = junc_hits / double ((*sam_rec_iter)->spliceway_vec.size() - 1);

			pair_rate = pair_rate / double ((*sam_rec_iter)->spliceway_vec.size() - 1);
		
			(*sam_rec_iter)->filter_score = sumscore;

			(*sam_rec_iter)->junc_anchor_len = junc_anchor_len;

			(*sam_rec_iter)->ave_intron_len = sum_intron_len;

			(*sam_rec_iter)->ave_junc_mis = sum_ave_mis;

			(*sam_rec_iter)->junc_hits = junc_hits;

			(*sam_rec_iter)->pair_rate = pair_rate;

			(*sam_rec_iter)->canon_rate = (double)(*sam_rec_iter)->canon_count / (double)((*sam_rec_iter)->canon_count + (*sam_rec_iter)->noncanon_count);
		}
	}

	return 0;
}

void AlignmentHandler::ProcessPEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			#ifdef DEBUG
			cout << "parse"<<endl;
			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_base_name;

			#ifdef DEBUG

			cout << "group same read alignments"<<endl;

			#endif

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_base_name)
				{
					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_base_name;

					previous_iter = sam_iter;
				}
			}

			#ifdef DEBUG

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			#endif

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			#ifdef DEBUG

			cout << "process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadFilterJuncMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			#ifdef DEBUG
			
			cout << "output to files"<<endl;

			#endif

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

				WriteAlignment(m_same_reads_alignments_paired_iter->first);

				WriteAlignment(m_same_reads_alignments_paired_iter->second);
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			#ifdef DEBUG

			cout << "remained lines:"<< index << endl;

			cout << "finish"<<endl;

			#endif
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{
		#ifdef DEBUG

		cout << "parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_base_name;

		#ifdef DEBUG

		cout << "group same read alignments"<<endl;

		#endif

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_base_name)
			{
				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_base_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		#ifdef DEBUG

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		#endif

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		#ifdef DEBUG

		cout << "process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadFilterJuncMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));

			WriteAlignment(m_same_reads_alignments_paired_iter->first);

			WriteAlignment(m_same_reads_alignments_paired_iter->second);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		#ifdef DEBUG

		cout << "finish"<<endl;

		#endif
	}

	#ifdef DEBUG

	cout << "processed_line:"<< processed_line << endl;

	#endif

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void AlignmentHandler::ProcessPERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			#ifdef DEBUG

			cout << "parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_base_name;

			#ifdef DEBUG

			cout << "group same read alignments"<<endl;

			#endif

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_base_name)
				{
					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_base_name;

					previous_iter = sam_iter;
				}
			}

			#ifdef DEBUG

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			#endif

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

			#ifdef DEBUG

			cout << "process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first), m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second), m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);

				WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			m_same_reads_alignments_paired_is_unspliced.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			#ifdef DEBUG

			cout << "remain line:"<< index << endl;

			cout << "finish"<<endl;

			#endif
			//keep the last read's alignments to be processed with next block of alignments
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}


	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	#ifdef DEBUG

	cout << "process the end block of file" << endl;

	#endif

	if (m_stored_sam_recs.size())
	{
		#ifdef DEBUG

		cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_base_name;

		#ifdef DEBUG

		cout << "group same read alignments"<<endl;

		#endif

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_base_name)
			{
				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_base_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		#ifdef DEBUG

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		#endif

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

		#ifdef DEBUG

		cout << "process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first), m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second), m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);

			WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		m_same_reads_alignments_paired_is_unspliced.clear();
	}

	#ifdef DEBUG

	cout << "processed_line:"<< processed_line << endl;

	#endif

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired_is_unspliced.clear();

	m_same_reads_alignments_paired.clear();

	m_prev_tag = "";

	index = 0;
}

void AlignmentHandler::ProcessPEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			#ifdef DEBUG

			cout << "parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_base_name;

			#ifdef DEBUG

			cout << "group same read alignments"<<endl;

			#endif

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_base_name)
				{
					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_base_name;

					previous_iter = sam_iter;
				}
			}

			#ifdef DEBUG

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			#endif

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			#ifdef DEBUG

			cout << "process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadSelectBestMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			cout << "remained lines:"<< index << endl;

			cout << "finish"<<endl;

			//keep the last read's alignments to be processed with next block of alignments
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{
		#ifdef DEBUG

		cout << "parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_base_name;

		#ifdef DEBUG

		cout << "group same read alignments"<<endl;

		#endif

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_base_name)
			{
				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_base_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		#ifdef DEBUG

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		#endif

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		#ifdef DEBUG

		cout << "process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessPEReadSelectBestMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->second);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->second));
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		#ifdef DEBUG

		cout << "finish"<<endl;

		#endif

		//keep the last read's alignments to be processed with next block of alignments
	}

	#ifdef DEBUG

	cout << "processed_line:"<< processed_line << endl;

	#endif

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void AlignmentHandler::ProcessSEReadFilterJunc(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			#ifdef DEBUG

			cout << "parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_name;

			#ifdef DEBUG

			cout << "group same read alignments"<<endl;

			#endif

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_name)
				{
					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_name;

					previous_iter = sam_iter;
				}
			}

			#ifdef DEBUG

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			#endif

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			#ifdef DEBUG

			cout << "process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadFilterJuncMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

				WriteAlignment(m_same_reads_alignments_paired_iter->first);
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			//cout << "remained lines:"<< index << endl;

			#ifdef DEBUG

			cout << "finish"<<endl;

			#endif

			//keep the last read's alignments to be processed with next block of alignments
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{

		#ifdef DEBUG

		cout << "parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_name;

		#ifdef DEBUG

		cout << "group same read alignments"<<endl;

		#endif

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_name)
			{

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		#ifdef DEBUG

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		#endif

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		#ifdef DEBUG

		cout << "process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadFilterJuncMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));

			WriteAlignment(m_same_reads_alignments_paired_iter->first);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		#ifdef DEBUG

		cout << "finish"<<endl;

		#endif

		//keep the last read's alignments to be processed with next block of alignments
	}

	#ifdef DEBUG

	cout << "processed_line:"<< processed_line << endl;

	#endif

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void AlignmentHandler::ProcessSERead(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			#ifdef DEBUG

			cout << "parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_name;

			#ifdef DEBUG

			cout << "group same read alignments"<<endl;

			#endif

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_name)
				{

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_name;

					previous_iter = sam_iter;
				}
			}

			#ifdef DEBUG

			cout << "m_stored_sam_recs size:" << m_stored_sam_recs.size() << endl;

			cout << "last same read index:" << previous_iter - m_stored_sam_recs.begin() << endl;

			//cout << previous_iter->tostring(0, 0) << endl;

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			#endif

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

			#ifdef DEBUG

			cout << "process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first), m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);

				WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			m_same_reads_alignments_paired_is_unspliced.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			#ifdef DEBUG

			cout << "remain line:"<< index << endl;

			cout << "finish"<<endl;

			#endif
			//keep the last read's alignments to be processed with next block of alignments
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}


	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	#ifdef DEBUG

	cout << "process the end block of file" << endl;

	#endif

	if (m_stored_sam_recs.size())
	{
		#ifdef DEBUG

		cout << m_stored_sam_recs.size() << endl;

		cout << "parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_name;

		#ifdef DEBUG

		cout << "group same read alignments"<<endl;

		#endif

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_name)
			{

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		#ifdef DEBUG

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		#endif

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		m_same_reads_alignments_paired_is_unspliced.resize(m_same_reads_alignments.size(), false);

		#ifdef DEBUG

		cout << "process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first), m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);

			WriteAlignmentUnspliced(m_same_reads_alignments_paired_iter, m_same_reads_alignments_paired_is_unspliced[m_same_reads_alignments_paired_iter -  m_same_reads_alignments_paired.begin()]);
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		m_same_reads_alignments_paired_is_unspliced.clear();
	}

	#ifdef DEBUG

	cout << "processed_line:"<< processed_line << endl;

	#endif

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	m_same_reads_alignments_paired_is_unspliced.clear();

	m_prev_tag = "";

	index = 0;
}

void AlignmentHandler::ProcessSEReadSelectBest(JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered)
{
	size_t index = 0;

#ifdef LINUX
	vector<pthread_t> iThreadIds(threads_number);

	vector<int> return_values(threads_number);
#endif

	size_t processed_line = 0;

	while (!m_ifs_alignments.eof())
	{
		if (index == buf_size)
		{
			#ifdef DEBUG

			cout << "parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait parse"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

			vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

			m_prev_tag = sam_iter->tag_name;

			#ifdef DEBUG

			cout << "group same read alignments"<<endl;

			#endif

			for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
			{
				if (m_prev_tag != sam_iter->tag_name)
				{

					m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

					processed_line += sam_iter - 1 - previous_iter + 1;

					m_prev_tag = sam_iter->tag_name;

					previous_iter = sam_iter;
				}
			}

			#ifdef DEBUG

			cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

			#endif

			m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

			#ifdef DEBUG

			cout << "process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadSelectBestMulti, NULL);
#endif
			}

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				if (return_values[i] < 0)
					cout << "error creating thread:"<< errno<<endl;
#endif
			}

			#ifdef DEBUG

			cout << "wait process same read alignments"<<endl;

			#endif

			for (size_t i = 0; i < threads_number; ++i)
			{
#ifdef LINUX
				pthread_join (iThreadIds[i], NULL);
#endif
			}

			vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

			for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
				m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
				++m_same_reads_alignments_paired_iter)
			{
				junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

				CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));
			}

			m_cur_read_alignments_index = 0;

			m_cur_line_index = 0;

			m_same_reads_alignments.clear();

			m_same_reads_alignments_paired.clear();

			size_t remained_lines = m_stored_sam_recs.end() - previous_iter;

			size_t processed_lines = m_stored_lines.size() - remained_lines;

			for (size_t i = 0; i < remained_lines; ++i)
			{
				m_stored_lines[i] = m_stored_lines[processed_lines + i];
			}

			index = remained_lines;

			#ifdef DEBUG

			cout << "finish"<<endl;

			#endif
		}

		getline(m_ifs_alignments, m_stored_lines[index++]);
	}

	//process the end block of file
	m_stored_lines.resize(index - 1);

	m_stored_sam_recs.resize(index - 1);

	if (m_stored_sam_recs.size())
	{

		#ifdef DEBUG

		cout << "parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ParseLines, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait parse"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<SamRec>::iterator sam_iter = m_stored_sam_recs.begin();

		vector<SamRec>::iterator previous_iter = m_stored_sam_recs.begin();

		m_prev_tag = sam_iter->tag_name;

		#ifdef DEBUG

		cout << "group same read alignments"<<endl;

		#endif

		for (; sam_iter !=  m_stored_sam_recs.end(); ++sam_iter)
		{
			if (m_prev_tag != sam_iter->tag_name)
			{

				m_same_reads_alignments.push_back(make_pair(previous_iter, sam_iter - 1));

				processed_line += sam_iter - 1 - previous_iter + 1;

				m_prev_tag = sam_iter->tag_name;

				previous_iter = sam_iter;
			}
		}

		m_same_reads_alignments.push_back(make_pair(previous_iter, m_stored_sam_recs.end() - 1));

		processed_line += m_stored_sam_recs.end() - 1 - previous_iter + 1;

		#ifdef DEBUG

		cout << "total grouped reads:"<< m_same_reads_alignments.size()<<endl;

		#endif

		m_same_reads_alignments_paired.resize(m_same_reads_alignments.size());

		#ifdef DEBUG

		cout << "process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			return_values[i] = pthread_create(&(iThreadIds[i]), NULL, &ProcessSEReadSelectBestMulti, NULL);
#endif
		}

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			if (return_values[i] < 0)
				cout << "error creating thread:"<< errno<<endl;
#endif
		}

		#ifdef DEBUG

		cout << "wait process same read alignments"<<endl;

		#endif

		for (size_t i = 0; i < threads_number; ++i)
		{
#ifdef LINUX
			pthread_join (iThreadIds[i], NULL);
#endif
		}

		vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator m_same_reads_alignments_paired_iter;

		for (m_same_reads_alignments_paired_iter = m_same_reads_alignments_paired.begin(); 
			m_same_reads_alignments_paired_iter != m_same_reads_alignments_paired.end();
			++m_same_reads_alignments_paired_iter)
		{
			junc_handler->SamRecVec2Junc(m_same_reads_alignments_paired_iter->first);

			CollectStatsStatic((void*)&(m_same_reads_alignments_paired_iter->first));
		}

		m_cur_read_alignments_index = 0;

		m_cur_line_index = 0;

		m_same_reads_alignments.clear();

		m_same_reads_alignments_paired.clear();

		#ifdef DEBUG

		cout << "finish"<<endl;

		#endif

		//keep the last read's alignments to be processed with next block of alignments
	}

	#ifdef DEBUG

	cout << "processed_line:"<< processed_line << endl;

	#endif

	m_stored_lines.resize(buf_size);

	m_stored_sam_recs.resize(buf_size);

	m_cur_read_alignments_index = 0;

	m_cur_line_index = 0;

	m_same_reads_alignments.clear();

	m_same_reads_alignments_paired.clear();

	index = 0;

	m_prev_tag = "";
}

void* AlignmentHandler::ParseLines(void * str)
{
	while(true) 
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_line_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_stored_lines.size())
		{
			break;
		}

		m_stored_sam_recs[read_index].Set(m_stored_lines[read_index]/*, 3*/);
	}

	return 0;
}

void* AlignmentHandler::ProcessPEReadMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		bool are_both_ends_unspliced = true;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if ((sam_rec_iter->isexonic == false && sam_rec_iter->is_fusion == false) || (sam_rec_iter->is_fusion == true))
					are_both_ends_unspliced = false;

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}
			}
			else
			{
				are_both_ends_unspliced = false;

				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					#ifdef DEBUG

					cout << "new format fusion alignment not together:"<<endl <<sam_rec_iter->tostring(0, 0 ) << endl;

					#endif

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}

				++sam_rec_iter;
			}
		}

		int tmp_do_filter = m_do_filter;

		if (are_both_ends_unspliced)
			tmp_do_filter = tmp_do_filter | 16;

		//cout << "remove duplication 1"<<endl;

		if (cur_same_read_alignments_paired.first.size())
			RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);

		//cout << "remove duplication 2"<<endl;
		if (cur_same_read_alignments_paired.second.size())
			RemoveDupStatic((void*)&cur_same_read_alignments_paired.second);

		SetSamRecUnpaired(cur_same_read_alignments_paired);

		vector<PairedSamRec> fusion_paired_reads;

		//cout << " EstablishPairingStatic " << endl;

		bool paired = false;

		EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads, tmp_do_filter);

		paired = EstablishPairingStatic((void *) &est_pairing);

		//cout << " FilterSingleMultiStatic " << endl;

		if (paired == false && are_both_ends_unspliced) //do filter multiple unpaired
		{
				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

				FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		if (are_both_ends_unspliced)
		{
			m_same_reads_alignments_paired_is_unspliced[read_index] = true;

			//Set expressed region
			SetHits(cur_same_read_alignments_paired.first);

			SetHits(cur_same_read_alignments_paired.second);
			//Set bits
		}

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			bool isend1mapped = false, isend2mapped = false;

			for (size_t j = 0; j < cur_same_read_alignments_paired.first.size(); ++j)
			{
				if (cur_same_read_alignments_paired.first[j]->isunmapped == false)
				{
					isend1mapped = true;

					break;
				}
			}

			for (size_t j = 0; j < cur_same_read_alignments_paired.second.size(); ++j)
			{
				if (cur_same_read_alignments_paired.second[j]->isunmapped == false)
				{
					isend2mapped = true;

					break;
				}
			}

			if (isend1mapped && isend2mapped/*cur_same_read_alignments_paired.second.size()*/)
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
			}
			else
			{
				if (isend1mapped)
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.first[i]->is_fusion)
							cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
				}
				else if (isend2mapped)
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.second[i]->is_fusion)
							cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
				}
				else
					;
			}
		}
		else
		{
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first), are_both_ends_unspliced);

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.second), are_both_ends_unspliced);

	}  

	return 0;
}

void* AlignmentHandler::ProcessSEReadMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		bool are_both_ends_unspliced = true;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if (sam_rec_iter->isexonic == false && sam_rec_iter->is_fusion == false)
					are_both_ends_unspliced = false;
			}

			cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
		}

		int tmp_do_filter = m_do_filter;

		if (are_both_ends_unspliced)
			tmp_do_filter = tmp_do_filter | 16;

		if (cur_same_read_alignments_paired.first.size() > 0)
			RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);

		bool paired = false;

		if (paired == false && are_both_ends_unspliced) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));
		}

		if (are_both_ends_unspliced)
		{
			m_same_reads_alignments_paired_is_unspliced[read_index] = true;
			//Set bits
		}

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
				}
				else
					;
			}
		}
		else
		{
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);
		}

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first), are_both_ends_unspliced);

	}  

	return 0;
}

void* AlignmentHandler::ProcessPEReadFilterJuncMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					#ifdef DEBUG

					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					#endif

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}

				++sam_rec_iter;
			}
		}

		SetSamRecUnpaired(cur_same_read_alignments_paired);

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.second));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		bool paired = false;

		vector<PairedSamRec> fusion_paired_reads;

		EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads, m_do_filter);

		paired = EstablishPairingStatic((void *) &est_pairing);

		if (paired == false) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));

			EstablishFusionPairingStatic((void *) &est_pairing);

			FindFusionJuncRegionVecStatic((void*)&fusion_paired_reads);
		}		

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			bool isend1mapped = false, isend2mapped = false;

			for (size_t j = 0; j < cur_same_read_alignments_paired.first.size(); ++j)
			{
				if (cur_same_read_alignments_paired.first[j]->isunmapped == false)
				{
					isend1mapped = true;

					break;
				}
			}

			for (size_t j = 0; j < cur_same_read_alignments_paired.second.size(); ++j)
			{
				if (cur_same_read_alignments_paired.second[j]->isunmapped == false)
				{
					isend2mapped = true;

					break;
				}
			}

			if (isend1mapped && isend2mapped)
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
			}
			else
			{
				if (isend1mapped)
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.first[i]->is_fusion)
						{
							cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
						}
					}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
				}
				else if (isend2mapped)
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.second[i]->is_fusion)
						{
							cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
						}
					}

					SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
				}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first));

		GenerateAlignment((void *)&(cur_same_read_alignments_paired.second));

	}  

	return 0;
}

void* AlignmentHandler::ProcessPEReadSelectBestMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					#ifdef DEBUG

					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					#endif

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				if (sam_rec_iter->end_id == 1)
				{
					cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
				}
				else if (sam_rec_iter->end_id == 2)
				{
					cur_same_read_alignments_paired.second.push_back(&(*sam_rec_iter));
				}

				++sam_rec_iter;
			}
		}

		SetSamRecUnpaired(cur_same_read_alignments_paired);

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.second));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		bool paired = false;

		vector<PairedSamRec> fusion_paired_reads;

		EstPairingStruct est_pairing(&cur_same_read_alignments_paired, &fusion_paired_reads, m_do_filter);

		paired = EstablishPairingStatic((void *) &est_pairing);

		if (paired == false) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.second));
		}

		//Set expressed region
		SetHits(cur_same_read_alignments_paired.first);

		SetHits(cur_same_read_alignments_paired.second);

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			bool isend1mapped = false, isend2mapped = false;

			for (size_t j = 0; j < cur_same_read_alignments_paired.first.size(); ++j)
			{
				if (cur_same_read_alignments_paired.first[j]->isunmapped == false)
				{
					isend1mapped = true;

					break;
				}
			}

			for (size_t j = 0; j < cur_same_read_alignments_paired.second.size(); ++j)
			{
				if (cur_same_read_alignments_paired.second[j]->isunmapped == false)
				{
					isend2mapped = true;

					break;
				}
			}

			if (isend1mapped && isend2mapped)
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);

				SetPairedType(cur_same_read_alignments_paired.second, FUSION_PAIRED);
			}
			else
			{
				if (isend1mapped)
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.first.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.first[i]->is_fusion)
							cur_same_read_alignments_paired.first[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
				}
				else if (isend2mapped)
				{
					for (size_t i = 0; i < cur_same_read_alignments_paired.second.size(); ++i)
					{
						if (!cur_same_read_alignments_paired.second[i]->is_fusion)
							cur_same_read_alignments_paired.second[i]->strand_t |= MATE_UNMAPPED;
					}

					SetPairedType(cur_same_read_alignments_paired.second, SINGLE);
				}
				else
					;
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);

			SetPairedType(cur_same_read_alignments_paired.second, NORMAL_PAIRED);
		}

		SetBitInfoStatic((void*)&(cur_same_read_alignments_paired));
	}  

	return 0;
}

void* AlignmentHandler::ProcessSEReadFilterJuncMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					#ifdef DEBUG

					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					#endif

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));

				++sam_rec_iter;
			}
		}

		RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));
		}

		bool paired = false;

		if (paired == false) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));
		}

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
				}
			}
		}
		else
		{
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);
		}
	
		GenerateAlignment((void *)&(cur_same_read_alignments_paired.first));
	}  

	return 0;
}

void* AlignmentHandler::ProcessSEReadSelectBestMulti(void * str)
{
	while (true)
	{
		size_t read_index = 0;

		#ifdef LINUX
		pthread_mutex_lock(&inc_num_threads);
		
		read_index = m_cur_read_alignments_index++;

		pthread_mutex_unlock(&inc_num_threads);
		#endif

		if (read_index >= m_same_reads_alignments.size())
		{
			break;
		}

		pair<vector<SamRec>::iterator, vector<SamRec>::iterator>& cur_same_read_alignments = m_same_reads_alignments[read_index];
		
		pair<vector<SamRec*>, vector<SamRec*> >& cur_same_read_alignments_paired = m_same_reads_alignments_paired[read_index];

		vector<SamRec>::iterator sam_rec_iter;

		for (sam_rec_iter = cur_same_read_alignments.first; sam_rec_iter <= cur_same_read_alignments.second; ++sam_rec_iter)
		{
			if (!sam_rec_iter->is_fusion_newfmt)
			{
				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));
			}
			else
			{
				if (sam_rec_iter + 1 > cur_same_read_alignments.second || (sam_rec_iter + 1)->is_fusion_newfmt == false)
				{
					#ifdef DEBUG

					cout << "new format fusion alignment not together:"<<endl<<sam_rec_iter->tostring(0, 0 ) << endl;

					#endif

					continue;
				}

				MergeStdFusion(sam_rec_iter, sam_rec_iter + 1);

				cur_same_read_alignments_paired.first.push_back(&(*sam_rec_iter));

				++sam_rec_iter;
			}
		}

		RemoveDupStatic((void*)&cur_same_read_alignments_paired.first);

		if (m_do_filter & 4)
		{
			MarkCanonNoncanonByReadsStatic((void*)&(cur_same_read_alignments_paired.first));

			FilterByFilteredJunctionStatic((void*)&(cur_same_read_alignments_paired.first));
		}

		bool paired = false;

		if (paired == false) //do filter multiple unpaired
		{
			FilterSingleMultiStatic((void*)&(cur_same_read_alignments_paired.first));
		}

		if (paired == false)
		{
			SetFusionBitInfoStatic((void*)&(cur_same_read_alignments_paired));

			if (cur_same_read_alignments_paired.first.size() && cur_same_read_alignments_paired.second.size())
			{
				SetPairedType(cur_same_read_alignments_paired.first, FUSION_PAIRED);
			}
			else
			{
				if (cur_same_read_alignments_paired.first.size())
				{
					SetPairedType(cur_same_read_alignments_paired.first, SINGLE);
				}
			}
		}
		else
		{
			//cerr <<"paired"<< endl;
			SetPairedType(cur_same_read_alignments_paired.first, NORMAL_PAIRED);
		}
	}  

	return 0;
}

void AlignmentHandler::ProcessFile(string alignment_file, JunctionHandler* junc_handler, JunctionHandler* junc_handler_filtered, int fileter_by_junc)
{
	m_junction_handler_prev_ptr = junc_handler_filtered;

	
	junc_handler->ClearJunction();

	if (fileter_by_junc == 0)
	{
		if (m_do_filter & 32)
		{
			if (!m_ofs_unspliced_fusion_paired.is_open() && !m_ofs_unspliced_single.is_open() && !m_ofs_unspliced_paired.is_open() && !m_ofs_onespliced.is_open())
			{
				string bothunspliced_fusion_paired = alignment_file; bothunspliced_fusion_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_fusion_paired.append(".bothunspliced.fusion_paired");

				m_ofs_unspliced_fusion_paired.open(bothunspliced_fusion_paired.c_str());

				string bothunspliced_single = alignment_file; bothunspliced_single.append(basename2(m_filtered_alignment_file)); bothunspliced_single.append(".bothunspliced.single");

				m_ofs_unspliced_single.open(bothunspliced_single.c_str());

				string bothunspliced_paired = alignment_file; bothunspliced_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_paired.append(".bothunspliced.paired");

				m_ofs_unspliced_paired.open(bothunspliced_paired.c_str());

				string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

				m_ofs_onespliced.open(onespliced.c_str());

				m_ifs_alignments.open(alignment_file.c_str());
			}
		}
		else
		{
			m_ofs_unspliced_fusion_paired.close();

			m_ofs_unspliced_single.close();

			m_ofs_unspliced_paired.close();

			m_ofs_onespliced.close();

			string bothunspliced_fusion_paired = alignment_file; bothunspliced_fusion_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_fusion_paired.append(".bothunspliced.fusion_paired");

			m_ofs_unspliced_fusion_paired.open(bothunspliced_fusion_paired.c_str());

			string bothunspliced_single = alignment_file; bothunspliced_single.append(basename2(m_filtered_alignment_file)); bothunspliced_single.append(".bothunspliced.single");

			m_ofs_unspliced_single.open(bothunspliced_single.c_str());

			string bothunspliced_paired = alignment_file; bothunspliced_paired.append(basename2(m_filtered_alignment_file)); bothunspliced_paired.append(".bothunspliced.paired");

			m_ofs_unspliced_paired.open(bothunspliced_paired.c_str());

			string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

			m_ofs_onespliced.open(onespliced.c_str());

			m_ifs_alignments.open(alignment_file.c_str());
		}
	}
	else
	{
		if (m_do_filter & 32)
		{
			if (!m_ofs_filtered_alignment.is_open() && 
				!m_ofs_fusion_paired.is_open() && 
				!m_ofs_single.is_open() && 
				!m_ofs_fusion_std.is_open() && 
				!m_ofs_paired.is_open())
			{
				string to_mapper_str = alignment_file; to_mapper_str.append(basename2(m_filtered_alignment_file));

				m_ofs_filtered_alignment.open(to_mapper_str.c_str());

				string fusion_paired_str = alignment_file; fusion_paired_str.append(basename2(m_filtered_alignment_file)); fusion_paired_str.append(".fusion_paired");

				m_ofs_fusion_paired.open(fusion_paired_str.c_str());

				string single_str = alignment_file; single_str.append(basename2(m_filtered_alignment_file)); single_str.append(".single");

				m_ofs_single.open(single_str.c_str());

				string to_fusion_str = alignment_file; to_fusion_str.append(basename2(m_filtered_alignment_file)); to_fusion_str.append(".fusion");

				m_ofs_fusion_std.open(to_fusion_str.c_str());

				string paired_str = alignment_file; paired_str.append(basename2(m_filtered_alignment_file)); paired_str.append(".paired");

				m_ofs_paired.open(paired_str.c_str());

				string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

				m_ifs_alignments.open(onespliced.c_str());
			}
			else
			{
				#ifdef DEBUG

				cout << "returned ProcessFile"<<endl;

				#endif

				return;
			}
		}
		else
		{
			m_ofs_filtered_alignment.close();

			m_ofs_fusion_std.close();

			m_ofs_fusion_paired.close();

			m_ofs_single.close();

			m_ofs_paired.close();

			string to_mapper_str = alignment_file; to_mapper_str.append(basename2(m_filtered_alignment_file));

			m_ofs_filtered_alignment.open(to_mapper_str.c_str());

			string fusion_paired_str = alignment_file; fusion_paired_str.append(basename2(m_filtered_alignment_file)); fusion_paired_str.append(".fusion_paired");

			m_ofs_fusion_paired.open(fusion_paired_str.c_str());

			string single_str = alignment_file; single_str.append(basename2(m_filtered_alignment_file)); single_str.append(".single");

			m_ofs_single.open(single_str.c_str());

			string to_fusion_str = alignment_file; to_fusion_str.append(basename2(m_filtered_alignment_file)); to_fusion_str.append(".fusion");

			m_ofs_fusion_std.open(to_fusion_str.c_str());

			string paired_str = alignment_file; paired_str.append(basename2(m_filtered_alignment_file)); paired_str.append(".paired");

			m_ofs_paired.open(paired_str.c_str());

			string onespliced = alignment_file; onespliced.append(basename2(m_filtered_alignment_file)); onespliced.append(".onespliced");

			m_ifs_alignments.open(onespliced.c_str());
		}
	}
	

	//parallel

	if (m_is_paired)
	{
		//	//process each read
		if (fileter_by_junc == 2)
			ProcessPEReadSelectBest(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc ==	1)
			ProcessPEReadFilterJunc(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc == 0)
			ProcessPERead(junc_handler, junc_handler_filtered);
	}
	else//SE reads
	{
		if (fileter_by_junc == 2)
			ProcessSEReadSelectBest(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc ==	1)
			ProcessSEReadFilterJunc(junc_handler, junc_handler_filtered);
		else if (fileter_by_junc == 0)
			ProcessSERead(junc_handler, junc_handler_filtered);
	}

	m_ifs_alignments.close();
}


void 
AlignmentHandler::FilterAlignment()
{
	if (!m_junction_file.empty())
	{
		vector<string> annotated_junction;

		annotated_junction.push_back(m_junction_file);

		#ifdef DEBUG

		cout << "read annotated junctions"<< endl;

		#endif

		m_junction_handler_annotated.ReadJunction(annotated_junction);
	}

	m_junction_handler.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch,  m_min_junc_anchor, m_min_fusion_coverage, m_do_filter);

	m_junction_handler_filtered.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor, m_min_fusion_coverage, m_do_filter);

	m_junction_handler_intermediate.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor, m_min_fusion_coverage, m_do_filter);

	doner_side_spanning_pairs_ofs = new ofstream; doner_side_spanning_pairs_ofs->open((m_filtered_alignment_file + ".doner_side_spanning_pairs").c_str());
	accetpr_side_spanning_pairs_ofs = new ofstream; accetpr_side_spanning_pairs_ofs->open((m_filtered_alignment_file + ".accetpr_side_spanning_pairs").c_str());
	single_spanning_ofs = new ofstream; single_spanning_ofs->open((m_filtered_alignment_file + ".single_spanning").c_str());
	spliceway_true_ofs = new ofstream; spliceway_true_ofs->open((m_filtered_alignment_file + ".spliceway_true").c_str());
	m_fusion_encompassing_reads_doner_ofs = new ofstream; m_fusion_encompassing_reads_doner_ofs->open((m_filtered_alignment_file + ".m_fusion_encompassing_reads_doner").c_str());
	m_fusion_encompassing_reads_acceptor_ofs = new ofstream; m_fusion_encompassing_reads_acceptor_ofs->open((m_filtered_alignment_file + ".m_fusion_encompassing_reads_acceptor").c_str());

	string unmapped_str = m_filtered_alignment_file; unmapped_str.append(".unmapped");

	m_ofs_unmapped.open(unmapped_str.c_str());

	#ifdef DEBUG

	cout << "Allocate memory for Hits" << endl;
	
	#endif

	m_junction_handler.ReadChromSize(m_chrom_size_file);
	
	AllocateHitsMemory();

	JunctionHandler* junction_handler_ptr = &m_junction_handler;

	JunctionHandler* junction_handler_filtered_ptr = &m_junction_handler_filtered;

	JunctionHandler* junction_handler_intermediate_ptr = &m_junction_handler_intermediate;

	string stat_file = m_filtered_alignment_file; stat_file.append(".stat");

	string final_stat_file = m_filtered_alignment_file; final_stat_file.append(".final_stat.txt");

	vector<string>::iterator sam_file_iter;

	clock_t t1, t2;

	time_t m1, m2;

	double m_any_time;

	t1=clock();

	m1 = time(NULL);

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_paired.close();

	m_ofs_unspliced_fusion_paired.close();

	m_ofs_unspliced_single.close();

	m_ofs_unspliced_paired.close();

	m_ofs_onespliced.close();

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_ptr, 0);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first process: " << m_any_time << endl;

	cout << "first process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;
	
	//cout << "load flank string 1

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first load flank string: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	WriteStats(stat_file, "original alignment stats");

	ClearStats();

	m_junction_handler.ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "original junction stats");

	string ori_junction_file = m_filtered_alignment_file; ori_junction_file.append(".ori.junc");

	string ori_junction_ins_file = m_filtered_alignment_file; ori_junction_ins_file.append(".ori.junc.ins");

	string ori_junction_del_file = m_filtered_alignment_file; ori_junction_del_file.append(".ori.junc.del");

	string ori_junction_fusion_file = m_filtered_alignment_file; ori_junction_fusion_file.append(".ori.junc.fusion");

	string ori_junction_filter_file = m_filtered_alignment_file; ori_junction_filter_file.append(".ori.junc.filter");

	m_junction_handler.WriteJunction(ori_junction_file, ori_junction_ins_file, ori_junction_del_file, ori_junction_fusion_file, ori_junction_filter_file);

	cout << "mark fitered junction "<< endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best
	//16 
	//32 if set, open one output file for all input sam files 
	//   if not, open an output file for each input sam file
	//64 if set, output all alignments to a file for 2nd and 3rd steps process
	//   if not set, output both unspliced to a file, the others to another file for 2nd and 3rd step process
	//128 if 64 is not set
	         // if 128 is not set, output all both unspliced to a file
	         // if 128 is set, output both unspliced and normal paired to a file, and single and fusion paired to another file for 2nd and 3rd step process	  
	//256 if set, not filter single reads if 128 is not set, for normal filtering with fusion step and paired end reads
	//    if not set, filter single end reads

	if (m_do_filter & 4)
	{
		m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);

		m_junction_handler.ClearStats();

		m_junction_handler.CollectStats();

		m_junction_handler.WriteStats(stat_file, "original filtered junction stats");

		string ori_fil_junction_file = m_filtered_alignment_file; ori_fil_junction_file.append(".ori.fil.junc");

		string ori_fil_junction_ins_file = m_filtered_alignment_file; ori_fil_junction_ins_file.append(".ori.fil.junc.ins");

		string ori_fil_junction_del_file = m_filtered_alignment_file; ori_fil_junction_del_file.append(".ori.fil.junc.del");

		string ori_fil_junction_fusion_file = m_filtered_alignment_file; ori_fil_junction_fusion_file.append(".ori.fil.junc.fusion");

		string ori_fil_junction_filter_file = m_filtered_alignment_file; ori_fil_junction_filter_file.append(".ori.fil.junc.filter");

		m_junction_handler.WriteJunction(ori_fil_junction_file, ori_fil_junction_ins_file, ori_fil_junction_del_file, ori_fil_junction_fusion_file, ori_fil_junction_filter_file);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "filtering: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	/***********add intermediate filtering here*****************/

	t1=clock();

	m1 = time(NULL);

	if (m_do_filter & 8)//do select best
		m_do_filter = m_do_filter | 16;

	cerr << "Second filtering"<< endl;

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_paired.close();

	m_ofs_unspliced_fusion_paired.close();

	m_ofs_unspliced_single.close();

	m_ofs_unspliced_paired.close();

	m_ofs_onespliced.close();

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_ptr, 2);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cout << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;


	cout << "Convert hits to expressed regions " << endl;

	cerr << "Convert hits to expressed regions " << endl;

	t1=clock();

	m1 = time(NULL);

	Hits2ExpressedRegions();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "Hits to expression regions: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	WriteStats(stat_file, "intermediate alignment stats");

	ClearStats();

	m_junction_handler.ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "intermediate junction stats");

	string int_junction_file = m_filtered_alignment_file; int_junction_file.append(".inter.junc");

	string int_junction_ins_file = m_filtered_alignment_file; int_junction_ins_file.append(".inter.junc.ins");

	string int_junction_del_file = m_filtered_alignment_file; int_junction_del_file.append(".inter.junc.del");

	string int_junction_fusion_file = m_filtered_alignment_file; int_junction_fusion_file.append(".inter.junc.fusion");

	string int_junction_filter_file = m_filtered_alignment_file; int_junction_filter_file.append(".inter.junc.filter");

	m_junction_handler.WriteJunction(int_junction_file, int_junction_ins_file, int_junction_del_file, int_junction_fusion_file, int_junction_filter_file);

	cout << "mark intermediate fitered junction "<< endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best
	//16 
	//32 if set, open one output file for all input sam files 
	//   if not, open an output file for each input sam file
	//64 if set, output all alignments to a file for 2nd and 3rd steps process
	//   if not set, output both unspliced to a file, the others to another file for 2nd and 3rd step process
	//128 if 64 is not set
	         // if 128 is not set, output all both unspliced to a file
	         // if 128 is set, output both unspliced and normal paired to a file, and single and fusion paired to another file for 2nd and 3rd step process	

	if (m_do_filter & 4)
	{
		m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);

		t1=clock();

		m1 = time(NULL);

		cout << "Union expression regions " << endl;

		cerr << "Union expression regions " << endl;

		UnionSets(&m_junction_handler, &m_expressed_regions);

		t2=clock();

		m2 = time(NULL);

		m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

		cerr << "Union sets: " << m_any_time << endl;

		cerr << "seconds passed " << m2 - m1 << endl;

		t1=clock();

		m1 = time(NULL);

		cout << "Filter fusion by isolated exons " << endl;

		cerr << "Filter fusion by isolated exons " << endl;

		t2=clock();

		m2 = time(NULL);

		m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

		cerr << "Filter fusion by isolated exons finished: " << m_any_time << endl;

		cerr << "seconds passed " << m2 - m1 << endl;

		m_junction_handler.ClearStats();

		m_junction_handler.CollectStats();

		m_junction_handler.WriteStats(stat_file, "intermediate filtered junction stats");

		string int_fil_junction_file = m_filtered_alignment_file; int_fil_junction_file.append(".inter.fil.junc");

		string int_fil_junction_ins_file = m_filtered_alignment_file; int_fil_junction_ins_file.append(".inter.fil.junc.ins");

		string int_fil_junction_del_file = m_filtered_alignment_file; int_fil_junction_del_file.append(".inter.fil.junc.del");

		string int_fil_junction_fusion_file = m_filtered_alignment_file; int_fil_junction_fusion_file.append(".inter.fil.junc.fusion");

		string int_fil_junction_filter_file = m_filtered_alignment_file; int_fil_junction_filter_file.append(".inter.fil.junc.filter");

		m_junction_handler.WriteJunction(int_fil_junction_file, int_fil_junction_ins_file, int_fil_junction_del_file, int_fil_junction_fusion_file, int_fil_junction_filter_file);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "filtering: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;


	/****************************/


	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFusionJuncToSortVec();

	m_junction_handler.SortFusionJunc();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process prepare: " << m_any_time << endl;

	cout << "third process prepare: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best
	//16 
	//32 if set, open one output file for all input sam files 
	//   if not, open an output file for each input sam file
	//64 if set, output all alignments to a file for 2nd and 3rd steps process
	//   if not set, output both unspliced to a file, the others to another file for 2nd and 3rd step process
	//128 if 64 is not set
	         // if 128 is not set, output all both unspliced to a file
	         // if 128 is set, output both unspliced and normal paired to a file, and single and fusion paired to another file for 2nd and 3rd step process	

	t1=clock();

	m1 = time(NULL);

	if (m_do_filter & 8)//do select best
		m_do_filter = m_do_filter | 16;

	if (m_do_filter & 128)//do select best
	{
		m_do_filter = m_do_filter | 512;

		global_do_filter = global_do_filter | 512;
	}

	cerr << "third filtering"<< endl;

	cout << "third filtering"<< endl;

	m_ofs_filtered_alignment.close();

	m_ofs_fusion_std.close();

	m_ofs_fusion_paired.close();

	m_ofs_single.close();

	m_ofs_paired.close();

	m_ofs_unspliced_fusion_paired.close();

	m_ofs_unspliced_single.close();

	m_ofs_unspliced_paired.close();

	m_ofs_onespliced.close();

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_ptr, 1);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	if (m_do_filter & 4)
	{
		//cerr << "final filtering"<<endl;

		//m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);
	}

	WriteStats(stat_file, "filtered alignment stats");

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "filtered junction");

	WriteStats(final_stat_file, "filtered alignment stats");

	m_junction_handler.WriteStats(final_stat_file, "filtered junction");

	string fil_junction_file = m_filtered_alignment_file; fil_junction_file.append(".fil.junc");

	string fil_junction_ins_file = m_filtered_alignment_file; fil_junction_ins_file.append(".fil.junc.ins");

	string fil_junction_del_file = m_filtered_alignment_file; fil_junction_del_file.append(".fil.junc.del");

	string fil_junction_fusion_file = m_filtered_alignment_file; fil_junction_fusion_file.append(".fil.junc.fusion");

	string fil_junction_filter_file = m_filtered_alignment_file; fil_junction_filter_file.append(".fil.junc.filter");

	m_junction_handler.WriteJunction(fil_junction_file, fil_junction_ins_file, fil_junction_del_file, fil_junction_fusion_file, fil_junction_filter_file);
	
	string ori_junction_fusion_encompass_file = m_filtered_alignment_file; ori_junction_fusion_encompass_file.append(".ori.junc.fusion.encompass");

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third write junction: " << m_any_time << endl;

	cout << "third write junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	cerr << m_do_filter << endl;

	cerr << "WriteFusionJunctionWithEncompassingAlignments filtering"<< endl;

	cout << "WriteFusionJunctionWithEncompassingAlignments filtering"<< endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.WriteFusionJunctionWithEncompassingAlignments(ori_junction_fusion_encompass_file, junction_handler_ptr);

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "write fusion junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;
}

void
AlignmentHandler::FilterAlignmentFiltered()
{
	if (!m_junction_file.empty())
	{
		vector<string> annotated_junction;

		annotated_junction.push_back(m_junction_file);

		cout << "read annotated junctions"<< endl;

		m_junction_handler_annotated.ReadJunction(annotated_junction);
	}

	m_junction_handler.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch,  m_min_junc_anchor, m_min_fusion_coverage, m_do_filter);

	m_junction_handler_filtered.Init(m_sam_files, m_max_read_width, m_min_ins, m_max_del, m_chrom_dir, m_min_anchor, m_min_mismatch, m_min_junc_anchor, m_min_fusion_coverage, m_do_filter);

	cout << "Allocate memory for Hits" << endl;

	JunctionHandler* junction_handler_ptr = &m_junction_handler;

	JunctionHandler* junction_handler_filtered_ptr = &m_junction_handler_filtered;

	JunctionHandler* junction_handler_intermediate_ptr = &m_junction_handler_intermediate;

	string stat_file = m_filtered_alignment_file; stat_file.append(".stat");

	vector<string>::iterator sam_file_iter;

	clock_t t1, t2;

	time_t m1, m2;

	double m_any_time;

	t1=clock();

	m1 = time(NULL);

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_ptr, junction_handler_filtered_ptr, 0);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "first load flank string: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	WriteStats(stat_file, "original alignment stats");

	ClearStats();

	m_junction_handler.ClearStats();

	m_junction_handler.CollectStats();

	m_junction_handler.WriteStats(stat_file, "original junction stats");

	string ori_junction_file = m_filtered_alignment_file; ori_junction_file.append(".ori.junc");

	string ori_junction_ins_file = m_filtered_alignment_file; ori_junction_ins_file.append(".ori.junc.ins");

	string ori_junction_del_file = m_filtered_alignment_file; ori_junction_del_file.append(".ori.junc.del");

	string ori_junction_fusion_file = m_filtered_alignment_file; ori_junction_fusion_file.append(".ori.junc.fusion");

	string ori_junction_filter_file = m_filtered_alignment_file; ori_junction_filter_file.append(".ori.junc.filter");

	m_junction_handler.WriteJunction(ori_junction_file, ori_junction_ins_file, ori_junction_del_file, ori_junction_fusion_file, ori_junction_filter_file);

	cout << "mark fitered junction "<< endl;

	//do_filter
	//1 pair
	//2 pairing filter
	//4 filter junc
	//8 select best

	if (m_do_filter & 4)
	{
		m_junction_handler.MarkFiltered(m_is_paired, m_junction_handler_annotated);

		m_junction_handler.ClearStats();

		m_junction_handler.CollectStats();

		m_junction_handler.WriteStats(stat_file, "original filtered junction stats");

		string ori_fil_junction_file = m_filtered_alignment_file; ori_fil_junction_file.append(".ori.fil.junc");

		string ori_fil_junction_ins_file = m_filtered_alignment_file; ori_fil_junction_ins_file.append(".ori.fil.junc.ins");

		string ori_fil_junction_del_file = m_filtered_alignment_file; ori_fil_junction_del_file.append(".ori.fil.junc.del");

		string ori_fil_junction_fusion_file = m_filtered_alignment_file; ori_fil_junction_fusion_file.append(".ori.fil.junc.fusion");

		string ori_fil_junction_filter_file = m_filtered_alignment_file; ori_fil_junction_filter_file.append(".ori.fil.junc.filter");

		m_junction_handler.WriteJunction(ori_fil_junction_file, ori_fil_junction_ins_file, ori_fil_junction_del_file, ori_fil_junction_fusion_file, ori_fil_junction_filter_file);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "filtering: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.LoadFusionJuncToSortVec();

	m_junction_handler.SortFusionJunc();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process prepare: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	////do_filter
	////1 pair
	////2 pairing filter
	////4 filter junc
	////8 select best

	t1=clock();

	m1 = time(NULL);

	if (m_do_filter & 8)//do select best
		m_do_filter = m_do_filter | 16;

	cerr << "third filtering"<< endl;

	for (sam_file_iter = m_sam_files.begin(); sam_file_iter != m_sam_files.end(); ++sam_file_iter)
	{
		ProcessFile(*sam_file_iter, junction_handler_filtered_ptr, junction_handler_ptr, 1);
	}

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "second process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler_filtered.LoadFlankString();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third process: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	WriteStats(stat_file, "filtered alignment stats");

	t1=clock();

	m1 = time(NULL);

	m_junction_handler_filtered.CollectStats();

	m_junction_handler_filtered.WriteStats(stat_file, "filtered junction");

	string fil_junction_file = m_filtered_alignment_file; fil_junction_file.append(".fil.junc");

	string fil_junction_ins_file = m_filtered_alignment_file; fil_junction_ins_file.append(".fil.junc.ins");

	string fil_junction_del_file = m_filtered_alignment_file; fil_junction_del_file.append(".fil.junc.del");

	string fil_junction_fusion_file = m_filtered_alignment_file; fil_junction_fusion_file.append(".fil.junc.fusion");

	string fil_junction_filter_file = m_filtered_alignment_file; fil_junction_filter_file.append(".fil.junc.filter");

	m_junction_handler_filtered.WriteJunction(fil_junction_file, fil_junction_ins_file, fil_junction_del_file, fil_junction_fusion_file, fil_junction_filter_file);
	
	string ori_junction_fusion_encompass_file = m_filtered_alignment_file; ori_junction_fusion_encompass_file.append(".ori.junc.fusion.encompass");

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "third write junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	//fusion encompassing 

	cout << "gene fuion struct" << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.GenerateFusionStruct();

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "generate structure: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;

	cout << "write fusion struct" << endl;

	t1=clock();

	m1 = time(NULL);

	m_junction_handler.WriteFusionJunctionWithEncompassingAlignments(ori_junction_fusion_encompass_file, junction_handler_filtered_ptr);

	t2=clock();

	m2 = time(NULL);

	m_any_time = (t2-t1)/(double)CLOCKS_PER_SEC;

	cerr << "write fusion junction: " << m_any_time << endl;

	cerr << "seconds passed " << m2 - m1 << endl;
}

void 
AlignmentHandler::WriteAlignment(vector<SamRec*>& sam_rec_ptr)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->isunmapped)
		{
			if (!m_ofs_unmapped.is_open())
				cout << "m_ofs_unmapped" <<" is not open"<<endl;

			m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 

			continue;
		}
		else if ((*sam_rec_ptr_iter)->is_fusion)
			m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;
		else
			m_ofs_filtered_alignment << (*sam_rec_ptr_iter)->output_line<<endl;

		if ((*sam_rec_ptr_iter)->is_fusion)
			m_ofs_fusion_std << (*sam_rec_ptr_iter)->output_line<<endl;

		if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
			m_ofs_single << (*sam_rec_ptr_iter)->output_line<<endl; 
		else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
			m_ofs_fusion_paired << (*sam_rec_ptr_iter)->output_line<<endl;
		else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
		{
			m_ofs_paired << (*sam_rec_ptr_iter)->output_line<<endl;
		}
		else
			cerr << (*sam_rec_ptr_iter)->paired_type << endl;
	}

}

void 
AlignmentHandler::WriteAlignmentUnspliced(vector<pair<vector<SamRec*>, vector<SamRec*> > >::iterator pair_samvec_iter, bool is_unspliced/*;vector<SamRec*>& sam_rec_ptr*/)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	if (is_unspliced && ((m_do_filter & 64) == 0))
	{
		if ((m_do_filter & 128) == 0)
		{
			for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end() && is_unspliced; ++sam_rec_ptr_iter)
			{

				if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
					m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_unspliced_single << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_unspliced_fusion_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
			}

			for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end() && is_unspliced; ++sam_rec_ptr_iter)
			{
				if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
					m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_unspliced_single << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_unspliced_fusion_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
			}
		}
		else
		{
			for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end() && is_unspliced; ++sam_rec_ptr_iter)
			{

				if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
					m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
			}

			for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end() && is_unspliced; ++sam_rec_ptr_iter)
			{
				if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
					m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == SINGLE)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl; 
				else if ((*sam_rec_ptr_iter)->paired_type == FUSION_PAIRED)
					m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line<<endl;
				else if ((*sam_rec_ptr_iter)->paired_type == NORMAL_PAIRED)
				{
					m_ofs_unspliced_paired << (*sam_rec_ptr_iter)->output_line<<endl;
				}
				else
					cerr << "not possible"<< endl << (*sam_rec_ptr_iter)->output_line << endl;
			}
		}
	}
	else
	{
		for (sam_rec_ptr_iter = pair_samvec_iter->first.begin(); sam_rec_ptr_iter != pair_samvec_iter->first.end(); ++sam_rec_ptr_iter)
		{
			if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
				m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 
			else
				m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line << endl;
		}

		for (sam_rec_ptr_iter = pair_samvec_iter->second.begin(); sam_rec_ptr_iter != pair_samvec_iter->second.end(); ++sam_rec_ptr_iter)
		{
			if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
				m_ofs_unmapped << (*sam_rec_ptr_iter)->output_line<<endl; 
			else
				m_ofs_onespliced << (*sam_rec_ptr_iter)->output_line << endl;
		}
	}
}

void* AlignmentHandler::SetBitInfoStatic(void * str)
{
	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *((pair<vector<SamRec*>, vector<SamRec*> >*)str);

	for (size_t i = 0; i < sam_rec_pe_ptr.first.size(); ++i)
	{
		sam_rec_pe_ptr.first[i]->strand_t |= IS_PAIRED;

		if (!sam_rec_pe_ptr.first[i]->is_fusion)
		{
			sam_rec_pe_ptr.first[i]->strand_t |= IS_FIRST_END;
		}
		else
			sam_rec_pe_ptr.first[i]->strand_t2 |= IS_PAIRED;
	}

	for (size_t i = 0; i < sam_rec_pe_ptr.second.size(); ++i)
	{
		sam_rec_pe_ptr.second[i]->strand_t |= IS_PAIRED;

		if (!sam_rec_pe_ptr.second[i]->is_fusion)
		{
			sam_rec_pe_ptr.second[i]->strand_t |= IS_SECOND_END;
		}
		else
			sam_rec_pe_ptr.second[i]->strand_t2 |= IS_PAIRED;
	}
	return 0;
}

void* AlignmentHandler::SetFusionBitInfoStatic(void * str)
{
	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *((pair<vector<SamRec*>, vector<SamRec*> >*)str);

	for (size_t i = 0; i < sam_rec_pe_ptr.first.size(); ++i)
	{
		if (sam_rec_pe_ptr.first[i]->isunmapped)
			continue;

		if (sam_rec_pe_ptr.first[i]->is_fusion)
		{
			sam_rec_pe_ptr.first[i]->setfusionbit();
		}
	}

	for (size_t i = 0; i < sam_rec_pe_ptr.second.size(); ++i)
	{
		if (sam_rec_pe_ptr.second[i]->isunmapped)
			continue;

		if (sam_rec_pe_ptr.second[i]->is_fusion)
		{
			sam_rec_pe_ptr.second[i]->setfusionbit();
		}
	}

	return 0;
}

void
AlignmentHandler::WriteAlignment(vector<SamRec*>& sam_rec_ptr, ofstream& cur_ofs)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		
		if ((*sam_rec_ptr_iter)->is_fusion)
			cur_ofs << (*sam_rec_ptr_iter)->tostandfusion(sam_rec_ptr.size(), sam_rec_ptr_iter - sam_rec_ptr.begin() + 1)<<endl;
		else
			cur_ofs << (*sam_rec_ptr_iter)->tostring(sam_rec_ptr.size(), sam_rec_ptr_iter - sam_rec_ptr.begin() + 1)<<endl;
	}
}

void AlignmentHandler::SetPairedType(vector<SamRec*>& sam_rec_ptr, PAIRED_TYPE paired_type)
{
	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->isunmapped)
			continue;

		if ((*sam_rec_ptr_iter)->paired_type == UNMAPPED)
			;
		else if ((*sam_rec_ptr_iter)->is_fusion)
			(*sam_rec_ptr_iter)->paired_type = paired_type;
		else
			(*sam_rec_ptr_iter)->paired_type = paired_type;
	}
}

void* AlignmentHandler::GenerateAlignment(void * str, bool fusion_std)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>* )str);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	if (sam_rec_ptr.empty())
		return 0;

	size_t is_primary_count = 0;

	//cout << "GenerateAlignment 1" << endl;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion)
		{
			if (fusion_std)
			{
				(*sam_rec_ptr_iter)->output_line = (*sam_rec_ptr_iter)->tostandfusion(sam_rec_ptr.size(), sam_rec_ptr_iter - sam_rec_ptr.begin() + 1);
			}
			else
			{
				(*sam_rec_ptr_iter)->output_line = (*sam_rec_ptr_iter)->tostring(0, 0);
			}
		}
		else
		{
			(*sam_rec_ptr_iter)->output_line = (*sam_rec_ptr_iter)->tostring(sam_rec_ptr.size(), sam_rec_ptr_iter - sam_rec_ptr.begin() + 1);
		}

		if (((*sam_rec_ptr_iter)->strand_t & IS_PRIMARY) == 0)
			++is_primary_count;

		if ((*sam_rec_ptr_iter)->is_fusion)
		{
			if (((*sam_rec_ptr_iter)->strand_t2 & IS_PRIMARY) == 0)
				++is_primary_count;
		}

	}

	//cout << "GenerateAlignment 2" << endl;

	if (is_primary_count != 1 && fusion_std)
	{
		cout << "is_primary_count:"<< is_primary_count << endl;

		cout <<sam_rec_ptr.front()->tag_name<<endl;
	}

	return 0;
}

void
AlignmentHandler::FilterCanonNonCanonByReads(vector<SamRec>& read_sam, vector<SamRec* >& filtered_read_sam_ptr)
{
	int count=0, tag_count = 0, unspliced = 0, filtered_canon = 0, filtered_noncanon = 0, filtered_canon_noncanon = 0;

	vector<SamRec>::iterator samrec_iter;

	bool canon = false, noncanon = false, is_insert = false, not_insert = false;

	for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
	{
		if (samrec_iter->canon_count)
			canon =true;

		if (samrec_iter->noncanon_count)
			noncanon =true;

		if (samrec_iter->issmallins)
			is_insert = true;
		else
			not_insert = true;
	}

	int insert_model;

	if (is_insert && not_insert)
		insert_model = 0;
	else if (!is_insert && not_insert)
		insert_model = 1;
	else if (is_insert && !not_insert)
		insert_model = -1;
	else
	{
		insert_model = -2;

		cout << "not possible"<<endl;
	}

	if ((canon && !noncanon) || (!canon && noncanon) || (!canon && !noncanon))
	{
		if (insert_model == 0)
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (samrec_iter->issmallins == false)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
				}					
			}
		}
		else
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				filtered_read_sam_ptr.push_back(&(*samrec_iter));
			}
		}

		++filtered_canon;
	}
	else if (canon && noncanon)
	{
		double max_canon_rate = 0;

		bool noncan_can = false, noncan = false, unspliced = false;

		if (insert_model == 0)
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (max_canon_rate < samrec_iter->canon_rate && samrec_iter->issmallins == false)
					max_canon_rate = samrec_iter->canon_rate;
			}

			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (samrec_iter->issmallins)
				{
				}
				else if (max_canon_rate == samrec_iter->canon_rate)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
				}
				else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
				{
					noncan_can = true;
				}
				else if (samrec_iter->noncanon_count)
				{
					noncan = true;
				}
				else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));

					unspliced = true;
				}

			}
		}
		else
		{
			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{
				if (max_canon_rate < samrec_iter->canon_rate)
					max_canon_rate = samrec_iter->canon_rate;
			}

			for (samrec_iter = read_sam.begin(); samrec_iter != read_sam.end(); ++samrec_iter)
			{

				if (max_canon_rate == samrec_iter->canon_rate)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));
				}
				else if (samrec_iter->canon_count && samrec_iter->noncanon_count)
				{
					noncan_can = true;
				}
				else if (samrec_iter->noncanon_count)
				{
					noncan = true;
				}
				else if (samrec_iter->noncanon_count == 0 && samrec_iter->canon_count == 0)
				{
					filtered_read_sam_ptr.push_back(&(*samrec_iter));

					unspliced = true;
				}
			}
		}
	}
}



void
AlignmentHandler::FilterPairedSamRec(vector<PairedSamRec>& pairedsamrec_vec, pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr)
{
	if (pairedsamrec_vec.size() == 0)
		return;

	vector<PairedSamRec>::iterator psm_iter;

	sort(pairedsamrec_vec.begin(), pairedsamrec_vec.end(), comp_dist);

	vector<PairedSamRec>::iterator vps_iter;

	if (m_do_filter & 16)
	{
		sam_rec_pe_ptr.first.clear();

		sam_rec_pe_ptr.second.clear();
	}

	vps_iter = pairedsamrec_vec.begin();

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

		vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;//static_cast<long> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.first->mate_match = "=";
		else
			vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

		if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
			vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
	}

	if (vps_iter->paired_sam_rec.first->is_fusion)
	{
		vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

		vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

		vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

		vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
	}

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

		vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist; //static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.second->mate_match = "=";
		else
			vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

		if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
			vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

	}

	if (vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

		vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

		vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

		vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
	}

	//output

	if (m_do_filter & 16)
	{
		sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

		sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
	}

	for (vps_iter = pairedsamrec_vec.begin() + 1; vps_iter != pairedsamrec_vec.end(); ++vps_iter)
	{
		if ((vps_iter - 1)->mappedlen == (vps_iter)->mappedlen && !comp_mate_dist(((vps_iter - 1)->mate_dist), ((vps_iter)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (vps_iter - 1)->total_mismatch == (vps_iter)->total_mismatch 
			&& (vps_iter - 1)->total_ave_mismatch == (vps_iter)->total_ave_mismatch
			&& !comp_intron_dist((vps_iter - 1)->intron_size, (vps_iter)->intron_size)
			&& (vps_iter - 1)->contiglen == (vps_iter)->contiglen
			&& (vps_iter - 1)->min_anchor == (vps_iter)->min_anchor)
		{

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

				vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.first->mate_match = "=";
				else
					vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

				if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.first->is_fusion)
			{
				vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

				vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

				vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
			}

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{

				vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

				vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist;//static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.second->mate_match = "=";
				else
					vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

				if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

				vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

				vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
			}

			if (m_do_filter & 16)
			{
				sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

				sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
			}
		}
		else
			break;
	}
}

void* 
AlignmentHandler::FilterPairedSamRecStatic(void * str)
{
	FilterPairingStruct* filter_pairing = (FilterPairingStruct*)str;

	vector<PairedSamRec>& pairedsamrec_vec = *(filter_pairing->pairedsamrec_vec);
	
	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(filter_pairing->sam_rec_pe_ptr);

	int do_filter = filter_pairing->do_filter;

	if (pairedsamrec_vec.size() == 0)
		return 0;

	vector<PairedSamRec>::iterator psm_iter;

	sort(pairedsamrec_vec.begin(), pairedsamrec_vec.end(), comp_dist_0);

	vector<PairedSamRec>::iterator same_vps_iter;

	for (same_vps_iter = pairedsamrec_vec.begin() + 1; same_vps_iter != pairedsamrec_vec.end(); ++same_vps_iter)
	{
		if ((same_vps_iter - 1)->mappedlen == (same_vps_iter)->mappedlen && !comp_mate_dist(((same_vps_iter - 1)->mate_dist), ((same_vps_iter)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (same_vps_iter - 1)->total_mismatch == (same_vps_iter)->total_mismatch 
			&& (same_vps_iter - 1)->total_canon_rate == (same_vps_iter)->total_canon_rate)
		{
		}
		else
			break;
	}

	if (same_vps_iter - pairedsamrec_vec.begin() > 1)
		sort(pairedsamrec_vec.begin(), same_vps_iter, comp_dist_1);

	vector<PairedSamRec>::iterator same_vps_iter2;

	for (same_vps_iter2 = pairedsamrec_vec.begin() + 1; same_vps_iter2 != pairedsamrec_vec.end(); ++same_vps_iter2)
	{
		if ((same_vps_iter2 - 1)->mappedlen == (same_vps_iter2)->mappedlen && !comp_mate_dist(((same_vps_iter2 - 1)->mate_dist), ((same_vps_iter2)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& !comp_intron_dist(((same_vps_iter2 - 1)->intron_size), ((same_vps_iter2)->intron_size))
			&& (same_vps_iter2 - 1)->contiglen == (same_vps_iter2)->contiglen 
			&& (same_vps_iter2 - 1)->total_mismatch == (same_vps_iter2)->total_mismatch 
			&& (same_vps_iter2 - 1)->total_canon_rate == (same_vps_iter2)->total_canon_rate)
		{
		}
		else
			break;
	}

	if (pairedsamrec_vec.begin()->intron_size && same_vps_iter2 - pairedsamrec_vec.begin() > 1)
		sort(pairedsamrec_vec.begin(), same_vps_iter2, comp_dist_2);

	vector<PairedSamRec>::iterator vps_iter;

	if (do_filter & 16)
	{
		sam_rec_pe_ptr.first.clear();

		sam_rec_pe_ptr.second.clear();
	}

	vps_iter = pairedsamrec_vec.begin();

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

		vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;//static_cast<long> (vps_iter->paired_sam_rec.first->mate_offset - vps_iter->paired_sam_rec.first->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.first->mate_match = "=";
		else
			vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

		if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
		vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
	}

	if (vps_iter->paired_sam_rec.first->is_fusion)
	{
		vps_iter->paired_sam_rec.first->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

		vps_iter->paired_sam_rec.first->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

		vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;


		vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

		if (do_filter & 512)
		{
			vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

			vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
		}

		vps_iter->paired_sam_rec.first->fusion_mate_ptr = vps_iter->paired_sam_rec.second;

		vps_iter->paired_sam_rec.first->fusion_mate_on_doner_side = vps_iter->fusion_mate_on_doner_side;
	}

	if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

		vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist; //static_cast<int> (vps_iter->paired_sam_rec.second->mate_offset - vps_iter->paired_sam_rec.second->start);

		if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
			vps_iter->paired_sam_rec.second->mate_match = "=";
		else
			vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

		if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
			vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

		vps_iter->paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

	}

	if (vps_iter->paired_sam_rec.second->is_fusion)
	{
		vps_iter->paired_sam_rec.second->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

		vps_iter->paired_sam_rec.second->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

		vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

		vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

		if (do_filter & 512)
		{
			vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

			vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
		}

		vps_iter->paired_sam_rec.second->fusion_mate_ptr = vps_iter->paired_sam_rec.first;

		vps_iter->paired_sam_rec.second->fusion_mate_on_doner_side = vps_iter->fusion_mate_on_doner_side;
	}

	if (vps_iter->paired_sam_rec.first->is_fusion)
	{
		if (vps_iter->paired_sam_rec.first->fusion_prefix_paired)
			vps_iter->paired_sam_rec.first->strand_t -= IS_PRIMARY;
		else if (vps_iter->paired_sam_rec.first->fusion_suffix_paired)
			vps_iter->paired_sam_rec.first->strand_t2 -= IS_PRIMARY;
	}
	else
		vps_iter->paired_sam_rec.first->strand_t -= IS_PRIMARY;

	if (vps_iter->paired_sam_rec.second->is_fusion)
	{
		if (vps_iter->paired_sam_rec.second->fusion_prefix_paired)
			vps_iter->paired_sam_rec.second->strand_t -= IS_PRIMARY;
		else if (vps_iter->paired_sam_rec.second->fusion_suffix_paired)
			vps_iter->paired_sam_rec.second->strand_t2 -= IS_PRIMARY;
	}
	else
		vps_iter->paired_sam_rec.second->strand_t -= IS_PRIMARY;

	if (do_filter & 16)
	{
		sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

		sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
	}

	for (vps_iter = pairedsamrec_vec.begin() + 1; vps_iter != pairedsamrec_vec.end(); ++vps_iter)
	{
		#ifdef DEBUG

		if ((vps_iter - 1)->paired_sam_rec.first->tag_name.find("R_15755/") != string::npos || (vps_iter - 1)->paired_sam_rec.first->tag_name.find("R_15755/") != string::npos ||
			(vps_iter - 1)->paired_sam_rec.first->tag_name.find("ILLUMINA-2F52BD_0038:2:114:18746:16774#0/") != string::npos || (vps_iter - 1)->paired_sam_rec.first->tag_name.find("seq.10012918/") != string::npos)
		{
			cout << (vps_iter - 1)->paired_sam_rec.first->tostring(1,1) << endl << (vps_iter - 1)->paired_sam_rec.second->tostring(1,1) << endl;

			cout << (vps_iter)->paired_sam_rec.first->tostring(1,1) << endl << (vps_iter)->paired_sam_rec.second->tostring(1,1) << endl;

			cout <<"mappedlen" << endl<< (vps_iter - 1)->mappedlen << endl <<  (vps_iter)->mappedlen << endl;

			cout <<"total_mismatch" << endl<< (vps_iter - 1)->total_mismatch << endl <<  (vps_iter)->total_mismatch << endl;

			cout <<"mate_dist" << endl<< (vps_iter - 1)->mate_dist << endl <<  (vps_iter)->mate_dist << endl;

			cout <<"total_ave_mismatch" << endl<< (vps_iter - 1)->total_ave_mismatch << endl <<  (vps_iter)->total_ave_mismatch << endl;

			cout <<"intron_size" << endl<< (vps_iter - 1)->intron_size << endl <<  (vps_iter)->intron_size << endl;

			cout <<"total_pairing_rate" << endl<< (vps_iter - 1)->total_pairing_rate << endl <<  (vps_iter)->total_pairing_rate << endl;

			cout <<"total_filter_score" << endl<< (vps_iter - 1)->total_filter_score << endl <<  (vps_iter)->total_filter_score << endl;

			cout <<"total_hits" << endl<< (vps_iter - 1)->total_hits << endl <<  (vps_iter)->total_hits << endl;

			cout <<"mappedlen" << endl<< (vps_iter - 1)->total_anchor_len << endl <<  (vps_iter)->total_anchor_len << endl;

			cout <<"min_anchor" << endl<< (vps_iter - 1)->min_anchor << endl <<  (vps_iter)->min_anchor << endl;

			cout <<"contiglen" << endl<< (vps_iter - 1)->contiglen << endl <<  (vps_iter)->contiglen << endl;
			
		}

		#endif

		if ((vps_iter - 1)->mappedlen == (vps_iter)->mappedlen && !comp_mate_dist(((vps_iter - 1)->mate_dist), ((vps_iter)->mate_dist)) //(vps_iter - 1)->mate_dist == (vps_iter)->mate_dist 
			&& (vps_iter - 1)->total_mismatch == (vps_iter)->total_mismatch 
			&& (vps_iter - 1)->intron_size == (vps_iter)->intron_size
			&& (vps_iter - 1)->total_canon_rate == (vps_iter)->total_canon_rate
			&& (vps_iter - 1)->total_hits == (vps_iter)->total_hits
			&& (vps_iter - 1)->total_filter_score == (vps_iter)->total_filter_score
			&& (vps_iter - 1)->contiglen == (vps_iter)->contiglen
			&& (vps_iter - 1)->min_anchor == (vps_iter)->min_anchor)
		{

			if (vps_iter->paired_sam_rec.first->mate_match != "*" || vps_iter->paired_sam_rec.second->mate_match != "*")
				continue;

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.first->mate_offset = vps_iter->paired_sam_rec.second->start;

				vps_iter->paired_sam_rec.first->mate_diff = -vps_iter->outter_dist;

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.first->mate_match = "=";
				else
					vps_iter->paired_sam_rec.first->mate_match = vps_iter->paired_sam_rec.second->chrom_name;

				if (vps_iter->paired_sam_rec.second->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.first->is_fusion)
			{
				vps_iter->paired_sam_rec.first->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

				vps_iter->paired_sam_rec.first->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

				vps_iter->paired_sam_rec.first->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.first->fusion_suffix_len = vps_iter->suffixlen;

				if (do_filter & 512)
				{
					vps_iter->paired_sam_rec.first->left_splice_ways = vps_iter->left_splice_ways;

					vps_iter->paired_sam_rec.first->right_splice_ways = vps_iter->right_splice_ways;
				}

				vps_iter->paired_sam_rec.first->fusion_mate_ptr = vps_iter->paired_sam_rec.second;

				vps_iter->paired_sam_rec.first->fusion_mate_on_doner_side = vps_iter->fusion_mate_on_doner_side;
			}

			if (!vps_iter->paired_sam_rec.first->is_fusion && !vps_iter->paired_sam_rec.second->is_fusion)
			{

				vps_iter->paired_sam_rec.second->mate_offset = vps_iter->paired_sam_rec.first->start;

				vps_iter->paired_sam_rec.second->mate_diff = vps_iter->outter_dist;

				if (vps_iter->paired_sam_rec.first->chrom_name == vps_iter->paired_sam_rec.second->chrom_name)
					vps_iter->paired_sam_rec.second->mate_match = "=";
				else
					vps_iter->paired_sam_rec.second->mate_match = vps_iter->paired_sam_rec.first->chrom_name;

				if (vps_iter->paired_sam_rec.first->strand_t & IS_REVERSE)
					vps_iter->paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

				vps_iter->paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;
			}

			if (vps_iter->paired_sam_rec.second->is_fusion)
			{
				vps_iter->paired_sam_rec.second->fusion_prefix_paired = vps_iter->fusion_prefix_matched;

				vps_iter->paired_sam_rec.second->fusion_suffix_paired = vps_iter->fusion_suffix_matched;

				vps_iter->paired_sam_rec.second->fusion_prefix_len = vps_iter->prefixlen;

				vps_iter->paired_sam_rec.second->fusion_suffix_len = vps_iter->suffixlen;

				if (do_filter & 512)
				{
					vps_iter->paired_sam_rec.second->left_splice_ways = vps_iter->left_splice_ways;

					vps_iter->paired_sam_rec.second->right_splice_ways = vps_iter->right_splice_ways;
				}

				vps_iter->paired_sam_rec.second->fusion_mate_ptr = vps_iter->paired_sam_rec.first;

				vps_iter->paired_sam_rec.second->fusion_mate_on_doner_side = vps_iter->fusion_mate_on_doner_side;
			}

			if (do_filter & 16)
			{
				sam_rec_pe_ptr.first.push_back(vps_iter->paired_sam_rec.first);

				sam_rec_pe_ptr.second.push_back(vps_iter->paired_sam_rec.second);
			}
		}
		else
			break;
	}

	return 0;
}

bool AlignmentHandler::EstablishPairingStatic(void * str)
{
	EstPairingStruct* est_pairing_ptr = (EstPairingStruct*) str;

	vector<PairedSamRec> paired_reads_ptr;

	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(est_pairing_ptr->sam_rec_pe_ptr);

	vector<PairedSamRec>& fusion_paired_reads_ptr = *(est_pairing_ptr->fusion_paired_reads_ptr);

	int do_filter = est_pairing_ptr->do_filter;
	
	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	vector<SamRec*>::iterator sam_rec_ptr_iter1, sam_rec_ptr_iter2;

	for (sam_rec_ptr_iter1 = sam_rec_pe_ptr.first.begin(); sam_rec_ptr_iter1 != sam_rec_pe_ptr.first.end(); ++sam_rec_ptr_iter1)
	{
		if ((*sam_rec_ptr_iter1)->isunmapped)
			continue;

		for (sam_rec_ptr_iter2 = sam_rec_pe_ptr.second.begin(); sam_rec_ptr_iter2 != sam_rec_pe_ptr.second.end(); ++sam_rec_ptr_iter2)
		{
			if ((*sam_rec_ptr_iter2)->isunmapped)
				continue;

			long s_1_st, s_1_end, s_2_st, s_2_end;

			if ((*sam_rec_ptr_iter1)->is_fusion)
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->fusion_prefix_st);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->fusion_suffix_end);
			}
			else
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->start);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->end);
			}


			if ((*sam_rec_ptr_iter2)->is_fusion)
			{
				s_2_st = static_cast <long>( (*sam_rec_ptr_iter2)->fusion_prefix_st);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->fusion_suffix_end);
			}
			else
			{
				s_2_st = static_cast <long>((*sam_rec_ptr_iter2)->start);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->end);
			}

			long mate_dist1 = 0;

			long mate_dist2 = 0;

			bool crossed = false;

			bool not_crossed = false;

			bool is_fusion_paired = false;

			bool fusion_prefix_matched = false, fusion_suffix_matched = false;

			bool fusion_mate_on_doner_side = false;

			size_t strand_1 = 0, strand_2 = 0;

			size_t prefix_len = 0, suffix_len = 0;

			vector<SpliceWay> left_splice_ways;

			vector<SpliceWay> right_splice_ways;

			string chr1 = "", chr2 = "";

			if (!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion &&            //normal paired
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec) && 
				!check_overlap((*sam_rec_ptr_iter2)->spliceway_vec, (*sam_rec_ptr_iter1)->spliceway_vec)
				)
			{
				not_crossed = true;

				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;
			}
			else if (m_do_filter & 16 &&
					!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion &&             //normal alignment fusion paired
					(*sam_rec_ptr_iter1)->chrom_name != (*sam_rec_ptr_iter2)->chrom_name)
			{
				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

			}
			else if ((*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion)   //end 1 fusion, end2 normal
			{
				bool t1 = false, t2 = false, t3 = false, t4 = false;

				bool tt1 = false, tt2 = false, tt3 = false, tt4 = false;

				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->end && abs(long(s_1_st - s_2_end)) < m_max_pair_dist)
				{
					t1 = true;
				}
				
				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist)
				{
					t2 = true;
				}
				
				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->start && abs(long(s_1_end - s_2_st)) < m_max_pair_dist)
				{
					t3 = true;
				}
				
				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist)
				{
					t4 = true;
				}

				if (t1 && t3)
				{
					if ( abs(long(s_1_st - s_2_end)) < abs(long(s_1_end - s_2_st)) )
						tt1 = true;
					else
						tt3 = true;
				}
				else if (t1 && t4)
				{
					if ( abs(long(s_1_st - s_2_end)) < abs(long(s_1_end - s_2_end)) )
						tt1 = true;
					else
						tt4 = true;
				}
				else if (t2 && t3)
				{
					if ( abs(long(s_1_st - s_2_st)) < abs(long(s_1_end - s_2_st)) )
						tt1 = true;
					else
						tt3 = true;
				}
				else if (t2 && t4)
				{
					if ( abs(long(s_1_st - s_2_st)) < abs(long(s_1_end - s_2_end)) )
						tt1 = true;
					else
						tt4 = true;
				}
				else if (t1)
				{
					tt1 = true;
				}
				else if (t2)
				{
					tt2 = true;
				}
				else if (t3)
				{
					tt3 = true;
				}
				else if (t4)
				{
					tt4 = true;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->end && abs(long(s_1_st - s_2_end)) < m_max_pair_dist && tt1)
				{
					fusion_mate_on_doner_side = true;
cout << "case1: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					is_fusion_paired = true; 

					fusion_prefix_matched = true;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1 + abs(mate_dist1) + (*sam_rec_ptr_iter2)->spliceway_vec.back().second;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end;
					
				}

				else if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist && tt2)
				{
					fusion_mate_on_doner_side = true;
cout << "case2: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1 + abs(mate_dist1) + (*sam_rec_ptr_iter2)->spliceway_vec.front().second;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->fusion_suffix_end - (*sam_rec_ptr_iter1)->fusion_prefix_st;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start;
				}

				else if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->start && abs(long(s_1_end - s_2_st)) < m_max_pair_dist && tt3)
				{
					fusion_mate_on_doner_side = false;
cout << "case3: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2 + abs(mate_dist2) + (*sam_rec_ptr_iter2)->spliceway_vec.front().second;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->start2 - (*sam_rec_ptr_iter2)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start2;

				}

				else if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name && (*sam_rec_ptr_iter1)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist && tt4)
				{
					fusion_mate_on_doner_side = false;
cout << "case4: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter1)->mappedlen2 + abs(mate_dist2) + (*sam_rec_ptr_iter2)->spliceway_vec.back().second;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name2), (*sam_rec_ptr_iter1)->start2, &((*sam_rec_ptr_iter1)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->fusion_prefix_st - (*sam_rec_ptr_iter1)->fusion_suffix_end;

					//fragment 2
					(*sam_rec_ptr_iter1)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_FIRST_END;

					(*sam_rec_ptr_iter1)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff2 = (*sam_rec_ptr_iter1)->end2 - (*sam_rec_ptr_iter2)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end2;
				}
			}
			else if (!(*sam_rec_ptr_iter1)->is_fusion && (*sam_rec_ptr_iter2)->is_fusion) //end 1 is normal, end 2 is fusion
			{
				size_t min_mate_dist = -1;

				bool t1 = false, t2 = false, t3 = false, t4 = false;

				bool tt1 = false, tt2 = false, tt3 = false, tt4 = false;

				if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_st)) < m_max_pair_dist)
				{
					t1 = true;
				}
				
				if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist)
				{
					t2 = true;
				}
				
				if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_end)) < m_max_pair_dist)
				{
					t3 = true;
				}
				
				if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist)
				{
					t4 = true;
				}

				if (t1 && t3)
				{
					if ( abs(long(s_1_end - s_2_st)) < abs(long(s_1_st - s_2_end)) )
						tt1 = true;
					else
						tt3 = true;
				}
				else if (t1 && t4)
				{
					if ( abs(long(s_1_end - s_2_st)) < abs(long(s_1_end - s_2_end)) )
						tt1 = true;
					else
						tt4 = true;
				}
				else if (t2 && t3)
				{
					if ( abs(long(s_1_st - s_2_st)) < abs(long(s_1_st - s_2_end)) )
						tt1 = true;
					else
						tt3 = true;
				}
				else if (t2 && t4)
				{
					if ( abs(long(s_1_st - s_2_st)) < abs(long(s_1_end - s_2_end)) )
						tt1 = true;
					else
						tt4 = true;
				}
				else if (t1)
				{
					tt1 = true;
				}
				else if (t2)
				{
					tt2 = true;
				}
				else if (t3)
				{
					tt3 = true;
				}
				else if (t4)
				{
					tt4 = true;
				}

				if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_st)) < m_max_pair_dist && tt1)
				{
					fusion_mate_on_doner_side = true;
cout << "case5: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->spliceway_vec.back().second + abs(mate_dist2) + (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->end - (*sam_rec_ptr_iter1)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;

				}

				else if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_prefix_st < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_st)) < m_max_pair_dist && tt2)
				{
					fusion_mate_on_doner_side = true;
cout << "case6: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					is_fusion_paired = true;

					fusion_prefix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;

					prefix_len = (*sam_rec_ptr_iter1)->spliceway_vec.front().second + abs(mate_dist1) + (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2;

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));


					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->start - (*sam_rec_ptr_iter1)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter2)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter2)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->fusion_suffix_end - (*sam_rec_ptr_iter2)->fusion_prefix_st;
				}

				else if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end < (*sam_rec_ptr_iter1)->start && abs(long(s_1_st - s_2_end)) < m_max_pair_dist && tt3)
				{
					fusion_mate_on_doner_side = false;
cout << "case7: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_end));

					mate_dist2 = (long(s_1_end - s_2_st));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;

					prefix_len = (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2 + (*sam_rec_ptr_iter1)->spliceway_vec.front().second + abs(mate_dist1);

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));


					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->end - (*sam_rec_ptr_iter2)->start2;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->start2 - (*sam_rec_ptr_iter1)->end;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter2)->fusion_suffix_end;
				}

				else if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter2)->fusion_suffix_end > (*sam_rec_ptr_iter1)->end && abs(long(s_1_end - s_2_end)) < m_max_pair_dist && tt4)
				{
					fusion_mate_on_doner_side = false;
cout << "case8: " << (*sam_rec_ptr_iter1)->tag_name << endl;
					not_crossed = true;

					is_fusion_paired = true;

					fusion_suffix_matched = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;

					prefix_len = (*sam_rec_ptr_iter2)->mappedlen1;

					suffix_len = (*sam_rec_ptr_iter2)->mappedlen2 + (*sam_rec_ptr_iter1)->spliceway_vec.back().second + abs(mate_dist2);

					left_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name), (*sam_rec_ptr_iter2)->start, &((*sam_rec_ptr_iter2)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter1)->chrom_name), (*sam_rec_ptr_iter1)->start, &((*sam_rec_ptr_iter1)->spliceway_vec)));

					right_splice_ways.push_back(SpliceWay(&((*sam_rec_ptr_iter2)->chrom_name2), (*sam_rec_ptr_iter2)->start2, &((*sam_rec_ptr_iter2)->spliceway_vec2)));

					//fragment 1
					(*sam_rec_ptr_iter1)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter1)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter1)->strand_t |= IS_FIRST_END;

					if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name2)
						(*sam_rec_ptr_iter1)->mate_match = "=";
					else
						(*sam_rec_ptr_iter1)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter1)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter1)->mate_diff = (*sam_rec_ptr_iter1)->start - (*sam_rec_ptr_iter2)->end2;

					//fragment 2
					(*sam_rec_ptr_iter2)->strand_t2 |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter1)->strand_t & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t2 |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t2 |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter1)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match2 = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match2 = (*sam_rec_ptr_iter1)->chrom_name;

					(*sam_rec_ptr_iter2)->mate_offset2 = (*sam_rec_ptr_iter1)->start;

					(*sam_rec_ptr_iter2)->mate_diff2 = (*sam_rec_ptr_iter2)->end2 - (*sam_rec_ptr_iter1)->start;

					//fragment 3
					(*sam_rec_ptr_iter2)->strand_t |= IS_PAIRED_MAPPED;

					if ((*sam_rec_ptr_iter2)->strand_t2 & IS_REVERSE)
						(*sam_rec_ptr_iter2)->strand_t |= IS_MATE_REVERSE;

					(*sam_rec_ptr_iter2)->strand_t |= IS_FIRST_END;

					(*sam_rec_ptr_iter2)->strand_t |= IS_SECOND_END;

					if ((*sam_rec_ptr_iter2)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name)
						(*sam_rec_ptr_iter2)->mate_match = "=";
					else
						(*sam_rec_ptr_iter2)->mate_match = (*sam_rec_ptr_iter2)->chrom_name2;

					(*sam_rec_ptr_iter2)->mate_offset = (*sam_rec_ptr_iter2)->start2;

					(*sam_rec_ptr_iter2)->mate_diff = (*sam_rec_ptr_iter2)->fusion_prefix_st - (*sam_rec_ptr_iter2)->fusion_suffix_end;
				}
			}
			else if ((*sam_rec_ptr_iter1)->is_fusion && (*sam_rec_ptr_iter2)->is_fusion)
			{
				continue;        //no fusion of both ends allowed at this time
				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
					(*sam_rec_ptr_iter1)->strand1 == '+' && (*sam_rec_ptr_iter2)->strand1 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st > (*sam_rec_ptr_iter2)->fusion_prefix_st)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
					(*sam_rec_ptr_iter1)->strand1 == '-' && (*sam_rec_ptr_iter2)->strand1 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_prefix_st < (*sam_rec_ptr_iter2)->fusion_prefix_st)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
					(*sam_rec_ptr_iter1)->strand2 == '+' && (*sam_rec_ptr_iter2)->strand2 == '-' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end < (*sam_rec_ptr_iter2)->fusion_suffix_end)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;
				}

				if ((*sam_rec_ptr_iter1)->chrom_name2 == (*sam_rec_ptr_iter2)->chrom_name2 && 
					(*sam_rec_ptr_iter1)->strand2 == '-' && (*sam_rec_ptr_iter2)->strand2 == '+' &&
					(*sam_rec_ptr_iter1)->fusion_suffix_end > (*sam_rec_ptr_iter2)->fusion_suffix_end)
				{
					not_crossed = true;

					mate_dist1 = (long(s_1_st - s_2_st));

					mate_dist2 = (long(s_1_end - s_2_end));

					strand_1 = (*sam_rec_ptr_iter1)->strand_t2;
				
					strand_2 = (*sam_rec_ptr_iter2)->strand_t2;

					chr1 = (*sam_rec_ptr_iter1)->chrom_name2;

					chr2 = (*sam_rec_ptr_iter2)->chrom_name2;
				}
			}

			if (!not_crossed)
				crossed = true;

			if (((((strand_1 ^ strand_2) & IS_REVERSE) || is_fusion_paired)&&
				!crossed &&
				(abs(mate_dist1) < m_max_pair_dist ||
				 abs(mate_dist2) < m_max_pair_dist)) || 
				 (m_do_filter & 2048 &&
				  (((strand_1 ^ strand_2) & IS_REVERSE) == 0 && !is_fusion_paired) &&
					!crossed &&
					(abs(mate_dist1) < m_min_mate_dist || abs(mate_dist2) < m_min_mate_dist)))
			{
				paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
					(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
					chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, 
					fusion_prefix_matched, fusion_suffix_matched, fusion_mate_on_doner_side));
				paired = true;
			}

		}
	}

	FilterPairingStruct filter_pairing(&paired_reads_ptr, &sam_rec_pe_ptr, do_filter);

	FilterPairedSamRecStatic( (void *) (&filter_pairing));

	return paired;

}

bool AlignmentHandler::EstablishFusionPairingStatic(void * str)
{
	EstPairingStruct* est_pairing_ptr = (EstPairingStruct*) str;

	vector<PairedSamRec> paired_reads_ptr;

	pair<vector<SamRec*>, vector<SamRec*> >& sam_rec_pe_ptr = *(est_pairing_ptr->sam_rec_pe_ptr);

	vector<PairedSamRec>& fusion_paired_reads_ptr = *(est_pairing_ptr->fusion_paired_reads_ptr);

	int do_filter = est_pairing_ptr->do_filter;
	
	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	vector<SamRec*>::iterator sam_rec_ptr_iter1, sam_rec_ptr_iter2;

	for (sam_rec_ptr_iter1 = sam_rec_pe_ptr.first.begin(); sam_rec_ptr_iter1 != sam_rec_pe_ptr.first.end(); ++sam_rec_ptr_iter1)
	{
		if ((*sam_rec_ptr_iter1)->isunmapped)
			continue;

		for (sam_rec_ptr_iter2 = sam_rec_pe_ptr.second.begin(); sam_rec_ptr_iter2 != sam_rec_pe_ptr.second.end(); ++sam_rec_ptr_iter2)
		{
			if ((*sam_rec_ptr_iter2)->isunmapped)
				continue;

			long s_1_st, s_1_end, s_2_st, s_2_end;

			if ((*sam_rec_ptr_iter1)->is_fusion)
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->fusion_prefix_st);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->fusion_suffix_end);
			}
			else
			{
				s_1_st = static_cast <long>((*sam_rec_ptr_iter1)->start);

				s_1_end = static_cast <long>((*sam_rec_ptr_iter1)->end);
			}


			if ((*sam_rec_ptr_iter2)->is_fusion)
			{
				s_2_st = static_cast <long>( (*sam_rec_ptr_iter2)->fusion_prefix_st);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->fusion_suffix_end);
			}
			else
			{
				s_2_st = static_cast <long>((*sam_rec_ptr_iter2)->start);

				s_2_end = static_cast <long>((*sam_rec_ptr_iter2)->end);
			}

			long mate_dist1 = 0; //= (long(s_1_st - s_2_end));

			long mate_dist2 = 0; //= (long(s_1_end - s_2_st));

			bool crossed = false;

			bool not_crossed = false;

			bool is_fusion_paired = false;

			bool fusion_prefix_matched = false, fusion_suffix_matched = false;

			size_t strand_1 = 0, strand_2 = 0;

			size_t prefix_len = 0, suffix_len = 0;

			vector<SpliceWay> left_splice_ways;

			vector<SpliceWay> right_splice_ways;

			string chr1 = "", chr2 = "";

			if (!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
				(*sam_rec_ptr_iter1)->chrom_name == (*sam_rec_ptr_iter2)->chrom_name && 
				!check_overlap((*sam_rec_ptr_iter1)->spliceway_vec, (*sam_rec_ptr_iter2)->spliceway_vec) && 
				!check_overlap((*sam_rec_ptr_iter2)->spliceway_vec, (*sam_rec_ptr_iter1)->spliceway_vec))
			{

				not_crossed = true;

				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				if (m_do_filter & 16 && (((strand_1 ^ strand_2) & IS_REVERSE) == 0) ||
					(abs(mate_dist1) >= m_max_pair_dist && abs(mate_dist2) >= m_max_pair_dist))
				{

					fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
						(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
						chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false, false));
				}
			}
			else if (m_do_filter & 16 &&
					!(*sam_rec_ptr_iter1)->is_fusion && !(*sam_rec_ptr_iter2)->is_fusion && 
					(*sam_rec_ptr_iter1)->chrom_name != (*sam_rec_ptr_iter2)->chrom_name)
			{
				mate_dist1 = (long(s_1_st - s_2_end));

				mate_dist2 = (long(s_1_end - s_2_st));

				strand_1 = (*sam_rec_ptr_iter1)->strand_t;
				
				strand_2 = (*sam_rec_ptr_iter2)->strand_t;

				chr1 = (*sam_rec_ptr_iter1)->chrom_name;

				chr2 = (*sam_rec_ptr_iter2)->chrom_name;

				fusion_paired_reads_ptr.push_back(PairedSamRec(mate_dist1, mate_dist2, (*sam_rec_ptr_iter1)->mis_match + (*sam_rec_ptr_iter2)->mis_match, 
					(*sam_rec_ptr_iter1)->intron_size + (*sam_rec_ptr_iter2)->intron_size, (*sam_rec_ptr_iter1)->mappedlen + (*sam_rec_ptr_iter2)->mappedlen, 
					chr1, chr2, strand_1, strand_2, (*sam_rec_ptr_iter1), (*sam_rec_ptr_iter2), prefix_len, suffix_len, left_splice_ways, right_splice_ways, false, false, false));
			}
			

		}
	}

	return paired;

}

void 
AlignmentHandler::FilterByMisMatch(vector<SamRec*>& sam_rec_ptr)
{
	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_mis);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->mis_match != (*(sam_rec_ptr_iter-1))->mis_match)
		{
			sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());
			break;
		}
	}
}

void* AlignmentHandler::FilterByMisMatchStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_mis);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->mis_match != (*(sam_rec_ptr_iter-1))->mis_match)
		{
			sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());
			break;
		}
	}
	return 0;
}


void* AlignmentHandler::RemoveDupStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	vector<SamRec*> filtered_sam_rec_ptr;

	//cout << "sort"<<endl;

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_dup);

	vector<SamRec*>::iterator sam_rec_iter;

	filtered_sam_rec_ptr.push_back(*sam_rec_ptr.begin());

	//cout << "remove"<<endl;

	for (sam_rec_iter = sam_rec_ptr.begin() + 1; sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->chrom_name != (*(sam_rec_iter - 1))->chrom_name ||
			(*sam_rec_iter)->start != (*(sam_rec_iter - 1))->start || 
			(*sam_rec_iter)->ori_splice_way != (*(sam_rec_iter - 1))->ori_splice_way ||
			(*sam_rec_iter)->chrom_name2 != (*(sam_rec_iter - 1))->chrom_name2 ||
			(*sam_rec_iter)->start2 != (*(sam_rec_iter - 1))->start2 || 
			(*sam_rec_iter)->ori_splice_way2 != (*(sam_rec_iter - 1))->ori_splice_way2)
			filtered_sam_rec_ptr.push_back(*sam_rec_iter);
		else
		{
		}
	}

	sam_rec_ptr = filtered_sam_rec_ptr;

	//cout << "finish"<<endl;

	return 0;
}

void 
AlignmentHandler::AllocateHitsMemory()
{
	ifstream ifs(m_chrom_size_file);

	if (ifs.is_open())
	{
		vector<bool> dumy;

		while (!ifs.eof() )
		{
			string line;

			getline(ifs, line);

			if (line.empty())
				continue;

			char chrom_name[1000];

			size_t chrom_size;

			sscanf(line.c_str(), "%s\t%llu", chrom_name, &chrom_size);

			m_mapped_pos.insert(make_pair(chrom_name, dumy));

			m_mapped_pos[chrom_name].resize(chrom_size, false);
		}
	}
	else
	{
		cerr <<"Can't open file: " << m_chrom_size_file << endl;
		exit(0);
	}

	ifs.close();
}

void
AlignmentHandler::Hits2ExpressedRegions()
{
	hash_map<string, vector<bool> >::iterator mapped_pos_iter;
	
	for (mapped_pos_iter = m_mapped_pos.begin(); mapped_pos_iter != m_mapped_pos.end(); ++mapped_pos_iter)
	{
		vector<pair<size_t, size_t> >& cur_islands = m_expressed_regions[mapped_pos_iter->first];

		vector<bool>& cur_hits = mapped_pos_iter->second;

		bool inisland = false;

		size_t tmpstart = 0, tmpend = 0;

		for (size_t i = 0; i < cur_hits.size(); ++i)
		{
			if (cur_hits[i] && !inisland)
			{
				inisland = true;
				tmpstart = i;
			}
			else if(!cur_hits[i] && inisland)
			{
				tmpend = i;

				inisland = false;

				if (cur_islands.size() && (tmpstart - cur_islands[cur_islands.size() - 1].second) <= boundary)
					cur_islands[cur_islands.size() - 1].second = tmpend/* + 45*/;
				else
				{
					cur_islands.push_back(make_pair(tmpstart/* - 45*/, tmpend/* + 45*/));
				}
			}
		}

		m_disjointset_map[mapped_pos_iter->first].Resize(cur_islands.size());
	}

	string islands_file = m_filtered_alignment_file; islands_file.append(".islands");

	ofstream test_islands(islands_file.c_str());

	hash_map<string, vector<pair<size_t, size_t> > >::iterator chrom_islands_iter;

	for (chrom_islands_iter = m_expressed_regions.begin(); chrom_islands_iter != m_expressed_regions.end(); ++chrom_islands_iter)
	{
		vector<pair<size_t, size_t> >::iterator island_iter;

		for (island_iter = chrom_islands_iter->second.begin(); island_iter != chrom_islands_iter->second.end(); ++island_iter)
		{
			test_islands << chrom_islands_iter->first << '\t' << island_iter->first << '\t' << island_iter->second << '\t' <<island_iter->second - island_iter->first<< endl;
		}
	}
}

bool
AlignmentHandler::FindRegion(pair<size_t, size_t>& cur_region, vector<pair<size_t, size_t> >& sorted_regions, size_t& find_region_idx)
{
	vector<pair<size_t, size_t> >::iterator express_region_iter;

	express_region_iter = lower_bound(sorted_regions.begin(), sorted_regions.end(), cur_region, compare_pair_region);

	if (express_region_iter != sorted_regions.begin())
		--express_region_iter;

	bool find_region = false;

	vector<pair<size_t, size_t> >::iterator cur_express_region_iter = express_region_iter;

	while (true)
	{
		if (cur_express_region_iter == sorted_regions.end() || cur_express_region_iter->first > cur_region.first)
		{
			break;
		}
		else if (cur_express_region_iter->first  <= cur_region.first + 3 && cur_express_region_iter->second + 3 >= cur_region.first)
		{
			find_region = true;
			break;
		}

		++cur_express_region_iter;
	}

	if (find_region)
		find_region_idx = cur_express_region_iter - sorted_regions.begin();

	return find_region;
}

void
AlignmentHandler::UnionSets(JunctionHandler* junc_handler, hash_map<string, vector<pair<size_t, size_t> > >* expressed_regions)
{
	junc_handler->LoadJuncToSortVec();

	junc_handler->SortJunc();

	hash_map<string, vector<pair<size_t, size_t> > >::iterator expressed_regions_chrom_iter;

	hash_map<size_t, UnionExpressedRegions > dumy;

	for (expressed_regions_chrom_iter = expressed_regions->begin(); expressed_regions_chrom_iter != expressed_regions->end(); ++expressed_regions_chrom_iter)
	{
		if (junc_handler->m_junc_sort.find(expressed_regions_chrom_iter->first) == junc_handler->m_junc_sort.end())
			continue;

		vector<JunctionSeed*>& sorted_juncs = (junc_handler->m_junc_sort.find(expressed_regions_chrom_iter->first))->second;

		vector<pair<size_t, size_t> >& sorted_regions = expressed_regions_chrom_iter->second;

		DisjointSet& cur_disjointSet = (m_disjointset_map.find(expressed_regions_chrom_iter->first))->second;

		vector<JunctionSeed*>::iterator junc_seed_iter;

		for (junc_seed_iter = sorted_juncs.begin(); junc_seed_iter != sorted_juncs.end(); ++junc_seed_iter)
		{
			if ((*junc_seed_iter)->m_filtered_type != NOT_FILTERED)
				continue;

			pair<size_t, size_t> cur_junc_region_start = make_pair((*junc_seed_iter)->m_start, (*junc_seed_iter)->m_start);

			pair<size_t, size_t> cur_junc_region_end = make_pair((*junc_seed_iter)->m_end, (*junc_seed_iter)->m_end);

			size_t start_region_idx = -1, end_region_idx = -1;

			bool find_start_region = FindRegion(cur_junc_region_start, sorted_regions, start_region_idx);

			bool find_end_region = FindRegion(cur_junc_region_end, sorted_regions, end_region_idx);

			if (find_start_region && find_end_region && start_region_idx != end_region_idx)
			{
				cur_disjointSet.Union(start_region_idx, end_region_idx);
			}
			else if (start_region_idx != end_region_idx)
			{
				cout << "Union sets, region not found: "<<expressed_regions_chrom_iter->first << '\t' << (*junc_seed_iter)->m_start << '\t' << (*junc_seed_iter)->m_end<< endl;//cout error information
				continue;
			}
		}

		m_unioned_expressed_regions_map.insert(make_pair(expressed_regions_chrom_iter->first, dumy));

		hash_map<size_t, UnionExpressedRegions >& cur_unioned_expressed_regions = m_unioned_expressed_regions_map.find(expressed_regions_chrom_iter->first)->second;

		//put each expressed regions into tag contig

		for (size_t i = 0; i < sorted_regions.size(); ++i)
		{
			size_t source_id = cur_disjointSet.FindSet(i);

			cur_unioned_expressed_regions[source_id].m_unioned_regions_ids.push_back(i);

			cur_unioned_expressed_regions[source_id].m_unioned_regions.push_back(sorted_regions[i]);

			cur_unioned_expressed_regions[source_id].m_unioned_region_length += sorted_regions[i].second - sorted_regions[i].first + 1;
		}
	}

	string expressed_regions_file = m_filtered_alignment_file; expressed_regions_file.append(".expressed_regions");

	ofstream expressed_regions_ofs(expressed_regions_file.c_str());

	hash_map<string, hash_map<size_t, UnionExpressedRegions > >::iterator chrom_expressed_regions_iter;

	hash_map<size_t, UnionExpressedRegions >::iterator source_expressed_regions_iter;

	for (chrom_expressed_regions_iter = m_unioned_expressed_regions_map.begin(); chrom_expressed_regions_iter != m_unioned_expressed_regions_map.end(); ++chrom_expressed_regions_iter)
	{
		for (source_expressed_regions_iter = chrom_expressed_regions_iter->second.begin(); source_expressed_regions_iter != chrom_expressed_regions_iter->second.end(); ++source_expressed_regions_iter)
		{
			expressed_regions_ofs << chrom_expressed_regions_iter->first << '\t' << source_expressed_regions_iter->first << '\t' << source_expressed_regions_iter->second.m_unioned_region_length << '\t';
			
			for (size_t i = 0; i < source_expressed_regions_iter->second.m_unioned_regions.size(); ++i)
				expressed_regions_ofs<<source_expressed_regions_iter->second.m_unioned_regions[i].first<<'\t'
				<<source_expressed_regions_iter->second.m_unioned_regions[i].second<< '\t';

			expressed_regions_ofs << endl;
		}
	}
}

void
AlignmentHandler::FilterFusionByIsolatedExons(JunctionHandler* junc_handler)
{
	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	JUNC_HASH_COMB::iterator offset_iter;

	hash_map<string, hash_map<size_t, int> > doner_not_filtered;

	hash_map<string, hash_map<size_t, int> > acceptor_not_filtered;

	for (chrom_iter = junc_handler->m_junc_hash.begin(); chrom_iter != junc_handler->m_junc_hash.end(); ++chrom_iter)
	{
		vector<pair<size_t, size_t> >* sorted_regions_ptr1;
		
		vector<pair<size_t, size_t> >* sorted_regions_ptr2;

		DisjointSet* cur_disjointSet1;
		
		DisjointSet* cur_disjointSet2;

		if (chrom_iter->first.find("~") != string::npos)
		{
			char chr1[1000], chr2[1000];

			sscanf(chrom_iter->first.c_str(), "%[^~]~%[^~]", chr1, chr2);

			if (m_expressed_regions.find(chr1) == m_expressed_regions.end())
			{
				cout << "m_expressed_region 1 not found:"<< chr1<<endl;

				continue;
			}

			sorted_regions_ptr1 = &(m_expressed_regions[chr1]);

			if (m_expressed_regions.find(chr2) == m_expressed_regions.end())
			{
				cout << "m_expressed_region 2 not found:"<< chr2<<endl;

				continue;
			}

			sorted_regions_ptr2 = &(m_expressed_regions[chr2]);


			//find corresponding disjoint set
			if (m_disjointset_map.find(chr1) == m_disjointset_map.end())
			{
				cout << "m_disjointset_map 1 not found:"<< chr1<<endl;

				continue;
			}

			cur_disjointSet1 = &(m_disjointset_map[chr1]);

			if (m_disjointset_map.find(chr2) == m_disjointset_map.end())
			{
				cout << "m_disjointset_map 2 not found:"<< chr2<<endl;

				continue;
			}

			cur_disjointSet2 = &(m_disjointset_map[chr2]);

			for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
			{
				if (offset_iter->second.m_is_fusion && offset_iter->second.m_filtered_type == NOT_FILTERED)
				{
					pair<size_t, size_t> cur_junc_region_start = make_pair(offset_iter->second.m_start, offset_iter->second.m_start);

					pair<size_t, size_t> cur_junc_region_end = make_pair(offset_iter->second.m_end, offset_iter->second.m_end);

					size_t start_region_idx = -1, end_region_idx = -1;

					bool find_start_region = FindRegion(cur_junc_region_start, *sorted_regions_ptr1, start_region_idx);

					bool find_end_region = FindRegion(cur_junc_region_end, *sorted_regions_ptr2, end_region_idx);
					
					if (find_start_region)
					{
						size_t parent_idx = (*cur_disjointSet1).FindSet(start_region_idx);

						int rank = (*cur_disjointSet1).Rank(parent_idx);

						size_t exon_length = (*sorted_regions_ptr1)[start_region_idx].second - (*sorted_regions_ptr1)[start_region_idx].first + 1;

						if (rank == 0)
						{
							offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_ISOLATED_EXON;

							if ((offset_iter->second.m_start == 59934577 || offset_iter->second.m_start == 60053472) &&
								(offset_iter->second.m_end == 59934577 || offset_iter->second.m_end == 60053472))
							{
								cerr << (*sorted_regions_ptr1)[start_region_idx].first<<'\t'<< (*sorted_regions_ptr1)[start_region_idx].second<<endl;
							}

							continue;
						}
					}
					else
					{
						cout << "Find fusion junction correpsonding exon, doner not found: "<<chrom_iter->first << '\t' << offset_iter->second.m_start << '\t' << offset_iter->second.m_end<< endl;//cout error information

						offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_NO_EXON;

						continue;
					}

					if (find_end_region)
					{
						size_t parent_idx = (*cur_disjointSet2).FindSet(end_region_idx);

						int rank = (*cur_disjointSet2).Rank(parent_idx);

						size_t exon_length = (*sorted_regions_ptr2)[end_region_idx].second - (*sorted_regions_ptr2)[end_region_idx].first + 1;

						if (rank == 0)
						{
							offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_ISOLATED_EXON;

							if ((offset_iter->second.m_start == 59934577 || offset_iter->second.m_start == 60053472) &&
								(offset_iter->second.m_end == 59934577 || offset_iter->second.m_end == 60053472))
							{
								cerr << (*sorted_regions_ptr1)[start_region_idx].first<<'\t'<< (*sorted_regions_ptr1)[start_region_idx].second<<endl;
							}

							continue;
						}
					}
					else
					{
						cout << "Find fusion junction correpsonding exon, acceptor not found: "<<chrom_iter->first << '\t' << offset_iter->second.m_start << '\t' << offset_iter->second.m_end<< endl;//cout error information

						offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_NO_EXON;

						continue;
					}
				}
			}
		}
		else
			continue;		
	}
}

bool comp_normal_pair(const pair<int, int>& lhs, const pair<int, int>& rhs)
{
	return lhs.first < rhs.first;
}

bool comp_normal_pair_intron(const pair<int, int>& lhs, const pair<int, int>& rhs)
{
	return lhs.second + rhs.first > rhs.second + lhs.first ;
}

bool comp_junc_seed_intron(const JunctionSeed* lhs, const JunctionSeed* rhs)
{
	return lhs->m_intronlen > rhs->m_intronlen;
}

void AlignmentHandler::EstimateFragmentLength(JunctionHandler* junc_handler)
{

	hash_map<string, vector<pair<int, int> > >::iterator chrom_iter;

	string ori_junc_file = m_filtered_alignment_file; ori_junc_file.append(".ori_junc.simp");

	ofstream ofs_ori_junc(ori_junc_file.c_str());

	for (chrom_iter = m_original_junc_locs.begin(); chrom_iter != m_original_junc_locs.end(); ++chrom_iter)
	{
		vector<pair<int, int> >::iterator junc_iter;

		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_normal_pair);

		for (junc_iter = chrom_iter->second.begin(); junc_iter != chrom_iter->second.end(); ++junc_iter)
		{
			ofs_ori_junc << chrom_iter->first << '\t' << junc_iter->first <<'\t' << junc_iter->second << endl;
		}
	}

	size_t sum_fragment_length = 0, fragment_count = 0, more_than_one_junc = 0;

	string more_than_one_junc_file = m_filtered_alignment_file; more_than_one_junc_file.append(".more_than_one_junc");

	ofstream ofs_junc(more_than_one_junc_file.c_str());

	string normal_pair_file = m_filtered_alignment_file; normal_pair_file.append(".normal_pair");

	ofstream ofs_normal_pair(normal_pair_file.c_str());

	for (chrom_iter = m_normal_paired_locs.begin(); chrom_iter != m_normal_paired_locs.end(); ++chrom_iter)
	{
		vector<pair<int, int> >& cur_normal_pairs = m_normal_paired_locs[chrom_iter->first];

		sort(cur_normal_pairs.begin(), cur_normal_pairs.end(), comp_normal_pair);

		if (m_original_junc_locs.find(chrom_iter->first) != m_original_junc_locs.end())
		{
			vector<pair<int, int> >& cur_sorted_juncs = m_original_junc_locs[chrom_iter->first];

			vector<pair<int, int> >::iterator normal_pair_iter;

			for (normal_pair_iter = cur_normal_pairs.begin(); normal_pair_iter != cur_normal_pairs.end(); ++normal_pair_iter)
			{
				ofs_normal_pair<< chrom_iter->first << '\t' << normal_pair_iter->first << '\t' << normal_pair_iter->second << endl;

				pair<int, int> cur_junc_seed(normal_pair_iter->first - m_max_read_width, normal_pair_iter->second + m_max_read_width);

				vector<pair<int, int> >::iterator lower_bound_iter;

				lower_bound_iter = lower_bound(cur_sorted_juncs.begin(), cur_sorted_juncs.end(), cur_junc_seed, comp_normal_pair);

				#ifdef DEBUG

				cout << chrom_iter->first << '\t' << normal_pair_iter->first << '\t' << normal_pair_iter->second << endl;

				#endif

				vector<pair<int, int> >::iterator junction_iter = lower_bound_iter;

				while(junction_iter != cur_sorted_juncs.end() && cur_junc_seed.second > (junction_iter)->first)
					++junction_iter;

				vector<pair<int, int> >::iterator junc_iter;

				vector<pair<int, int> > inner_juncs;

				for (junc_iter = lower_bound_iter; junc_iter != junction_iter; ++junc_iter)
				{
					if (cur_junc_seed.first < (junc_iter)->first && cur_junc_seed.second > (junc_iter)->second)
						inner_juncs.push_back(*junc_iter);
				}

				#ifdef DEBUG

				if (junction_iter != cur_sorted_juncs.end())
					cout << chrom_iter->first<< '\t'<<junction_iter->first <<'\t' << junction_iter->second << endl;

				#endif

				if (inner_juncs.size())
				{
					#ifdef DEBUG

					cout << "branch 1"<< endl;

					cout << normal_pair_iter->second << '\t' <<  normal_pair_iter->first<< endl;

					#endif

					sort(inner_juncs.begin(), inner_juncs.end(), comp_normal_pair_intron);

					#ifdef DEBUG

					cout << inner_juncs.front().second - inner_juncs.front().first << endl;

					cout << m_max_read_width<< endl;

					#endif

					normal_pair_iter->second = normal_pair_iter->second - normal_pair_iter->first - 1 - (inner_juncs.front().second - inner_juncs.front().first) + m_max_read_width + m_max_read_width;

					sum_fragment_length += normal_pair_iter->second;

					++fragment_count;

					if (inner_juncs.size() > 1)
					{
						++more_than_one_junc;

						ofs_junc << normal_pair_iter->first << '\t' << normal_pair_iter->second << endl;

						ofs_junc << inner_juncs.size() << " junctions"<<endl;

						for (size_t i = 0; i < inner_juncs.size(); ++i)
						{
							ofs_junc <<chrom_iter->first << '\t' << inner_juncs[i].first <<'\t' <<inner_juncs[i].second <<endl; 
						}
					}
				}
				else
				{
					#ifdef DEBUG

					cout << "branch 2"<< endl;

					cout << normal_pair_iter->second << '\t' <<  normal_pair_iter->first<< endl;
					cout << m_max_read_width<< endl;

					#endif

					normal_pair_iter->second = normal_pair_iter->second - normal_pair_iter->first - 1 + m_max_read_width + m_max_read_width;

					sum_fragment_length += normal_pair_iter->second;

					++fragment_count;
				}
			}
		}
		else
		{
			vector<pair<int, int> >::iterator normal_pair_iter;

			for (normal_pair_iter = cur_normal_pairs.begin(); normal_pair_iter != cur_normal_pairs.end(); ++normal_pair_iter)
			{
				#ifdef DEBUG

				cout <<"branch 3"<<endl;

				#endif

				ofs_normal_pair<< chrom_iter->first << '\t' << normal_pair_iter->first << '\t' << normal_pair_iter->second << endl;

				#ifdef DEBUG

				cout << normal_pair_iter->second << '\t' <<  normal_pair_iter->first<< endl;

				cout << m_max_read_width<< endl;

				#endif

				normal_pair_iter->second = normal_pair_iter->second - normal_pair_iter->first - 1 + m_max_read_width + m_max_read_width;

				sum_fragment_length += normal_pair_iter->second;

				++fragment_count;
			}
			
		}
	}

	double average_fragment = (double) sum_fragment_length / (double) fragment_count;

	double variance_sum = 0;

	string normal_fragment_lengths_file = m_filtered_alignment_file; normal_fragment_lengths_file.append(".fragments");

	ofstream ofs(normal_fragment_lengths_file.c_str());

	for (chrom_iter = m_normal_paired_locs.begin(); chrom_iter != m_normal_paired_locs.end(); ++chrom_iter)
	{
		vector<pair<int, int> >& cur_normal_pairs = m_normal_paired_locs[chrom_iter->first];

		vector<pair<int, int> >::iterator normal_pair_iter;

		for (normal_pair_iter = cur_normal_pairs.begin(); normal_pair_iter != cur_normal_pairs.end(); ++normal_pair_iter)
		{
			ofs << normal_pair_iter->second << endl;

			double diff = (double) normal_pair_iter->second - average_fragment;

			variance_sum += diff * diff;
		}
	}

	double variance = variance_sum / (double) fragment_count;

	cerr << "average_fragment:" << average_fragment << endl;

	cerr << "variance:"<<variance << endl;
}

void*
AlignmentHandler::FindFusionJuncRegionVecStatic(void * str)
{
	JunctionHandler* junc_handler = m_junction_handler_prev_ptr;
	
	vector<PairedSamRec>& fusion_paired_reads_ptr = *((vector<PairedSamRec>*)str);

	vector<PairedSamRec>::iterator paired_sam_rec_iter;

	for (paired_sam_rec_iter = fusion_paired_reads_ptr.begin(); paired_sam_rec_iter != fusion_paired_reads_ptr.end(); ++paired_sam_rec_iter)
	{
		FusionJuncRegion cur_fusion_junc_region;
		
		string chrom_name;

		//not reverse alignment 3
		if (paired_sam_rec_iter->paired_sam_rec.first->chrom_name < paired_sam_rec_iter->paired_sam_rec.second->chrom_name ||
		    (paired_sam_rec_iter->paired_sam_rec.first->chrom_name == paired_sam_rec_iter->paired_sam_rec.second->chrom_name &&
			 paired_sam_rec_iter->paired_sam_rec.first->start <= paired_sam_rec_iter->paired_sam_rec.second->start)
			 )
		{
			chrom_name = paired_sam_rec_iter->paired_sam_rec.first->chrom_name + "~" + paired_sam_rec_iter->paired_sam_rec.second->chrom_name;

			cur_fusion_junc_region.set(paired_sam_rec_iter->paired_sam_rec.first->start, paired_sam_rec_iter->paired_sam_rec.first->end,
				paired_sam_rec_iter->paired_sam_rec.second->start, paired_sam_rec_iter->paired_sam_rec.second->end);
		}
		else
		{
			chrom_name = paired_sam_rec_iter->paired_sam_rec.second->chrom_name + "~" + paired_sam_rec_iter->paired_sam_rec.first->chrom_name;

			cur_fusion_junc_region.set(paired_sam_rec_iter->paired_sam_rec.second->start, paired_sam_rec_iter->paired_sam_rec.second->end,
				paired_sam_rec_iter->paired_sam_rec.first->start, paired_sam_rec_iter->paired_sam_rec.first->end);
		}

		if (junc_handler->m_fuson_junc_sort.find(chrom_name) != junc_handler->m_fuson_junc_sort.end())
		{
			vector<FusionJuncRegion>& cur_chrom_fusion_junc_regions = (junc_handler->m_fuson_junc_sort.find(chrom_name))->second;

			FindFusionJuncRegionStruct ffr(&cur_fusion_junc_region, &(*paired_sam_rec_iter), &cur_chrom_fusion_junc_regions);

			FindFusionJuncRegionStatic((void*)&ffr);
		}
	}
	
	return 0;
}

void*
AlignmentHandler::FindFusionJuncRegionStatic(void * str)
{
	FindFusionJuncRegionStruct* ffrs = (FindFusionJuncRegionStruct*) str;

	FusionJuncRegion& cur_fusion_region = *(ffrs->cur_fusion_region);

	PairedSamRec& cur_paired_sam_rec = *(ffrs->cur_paired_sam_rec);

	vector<FusionJuncRegion>& sorted_fusion_regions = *(ffrs->sorted_fusion_regions);

	vector<FusionJuncRegion>::iterator fusion_region_iter;

	fusion_region_iter = lower_bound(sorted_fusion_regions.begin(), sorted_fusion_regions.end(), cur_fusion_region, comp_fusion_junc_sort);

	if (fusion_region_iter != sorted_fusion_regions.begin())
		--fusion_region_iter;

	while(fusion_region_iter != sorted_fusion_regions.begin() && fusion_region_iter->m_doner_end > cur_fusion_region.m_doner_end)
		--fusion_region_iter;

	bool find_region = false;

	size_t regions_count = 0;

	vector<FusionJuncRegion>::iterator cur_fusion_region_iter = fusion_region_iter;

	while (true)
	{
		if (cur_fusion_region_iter == sorted_fusion_regions.end() || cur_fusion_region_iter->m_doner_st > cur_fusion_region.m_doner_end)
		{
			break;
		}
		else if ( ((cur_fusion_region_iter->m_doner_st <= cur_fusion_region.m_doner_st && cur_fusion_region_iter->m_doner_end >= cur_fusion_region.m_doner_end)) && 
				   ((cur_fusion_region_iter->m_acceptor_st <= cur_fusion_region.m_acceptor_st && cur_fusion_region_iter->m_acceptor_end >= cur_fusion_region.m_acceptor_end)))
		{
			++regions_count;

			//fragment 1
			if (cur_paired_sam_rec.paired_sam_rec.first->mate_match == "*" && cur_paired_sam_rec.paired_sam_rec.second->mate_match == "*")
			{
				cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_PAIRED_MAPPED;

				if (cur_paired_sam_rec.paired_sam_rec.second->strand_t & IS_REVERSE)
					cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_MATE_REVERSE;

				cur_paired_sam_rec.paired_sam_rec.first->strand_t |= IS_FIRST_END;

				if (cur_paired_sam_rec.paired_sam_rec.first->chrom_name == cur_paired_sam_rec.paired_sam_rec.second->chrom_name)
					cur_paired_sam_rec.paired_sam_rec.first->mate_match = "=";
				else
					cur_paired_sam_rec.paired_sam_rec.first->mate_match = cur_paired_sam_rec.paired_sam_rec.second->chrom_name;

				cur_paired_sam_rec.paired_sam_rec.first->mate_offset = cur_paired_sam_rec.paired_sam_rec.second->start;

				cur_paired_sam_rec.paired_sam_rec.first->mate_diff = static_cast<long> (cur_paired_sam_rec.paired_sam_rec.first->end) - static_cast<long> (cur_paired_sam_rec.paired_sam_rec.second->start);

				//fragment 2
				cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_PAIRED_MAPPED;

				if (cur_paired_sam_rec.paired_sam_rec.first->strand_t & IS_REVERSE)
					cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_MATE_REVERSE;

				cur_paired_sam_rec.paired_sam_rec.second->strand_t |= IS_SECOND_END;

				if (cur_paired_sam_rec.paired_sam_rec.first->chrom_name == cur_paired_sam_rec.paired_sam_rec.second->chrom_name)
					cur_paired_sam_rec.paired_sam_rec.second->mate_match = "=";
				else
					cur_paired_sam_rec.paired_sam_rec.second->mate_match = cur_paired_sam_rec.paired_sam_rec.first->chrom_name;

				cur_paired_sam_rec.paired_sam_rec.second->mate_offset = cur_paired_sam_rec.paired_sam_rec.first->start;

				cur_paired_sam_rec.paired_sam_rec.second->mate_diff = static_cast<long> (cur_paired_sam_rec.paired_sam_rec.second->end) - static_cast<long> (cur_paired_sam_rec.paired_sam_rec.first->start);
			}

			//add encompass read to junction 
			if (cur_fusion_region.m_doner_st == cur_paired_sam_rec.paired_sam_rec.first->start &&
				cur_fusion_region.m_doner_end == cur_paired_sam_rec.paired_sam_rec.first->end &&
				cur_fusion_region.m_acceptor_st == cur_paired_sam_rec.paired_sam_rec.second->start &&
				cur_fusion_region.m_acceptor_end == cur_paired_sam_rec.paired_sam_rec.second->end)
			{
				cur_paired_sam_rec.paired_sam_rec.first->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, true));

				cur_paired_sam_rec.paired_sam_rec.second->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, false));
			}
			else if (cur_fusion_region.m_doner_st == cur_paired_sam_rec.paired_sam_rec.second->start &&
					 cur_fusion_region.m_doner_end == cur_paired_sam_rec.paired_sam_rec.second->end &&
					 cur_fusion_region.m_acceptor_st == cur_paired_sam_rec.paired_sam_rec.first->start &&
					 cur_fusion_region.m_acceptor_end == cur_paired_sam_rec.paired_sam_rec.first->end)
			{
				cur_paired_sam_rec.paired_sam_rec.first->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, false));

				cur_paired_sam_rec.paired_sam_rec.second->encompassed_juncs.push_back(make_pair(cur_fusion_region_iter->m_junc_seed_ptr, true));
			}
			else
			{
				cerr <<"no matched alignments?"<<endl;
			}
		}

		++cur_fusion_region_iter;
	}

	return 0;
}

void 
AlignmentHandler::SetHits(vector<SamRec*>& samrecs)
{
	vector<SamRec*>::iterator sam_rec_vec_iter;

	for (sam_rec_vec_iter = samrecs.begin(); sam_rec_vec_iter != samrecs.end(); ++sam_rec_vec_iter)
	{
		if ((*sam_rec_vec_iter)->isunmapped)
			continue;

		SamRec2Hits(*sam_rec_vec_iter);
	}
}

void
AlignmentHandler::SamRec2Hits(SamRec* sam_rec_ptr)
{
	vector<pair<size_t, int> >::iterator splice_way_iter;

	hash_map<string, vector<bool> >::iterator mapped_pos_iter;
	
	mapped_pos_iter = m_mapped_pos.find(sam_rec_ptr->chrom_name);

	if (mapped_pos_iter == m_mapped_pos.end())
		return;

	vector<bool>& mapped_pos_ref = mapped_pos_iter->second;

	for (splice_way_iter = sam_rec_ptr->spliceway_vec.begin(); splice_way_iter != sam_rec_ptr->spliceway_vec.end(); ++splice_way_iter)
	{
		if (splice_way_iter->second > 0)
		{
			for (size_t offset = splice_way_iter->first; offset < splice_way_iter->first + splice_way_iter->second; ++offset)
			{
				mapped_pos_ref[offset] = true;
			}
		}
	}

	if (sam_rec_ptr->is_fusion)
	{
		mapped_pos_iter = m_mapped_pos.find(sam_rec_ptr->chrom_name2);

		if (mapped_pos_iter == m_mapped_pos.end())
			return;

		vector<bool>& mapped_pos_ref2 = mapped_pos_iter->second;

		for (splice_way_iter = sam_rec_ptr->spliceway_vec2.begin(); splice_way_iter != sam_rec_ptr->spliceway_vec2.end(); ++splice_way_iter)
		{
			if (splice_way_iter->second > 0)
			{
				for (size_t offset = splice_way_iter->first; offset < splice_way_iter->first + splice_way_iter->second; ++offset)
				{
					mapped_pos_ref2[offset] = true;
				}
			}
		}
	}
}

void
AlignmentHandler::FindSamVecContigLen(vector<SamRec*>& sam_rec_ptr)
{
	vector<SamRec*>::iterator sam_rec_iter;

	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		hash_map<string, vector<pair<size_t, size_t> > >::iterator chrom_expressed_regions_iter;

		chrom_expressed_regions_iter = m_expressed_regions.find((*sam_rec_iter)->chrom_name);

		if (chrom_expressed_regions_iter == m_expressed_regions.end())
			continue;

		vector<pair<size_t, size_t> >& sorted_regions = chrom_expressed_regions_iter->second;

		if (m_unioned_expressed_regions_map.find((*sam_rec_iter)->chrom_name) == m_unioned_expressed_regions_map.end())
			continue;

		hash_map<size_t, UnionExpressedRegions >& cur_unioned_expressed_regions = m_unioned_expressed_regions_map.find((*sam_rec_iter)->chrom_name)->second;

		pair<size_t, size_t> cur_junc_region_start = make_pair((*sam_rec_iter)->start, (*sam_rec_iter)->start);
		
		size_t start_region_idx = -1;

		bool find_start_region = FindRegion(cur_junc_region_start, sorted_regions, start_region_idx);

		DisjointSet& cur_disjointSet = (m_disjointset_map.find((*sam_rec_iter)->chrom_name))->second;

		if (find_start_region)
		{
			size_t source_idx = cur_disjointSet.FindSet(start_region_idx);

			(*sam_rec_iter)->m_contig_len = cur_unioned_expressed_regions[source_idx].m_unioned_region_length;
		}
	}
}


void* AlignmentHandler::FilterByFilteredJunctionStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);;
	
	JunctionHandler* junc_handler = m_junction_handler_prev_ptr;

	vector<SamRec*> filtered_sam_rec_ptr;

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->isunmapped)
			continue;

		vector<JunctionSeed*>::iterator junc_iter;

		bool filtered = false;

		if ((*sam_rec_ptr_iter)->filter_type != NOT_FILTERED) 
		{
			//cout << "corresponding junction is filtered by small anchor" << endl;
			//cerr << "corresponding junction is filtered by small anchor" << endl;
			filtered = true;
		}

		bool trimmed = false;

		if (filtered)
		{
			if (!(*sam_rec_ptr_iter)->is_fusion)
				trimmed = (*sam_rec_ptr_iter)->clip_by_small_anchor(m_add_S);
			else
			{
				trimmed = false;
			}
		}

		if (!filtered || trimmed)
		{
			filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);
		}
	}

	if (filtered_sam_rec_ptr.empty() && !sam_rec_ptr.empty())
	{
		filtered_sam_rec_ptr.push_back(sam_rec_ptr.front());

		filtered_sam_rec_ptr.front()->isunmapped = true;

		filtered_sam_rec_ptr.front()->paired_type = UNMAPPED;

		#ifdef DEBUG

		cout << "mark unmapped"<<endl;

		cout << filtered_sam_rec_ptr.front()->tostring(0, 0) << endl;

		if (filtered_sam_rec_ptr.front()->tag_name.find("HWI-EAS217:4:55:1764:913#0/1") != string::npos)
		{
			cout << "set unmapped in filter single multi"<<endl;

			cout << filtered_sam_rec_ptr.front()->tostring(0, 0) << endl;
		}

		#endif
	}

	sam_rec_ptr = filtered_sam_rec_ptr;

	if (sam_rec_ptr.size() > 1)
		RemoveDupStatic((void*)&sam_rec_ptr);
	return 0;
}

void
AlignmentHandler::FilterSingleMulti(vector<SamRec*>& sam_rec_ptr)
{
	if (sam_rec_ptr.size() <= 1)
	{
		return;
	}

	FilterByMisMatch(sam_rec_ptr);

	sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_filterscore);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*(sam_rec_ptr_iter - 1))->isexonic && (*(sam_rec_ptr_iter))->isexonic &&
			(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen)
			;
		else if (!(*(sam_rec_ptr_iter - 1))->isexonic && !(*(sam_rec_ptr_iter))->isexonic && 
			(*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
			(*(sam_rec_ptr_iter - 1))->ave_intron_len == (*(sam_rec_ptr_iter))->ave_intron_len &&
			(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
			(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
			&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
			)
			;
		else
		{
			sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

			break;
		}
	}
}

void* AlignmentHandler::FilterSingleMultiStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion)
		{
			//fragment 1
			(*sam_rec_ptr_iter)->strand_t2 |= IS_PAIRED_MAPPED;

			if ((*sam_rec_ptr_iter)->strand_t & IS_REVERSE)
				(*sam_rec_ptr_iter)->strand_t2 |= IS_MATE_REVERSE;

			(*sam_rec_ptr_iter)->strand_t2 |= IS_FIRST_END;

			(*sam_rec_ptr_iter)->strand_t2 |= IS_SECOND_END;

			if ((*sam_rec_ptr_iter)->chrom_name2 == (*sam_rec_ptr_iter)->chrom_name)
				(*sam_rec_ptr_iter)->mate_match2 = "=";
			else
				(*sam_rec_ptr_iter)->mate_match2 = (*sam_rec_ptr_iter)->chrom_name;

			(*sam_rec_ptr_iter)->mate_offset2 = (*sam_rec_ptr_iter)->start;

			(*sam_rec_ptr_iter)->mate_diff2 = (*sam_rec_ptr_iter)->fusion_suffix_end - (*sam_rec_ptr_iter)->fusion_prefix_st;

			//fragment 2
			(*sam_rec_ptr_iter)->strand_t |= IS_PAIRED_MAPPED;

			if ((*sam_rec_ptr_iter)->strand_t2 & IS_REVERSE)
				(*sam_rec_ptr_iter)->strand_t |= IS_MATE_REVERSE;

			//(*sam_rec_ptr_iter)->strand_t |= IS_FIRST_END;

			(*sam_rec_ptr_iter)->strand_t |= ((*sam_rec_ptr_iter)->end_id == 1 ? IS_FIRST_END : IS_SECOND_END);

			if ((*sam_rec_ptr_iter)->chrom_name == (*sam_rec_ptr_iter)->chrom_name2)
				(*sam_rec_ptr_iter)->mate_match = "=";
			else
				(*sam_rec_ptr_iter)->mate_match = (*sam_rec_ptr_iter)->chrom_name2;

			(*sam_rec_ptr_iter)->mate_offset = (*sam_rec_ptr_iter)->start2;

			(*sam_rec_ptr_iter)->mate_diff = (*sam_rec_ptr_iter)->fusion_prefix_st - (*sam_rec_ptr_iter)->fusion_suffix_end;
		}
	}

	if (sam_rec_ptr.size() == 1)
	{
		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;
	}

	if (sam_rec_ptr.size() <= 1)
	{
		return 0;
	}

	bool is_fusion = false, is_not_fusion = true;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion == true)
		{
			is_fusion = true;
		}
		else
			is_not_fusion = false;
	}

	if (is_fusion && is_not_fusion == false)  //contains both fusion and normal
	{
		FilterByMisMatchStatic((void*)&sam_rec_ptr);

		sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_fusion_filterscore);

		vector<SamRec*>::iterator sam_rec_ptr_iter;

		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;

		if (!(m_do_filter & 16))
		{
			return 0;
		}

		for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
		{
			if ((*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
				(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
				&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
				)
				;
			else
			{
				if ((m_do_filter & 256) == 0)
					sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());
				break;
			}
		}
	}
	else if (is_fusion && is_not_fusion)  // all fusion alignments
	{
		FilterByMisMatchStatic((void*)&sam_rec_ptr);

		sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_fusion_filterscore);

		vector<SamRec*>::iterator sam_rec_ptr_iter;

		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;

		if (!(m_do_filter & 16))
		{
			return 0;
		}

		for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
		{
			if ((*(sam_rec_ptr_iter - 1))->filter_score == (*(sam_rec_ptr_iter))->filter_score &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
				(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
				&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
				)
				;
			else
			{
				if ((m_do_filter & 256) == 0)
					sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

				break;
			}
		}
	}
	else  //all normal alignment
	{
		FilterByMisMatchStatic((void*)&sam_rec_ptr);

		sort(sam_rec_ptr.begin(), sam_rec_ptr.end(), comp_filterscore);

		if (sam_rec_ptr.front()->isexonic)
		{
		}
		else
		{
			vector<SamRec*>::iterator same_vps_iter;

			for (same_vps_iter = sam_rec_ptr.begin() + 1; same_vps_iter != sam_rec_ptr.end(); ++same_vps_iter)
			{
				if ((*(same_vps_iter - 1))->mappedlen == (*same_vps_iter)->mappedlen && !comp_filter_score_range((*(same_vps_iter - 1))->filter_score, (*same_vps_iter)->filter_score))
				{
				}
				else
					break;
			}

			if (same_vps_iter - sam_rec_ptr.begin() > 1)
				sort(sam_rec_ptr.begin(), same_vps_iter, comp_filterscore1);

			for (same_vps_iter = sam_rec_ptr.begin() + 1; same_vps_iter != sam_rec_ptr.end(); ++same_vps_iter)
			{
				if ((*(same_vps_iter - 1))->mappedlen == (*same_vps_iter)->mappedlen && !comp_filter_score_range((*(same_vps_iter - 1))->filter_score, (*same_vps_iter)->filter_score)
					&&!comp_intron_dist((*(same_vps_iter - 1))->ave_intron_len, (*same_vps_iter)->ave_intron_len))
				{
				}
				else
					break;
			}

			if (same_vps_iter - sam_rec_ptr.begin() > 1)
				sort(sam_rec_ptr.begin(), same_vps_iter, comp_filterscore2);
		}
		vector<SamRec*>::iterator sam_rec_ptr_iter;

		sam_rec_ptr.front()->strand_t -= IS_PRIMARY;

		if (!(m_do_filter & 16))
		{
			return 0;
		}

		for (sam_rec_ptr_iter = sam_rec_ptr.begin() + 1; sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
		{
			if ((*(sam_rec_ptr_iter - 1))->isexonic && (*(sam_rec_ptr_iter))->isexonic &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen)
				;
			else if (!(*(sam_rec_ptr_iter - 1))->isexonic && !(*(sam_rec_ptr_iter))->isexonic && 
				!comp_filter_score_range((*(sam_rec_ptr_iter - 1))->filter_score, (*(sam_rec_ptr_iter))->filter_score) &&
				!comp_intron_dist((*(sam_rec_ptr_iter - 1))->ave_intron_len, (*(sam_rec_ptr_iter))->ave_intron_len) &&
				(*(sam_rec_ptr_iter - 1))->mappedlen == (*(sam_rec_ptr_iter))->mappedlen &&
				(*(sam_rec_ptr_iter - 1))->ave_junc_mis == (*(sam_rec_ptr_iter))->ave_junc_mis
				&&(*(sam_rec_ptr_iter - 1))->m_contig_len == (*(sam_rec_ptr_iter))->m_contig_len
				)
				;
			else
			{
				if ((m_do_filter & 256) == 0)
					sam_rec_ptr.resize(sam_rec_ptr_iter - sam_rec_ptr.begin());

				break;
			}
		}
	}

	if (sam_rec_ptr.size() > m_max_hits)
	{
		sam_rec_ptr.resize(1);

		sam_rec_ptr.front()->isunmapped = true;

		sam_rec_ptr.front()->paired_type = UNMAPPED;
	}

	return 0;
}

void* AlignmentHandler::RemoveFusionUnPairedStatic(void * str)
{
	vector<SamRec*>& sam_rec_ptr = *((vector<SamRec*>*)str);

	vector<SamRec*>::iterator sam_rec_ptr_iter;

	if (!(m_do_filter & 16))
	{
		return 0;
	}

	bool is_fusion = false, is_not_fusion = false;

	for (sam_rec_ptr_iter = sam_rec_ptr.begin(); sam_rec_ptr_iter != sam_rec_ptr.end(); ++sam_rec_ptr_iter)
	{
		if ((*sam_rec_ptr_iter)->is_fusion == true)
			is_fusion = true;
		else
			is_not_fusion = false;
	}

	if (is_fusion && is_not_fusion)
	{
		cout << "is fusion and is not fusion:"<<(*sam_rec_ptr_iter)->tag_name<<endl;
	}
	else if (is_fusion)
	{
		sam_rec_ptr.clear();
	}
}

void AlignmentHandler::SetSamRecUnpaired(pair<vector<SamRec>, vector<SamRec> >& m_sam_rec_pe)
{
	vector<SamRec>::iterator sam_rec_iter;

	for (sam_rec_iter = m_sam_rec_pe.first.begin(); sam_rec_iter != m_sam_rec_pe.first.end(); ++sam_rec_iter)
	{
		sam_rec_iter->set_unmapped();
	}

	for (sam_rec_iter = m_sam_rec_pe.second.begin(); sam_rec_iter != m_sam_rec_pe.second.end(); ++sam_rec_iter)
	{
		sam_rec_iter->set_unmapped();
	}
}

void AlignmentHandler::SetSamRecUnpaired(pair<vector<SamRec*>, vector<SamRec*> >& m_sam_rec_pe_ptr)
{
	vector<SamRec*>::iterator sam_rec_iter;

	for (sam_rec_iter = m_sam_rec_pe_ptr.first.begin(); sam_rec_iter != m_sam_rec_pe_ptr.first.end(); ++sam_rec_iter)
	{
		(*sam_rec_iter)->set_unmapped();
	}

	for (sam_rec_iter = m_sam_rec_pe_ptr.second.begin(); sam_rec_iter != m_sam_rec_pe_ptr.second.end(); ++sam_rec_iter)
	{
		(*sam_rec_iter)->set_unmapped();
	}
}

bool AlignmentHandler::MergeStdFusion(vector<SamRec>::iterator sam_rec_ptr1, vector<SamRec>::iterator sam_rec_ptr2)
{
	//fusion 1
	sam_rec_ptr1->splice_way2 = sam_rec_ptr2->splice_way;

	sam_rec_ptr1->ori_splice_way2 = sam_rec_ptr2->ori_splice_way;

	sam_rec_ptr1->chrom_name2 = sam_rec_ptr2->chrom_name;

	sam_rec_ptr1->start2 = sam_rec_ptr2->start;

	sam_rec_ptr1->end2 = sam_rec_ptr2->end;

	sam_rec_ptr1->strand_t2 = sam_rec_ptr2->strand_t;

	sam_rec_ptr1->mappedlen1 = sam_rec_ptr1->mappedlen;

	sam_rec_ptr1->insertlen1 = sam_rec_ptr1->insertlen1;

	sam_rec_ptr1->mappedlen2 = sam_rec_ptr2->mappedlen;

	sam_rec_ptr1->insertlen2 = sam_rec_ptr2->insertlen1;

	sam_rec_ptr1->mappedlen = sam_rec_ptr1->mappedlen1 + sam_rec_ptr1->mappedlen2;

	sam_rec_ptr1->fusion_prefix_len = sam_rec_ptr1->mappedlen1;

	sam_rec_ptr1->fusion_suffix_len = sam_rec_ptr1->mappedlen2;

	sam_rec_ptr1->intron_size = sam_rec_ptr1->intron_size + sam_rec_ptr2->intron_size;

	sam_rec_ptr1->spliceway_vec2 = sam_rec_ptr2->spliceway_vec;

	sam_rec_ptr1->is_fusion = true;

	if (sam_rec_ptr1->strand_t & IS_REVERSE)
	{
		sam_rec_ptr1->strand1 = '-';
	}
	else
	{
		sam_rec_ptr1->strand1 = '+';
	}

	if (sam_rec_ptr1->strand_t2 & IS_REVERSE)
	{
		sam_rec_ptr1->strand2 = '-';
	}
	else
	{
		sam_rec_ptr1->strand2 = '+';
	}

	if (sam_rec_ptr1->strand1 == '+')
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->start;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->end;
	}
	else
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->end;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->start;
	}

	if (sam_rec_ptr1->strand2 == '+')
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->start2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->end2;
	}
	else
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->end2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->start2;
	}

	if (sam_rec_ptr1->min_anchor > sam_rec_ptr2->min_anchor)
	{
		sam_rec_ptr1->min_anchor = sam_rec_ptr2->min_anchor;
	}

	if (sam_rec_ptr1->max_anchor < sam_rec_ptr2->max_anchor)
	{
		sam_rec_ptr1->max_anchor = sam_rec_ptr2->max_anchor;
	}

	if (true|| sam_rec_ptr1->chrom_name < sam_rec_ptr1->chrom_name2 || (sam_rec_ptr1->chrom_name == sam_rec_ptr1->chrom_name2 && sam_rec_ptr1->start <= sam_rec_ptr1->start2))
	{
		if (!(sam_rec_ptr1->chrom_name < sam_rec_ptr1->chrom_name2 || (sam_rec_ptr1->chrom_name == sam_rec_ptr1->chrom_name2 && sam_rec_ptr1->start <= sam_rec_ptr1->start2)))
			sam_rec_ptr1->need_swap = true;
	}
	else
	{
		sam_rec_ptr1->Swap();
	}
	return 0;
}


void AlignmentHandler::CollectStats(const vector<SamRec*>& sam_rec_ptr)
{
	if (sam_rec_ptr.empty())
		return;

	vector<SamRec*>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;
	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((*sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((*sam_rec_iter)->isexonic)
			isexonic = true;

		if ((*sam_rec_iter)->isspliced)
			isspliced = true;

		if ((*sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((*sam_rec_iter)->issmallins)
			issmallins = true;

		if ((*sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((*sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((*sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((*sam_rec_iter)->isclipped)
			isclipped = true;

		if ((*sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}

	if (isexonic)
		++m_unspliced;

	if (isspliced)
		++m_spliced;

	if (issmalldel)
		++m_deletion;

	if (issmallins)
		++m_insertion;

	if (ispaired)
		++m_paired;
	else if (!isunmapped)
		++m_single;

	if (iscanon)
		++m_canoical;

	if (issemicanon)
		++m_semi_canonical;

	if (isnoncanon)
		++m_non_canonical;

	if (isclipped)
		++m_clipped;

	if (isfusion)
		++m_fusion;

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
		++m_multiple;
	else if (sam_rec_ptr.size() == 1)
		++m_unique;
}

void* AlignmentHandler::CollectStatsStatic(void * str, bool both_unspliced)
{
	const vector<SamRec*>& sam_rec_ptr = (*(vector<SamRec*>*)str);

	if (sam_rec_ptr.empty())
		return 0;

	vector<SamRec*>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;

	bool paired_type_single = false, paired_type_fusion_paired = false, paired_type_normal_paired = false;

	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((*sam_rec_iter)->paired_type == SINGLE)
			paired_type_single = true;
		else if ((*sam_rec_iter)->paired_type == FUSION_PAIRED)
			paired_type_fusion_paired = true;
		else if ((*sam_rec_iter)->paired_type == NORMAL_PAIRED)
			paired_type_normal_paired = true;

		if ((*sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((*sam_rec_iter)->isexonic && !(*sam_rec_iter)->is_fusion)
			isexonic = true;

		if ((*sam_rec_iter)->isspliced)
			isspliced = true;

		if ((*sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((*sam_rec_iter)->issmallins)
			issmallins = true;

		if ((*sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((*sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((*sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((*sam_rec_iter)->isclipped)
			isclipped = true;

		if ((*sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}

	if (isexonic)
	{
		++m_unspliced;

		if (paired_type_normal_paired)
			++m_unspliced_paired;
		else if (paired_type_single)
			++m_unspliced_single;
		else if (paired_type_fusion_paired)
			++m_unspliced_fusion_paired;
	}

	if (isspliced)
	{
		++m_spliced;

		if (paired_type_normal_paired)
			++m_spliced_paired;
		else if (paired_type_single)
			++m_spliced_single;
		else if (paired_type_fusion_paired)
			++m_spliced_fusion_paired;
	}

	if (issmalldel)
	{
		++m_deletion;

		if (paired_type_normal_paired)
			++m_deletion_paired;
		else if (paired_type_single)
			++m_deletion_single;
		else if (paired_type_fusion_paired)
			++m_deletion_fusion_paired;
	}

	if (issmallins)
	{
		++m_insertion;

		if (paired_type_normal_paired)
			++m_insertion_paired;
		else if (paired_type_single)
			++m_insertion_single;
		else if (paired_type_fusion_paired)
			++m_insertion_fusion_paired;
	}

	if (paired_type_normal_paired)
		++m_paired;
	else if (paired_type_single)
		++m_single;
	else if (paired_type_fusion_paired)
		++m_fusion_paired;

	if (iscanon && !isfusion && !issmalldel)
	{
		++m_canoical;

		if (paired_type_normal_paired)
			++m_canoical_paired;
		else if (paired_type_single)
			++m_canoical_single;
		else if (paired_type_fusion_paired)
			++m_canoical_fusion_paired;
	}

	if (issemicanon && !isfusion && !issmalldel)
	{
		++m_semi_canonical;

		if (paired_type_normal_paired)
			++m_semi_canonical_paired;
		else if (paired_type_single)
			++m_semi_canonical_single;
		else if (paired_type_fusion_paired)
			++m_semi_canonical_fusion_paired;
	}

	if (isnoncanon && !isfusion && !issmalldel)
	{
		++m_non_canonical;

		if (paired_type_normal_paired)
			++m_non_canonical_paired;
		else if (paired_type_single)
			++m_non_canonical_single;
		else if (paired_type_fusion_paired)
			++m_non_canonical_fusion_paired;
	}

	if (isclipped)
	{
		++m_clipped;

		if (paired_type_normal_paired)
			++m_clipped_paired;
		else if (paired_type_single)
			++m_clipped_single;
		else if (paired_type_fusion_paired)
			++m_clipped_fusion_paired;
	}

	if (isfusion)
	{
		++m_fusion;

		if (iscanon)
			++m_fusion_canonical;

		if (issemicanon)
			++m_fusion_semi_canonical;

		if (isnoncanon)
			++m_fusion_non_canonical;

		if (sam_rec_ptr.size() > 1)
			++m_fusion_multiple;
		else if (sam_rec_ptr.size() == 1)
			++m_fusion_unique;
	}

	if (both_unspliced)
	{
		if (isunmapped)
			;
		else if (sam_rec_ptr.size() == 1)
		{
			++m_both_unspliced_unique;

			if (paired_type_normal_paired)
				++m_both_unspliced_paired_unique;
			else if (paired_type_single)
				++m_both_unspliced_single_unique;
			else if (paired_type_fusion_paired)
				++m_both_unspliced_fusion_paired_unique;
		}
		else if (sam_rec_ptr.size() > 1)
		{
			++m_both_unspliced_multiple;

			if (paired_type_normal_paired)
				++m_both_unspliced_paired_multiple;
			else if (paired_type_single)
				++m_both_unspliced_single_multiple;
			else if (paired_type_fusion_paired)
				++m_both_unspliced_fusion_paired_multiple;
		}		

	}

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
	{
		++m_multiple;

		if (paired_type_normal_paired)
			++m_paired_multiple;
		else if (paired_type_single)
			++m_single_multiple;
		else if (paired_type_fusion_paired)
			++m_fusion_paired_multiple;
	}
	else if (sam_rec_ptr.size() == 1)
	{
		++m_unique;

		if (paired_type_normal_paired)
			++m_paired_unique;
		else if (paired_type_single)
			++m_single_unique;
		else if (paired_type_fusion_paired)
			++m_fusion_paired_unique;
	}

	return 0;
}

void* AlignmentHandler::CollectStatsStatic(void * str)
{
	const vector<SamRec*>& sam_rec_ptr = (*(vector<SamRec*>*)str);

	if (sam_rec_ptr.empty())
		return 0;

	vector<SamRec*>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;

	bool paired_type_single = false, paired_type_fusion_paired = false, paired_type_normal_paired = false;

	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((*sam_rec_iter)->paired_type == SINGLE)
			paired_type_single = true;
		else if ((*sam_rec_iter)->paired_type == FUSION_PAIRED)
			paired_type_fusion_paired = true;
		else if ((*sam_rec_iter)->paired_type == NORMAL_PAIRED)
			paired_type_normal_paired = true;

		if ((*sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((*sam_rec_iter)->isexonic && !(*sam_rec_iter)->is_fusion)
			isexonic = true;

		if ((*sam_rec_iter)->isspliced)
			isspliced = true;

		if ((*sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((*sam_rec_iter)->issmallins)
			issmallins = true;

		if ((*sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((*sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((*sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((*sam_rec_iter)->isclipped)
			isclipped = true;

		if ((*sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}

	if (isexonic)
	{
		++m_unspliced;

		if (paired_type_normal_paired)
			++m_unspliced_paired;
		else if (paired_type_single)
			++m_unspliced_single;
		else if (paired_type_fusion_paired)
			++m_unspliced_fusion_paired;
	}

	if (isspliced)
	{
		++m_spliced;

		if (paired_type_normal_paired)
			++m_spliced_paired;
		else if (paired_type_single)
			++m_spliced_single;
		else if (paired_type_fusion_paired)
			++m_spliced_fusion_paired;
	}

	if (issmalldel)
	{
		++m_deletion;

		if (paired_type_normal_paired)
			++m_deletion_paired;
		else if (paired_type_single)
			++m_deletion_single;
		else if (paired_type_fusion_paired)
			++m_deletion_fusion_paired;
	}

	if (issmallins)
	{
		++m_insertion;

		if (paired_type_normal_paired)
			++m_insertion_paired;
		else if (paired_type_single)
			++m_insertion_single;
		else if (paired_type_fusion_paired)
			++m_insertion_fusion_paired;
	}

	if (paired_type_normal_paired)
		++m_paired;
	else if (paired_type_single)
		++m_single;
	else if (paired_type_fusion_paired)
		++m_fusion_paired;

	if (iscanon && !isfusion && !issmalldel)
	{
		++m_canoical;

		if (paired_type_normal_paired)
			++m_canoical_paired;
		else if (paired_type_single)
			++m_canoical_single;
		else if (paired_type_fusion_paired)
			++m_canoical_fusion_paired;
	}

	if (issemicanon && !isfusion && !issmalldel)
	{
		++m_semi_canonical;

		if (paired_type_normal_paired)
			++m_semi_canonical_paired;
		else if (paired_type_single)
			++m_semi_canonical_single;
		else if (paired_type_fusion_paired)
			++m_semi_canonical_fusion_paired;
	}

	if (isnoncanon && !isfusion && !issmalldel)
	{
		++m_non_canonical;

		if (paired_type_normal_paired)
			++m_non_canonical_paired;
		else if (paired_type_single)
			++m_non_canonical_single;
		else if (paired_type_fusion_paired)
			++m_non_canonical_fusion_paired;
	}

	if (isclipped)
	{
		++m_clipped;

		if (paired_type_normal_paired)
			++m_clipped_paired;
		else if (paired_type_single)
			++m_clipped_single;
		else if (paired_type_fusion_paired)
			++m_clipped_fusion_paired;
	}

	if (isfusion)
	{
		++m_fusion;

		if (iscanon)
			++m_fusion_canonical;

		if (issemicanon)
			++m_fusion_semi_canonical;

		if (isnoncanon)
			++m_fusion_non_canonical;

		if (sam_rec_ptr.size() > 1)
			++m_fusion_multiple;
		else if (sam_rec_ptr.size() == 1)
			++m_fusion_unique;
	}

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
	{
		++m_multiple;

		if (paired_type_normal_paired)
			++m_paired_multiple;
		else if (paired_type_single)
			++m_single_multiple;
		else if (paired_type_fusion_paired)
			++m_fusion_paired_multiple;
	}
	else if (sam_rec_ptr.size() == 1)
	{
		++m_unique;

		if (paired_type_normal_paired)
			++m_paired_unique;
		else if (paired_type_single)
			++m_single_unique;
		else if (paired_type_fusion_paired)
			++m_fusion_paired_unique;
	}

	return 0;
}

void AlignmentHandler::CollectStats(const vector<SamRec>& sam_rec_ptr)
{
	if (sam_rec_ptr.empty())
		return;

	vector<SamRec>::const_iterator sam_rec_iter;

	bool isexonic = false, isfusion = false, isspliced = false, issmalldel = false, issmallins = false, ispaired = false, iscanon = false, issemicanon = false, isnoncanon = false, isunmapped = false, isclipped = false;
	for (sam_rec_iter = sam_rec_ptr.begin(); sam_rec_iter != sam_rec_ptr.end(); ++sam_rec_iter)
	{
		if ((sam_rec_iter)->isunmapped)
		{
			isunmapped = true;
			continue;
		}

		if ((sam_rec_iter)->is_fusion)
			isfusion = true;

		if ((sam_rec_iter)->isexonic)
			isexonic = true;

		if ((sam_rec_iter)->isspliced)
			isspliced = true;

		if ((sam_rec_iter)->issmalldel)
			issmalldel = true;

		if ((sam_rec_iter)->issmallins)
			issmallins = true;

		if ((sam_rec_iter)->iscanonical)
			iscanon = true;

		if ((sam_rec_iter)->issemicanonical)
			issemicanon = true;

		if ((sam_rec_iter)->isnoncanoical)
			isnoncanon = true;

		if ((sam_rec_iter)->isclipped)
			isclipped = true;

		if ((sam_rec_iter)->mate_match != "*" )
			ispaired = true;
	}
	

	if (isexonic)
		++m_unspliced;

	if (isspliced)
		++m_spliced;

	if (issmalldel)
		++m_deletion;

	if (issmallins)
		++m_insertion;

	if (ispaired)
		++m_paired;
	else if (!isunmapped)
		++m_single;

	if (iscanon)
		++m_canoical;

	if (issemicanon)
		++m_semi_canonical;

	if (isnoncanon)
		++m_non_canonical;

	if (isclipped)
		++m_clipped;

	if (isfusion)
		++m_fusion;

	if (isunmapped)
		++m_unmapped;
	else if (sam_rec_ptr.size() > 1)
		++m_multiple;
	else if (sam_rec_ptr.size() == 1)
		++m_unique;
}

void AlignmentHandler::WriteStats(string stat_file, string headline)
{
	cout << "write stat"<<endl;

	ofstream ofs(stat_file.c_str(), ofstream::app);

	ofs << headline << endl;

	ofs <<"unmapped\t"<<m_unmapped<<endl;

	ofs <<"unspliced\t"<<m_unspliced<<endl;

	ofs <<"spliced\t"<<m_spliced<<endl;

	ofs <<"deletion\t"<<m_deletion<<endl;

	ofs <<"insertion\t"<<m_insertion<<endl;

	ofs <<"fusion_reads\t"<<m_fusion<<endl;

	ofs <<"paired\t"<<m_paired<<endl;

	ofs <<"single\t"<<m_single<<endl;

	ofs <<"fusion_pair\t"<<m_fusion_paired<<endl;

	ofs <<"canoical\t"<<m_canoical<<endl;

	ofs <<"semi_canonical\t"<<m_semi_canonical<<endl;

	ofs <<"non_canonical\t"<<m_non_canonical<<endl;

	ofs <<"multiple\t"<<m_multiple<<endl;

	ofs <<"unique\t"<<m_unique<<endl;

	ofs <<"clipped\t"<<m_clipped<<endl;

	ofs <<"both_unspliced_multiple\t"<<m_both_unspliced_multiple<<endl;

	ofs <<"both_unspliced_unique\t"<<m_both_unspliced_unique<<endl;

	ofs <<"both_unspliced_paired_multiple\t"<<m_both_unspliced_paired_multiple<<endl;

	ofs <<"both_unspliced_paired_unique\t"<<m_both_unspliced_paired_unique<<endl;

	ofs <<"both_unspliced_single_multiple\t"<<m_both_unspliced_single_multiple<<endl;

	ofs <<"both_unspliced_single_unique\t"<<m_both_unspliced_single_unique<<endl;

	ofs <<"both_unspliced_fusion_paired_multiple\t"<<m_both_unspliced_fusion_paired_multiple<<endl;

	ofs <<"both_unspliced_fusion_paired_unique\t"<<m_both_unspliced_fusion_paired_unique<<endl;

	///////////

	ofs <<"paired_multiple\t"<<m_paired_multiple<<endl;

	ofs <<"paired_unique\t"<<m_paired_unique<<endl;

	ofs <<"single_multiple\t"<<m_single_multiple<<endl;

	ofs <<"single_unique\t"<<m_single_unique<<endl;

	ofs <<"fusion_paired_multiple\t"<<m_fusion_paired_multiple<<endl;

	ofs <<"fusion_paired_unique\t"<<m_fusion_paired_unique<<endl;

	/////////////

	ofs <<"fusion_canonical\t"<<m_fusion_canonical<<endl;

	ofs <<"fusion_semi_canonical\t"<<m_fusion_semi_canonical<<endl;

	ofs <<"fusion_non_canonical\t"<<m_fusion_non_canonical<<endl;

	ofs <<"fusion_multiple\t"<<m_fusion_multiple<<endl;

	ofs <<"fusion_unique\t"<<m_fusion_unique<<endl;

	///////////////

	ofs <<"spliced_paired\t"<<m_spliced_paired<<endl;

	ofs <<"spliced_single\t"<<m_spliced_single<<endl;

	ofs <<"spliced_fusion_paired\t"<<m_spliced_fusion_paired<<endl;

	ofs <<"unspliced_paired\t"<<m_unspliced_paired<<endl;

	ofs <<"unspliced_single\t"<<m_unspliced_single<<endl;

	ofs <<"unspliced_fusion_paired\t"<<m_unspliced_fusion_paired<<endl;

	ofs <<"insertion_paired\t"<<m_insertion_paired<<endl;

	ofs <<"insertion_single\t"<<m_insertion_single<<endl;

	ofs <<"insertion_fusion_paired\t"<<m_insertion_fusion_paired<<endl;

	ofs <<"deletion_paired\t"<<m_deletion_paired<<endl;

	ofs <<"deletion_single\t"<<m_deletion_single<<endl;

	ofs <<"deletion_fusion_paired\t"<<m_deletion_fusion_paired<<endl;

	ofs <<"canoical_paired\t"<<m_canoical_paired<<endl;

	ofs <<"canoical_single\t"<<m_canoical_single<<endl;

	ofs <<"canoical_fusion_paired\t"<<m_canoical_fusion_paired<<endl;

	ofs <<"semi_canonical_paired\t"<<m_semi_canonical_paired<<endl;

	ofs <<"semi_canonical_single\t"<<m_semi_canonical_single<<endl;

	ofs <<"semi_canonical_fusion_paired\t"<<m_semi_canonical_fusion_paired<<endl;

	ofs <<"non_canonical_paired\t"<<m_non_canonical_paired<<endl;

	ofs <<"non_canonical_single\t"<<m_non_canonical_single<<endl;

	ofs <<"non_canonical_fusion_paired\t"<<m_non_canonical_fusion_paired<<endl;

	ofs <<"clipped_paired\t"<<m_clipped_paired<<endl;

	ofs <<"clipped_single\t"<<m_clipped_single<<endl;

	ofs <<"clipped_fusion_paired\t"<<m_clipped_fusion_paired<<endl;

	//m_junction_handler.WriteStats(stat_file);
}

void AlignmentHandler::WriteFinalStats(string stat_file, string headline)
{
	cout << "write stat"<<endl;

	ofstream ofs(stat_file.c_str(), ofstream::app);

	ofs << headline << endl;

	ofs <<"unmapped\t"<<m_unmapped<<endl;

	ofs <<"unspliced\t"<<m_unspliced<<endl;

	ofs <<"spliced\t"<<m_spliced<<endl;

	ofs <<"deletion\t"<<m_deletion<<endl;

	ofs <<"insertion\t"<<m_insertion<<endl;

	ofs <<"fusion_reads\t"<<m_fusion<<endl;

	ofs <<"paired\t"<<m_paired<<endl;

	ofs <<"single\t"<<m_single<<endl;

	ofs <<"fusion_paired\t"<<m_fusion_paired<<endl;

	ofs <<"canoical\t"<<m_canoical<<endl;

	ofs <<"semi_canonical\t"<<m_semi_canonical<<endl;

	ofs <<"non_canonical\t"<<m_non_canonical<<endl;

	ofs <<"multiple\t"<<m_multiple<<endl;

	ofs <<"unique\t"<<m_unique<<endl;

	ofs <<"clipped\t"<<m_clipped<<endl;

	ofs <<"both_unspliced_multiple\t"<<m_both_unspliced_multiple<<endl;

	ofs <<"both_unspliced_unique\t"<<m_both_unspliced_unique<<endl;

	ofs <<"both_unspliced_paired_multiple\t"<<m_both_unspliced_paired_multiple<<endl;

	ofs <<"both_unspliced_paired_unique\t"<<m_both_unspliced_paired_unique<<endl;

	ofs <<"both_unspliced_single_multiple\t"<<m_both_unspliced_single_multiple<<endl;

	ofs <<"both_unspliced_single_unique\t"<<m_both_unspliced_single_unique<<endl;

	ofs <<"both_unspliced_fusion_paired_multiple\t"<<m_both_unspliced_fusion_paired_multiple<<endl;

	ofs <<"both_unspliced_fusion_paired_unique\t"<<m_both_unspliced_fusion_paired_unique<<endl;

	///////////

	ofs <<"paired_multiple\t"<<m_paired_multiple<<endl;

	ofs <<"paired_unique\t"<<m_paired_unique<<endl;

	ofs <<"single_multiple\t"<<m_single_multiple<<endl;

	ofs <<"single_unique\t"<<m_single_unique<<endl;

	ofs <<"fusion_paired_multiple\t"<<m_fusion_paired_multiple<<endl;

	ofs <<"fusion_paired_unique\t"<<m_fusion_paired_unique<<endl;

	/////////////

	ofs <<"fusion_canonical\t"<<m_fusion_canonical<<endl;

	ofs <<"fusion_semi_canonical\t"<<m_fusion_semi_canonical<<endl;

	ofs <<"fusion_non_canonical\t"<<m_fusion_non_canonical<<endl;

	ofs <<"fusion_multiple\t"<<m_fusion_multiple<<endl;

	ofs <<"fusion_unique\t"<<m_fusion_unique<<endl;

	///////////////

	ofs <<"spliced_paired\t"<<m_spliced_paired<<endl;

	ofs <<"spliced_single\t"<<m_spliced_single<<endl;

	ofs <<"spliced_fusion_paired\t"<<m_spliced_fusion_paired<<endl;

	ofs <<"unspliced_paired\t"<<m_unspliced_paired<<endl;

	ofs <<"unspliced_single\t"<<m_unspliced_single<<endl;

	ofs <<"unspliced_fusion_paired\t"<<m_unspliced_fusion_paired<<endl;

	ofs <<"insertion_paired\t"<<m_insertion_paired<<endl;

	ofs <<"insertion_single\t"<<m_insertion_single<<endl;

	ofs <<"insertion_fusion_paired\t"<<m_insertion_fusion_paired<<endl;

	ofs <<"deletion_paired\t"<<m_deletion_paired<<endl;

	ofs <<"deletion_single\t"<<m_deletion_single<<endl;

	ofs <<"deletion_fusion_paired\t"<<m_deletion_fusion_paired<<endl;

	ofs <<"canoical_paired\t"<<m_canoical_paired<<endl;

	ofs <<"canoical_single\t"<<m_canoical_single<<endl;

	ofs <<"canoical_fusion_paired\t"<<m_canoical_fusion_paired<<endl;

	ofs <<"semi_canonical_paired\t"<<m_semi_canonical_paired<<endl;

	ofs <<"semi_canonical_single\t"<<m_semi_canonical_single<<endl;

	ofs <<"semi_canonical_fusion_paired\t"<<m_semi_canonical_fusion_paired<<endl;

	ofs <<"non_canonical_paired\t"<<m_non_canonical_paired<<endl;

	ofs <<"non_canonical_single\t"<<m_non_canonical_single<<endl;

	ofs <<"non_canonical_fusion_paired\t"<<m_non_canonical_fusion_paired<<endl;

	ofs <<"clipped_paired\t"<<m_clipped_paired<<endl;

	ofs <<"clipped_single\t"<<m_clipped_single<<endl;

	ofs <<"clipped_fusion_paired\t"<<m_clipped_fusion_paired<<endl;
	
	//m_junction_handler.WriteStats(stat_file);
}

