#ifndef SAM_INFO_H
#define SAM_INFO_H

#include <string.h>
using namespace std;

class Sam_Info
{
public:
	Sam_Info()
	{
		start_pos=-1;
		printed = false;
		forward_strand_junc = false;
		reverse_strand_junc = false;
		unknown_strand_junc = false;
		deletion = false;
		imcomplete_single_anchor = false;
		num_junc_penalty = 0;
	}

	~Sam_Info()
	{
		mismatchs.clear();
		jump_code.clear();
		junc_flank_seq.clear();
		fusion_internal.clear();
		fusion_head.clear();
		fusion_tail.clear();
	}

	Sam_Info(const Sam_Info& my_sam_info)
	{
		start_pos = my_sam_info.start_pos;
		end_pos = my_sam_info.end_pos;
		mapped_len = my_sam_info.mapped_len;
		intron_len = my_sam_info.intron_len;
		insert_len = my_sam_info.insert_len;
		for(size_t i=0;i<my_sam_info.mismatchs.size();i++)
		{
			mismatchs.push_back(my_sam_info.mismatchs[i]);
		}
		for(size_t i=0;i<my_sam_info.jump_code.size();i++)
		{
			jump_code.push_back(my_sam_info.jump_code[i]);
		}
		start_contig_index=my_sam_info.start_contig_index;
		start_seg_no=my_sam_info.start_seg_no;
		end_seg_no=my_sam_info.end_seg_no;
		chrom=my_sam_info.chrom;
		pair_no = my_sam_info.pair_no;
		strand = my_sam_info.strand;
		printed = false;
		forward_strand_junc = my_sam_info.forward_strand_junc;
		reverse_strand_junc = my_sam_info.reverse_strand_junc;
		unknown_strand_junc = my_sam_info.unknown_strand_junc;
		num_junc_penalty = my_sam_info.num_junc_penalty;
		deletion = my_sam_info.deletion;
		imcomplete_single_anchor = my_sam_info.imcomplete_single_anchor;
		for(size_t i = 0; i< my_sam_info.junc_flank_seq.size(); i++)
		{
			junc_flank_seq.push_back(my_sam_info.junc_flank_seq[i]);	
		}
		for(size_t i = 0; i < my_sam_info.fusion_internal.size(); i++)
			fusion_internal.push_back(my_sam_info.fusion_internal[i]);
		for(size_t i = 0; i < my_sam_info.fusion_head.size(); i++)
			fusion_head.push_back(my_sam_info.fusion_head[i]);
		for(size_t i = 0; i < my_sam_info.fusion_tail.size(); i++)
			fusion_tail.push_back(my_sam_info.fusion_tail[i]);
	}

	Sam_Info(Bwtmap_Info& my_bwt_info, size_t current_index, Read_Block* r_block, int segment_number)
	{
		start_pos = my_bwt_info.start;
		int new_mapped_len = my_bwt_info.end - my_bwt_info.start + 1; 
		Jump_Code new_jumpcode(new_mapped_len, "M");
		jump_code.push_back(new_jumpcode);
		start_contig_index=current_index;
		start_seg_no=my_bwt_info.start_seg_no;
		end_seg_no=my_bwt_info.end_seg_no;
		chrom=my_bwt_info.chrom;
		pair_no = my_bwt_info.pair_no;
		strand = my_bwt_info.strand;
		printed = false;
		forward_strand_junc = false;
		reverse_strand_junc = false;
		unknown_strand_junc = false;
		num_junc_penalty = 0;
		deletion = false;
		imcomplete_single_anchor = false;
	}

	void Get_Extra_Info(string& chrom_seq, int max_del, Read_Block* r_block)
	{
		end_pos = start_pos - 1;
		mapped_len = 0;
		intron_len = 0;
		insert_len = 0;
		string map_seq;
		string read_seq;
		if(strand == "+")     //output sequence
		{
			for(int i = start_seg_no; i <= end_seg_no; i++)
			{
				read_seq.append( r_block->get_seg_seq(i) );
			}
		}
		else                   //out put - sequence
		{
			for(int i = start_seg_no; i >= end_seg_no; i--)
			{
				read_seq.append( r_block->get_revcom_seg_seq(i) );
			}
		}
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
				int junc_start = end_pos + 1;
				end_pos += jump_code[i].len;
				int junc_end = end_pos;
				intron_len += jump_code[i].len;
				string flank_seq = chrom_seq.substr(junc_start, 2) + chrom_seq.substr(junc_end - 1, 2);
				junc_flank_seq.push_back(flank_seq);
				if(jump_code[i].len <= max_del)
					deletion = true;
				else if(flank_seq == "GTAG" || flank_seq == "ATAC" || flank_seq == "GCAG")
					forward_strand_junc = true;
				else if(flank_seq == "CTAC" || flank_seq == "CTGC" || flank_seq == "GTAT")
					reverse_strand_junc = true;
				else
					unknown_strand_junc = true;
				num_junc_penalty ++;
			}
			else   // "I"
			{
				map_seq.append(read_seq.substr(mapped_len, jump_code[i].len));
				mapped_len += jump_code[i].len;
				insert_len += jump_code[i].len;
			}
		}
		for(size_t i = 0; i < read_seq.length(); i++)
		{
			if(map_seq[i] != read_seq[i])
			{
				Mismatch new_mis(1, i, map_seq[i], read_seq[i]);
				mismatchs.push_back(new_mis);		
			}
		}
		/*cout << r_block->get_read_id() << "\t" << chrom << "\t" << start_pos << "\t" << end_pos << "\t" << read_seq << "\t" << map_seq << "\t";
		for(size_t i = 0; i < jump_code.size(); i++)
			cout << jump_code[i].len << jump_code[i].type;
		cout << endl;*/
	}
	
	int count_gap()
	{
		int num_gap = 0;
		for(size_t i = 0; i < jump_code.size(); i++)
	  {
	  	if(jump_code[i].type != "M")
	  		num_gap ++;
	  }
	  return num_gap;
	}
	string chrom;
	string strand;
	int start_pos;
	int end_pos;
	int intron_len;
	int insert_len;
	int mapped_len;
	int num_junc_penalty;
	vector<Mismatch> mismatchs;
	vector<Jump_Code> jump_code;
	vector<string> junc_flank_seq;
	bool forward_strand_junc;
	bool reverse_strand_junc;
	bool unknown_strand_junc;
	bool deletion;
	bool imcomplete_single_anchor;
	size_t start_contig_index;
	int start_seg_no;
	int end_seg_no;
	int pair_no;
	bool printed;
	vector<Fusion_Splice> fusion_internal;
	vector<Fusion_Splice> fusion_head;
	vector<Fusion_Splice> fusion_tail;
};

class Paired_Sam_Info 
{
public:
	
	Paired_Sam_Info(Sam_Info* si1, Sam_Info* si2, int expected_mate_dist)
	{
		mate_dist = (si1->strand == "+") ? (si2->start_pos - si1->end_pos) : (si1->start_pos - si2->end_pos);
		mate_dist_diff = abs(mate_dist - expected_mate_dist);
		mappedlen = si1->mapped_len + si2->mapped_len;
		intron_size = si1->intron_len + si2->intron_len;
		insert_size = si1->insert_len + si2->insert_len;
		total_mismatch = si1->mismatchs.size() + si2->mismatchs.size() + si1->num_junc_penalty + si2->num_junc_penalty;
		paired_sam = make_pair(si1, si2);
	}

	~Paired_Sam_Info()
	{
	}

	int mate_dist;
	int mate_dist_diff;
	int intron_size;
	int insert_size;
	int mappedlen;
	size_t total_mismatch;
	pair<Sam_Info*, Sam_Info*> paired_sam;
};

#endif

