#include "SamRec.h"
#include "ReadNextTagAlignHandler.h"

size_t convert_point_coordinate(size_t cur_point, size_t x, size_t y, size_t x_prime, size_t y_prime, char strand)
{
	size_t cur_point_prime;

	if (strand == '+')
	{
		cur_point_prime = cur_point - x + x_prime;
	}
	else if (strand == '-')
	{
		cur_point_prime = y - cur_point + x_prime;
	}
	else
	{
		cerr << "strand?"<<endl;

		exit(0);
	}

	return cur_point_prime;
}


bool MergeStdFusion(vector<SamRec>::iterator sam_rec_ptr1, vector<SamRec>::iterator sam_rec_ptr2)
{
	//fusion 1
	sam_rec_ptr1->splice_way2 = sam_rec_ptr2->splice_way;

	sam_rec_ptr1->ori_splice_way2 = sam_rec_ptr2->ori_splice_way;

	sam_rec_ptr1->chrom_name2 = sam_rec_ptr2->chrom_name;

	sam_rec_ptr1->start2 = sam_rec_ptr2->start;

	sam_rec_ptr1->end2 = sam_rec_ptr2->end;

	sam_rec_ptr1->strand_t2 = sam_rec_ptr2->strand_t;

	sam_rec_ptr1->mappedlen1 = sam_rec_ptr1->mappedlen;

	sam_rec_ptr1->mappedlen2 = sam_rec_ptr2->mappedlen;

	sam_rec_ptr1->mappedlen = sam_rec_ptr1->mappedlen1 + sam_rec_ptr1->mappedlen2;

	sam_rec_ptr1->insertlen1 = sam_rec_ptr1->insertlen;

	sam_rec_ptr1->insertlen2 = sam_rec_ptr2->insertlen;

	sam_rec_ptr1->insertlen = sam_rec_ptr1->insertlen1 + sam_rec_ptr1->insertlen2;

	sam_rec_ptr1->fusion_prefix_len = sam_rec_ptr1->mappedlen1;

	sam_rec_ptr1->fusion_suffix_len = sam_rec_ptr1->mappedlen2;

	//cout << "sam_rec_ptr1->fusion_prefix_len:"<<sam_rec_ptr1->fusion_prefix_len<<endl;

	//cout << "sam_rec_ptr1->fusion_suffix_len:"<<sam_rec_ptr1->fusion_suffix_len<<endl;

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

	if (true || sam_rec_ptr1->chrom_name < sam_rec_ptr1->chrom_name2 || (sam_rec_ptr1->chrom_name == sam_rec_ptr1->chrom_name2 && sam_rec_ptr1->start <= sam_rec_ptr1->start2))
	{
	}
	else
	{
		//sam_rec_ptr1->Swap();
	}

	return true;
}

bool comp_sam_rec(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.zfzfusstr < rhs.zfzfusstr;

}

void fusiostd2oneline(pair<vector<SamRec>, vector<SamRec> >& stored_tag_aligns_pe, pair<vector<SamRec*>, vector<SamRec*> >& merged_tag_aligns_pe)
{
	vector<SamRec>::iterator align_iter;

	if (stored_tag_aligns_pe.first.size() > 2)
	{
		for (size_t i = 0; i < stored_tag_aligns_pe.first.size(); ++i)
		{
			cout << stored_tag_aligns_pe.first[i].tostring(1, 1)<<endl;
		}

		stored_tag_aligns_pe.first.clear();
	}
		//sort(stored_tag_aligns_pe.first.begin(), stored_tag_aligns_pe.first.end(), comp_sam_rec);

	if (stored_tag_aligns_pe.second.size() > 2)
	{
		for (size_t i = 0; i < stored_tag_aligns_pe.second.size(); ++i)
		{
			cout << stored_tag_aligns_pe.second[i].tostring(1, 1)<<endl;
		}

		stored_tag_aligns_pe.second.clear();
	}
		//sort(stored_tag_aligns_pe.second.begin(), stored_tag_aligns_pe.second.end(), comp_sam_rec);

	for (align_iter = stored_tag_aligns_pe.first.begin(); align_iter != stored_tag_aligns_pe.first.end(); ++align_iter)
	{
		if (align_iter->is_fusion_newfmt)
		{
			if (align_iter + 1 == stored_tag_aligns_pe.first.end() || (align_iter + 1)->is_fusion_newfmt == false)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			if (align_iter->isdoner == true && (align_iter + 1)->isacceptor == true)
			{
				MergeStdFusion(align_iter, align_iter + 1);

				merged_tag_aligns_pe.first.push_back(&(*align_iter));
			}
			else if (align_iter->isacceptor == true && (align_iter + 1)->isdoner == true)
			{
				MergeStdFusion(align_iter + 1, align_iter);

				merged_tag_aligns_pe.first.push_back(&(*(align_iter+1)));
			}
			else
			{
				cerr << "no matched doenr patten"<<endl;

				cerr << align_iter->tag_name<<endl;
			}

			++align_iter;
		}
		else
		{
			merged_tag_aligns_pe.first.push_back(&(*align_iter));
		}
	}

	for (align_iter = stored_tag_aligns_pe.second.begin(); align_iter != stored_tag_aligns_pe.second.end(); ++align_iter)
	{
		if (align_iter->is_fusion_newfmt)
		{
			if (align_iter + 1 == stored_tag_aligns_pe.second.end() || (align_iter + 1)->is_fusion_newfmt == false)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			if (align_iter->isdoner == true && (align_iter + 1)->isacceptor == true)
			{
				MergeStdFusion(align_iter, align_iter + 1);

				merged_tag_aligns_pe.second.push_back(&(*align_iter));
			}
			else if (align_iter->isacceptor == true && (align_iter + 1)->isdoner == true)
			{
				MergeStdFusion(align_iter + 1, align_iter);

				merged_tag_aligns_pe.second.push_back(&(*(align_iter + 1)));
			}
			else
			{
				cerr << "no matched doenr patten"<<endl;

				cerr << align_iter->tag_name<<endl;
			}


			++align_iter;
		}
		else
		{
			merged_tag_aligns_pe.second.push_back(&(*align_iter));
		}
	}
}

string reverse_jump_code(vector<pair<size_t, int> >& spliceway_vec)
{
	string jump_code;

	for (int i = spliceway_vec.size() - 1; i >= 0; --i)
	{
		char jump_code_chr[100];

		if (spliceway_vec[i].second > 0)
		{
			if ((i != spliceway_vec.size() - 1) && spliceway_vec[i + 1].second > 0)
			{
				int intron = spliceway_vec[i + 1].first - spliceway_vec[i].first - spliceway_vec[i].second;

				sprintf(jump_code_chr, "%dN%dM", intron, spliceway_vec[i].second);
			}
			else
			{
				sprintf(jump_code_chr, "%dM", spliceway_vec[i].second);
			}
		}
		else
		{
			int insertlen = -spliceway_vec[i].second;

			sprintf(jump_code_chr, "%dI", insertlen);
		}

		jump_code.append(jump_code_chr);
	}

	return jump_code;
}

void set_mate_pair_distance(SamRec& samrec1, SamRec& samrec2)
{
	int st1 = samrec1.start, end1 = samrec1.end, st2 = samrec2.start, end2 = samrec2.end;
	
	int mate_dist1, mate_dist2;

	mate_dist1 = ((st1 - end2));

	mate_dist2 = ((end1 - st2));

	int outter_dist;

	if (abs(mate_dist1) < abs(mate_dist2))
	{
		if (mate_dist2 > 0)
			outter_dist = mate_dist2 + 1;
		else
			outter_dist = mate_dist2 - 1;
	}
	else
	{
		if (mate_dist1 > 0)
			outter_dist = mate_dist1 + 1;
		else
			outter_dist = mate_dist1 - 1;
	}
	
	samrec1.mate_offset = samrec2.start;

	samrec1.mate_diff = -outter_dist;

	samrec1.mate_match = "=";

	if (samrec2.strand_t & IS_REVERSE)
		samrec1.strand_t |= IS_MATE_REVERSE;

	samrec1.strand_t |= IS_PAIRED_MAPPED;

	samrec1.strand_t |= IS_PAIRED;

	samrec1.strand_t |= IS_FIRST_END;
	//

	samrec2.mate_offset = samrec1.start;

	samrec2.mate_diff = outter_dist;

	samrec2.mate_match = "=";

	if (samrec1.strand_t & IS_REVERSE)
		samrec2.strand_t |= IS_MATE_REVERSE;

	samrec2.strand_t |= IS_PAIRED_MAPPED;

	samrec2.strand_t |= IS_PAIRED;

	samrec2.strand_t |= IS_FIRST_END;
}

void extract_and_convert_sam(char* fusion_alignment_file, char* ordered_fusion_alignment_file, int stdfusion)
{
	ofstream ofs(ordered_fusion_alignment_file);

	//write header of sam file

	//@SQ     SN:chr20        LN:63025520

	vector<pair<string, size_t> >::iterator chriter;

	ReadNextTagAlignHandler<SamRec> readnextalignmenthandler(fusion_alignment_file, 5);

	pair<vector<SamRec>, vector<SamRec> > stored_tag_aligns_pe;

	size_t cur_read_id;

	while(readnextalignmenthandler.ReadNextTagAlignPE(stored_tag_aligns_pe, cur_read_id))
	{
		pair<vector<SamRec*>, vector<SamRec*> > merged_tag_aligns_pe;

		fusiostd2oneline(stored_tag_aligns_pe, merged_tag_aligns_pe);

		if (merged_tag_aligns_pe.first.size() > 1 || merged_tag_aligns_pe.second.size() > 1)//multiple alignments
		{
			cout << "multi alignments"<<endl;

			for (size_t i = 0; i < merged_tag_aligns_pe.first.size(); ++i)
			{
				cout << merged_tag_aligns_pe.first[i]->tostring(merged_tag_aligns_pe.first.size(), i+1) << endl;

				if (stdfusion && merged_tag_aligns_pe.first[i]->is_fusion)
					ofs << merged_tag_aligns_pe.first[i]->tostandfusion() << endl;
				else
					ofs << merged_tag_aligns_pe.first[i]->tostring(1, 1) << endl;
			}

			for (size_t i = 0; i < merged_tag_aligns_pe.second.size(); ++i)
			{
				cout << merged_tag_aligns_pe.second[i]->tostring(merged_tag_aligns_pe.second.size(), i+1) << endl;

				if (stdfusion && merged_tag_aligns_pe.second[i]->is_fusion)
					ofs << merged_tag_aligns_pe.second[i]->tostandfusion() << endl;
				else
					ofs << merged_tag_aligns_pe.second[i]->tostring(1, 1) << endl;
			}
			//if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			//{
			//	cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:multiple"<<endl;
			//}

			//cout << "multiple "<<endl;

			//if (merged_tag_aligns_pe.first.size() > 1)
			//	cout << merged_tag_aligns_pe.first.front()->tag_name << endl;
			//else if (merged_tag_aligns_pe.second.size() > 1)
			//	cout << merged_tag_aligns_pe.second.front()->tag_name << endl;
		}
		else if (merged_tag_aligns_pe.first.size() == 1 && merged_tag_aligns_pe.second.size() && 1)
		{
			if (stdfusion && merged_tag_aligns_pe.first[0]->is_fusion)
				ofs << merged_tag_aligns_pe.first[0]->tostandfusion() << endl;
			else
				ofs << merged_tag_aligns_pe.first[0]->tostring(1, 1) << endl;

			if (stdfusion && merged_tag_aligns_pe.second[0]->is_fusion)
				ofs << merged_tag_aligns_pe.second[0]->tostandfusion() << endl;
			else
				ofs << merged_tag_aligns_pe.second[0]->tostring(1, 1) << endl;

			//if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			//{
			//	cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:unique"<<endl;
			//}

			//vector<GeneGenomeLoc*> fitted_genome_locs1;

			//vector<GeneGenomeLoc*> fitted_genome_locs2;
			//
			//bool find_loc1 = find_corresponding_fusion_gene(merged_tag_aligns_pe.first.front(), genes_genome_locs, max_gene_len, fitted_genome_locs1);

			//bool find_loc2 = find_corresponding_fusion_gene(merged_tag_aligns_pe.second.front(), genes_genome_locs, max_gene_len, fitted_genome_locs2);

			//vector<GeneGenomeLoc*>::iterator genome_loc_iter1, genome_loc_iter2;

			//if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			//{
			//	cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:"<<find_loc1 <<'\t' <<  find_loc2<<endl;
			//}

			//if (find_loc1 && find_loc2)
			//{
			//	vector<pair<GeneGenomeLoc*, GeneGenomeLoc*> > valid_genome_loc_pairs;

			//	for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
			//	{
			//		for (genome_loc_iter2 = fitted_genome_locs2.begin(); genome_loc_iter2 != fitted_genome_locs2.end(); ++genome_loc_iter2)
			//		{
			//			if ((*genome_loc_iter1)->fusiongene_ptr == (*genome_loc_iter2)->fusiongene_ptr)//both ends fall in same fusion gene!
			//			{
			//				if (merged_tag_aligns_pe.first.size() && merged_tag_aligns_pe.first.front()->tag_name == "HWUSI-EAS1531_0001:4:100:1066:348#0/1")
			//				{
			//					cout << "HWUSI-EAS1531_0001:4:100:1066:348#0/1:paired matched:"<<find_loc1 <<'\t' <<  find_loc2<<endl;
			//				}

			//				valid_genome_loc_pairs.push_back(make_pair(*genome_loc_iter1, *genome_loc_iter2));
			//			}
			//		}
			//	}

			//	if (!valid_genome_loc_pairs.empty())
			//	{
			//		vector<pair<GeneGenomeLoc*, GeneGenomeLoc*> >::iterator genome_loc_pair_iter;

			//		for (genome_loc_pair_iter = valid_genome_loc_pairs.begin(); genome_loc_pair_iter != valid_genome_loc_pairs.end(); ++genome_loc_pair_iter)
			//		{
			//			SamRec samrecdup1 = *(merged_tag_aligns_pe.first.front());

			//			convert_alignment_coordinate(genome_loc_pair_iter->first, samrecdup1);

			//			SamRec samrecdup2 = *(merged_tag_aligns_pe.second.front());

			//			convert_alignment_coordinate(genome_loc_pair_iter->second, samrecdup2);

			//			set_mate_pair_distance(samrecdup1, samrecdup2);

			//			ofs << samrecdup1.tostring(1, 1) << endl;

			//			ofs << samrecdup2.tostring(1, 1) << endl;
			//		}
			//	}
			//	else
			//	{
			//	}
			//}
			//else if (find_loc1)
			//{
			//	for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
			//	{
			//		SamRec samrecdup1 = *(merged_tag_aligns_pe.first.front());

			//		convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

			//		ofs << samrecdup1.tostring(1, 1) << endl;
			//	}				
			//}
			//else if (find_loc2)
			//{
			//	for (genome_loc_iter1 = fitted_genome_locs2.begin(); genome_loc_iter1 != fitted_genome_locs2.end(); ++genome_loc_iter1)
			//	{
			//		SamRec samrecdup1 = *(merged_tag_aligns_pe.second.front());

			//		convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

			//		ofs << samrecdup1.tostring(1, 1) << endl;
			//	}
			//}
			//else
			//{
			//}
		}
		else if (merged_tag_aligns_pe.first.size() == 1)
		{
			if (stdfusion && merged_tag_aligns_pe.first[0]->is_fusion)
				ofs << merged_tag_aligns_pe.first[0]->tostandfusion() << endl;
			else
				ofs << merged_tag_aligns_pe.first[0]->tostring(1, 1) << endl;
			//vector<GeneGenomeLoc*> fitted_genome_locs1;
			//
			//bool find_loc1 = find_corresponding_fusion_gene(merged_tag_aligns_pe.first.front(), genes_genome_locs, max_gene_len, fitted_genome_locs1);

			//vector<GeneGenomeLoc*>::iterator genome_loc_iter1;

			//if (find_loc1)
			//{
			//	for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
			//	{
			//		SamRec samrecdup1 = *(merged_tag_aligns_pe.first.front());

			//		convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

			//		ofs << samrecdup1.tostring(1, 1) << endl;
			//	}
			//}
		}
		else if (merged_tag_aligns_pe.second.size() == 1)
		{
			if (stdfusion && merged_tag_aligns_pe.second[0]->is_fusion)
				ofs << merged_tag_aligns_pe.second[0]->tostandfusion() << endl;
			else
				ofs << merged_tag_aligns_pe.second[0]->tostring(1, 1) << endl;
			//vector<GeneGenomeLoc*> fitted_genome_locs1;
			//
			//bool find_loc1 = find_corresponding_fusion_gene(merged_tag_aligns_pe.second.front(), genes_genome_locs, max_gene_len, fitted_genome_locs1);

			//vector<GeneGenomeLoc*>::iterator genome_loc_iter1;

			//if (find_loc1)
			//{
			//	for (genome_loc_iter1 = fitted_genome_locs1.begin(); genome_loc_iter1 != fitted_genome_locs1.end(); ++genome_loc_iter1)
			//	{
			//		SamRec samrecdup1 = *(merged_tag_aligns_pe.second.front());

			//		convert_alignment_coordinate((*genome_loc_iter1), samrecdup1);

			//		ofs << samrecdup1.tostring(1, 1) << endl;
			//	}
			//}
		}
		else
		{
			//cout << "what"<<endl;
		}
	}
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "fusion_alignment_file converted_fusion_alignment_file" <<endl;
		exit(0);
	}

	char* fusion_alignment_file = argv[1];

	char* converted_fusion_alignment_file = argv[2];

	int stdfusion = atoi(argv[3]);
	
	extract_and_convert_sam(fusion_alignment_file, converted_fusion_alignment_file, stdfusion);
}