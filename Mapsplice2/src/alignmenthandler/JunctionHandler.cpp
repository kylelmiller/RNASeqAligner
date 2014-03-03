#include "JunctionHandler.h"
#include "ReadNextTagAlignHandler.h"

bool comp_junc_sort(const JunctionSeed* lhs, const JunctionSeed* rhs)
{
	if (lhs->m_start == rhs->m_start)
		return lhs->m_end < rhs->m_end;
	else
		return lhs->m_start < rhs->m_start;
}

bool comp_fusion_junc_sort(const FusionJuncRegion& lhs, const FusionJuncRegion& rhs)
{
	if (lhs.m_doner_st == rhs.m_doner_st)
	{
		if (lhs.m_doner_end == rhs.m_doner_end)
		{
			if (lhs.m_acceptor_st == rhs.m_acceptor_st)
			{
				return lhs.m_acceptor_end == rhs.m_acceptor_end;
			}
			else
				return lhs.m_acceptor_st < rhs.m_acceptor_st;
		}
		else
			return lhs.m_doner_end < rhs.m_doner_end;
	}
	else
		return lhs.m_doner_st < rhs.m_doner_st;
}

bool comp_junc_sort_by_end(const JunctionSeed* lhs, const JunctionSeed* rhs)
{
	if (lhs->m_end == rhs->m_end)
		return lhs->m_start < rhs->m_start;
	else
		return lhs->m_end < rhs->m_end;	
}

JunctionHandler::JunctionHandler(): m_canonical(0), m_semi_canonical(0), m_non_canonical(0), m_small_del(0), m_small_ins(0), m_min_mismatch(0), m_min_anchor(0)
{
}

bool JunctionHandler::Init(const vector<string>& alignment_files, size_t max_read_width, int insert_len, int delete_len, 
	const string& chrom_dir, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor)
{
	Clear();

	m_alignment_files = alignment_files;

	m_max_read_width = max_read_width;

	m_insert_len = insert_len;

	m_delete_len = delete_len;

	m_chrom_dir = chrom_dir;

	m_min_anchor = min_anchor;

	m_min_mismatch = min_mismatch;

	m_min_junc_anchor = min_junc_anchor;

	return true;
}

bool JunctionHandler::Clear()
{
	m_alignment_files.clear();

	m_head_line.clear();

	m_junc_files.clear();

	m_junc_hash.clear();

	m_max_read_width = 0;

	m_chrom_dir.clear();

	m_junc_sort.clear();

	m_canonical = 0;
	
	m_semi_canonical = 0;
	
	m_non_canonical = 0;
	
	m_small_del = 0;
	
	m_small_ins = 0;

	m_min_anchor = 0;

	m_min_mismatch = 0;

	m_min_junc_anchor = 0;

	return true;
}

size_t JunctionHandler::ReadJunction(vector<string>& junc_files)
{
	vector<string>::iterator junc_file_iter;

	size_t count = 0;

	for (junc_file_iter = junc_files.begin(); junc_file_iter != junc_files.end(); ++junc_file_iter)
	{
		m_junc_files.push_back(*junc_file_iter);

		ifstream ifs(junc_file_iter->c_str());

		if (ifs.is_open())
		{
			getline(ifs, m_head_line);

			while (!ifs.eof() )
			{
				string line;

				getline(ifs,line);

				if (line == "")
					continue;

				char chromname[100], juncname[100], strand, rgb[100], flankchr[10];

				unsigned int hits;

				unsigned short kinds, flankcase;

				size_t juncst, juncend, prefixend, suffixst;

				double rank, lpq, il;

				unsigned short min_mis, max_mis;

				double ave_mis;

				size_t max_prefix_len, max_suffix_len;

				size_t start_block_offset, end_block_offset;

				unsigned int unique_count, multi_count;

				unsigned int paired_count, single_count;

				unsigned int left_paired_count, right_paired_count;

				unsigned int paired_mutiple_count, paired_unique_count;

				unsigned short min_anchor_diff;

				sscanf(line.c_str(), "%s\t%llu\t%llu\t%s\t%u\t%c\t%llu\t%llu\t%s\t%hu\t%llu,%llu,\t%llu,%llu,\t%lf\t%hu\t%s\t%lf\t%lf\t%hu\t%hu\t%lf\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%hu", chromname,  &prefixend, &suffixst, 
					juncname, &hits, &strand, &juncst, &juncend, rgb, &kinds, &max_prefix_len, &max_suffix_len, &start_block_offset, &end_block_offset, 
					&rank, &flankcase, flankchr, &il, &lpq, &min_mis, &max_mis, &ave_mis, &unique_count, &multi_count, &paired_count, &left_paired_count, &right_paired_count,
					&paired_mutiple_count, &paired_unique_count, &single_count, &min_anchor_diff);

				size_t comb_offset = (prefixend << THIRTY_TWO) + suffixst;

				CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = m_junc_hash.find(chromname);

				if (chrom_junc_hash_iter == m_junc_hash.end())
				{
					JUNC_HASH_COMB junc_hash_comb;

					chrom_junc_hash_iter = (m_junc_hash.insert(CHROM_JUNC_HASH_COMB::value_type(chromname, junc_hash_comb))).first;							
				}

				JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

				JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(comb_offset);

				if (junc_hash_comb_iter != junc_hash_comb.end())
				{
					junc_hash_comb_iter->second.m_hits += hits;

					if (max_prefix_len > junc_hash_comb_iter->second.m_max_prefix_len) 
						junc_hash_comb_iter->second.m_max_prefix_len = max_prefix_len;

					if (max_suffix_len > junc_hash_comb_iter->second.m_max_suffix_len) 
						junc_hash_comb_iter->second.m_max_suffix_len = max_suffix_len;
				}
				else
				{
					chrom_junc_hash_iter->second.insert(JUNC_HASH_COMB::value_type(comb_offset, JunctionSeed(juncname, hits, 
						strand, kinds, max_prefix_len, max_suffix_len, start_block_offset, end_block_offset, rank, flankcase, flankchr, il, lpq, min_mis, max_mis, ave_mis, prefixend, 
						suffixst, chromname, unique_count, multi_count, paired_count, left_paired_count, right_paired_count, paired_mutiple_count, paired_unique_count, single_count, min_anchor_diff)));
				}

				++count;
			}

			ifs.close();
		}
		else
			cerr<<"can't open file :"<< *junc_file_iter << endl;

	}

	return count;
}

void JunctionHandler::UpdateJuncTable(const string& chrom_name, size_t combined_offset, size_t prefixlen, size_t suffixlen, size_t start, size_t end, 
	const string& ins_str, const SamRec& samrec, size_t sam_count, bool is_fusion)
{
	if ((start == 98658330 || end == 98658330) && chrom_name.find("~") != string::npos)
		cout << samrec.tostring(0, 0) << endl;

	//string chrom_name;

	//if (is_fusion)
	//	chrom_name = samrec.chrom_name + "~" + samrec.chrom_name2;
	//else
	//	chrom_name = samrec.chrom_name;

	CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = m_junc_hash.find(chrom_name);

	if (chrom_junc_hash_iter == m_junc_hash.end())
	{
		JUNC_HASH_COMB junc_hash_comb;

		chrom_junc_hash_iter = (m_junc_hash.insert(CHROM_JUNC_HASH_COMB::value_type(chrom_name, junc_hash_comb))).first;							
	}

	JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

	JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

	if (junc_hash_comb_iter != junc_hash_comb.end())
	{
		if (samrec.mate_match != "*" && samrec.mate_diff == 0)
			cout << samrec.tostring(1, 1) << endl;
		junc_hash_comb_iter->second.inc_hits(prefixlen, suffixlen, samrec.tagidx, samrec.mis_match, samrec.strand_t, sam_count, samrec.mate_match, 
			samrec.mate_diff, samrec.left_splice_ways, samrec.right_splice_ways,  ins_str);
	}
	else
	{
		chrom_junc_hash_iter->second.insert(JUNC_HASH_COMB::value_type(combined_offset, JunctionSeed(/*prim, *//*flankseq, */
			prefixlen, suffixlen, m_max_read_width, samrec.tagidx, samrec.mis_match, samrec.strand_t, samrec.strand_t2, start, end, samrec.chrom_name, samrec.chrom_name2, 
			sam_count, samrec.mate_match, samrec.mate_diff, samrec.left_splice_ways, samrec.right_splice_ways, is_fusion, ins_str)));
	}

	//size_t junc_seed_size = sizeof(unsigned short) * junc_hash_comb.find(combined_offset)->second.m_prefix_count.size();
}


void JunctionHandler::SamRec2Junc(SamRec& samrec, size_t sam_count, hash_map<string, hash_map<size_t, JuncSeedSp> >& cur_read_junc)
{
	vector<pair<size_t, int> >& spliceway_vec = samrec.spliceway_vec;

	vector<pair<size_t, int> >& spliceway_vec2 = samrec.spliceway_vec2;

	const string& readstr = samrec.mapped_seq;

	size_t mappedlen = 0;

	if (spliceway_vec.size() > 1)
	{
		vector<pair<size_t, int> >::iterator vp_iter;

		for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
		{
			size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

			string ins_str;

			if (vp_iter->second < 0)
			{
				if (vp_iter == spliceway_vec.begin() || (vp_iter - 1)->second  < 0)
					prefixlen = 1;
				else
					prefixlen = (vp_iter - 1)->second + 1;

				if (vp_iter == spliceway_vec.end() - 1 || (vp_iter + 1)->second  < 0)
					suffixlen = 1;
				else
					suffixlen = (vp_iter + 1)->second + 1;

				ins_str = readstr.substr(mappedlen, -(vp_iter->second));

				prefixend = vp_iter->first - 1;

				suffixst = prefixend;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;
			}
			else if (vp_iter == spliceway_vec.end() - 1 || (vp_iter + 1)->second < 0)
			{
				mappedlen += abs(vp_iter->second);

				continue;
			}
			else
			{
				prefixst = vp_iter->first;

				prefixend = vp_iter->first + vp_iter->second - 1;

				suffixst = (vp_iter + 1)->first;

				prefixlen = vp_iter->second;

				suffixlen = (vp_iter + 1)->second;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;
			}

			mappedlen += abs(vp_iter->second);

			//string chrom_name;

			//if (is_fusion)
			//	chrom_name = samrec.chrom_name + samrec.chrom_name2;
			//else
			//	chrom_name = samrec.chrom_name;

			if (prefixlen >= m_min_anchor && suffixlen >= m_min_anchor && samrec.mis_match < m_min_mismatch)
			{
				hash_map<string, hash_map<size_t, JuncSeedSp> >::iterator chrom_junc_seed_iter = cur_read_junc.find(samrec.chrom_name);

				if (chrom_junc_seed_iter == cur_read_junc.end())
				{
					hash_map<size_t, JuncSeedSp> junc_seed_comb;

					chrom_junc_seed_iter = (cur_read_junc.insert(hash_map<string, hash_map<size_t, JuncSeedSp> >::value_type(samrec.chrom_name, junc_seed_comb))).first;							
				}

				hash_map<size_t, JuncSeedSp>& junc_hash_comb = chrom_junc_seed_iter->second;

				hash_map<size_t, JuncSeedSp>::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

				if (junc_hash_comb_iter != junc_hash_comb.end())
				{
					//remove same read map to same junction multiple times
					if (junc_hash_comb_iter->second.m_sam_rec_ptr->mate_match == "*" && samrec.mate_match != "*")
						junc_hash_comb_iter->second.set(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false);
					else if (junc_hash_comb_iter->second.m_sam_rec_ptr->mate_match == samrec.mate_match)
					{
						if (junc_hash_comb_iter->second.m_sam_rec_ptr->mis_match > samrec.mis_match)
							junc_hash_comb_iter->second.set(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false);
					}
				}
				else
				{
					junc_hash_comb.insert(hash_map<size_t, JuncSeedSp>::value_type(combined_offset, JuncSeedSp(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false)));

				}
			}


			//if (cur_read_junc.find(combined_offset) == cur_read_junc.end())
			//	cur_read_junc.

			//if (prefixlen >= m_min_anchor && suffixlen >= m_min_anchor && samrec.mis_match < m_min_mismatch)
			//	UpdateJuncTable(combined_offset, prefixlen, suffixlen, prefixend, suffixst, ins_str, samrec, sam_count);
		}
	}

	if (spliceway_vec2.size() > 1)
	{
		vector<pair<size_t, int> >::iterator vp_iter;

		for (vp_iter =  spliceway_vec2.begin(); vp_iter != spliceway_vec2.end(); ++vp_iter)
		{
			size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

			string ins_str;

			if (vp_iter->second < 0)
			{
				if (vp_iter == spliceway_vec2.begin() || (vp_iter - 1)->second  < 0)
					prefixlen = 1;
				else
					prefixlen = (vp_iter - 1)->second + 1;

				if (vp_iter == spliceway_vec2.end() - 1 || (vp_iter + 1)->second  < 0)
					suffixlen = 1;
				else
					suffixlen = (vp_iter + 1)->second + 1;

				ins_str = readstr.substr(mappedlen, -(vp_iter->second));

				prefixend = vp_iter->first - 1;

				suffixst = prefixend;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;
			}
			else if (vp_iter == spliceway_vec2.end() - 1 || (vp_iter + 1)->second < 0)
			{
				mappedlen += abs(vp_iter->second);

				continue;
			}
			else
			{
				prefixst = vp_iter->first;

				prefixend = vp_iter->first + vp_iter->second - 1;

				suffixst = (vp_iter + 1)->first;

				prefixlen = vp_iter->second;

				suffixlen = (vp_iter + 1)->second;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;
			}

			mappedlen += abs(vp_iter->second);

			if (prefixlen >= m_min_anchor && suffixlen >= m_min_anchor && samrec.mis_match < m_min_mismatch)
			{
				hash_map<string, hash_map<size_t, JuncSeedSp> >::iterator chrom_junc_seed_iter = cur_read_junc.find(samrec.chrom_name2);

				if (chrom_junc_seed_iter == cur_read_junc.end())
				{
					hash_map<size_t, JuncSeedSp> junc_seed_comb;

					chrom_junc_seed_iter = (cur_read_junc.insert(hash_map<string, hash_map<size_t, JuncSeedSp> >::value_type(samrec.chrom_name2, junc_seed_comb))).first;							
				}

				hash_map<size_t, JuncSeedSp>& junc_hash_comb = chrom_junc_seed_iter->second;

				hash_map<size_t, JuncSeedSp>::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

				if (junc_hash_comb_iter != junc_hash_comb.end())
				{
					//remove same read map to same junction multiple times
					if (junc_hash_comb_iter->second.m_sam_rec_ptr->mate_match == "*" && samrec.mate_match != "*")
						junc_hash_comb_iter->second.set(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false);
					else if (junc_hash_comb_iter->second.m_sam_rec_ptr->mate_match == samrec.mate_match)
					{
						if (junc_hash_comb_iter->second.m_sam_rec_ptr->mis_match > samrec.mis_match)
							junc_hash_comb_iter->second.set(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false);
					}
				}
				else
				{
					junc_hash_comb.insert(hash_map<size_t, JuncSeedSp>::value_type(combined_offset, JuncSeedSp(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false)));

				}
			}
				//UpdateJuncTable(combined_offset, prefixlen, suffixlen, prefixend, suffixst, ins_str, samrec, sam_count);
		}
	}

	if (samrec.is_fusion)
	{
		size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

		string ins_str;

		prefixst = samrec.fusion_prefix_st;

		prefixend = samrec.fusion_prefix_end;

		prefixlen = samrec.fusion_prefix_len;//spliceway_vec.back().second;

		suffixst = samrec.fusion_suffix_st;

		suffixlen = samrec.fusion_suffix_len;//spliceway_vec2.front().second;

		if (prefixend == 98658330 || suffixst == 98658330)
			cout << samrec.tostring(0, 0) << endl;

		combined_offset = (prefixend << THIRTY_TWO) + suffixst;

		if (prefixlen >= m_min_anchor && suffixlen >= m_min_anchor && samrec.mis_match < m_min_mismatch)
		{
			string chrom_name;

			chrom_name = samrec.chrom_name + "~" + samrec.chrom_name2;

			hash_map<string, hash_map<size_t, JuncSeedSp> >::iterator chrom_junc_seed_iter = cur_read_junc.find(chrom_name);

			if (chrom_junc_seed_iter == cur_read_junc.end())
			{
				hash_map<size_t, JuncSeedSp> junc_seed_comb;

				chrom_junc_seed_iter = (cur_read_junc.insert(hash_map<string, hash_map<size_t, JuncSeedSp> >::value_type(chrom_name, junc_seed_comb))).first;							
			}

			hash_map<size_t, JuncSeedSp>& junc_hash_comb = chrom_junc_seed_iter->second;

			hash_map<size_t, JuncSeedSp>::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

			if (junc_hash_comb_iter != junc_hash_comb.end())
			{
				//remove same read map to same junction multiple times
				if (junc_hash_comb_iter->second.m_sam_rec_ptr->mate_match == "*" && samrec.mate_match != "*")
					junc_hash_comb_iter->second.set(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, true);
				else if (junc_hash_comb_iter->second.m_sam_rec_ptr->mate_match == samrec.mate_match)
				{
					if (junc_hash_comb_iter->second.m_sam_rec_ptr->mis_match > samrec.mis_match)
						junc_hash_comb_iter->second.set(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, true);
				}
			}
			else
			{
				junc_hash_comb.insert(hash_map<size_t, JuncSeedSp>::value_type(combined_offset, JuncSeedSp(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, true)));
			}
		}
			//UpdateJuncTable(combined_offset, prefixlen, suffixlen, prefixend, suffixst, ins_str, samrec, sam_count, true);
	}
}

void JunctionHandler::SamRecVec2Junc(vector<SamRec*>& samrecs)
{
	size_t sam_count = samrecs.size();

	vector<SamRec*>::iterator sam_rec_iter;

	size_t small_anchor_count = 0;

	for (sam_rec_iter = samrecs.begin(); sam_rec_iter != samrecs.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->min_anchor < m_min_anchor)
			++small_anchor_count;
	}

	size_t real_sam_count = sam_count - small_anchor_count;

	//vector<SamRec>::iterator sam_rec_iter;

	hash_map<string, hash_map<size_t, JuncSeedSp> > cur_read_junc;

	for (sam_rec_iter = samrecs.begin(); sam_rec_iter != samrecs.end(); ++sam_rec_iter)
	{
		if (((*sam_rec_iter)->isspliced || (*sam_rec_iter)->issmallins || (*sam_rec_iter)->is_fusion) && real_sam_count)
			SamRec2Junc(**sam_rec_iter, real_sam_count, cur_read_junc);
	}

	hash_map<string, hash_map<size_t, JuncSeedSp> >::iterator chrom_iter;

	for (chrom_iter = cur_read_junc.begin(); chrom_iter != cur_read_junc.end(); ++chrom_iter)
	{
		hash_map<size_t, JuncSeedSp>::iterator offset_iter;

		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{
			UpdateJuncTable(chrom_iter->first, offset_iter->first, offset_iter->second.m_max_prefix_len, offset_iter->second.m_max_suffix_len, offset_iter->second.m_start,
				offset_iter->second.m_end, offset_iter->second.m_ins_str, *(offset_iter->second.m_sam_rec_ptr), sam_count, offset_iter->second.m_is_fusion);
				//prefixlen, suffixlen, prefixend, suffixst, ins_str, samrec, sam_count);
		}
	}
}


void JunctionHandler::ReadSam()
{
	vector<string>::iterator alignment_file_iter;

	for (alignment_file_iter = m_alignment_files.begin(); alignment_file_iter != m_alignment_files.end(); ++alignment_file_iter)
	{
		vector<SamRec> sam_rec;

		size_t read_id;

		ReadNextTagAlignHandler<SamRec> sam_file_handler(alignment_file_iter->c_str(), m_insert_len);

		size_t sam_count = sam_file_handler.ReadNextTagAlign(sam_rec, read_id);

		vector<SamRec>::iterator sam_rec_iter;

		size_t small_anchor_count = 0;

		for (sam_rec_iter = sam_rec.begin(); sam_rec_iter != sam_rec.end(); ++sam_rec_iter)
		{
			if (sam_rec_iter->min_anchor < m_min_anchor)
				++small_anchor_count;
		}

		size_t real_sam_count = sam_count - small_anchor_count;

		while (true)
		{
			// if it is end of file then stop
			if (sam_count == 0)
				break;

			vector<SamRec>::iterator sam_rec_iter;

			hash_map<string, hash_map<size_t, JuncSeedSp> > cur_read_junc;

			for (sam_rec_iter = sam_rec.begin(); sam_rec_iter != sam_rec.end(); ++sam_rec_iter)
			{
				if ((sam_rec_iter->isspliced || sam_rec_iter->issmallins) && real_sam_count)
					SamRec2Junc(*sam_rec_iter, real_sam_count, cur_read_junc);
			}

			sam_count = sam_file_handler.ReadNextTagAlign(sam_rec, read_id);

			small_anchor_count = 0;

			for (sam_rec_iter = sam_rec.begin(); sam_rec_iter != sam_rec.end(); ++sam_rec_iter)
			{
				if (sam_rec_iter->min_anchor < m_min_anchor)
					++small_anchor_count;
			}

			real_sam_count = sam_count - small_anchor_count;
		}
	}

	LoadFlankString();
}

void JunctionHandler::LoadFlankString()
{
	CHROM_JUNC_HASH_COMB::iterator chm_iter;

	for (chm_iter = m_junc_hash.begin(); chm_iter != m_junc_hash.end(); ++chm_iter)
	{
		//if (m_junc_sort.find(chm_iter->first) == m_junc_sort.end())
		//	m_junc_sort[chm_iter->first];

		//vector<JunctionSeed* >& junc_sort_vec = m_junc_sort.find(chm_iter->first)->second;

		cout << chm_iter->first << endl;

		if (chm_iter->first.find("~") == string::npos)
		{
			string chromfile = m_chrom_dir + chm_iter->first;

			chromfile.append(".fa");

			string chromseq;

			readchrom(chromfile.c_str(), chromseq);

			if (chromseq.empty())
			{
				cout <<"empty chrom: "<<chromfile<<endl;
				exit(1);
			}

			size_t chrom_size = chromseq.size() - 1;

			JUNC_HASH_COMB::iterator iter_conj;

			for (iter_conj = chm_iter->second.begin(); iter_conj != chm_iter->second.end(); ++iter_conj)
			{
				size_t comb_offset = iter_conj->first;

				size_t prefix_end = comb_offset >> THIRTY_TWO;

				size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

				cout << prefix_end << '\t' << suffix_st << endl;

				iter_conj->second.set_coverage();

				iter_conj->second.set_entropy();

				string flankstr = chromseq.substr(prefix_end, 2) + chromseq.substr(suffix_st - 3, 2);

				for (size_t i = 0; i < flankstr.length(); ++i)
				{
					if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
						flankstr[i] = flankstr[i] + 'A' - 'a';
				}

				iter_conj->second.set_flankstring(flankstr);

				size_t intron_len = suffix_st - prefix_end - 1;

				iter_conj->second.set_pq_score(intron_len, chrom_size);

				iter_conj->second.set_il_score(prefix_end, suffix_st);

				iter_conj->second.set_ave_mis();

				iter_conj->second.set_block_offset(prefix_end, suffix_st);

				//junc_sort_vec.push_back(&iter_conj->second);

			}
		}
		else
		{
			char chr1[1000], chr2[1000];

			sscanf(chm_iter->first.c_str(), "%[^~]~%[^~]", chr1, chr2);

			string chromfile1 = m_chrom_dir + chr1;

			chromfile1.append(".fa");

			string chromseq1;

			readchrom(chromfile1.c_str(), chromseq1);

			if (chromseq1.empty())
			{
				cout <<"empty chrom: "<<chromfile1<<endl;
				exit(1);
			}

			string chromfile2 = m_chrom_dir + chr2;

			chromfile2.append(".fa");

			string chromseq2;

			readchrom(chromfile2.c_str(), chromseq2);

			if (chromseq2.empty())
			{
				cout <<"empty chrom: "<<chromfile2<<endl;
				exit(1);
			}

			size_t chrom_size1 = chromseq1.size() - 1;

			size_t chrom_size2 = chromseq2.size() - 1;

			JUNC_HASH_COMB::iterator iter_conj;

			for (iter_conj = chm_iter->second.begin(); iter_conj != chm_iter->second.end(); ++iter_conj)
			{
				size_t comb_offset = iter_conj->first;

				size_t prefix_end = comb_offset >> THIRTY_TWO;

				size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

				iter_conj->second.set_coverage();

				iter_conj->second.set_entropy();

				cout << prefix_end << '\t' << suffix_st << endl;
				//if (

				string flankstr1, flankstr2;

				if (iter_conj->second.m_strand1 == '+')
					flankstr1 = chromseq1.substr(prefix_end/* + 1*/, 2);
				else
				{
					flankstr1 = chromseq1.substr(prefix_end - 3/* - 2*/, 2);
					flankstr1 = revcomp(flankstr1);
				}

				if (iter_conj->second.m_strand2 == '+')
					flankstr2 = chromseq2.substr(suffix_st - 3/* - 2*/, 2);
				else
				{
					flankstr2 = chromseq2.substr(suffix_st/* + 1*/, 2);
					flankstr2 = revcomp(flankstr2);
				}


				string flankstr = flankstr1 + flankstr2; //chromseq1.substr(prefix_end, 2) + chromseq2.substr(suffix_st - 3, 2);

				for (size_t i = 0; i < flankstr.length(); ++i)
				{
					if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
						flankstr[i] = flankstr[i] + 'A' - 'a';
				}

				iter_conj->second.set_flankstring(flankstr);

				//size_t intron_len = suffix_st - prefix_end - 1;

				//iter_conj->second.set_pq_score(intron_len, chrom_size);

				//iter_conj->second.set_il_score(prefix_end, suffix_st);

				iter_conj->second.set_ave_mis();

				iter_conj->second.set_block_offset(prefix_end, suffix_st);

				//junc_sort_vec.push_back(&iter_conj->second);

			}

			//char chrom1[1000], chrom2[1000];

			//sscanf(chm_iter->first.c_str(), "%[^_]~%[^_]", chrom1, chrom2);
		}
	}
}

void JunctionHandler::SortJunc()
{
	map<string, vector<JunctionSeed*> >::iterator chrom_iter;

	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
	}

}

void JunctionHandler::SortFusionJunc()
{
	map<string, vector<FusionJuncRegion> >::iterator chrom_iter;

	for (chrom_iter = m_fuson_junc_sort.begin(); chrom_iter != m_fuson_junc_sort.end(); ++chrom_iter)
	{
		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_fusion_junc_sort);
	}
}

void JunctionHandler::SortJuncByEnd()
{
	map<string, vector<JunctionSeed*> >::iterator chrom_iter;
	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort_by_end);
	}
}

void JunctionHandler::LoadJuncToSortVec()
{
	m_junc_sort.clear();

	CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter;

	JUNC_HASH_COMB::iterator junc_hash_comb_iter;

	for (chrom_junc_hash_iter = m_junc_hash.begin(); chrom_junc_hash_iter != m_junc_hash.end(); ++chrom_junc_hash_iter)
	{
		for (junc_hash_comb_iter = chrom_junc_hash_iter->second.begin(); junc_hash_comb_iter != chrom_junc_hash_iter->second.end(); ++junc_hash_comb_iter)
		{
			m_junc_sort[chrom_junc_hash_iter->first].push_back(&(junc_hash_comb_iter->second));
		}
	}
}

void
JunctionHandler::LoadFusionJuncToSortVec()
{
	m_fuson_junc_sort.clear();

	CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter;

	JUNC_HASH_COMB::iterator junc_hash_comb_iter;

	for (chrom_junc_hash_iter = m_junc_hash.begin(); chrom_junc_hash_iter != m_junc_hash.end(); ++chrom_junc_hash_iter)
	{
		if (!chrom_junc_hash_iter->first.find("~") != string::npos)//.m_chrom2.empty())
		{
			for (junc_hash_comb_iter = chrom_junc_hash_iter->second.begin(); junc_hash_comb_iter != chrom_junc_hash_iter->second.end(); ++junc_hash_comb_iter)
			{
				m_fuson_junc_sort[chrom_junc_hash_iter->first].push_back(FusionJuncRegion(&(junc_hash_comb_iter->second)));
			}
		}
	}
}

void JunctionHandler::WriteJunction(string junction_file, string junction_ins_file, string junction_del_file, string junction_fusion_file, string filtered_junc)
{
	ofstream ofs_normal(junction_file.c_str());

	ofstream ofs_insert(junction_ins_file.c_str());

	ofstream ofs_deltion(junction_del_file.c_str());

	ofstream ofs_fusion(junction_fusion_file.c_str());

	ofstream ofs_filtered(filtered_junc.c_str());

	map<string, vector<JunctionSeed*> >::iterator chrom_iter;

	size_t normal_count = 0, ins_count = 0, del_count = 0, filtered_count = 0, fusion_count = 0;

	//if (m_junc_sort
	LoadJuncToSortVec();

	SortJunc();	

	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		//sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
		vector<JunctionSeed*>::iterator junc_sort_iter;

		string chrom = chrom_iter->first;

		vector<JunctionSeed*>& junc_sort_vec = chrom_iter->second;

		for (junc_sort_iter = junc_sort_vec.begin(); junc_sort_iter != junc_sort_vec.end(); ++junc_sort_iter)
		{
			if ((*junc_sort_iter)->m_is_fusion)
			{
				++fusion_count;

				ofs_fusion << (*junc_sort_iter)->to_normal_junction(fusion_count)<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_SMALL_ANCHOR)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_SMALL_ANCHOR"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_SMALL_DELETION)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_SMALL_DELETION"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_LARGE_MULTIPLE_PAIRED)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_LARGE_MULTIPLE_PAIRED"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_LARGE_MIN_ANCHOR_DIFF)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_LARGE_MIN_ANCHOR_DIFF"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_NOPAIRED)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_NOPAIRED"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_INSERTION)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_INSERTION"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_LARGE_MISMATCH)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_LARGE_MISMATCH"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_NONCAN_LEFT_RIGHT_PAIR)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_NONCAN_LEFT_RIGHT_PAIR"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_NONCAN_ENTROPY)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_NONCAN_ENTROPY"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_NONCAN_MULTI)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_NONCAN_MULTI"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_NONCAN_ERROR)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_NONCAN_ERROR"<<endl;
			}
			else if ((*junc_sort_iter)->m_start == (*junc_sort_iter)->m_end)
			{
				++ins_count;

				ofs_insert << (*junc_sort_iter)->to_insert_junction(ins_count)<<endl;

			}
			else if ( (*junc_sort_iter)->m_end - (*junc_sort_iter)->m_start < m_delete_len) 
			{
				++del_count;

				ofs_deltion << (*junc_sort_iter)->to_insert_junction(del_count)<<endl;
			}
			else
			{
				++normal_count;

				ofs_normal << (*junc_sort_iter)->to_normal_junction(normal_count)<<endl;
			}
		}
	}

	cout << "insert: "<< ins_count << endl;

	cout << "delete: "<< del_count << endl;

	cout << "normal: "<< normal_count << endl;
}


bool less_than(int lhs, int rhs)
{
	return abs(lhs) < abs(rhs);

}
void JunctionHandler::WriteFusionJunctionWithEncompassingAlignments(string junction_fusion_file)
{
	ofstream ofs_fusion(junction_fusion_file.c_str());

	string junction_fusion_file_no_compass = junction_fusion_file; junction_fusion_file_no_compass.append(".no.compass");

	ofstream ofs_fusion_no_compass(junction_fusion_file_no_compass.c_str());

	map<string, vector<JunctionSeed*> >::iterator chrom_iter;

	size_t fusion_count = 0;

	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		//sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
		vector<JunctionSeed*>::iterator junc_sort_iter;

		string chrom = chrom_iter->first;

		vector<JunctionSeed*>& junc_sort_vec = chrom_iter->second;

		for (junc_sort_iter = junc_sort_vec.begin(); junc_sort_iter != junc_sort_vec.end(); ++junc_sort_iter)
		{
			if ((*junc_sort_iter)->m_is_fusion)
			{
				++fusion_count;

				//if ((*junc_sort_iter)->m_fusion_encompassing_reads.size())
				//{
				//	(*junc_sort_iter)->m_paired_count += (*junc_sort_iter)->m_fusion_encompassing_reads.size();
				//}

				ofs_fusion << (*junc_sort_iter)->to_normal_junction(fusion_count)<<endl;

				ofs_fusion_no_compass << (*junc_sort_iter)->to_normal_junction(fusion_count)<<'\t';//<<endl;

				
				for (size_t i = 0; i < (*junc_sort_iter)->left_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					vector<int> cur_path = (*junc_sort_iter)->left_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					ofs_fusion_no_compass << (*junc_sort_iter)->left_exons[abs(cur_path[0]) - 1].first <<",";

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						ofs_fusion_no_compass << (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1<<'M';
						
						if (j + 1 < cur_path.size())
						{
							if (cur_path[j] < 0 && cur_path[j+1] < 0)
								ofs_fusion_no_compass <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'P';
							else
								ofs_fusion_no_compass <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'N';
						}

					}

					ofs_fusion_no_compass << '|';
				}

				ofs_fusion_no_compass << '\t';

				for (size_t i = 0; i < (*junc_sort_iter)->right_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					//for (size_t j = 0; j < (*junc_sort_iter)->right_paths[i].size(); ++j)
					//	ofs_fusion_no_compass << (*junc_sort_iter)->right_exons[(*junc_sort_iter)->right_paths[i][j]].first<<':' << 
					//	(*junc_sort_iter)->right_exons[(*junc_sort_iter)->right_paths[i][j]].second << '=';

					//ofs_fusion_no_compass << '~';


					vector<int> cur_path = (*junc_sort_iter)->right_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					ofs_fusion_no_compass << (*junc_sort_iter)->right_exons[abs(cur_path[0]) - 1].first <<",";

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						ofs_fusion_no_compass << (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1<<'M';
						
						if (j + 1 < cur_path.size())
						{
							if (cur_path[j] < 0 && cur_path[j+1] < 0)
								ofs_fusion_no_compass <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'P';
							else
								ofs_fusion_no_compass <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'N';
						}

					}

					ofs_fusion_no_compass << '|';
				}

				ofs_fusion_no_compass << endl;;


				for (size_t i = 0; i < (*junc_sort_iter)->m_fusion_encompassing_reads.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->m_fusion_encompassing_reads[i]<<endl;
				}

				ofs_fusion << "left splice ways" << endl;

				for (size_t i = 0; i < (*junc_sort_iter)->left_splice_ways.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < (*junc_sort_iter)->left_splice_ways[i].spliceway_vec.size(); ++j)
						ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].spliceway_vec[j].first<<':' << (*junc_sort_iter)->left_splice_ways[i].spliceway_vec[j].second << '\t';

					ofs_fusion << endl;
				}

				ofs_fusion << "left exons" << endl;

				for (size_t i = 0; i < (*junc_sort_iter)->left_exons.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					ofs_fusion << (*junc_sort_iter)->left_exons[i].first<<':' << 
						(*junc_sort_iter)->left_exons[i].second << '\t';

				}

				ofs_fusion << endl;

				ofs_fusion << "left structure" << endl;

				for (size_t i = 0; i < (*junc_sort_iter)->left_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < (*junc_sort_iter)->left_paths[i].size(); ++j)
						ofs_fusion << (*junc_sort_iter)->left_exons[abs((*junc_sort_iter)->left_paths[i][j]) - 1].first<<':' << 
						(*junc_sort_iter)->left_exons[abs((*junc_sort_iter)->left_paths[i][j]) - 1].second << '\t';

					ofs_fusion << endl;
				}

				ofs_fusion << "right splice ways" << endl;

				for (size_t i = 0; i < (*junc_sort_iter)->right_splice_ways.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->right_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->right_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < (*junc_sort_iter)->right_splice_ways[i].spliceway_vec.size(); ++j)
						ofs_fusion << (*junc_sort_iter)->right_splice_ways[i].spliceway_vec[j].first<<':' << (*junc_sort_iter)->right_splice_ways[i].spliceway_vec[j].second << '\t';

					ofs_fusion << endl;
				}

				ofs_fusion << "right exons" << endl;

				for (size_t i = 0; i < (*junc_sort_iter)->right_exons.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					ofs_fusion << (*junc_sort_iter)->right_exons[i].first<<':' << 
						(*junc_sort_iter)->right_exons[i].second << '\t';

				}

				ofs_fusion << endl;

				ofs_fusion << "right structure" << endl;

				for (size_t i = 0; i < (*junc_sort_iter)->right_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < (*junc_sort_iter)->right_paths[i].size(); ++j)
						ofs_fusion << (*junc_sort_iter)->right_exons[abs((*junc_sort_iter)->right_paths[i][j]) - 1].first<<':' << 
						(*junc_sort_iter)->right_exons[abs((*junc_sort_iter)->right_paths[i][j]) - 1].second << '\t';

					ofs_fusion << endl;
				}
			}
		}
	}
}

vector<vector<int> >* graph_ptr;
vector<size_t> DFS_stack;
vector<int> DFS_in_stack;
vector<vector<int> >* stored_path_ptr;

extern void DFS_VISIT(size_t u, int path_type);

void DFS(vector<vector<int> >& graph, vector<vector<int> >& stored_path, size_t u)
{
	DFS_stack.clear();
	DFS_in_stack.clear();
	graph_ptr = &graph;
	stored_path_ptr = &stored_path;

	DFS_in_stack.resize(graph.size(), 0);

	DFS_VISIT(u, 1);
	//for (size_t i = 0; i < (*graph_ptr)[u].size(); ++i)
	//{
	//	//if ((*graph_ptr)[u][i] != 0)
	//}
}

void DFS_VISIT(size_t u, int path_type)
{
	DFS_in_stack[u] = path_type;

	DFS_stack.push_back(u);

	bool found_new_node = false;

	for (size_t i = 0; i < (*graph_ptr)[u].size(); ++i)
	{
		//adjacent node
		if ((*graph_ptr)[u][i] > 0)
		{
			if (DFS_in_stack[i] == 0)
			{
				found_new_node = true;

				DFS_VISIT(i, (*graph_ptr)[u][i]);
			}
		}
	}

	if (!found_new_node)//node can't find new node, print path
	{
		vector<int> cur_path;

		for(size_t i = 0; i < DFS_stack.size(); ++i)
		{
			if (DFS_in_stack[DFS_stack[i]] == 1)
				cur_path.push_back(static_cast<int> (DFS_stack[i] + 1));
			else if (DFS_in_stack[DFS_stack[i]] == 2)
			{
				if (cur_path.back() > 0)
					cur_path[cur_path.size() - 1] = -cur_path.back();

				cur_path.push_back(- (static_cast<int> (DFS_stack[i] + 1)));
			}
			else
				cout << "graph abnormal"<<endl;
			//cout << DFS_stack[i] <<'\t';
		}

		(*stored_path_ptr).push_back(cur_path);
		//cout << endl;
	}

	//pop stack
	DFS_in_stack[u] = 0;

	DFS_stack.pop_back();
}

void JunctionHandler::GenerateFusionStruct()
{
	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	JUNC_HASH_COMB::iterator offset_iter;

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		cout << 1 << endl;

		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{

			cout << 2 << endl;

			if (offset_iter->second.m_is_fusion)
			{
				cout << "junction "<<endl;

				cout << offset_iter->second.to_normal_junction(0)<<endl;

				cout << "left splice ways" << endl;

				for (size_t i = 0; i < offset_iter->second.left_splice_ways.size(); ++i)
				{
					cout << offset_iter->second.left_splice_ways[i].chrom_name<<'\t'<< offset_iter->second.left_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < offset_iter->second.left_splice_ways[i].spliceway_vec.size(); ++j)
						cout << offset_iter->second.left_splice_ways[i].spliceway_vec[j].first<<':' << offset_iter->second.left_splice_ways[i].spliceway_vec[j].second << '\t';

					cout << endl;
				}

				//cout << 
				if (offset_iter->second.left_splice_ways.size())
				{
					offset_iter->second.spliceways2exons_doner();

					vector<int> doner_array(offset_iter->second.left_exons.size(), -1);

					vector<vector<int> > doner_graph(offset_iter->second.left_exons.size(), doner_array);

					if (offset_iter->second.m_strand1 == '+')
					{
						offset_iter->second.construct_graph(offset_iter->second.left_exons, offset_iter->second.left_splice_ways, doner_graph, '-');

						DFS(doner_graph, offset_iter->second.left_paths, offset_iter->second.left_exons.size() - 1);
					}
					else
					{
						offset_iter->second.construct_graph(offset_iter->second.left_exons, offset_iter->second.left_splice_ways, doner_graph, '+');

						DFS(doner_graph, offset_iter->second.left_paths, 0);
					}
				}

				cout << "right splice ways" << endl;

				for (size_t i = 0; i < offset_iter->second.right_splice_ways.size(); ++i)
				{
					cout << offset_iter->second.right_splice_ways[i].chrom_name<<'\t'<< offset_iter->second.right_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < offset_iter->second.right_splice_ways[i].spliceway_vec.size(); ++j)
						cout << offset_iter->second.right_splice_ways[i].spliceway_vec[j].first<<':' << offset_iter->second.right_splice_ways[i].spliceway_vec[j].second << '\t';

					cout << endl;
				}

				if (offset_iter->second.right_splice_ways.size())
				{
					offset_iter->second.spliceways2exons_acceptor();

					vector<int> acceptor_array(offset_iter->second.right_exons.size(), -1);

					vector<vector<int> > acceptor_graph(offset_iter->second.right_exons.size(), acceptor_array);

					if (offset_iter->second.m_strand2 == '+')
					{
						offset_iter->second.construct_graph(offset_iter->second.right_exons, offset_iter->second.right_splice_ways, acceptor_graph, '+');

						DFS(acceptor_graph, offset_iter->second.right_paths, 0);
					}
					else
					{
						offset_iter->second.construct_graph(offset_iter->second.right_exons, offset_iter->second.right_splice_ways, acceptor_graph, '-');

						DFS(acceptor_graph, offset_iter->second.right_paths, offset_iter->second.right_exons.size() - 1);
					}
				}
			}
		}
	}
}



void JunctionHandler::CollectStats()
{
	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	JUNC_HASH_COMB::iterator offset_iter;

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{
			if (offset_iter->second.m_filtered_type == NOT_FILTERED)
			{
				if (offset_iter->second.m_start == offset_iter->second.m_end)
					++m_small_ins;
				else if (offset_iter->second.m_end - offset_iter->second.m_start <= (size_t) m_delete_len )
					++m_small_del;
				else if (offset_iter->second.m_flankcase >= 5)
					++m_canonical;
				else if (offset_iter->second.m_flankcase >= 1)
					++m_semi_canonical;
				else if (offset_iter->second.m_flankcase == 0)
					++m_non_canonical;
				else
					cout << "what?"<<endl;
			}
		}
	}
}

void JunctionHandler::ClearStats()
{
	m_small_ins = 0;

	m_small_del = 0;

	m_canonical = 0;

	m_semi_canonical = 0;

	m_non_canonical = 0;
}

void JunctionHandler::WriteStats(string stat_file, string headline)
{
	ofstream ofs(stat_file.c_str(), ofstream::app);

	ofs << headline << endl;

	ofs << "small_insertion_junction\t"<< m_small_ins << endl;

	ofs << "small_deletion_junction\t"<< m_small_del << endl;

	ofs << "canonical_junction\t"<< m_canonical << endl;

	ofs << "semi_canonical_junction\t"<< m_semi_canonical << endl;

	ofs << "nonical_junction\t"<< m_non_canonical << endl;
}

void JunctionHandler::MarkFiltered(bool paired)
{
	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	JUNC_HASH_COMB::iterator offset_iter;

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{
			if (offset_iter->second.m_is_fusion)
				continue;

			bool isindel = false;

			if (offset_iter->second.m_start == offset_iter->second.m_end || offset_iter->second.m_end - offset_iter->second.m_start <= (size_t) m_delete_len)
				isindel = true;


/*			if (offset_iter->second.m_end == offset_iter->second.m_start)
				offset_iter->second.m_filtered_type = FILTERED_BY_INSERTION;
			else */if ((offset_iter->second.m_max_prefix_len < m_min_junc_anchor || offset_iter->second.m_max_suffix_len < m_min_junc_anchor) && 
				(offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) )
				offset_iter->second.m_filtered_type = FILTERED_BY_SMALL_ANCHOR;
			//else if (offset_iter->second.m_end - offset_iter->second.m_start <= (size_t) m_delete_len && 
			//	offset_iter->second.m_end - offset_iter->second.m_start > 0/* && offset_iter->second.m_flankcase < 5*/)
			//	offset_iter->second.m_filtered_type = FILTERED_BY_SMALL_DELETION;
			else if (offset_iter->second.m_min_anchor_difference > max_anchor_diff && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
				offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MIN_ANCHOR_DIFF;
//			else if (paired && offset_iter->second.m_paired_count && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) )
//			{
///*				if ((double)offset_iter->second.m_paired_unique_count / (double) offset_iter->second.m_paired_count < 0.6)
//					offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MULTIPLE_PAIRED;			
//				else*/ //if (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)
//					offset_iter->second.m_filtered_type = FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR;
//			}
			else if (offset_iter->second.m_paired_count == 0 && paired)
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_NOPAIRED;
			}
			else if (offset_iter->second.m_ave_mismatch > 2 && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
				offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MISMATCH;
			//else if (!isindel && paired && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_left_paired_count < 2 || offset_iter->second.m_right_paired_count < 2))
			//{
			//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_LEFT_RIGHT_PAIR;
			//}
			else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_entropy < 2.5) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ENTROPY;
			}
			else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_multi_count > 0) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_MULTI;
			}
			else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_ave_mismatch > 1) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ERROR;
			}

			offset_iter->second.m_five_prime_known_id.clear();

			offset_iter->second.m_mapped_idx.clear();

			offset_iter->second.m_prefix_count.clear();

			offset_iter->second.m_three_prime_known_id.clear();
		}
	}
}

void JunctionHandler::FindCorrespondingJunction(SamRec* sam_rec)
{
	if (sam_rec->spliceway_vec.size() <= 1)
		return;

	vector<pair<size_t, int> >::iterator oft_mpl_iter;

	for (oft_mpl_iter = (sam_rec)->spliceway_vec.begin(); oft_mpl_iter != (sam_rec)->spliceway_vec.end() - 1; ++oft_mpl_iter)
	{
		if (oft_mpl_iter->second < 0 || (oft_mpl_iter + 1)->second < 0)
		{
			continue;
		}

		size_t comb_offset = ((oft_mpl_iter->first + oft_mpl_iter->second - 1) << THIRTY_TWO) + (oft_mpl_iter + 1)->first;

		CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = m_junc_hash.find(sam_rec->chrom_name);

		if (chrom_junc_hash_iter == m_junc_hash.end())
		{
			//cout << "corresponding junction not found" << endl;

			//cout << sam_rec->tostring((sam_rec)->spliceway_vec.size(), oft_mpl_iter - (sam_rec)->spliceway_vec.begin() + 1)<<endl;

			//cout << sam_rec->chrom_name << '\t' << oft_mpl_iter->first + oft_mpl_iter->second - 1 << '\t' << (oft_mpl_iter + 1)->first << endl;

			continue;
		}

		JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

		JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(comb_offset);

		if (junc_hash_comb_iter == junc_hash_comb.end())
		{
			//cout << "corresponding junction not found" << endl;

			//cout << sam_rec->tostring((sam_rec)->spliceway_vec.size(), oft_mpl_iter - (sam_rec)->spliceway_vec.begin() + 1)<<endl;

			//cout << sam_rec->chrom_name << '\t' << oft_mpl_iter->first + oft_mpl_iter->second - 1 << '\t' << (oft_mpl_iter + 1)->first << endl;

			continue;
		}

		sam_rec->filter_type.insert(junc_hash_comb_iter->second.m_filtered_type);

		sam_rec->corresponding_juncs.push_back(&junc_hash_comb_iter->second);

		if (junc_hash_comb_iter->second.m_flankcase >= 5)
		{
			(sam_rec)->canon_count++;
			(sam_rec)->iscanonical = true;
		}
		else if (junc_hash_comb_iter->second.m_flankcase >= 1)
		{
			(sam_rec)->canon_count++;

			(sam_rec)->noncanon_count++;

			(sam_rec)->issemicanonical = true;
		}
		else
		{
			(sam_rec)->noncanon_count++;

			(sam_rec)->isnoncanoical = true;
		}
	}
}


void JunctionHandler::FindMinPreviousExon(vector<JunctionSeed*>& junction_set, size_t curr_junc)
{
	for(int i=(int)curr_junc - 1;i>=0; i--) 
	{
		if(junction_set[curr_junc]->m_flankcase!= 0 && junction_set[i]->m_flankcase!= 0 && junction_set[curr_junc]->m_strand != junction_set[i]->m_strand)
			continue;
		if(junction_set[curr_junc]->m_start <= junction_set[i]->m_end)
			continue;
		junction_set[curr_junc]->m_left_exon = (unsigned int)(junction_set[curr_junc]->m_start - junction_set[i]->m_end + 1);
		break;
	}
}	

void JunctionHandler::FindMinNextExon(vector<JunctionSeed*>& junction_set, size_t curr_junc)
{
	for(int i=(int)curr_junc+1; i<(int)junction_set.size(); i++)
	{
		if(junction_set[curr_junc]->m_flankcase!= 0 && junction_set[i]->m_flankcase != 0 && junction_set[curr_junc]->m_strand != junction_set[i]->m_strand)
			continue;
		if(junction_set[i]->m_start <= junction_set[curr_junc]->m_end)
			continue;
		junction_set[curr_junc]->m_right_exon = (unsigned int)(junction_set[i]->m_start - junction_set[curr_junc]->m_end + 1 );
		break;
	}	
}

void JunctionHandler::GetMinimumExon()
{		
	SortJuncByEnd();
	map<string, vector<JunctionSeed*> >::iterator chrom_iter;
	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		for(size_t i =0; i< chrom_iter->second.size(); i++)
			FindMinPreviousExon(chrom_iter->second, i);
	}
	SortJunc();
	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		for(size_t i =0; i< chrom_iter->second.size(); i++)
			FindMinNextExon(chrom_iter->second, i);
	}
	//m_junc_sort.clear();
}

Jump_Code::Jump_Code(size_t _len, Jump_Type _type)
{
	len = _len;
	type = _type;
}

Junc_Seq::Junc_Seq()
{
	start = 0;
	total_len = 0;
}

Junc_Seq::~Junc_Seq()
{}

Junc_Seq::Junc_Seq(const Junc_Seq& my_junc_seq)
{
	start = my_junc_seq.start;
	jump_code = my_junc_seq.jump_code;
	sequence = my_junc_seq.sequence;
	strand = my_junc_seq.strand;
	flank_case = my_junc_seq.flank_case;
	junc_id= my_junc_seq.junc_id;
	total_len = my_junc_seq.total_len;
}

Junc_Seq::Junc_Seq(char _strand, 	unsigned short _flank_case, size_t _junc_id, Jump_Code& _new_code)
{
	strand = _strand;
	flank_case = _flank_case;
	junc_id = _junc_id;
	start = 0;
	total_len = 0;
	jump_code.push_back(_new_code);
}

Junc_Seq::Junc_Seq(char _strand, 	unsigned short _flank_case, size_t _junc_id, size_t _start)
{
	strand = _strand;
	flank_case = _flank_case;
	junc_id = _junc_id;
	start = _start;
	total_len = 0;
}

size_t JunctionHandler::CheckJuncSeqNum(map<size_t, vector<Junc_Seq> >::iterator& it_head, map<size_t, vector<Junc_Seq> >::iterator& it_tail)
{
	size_t total_num = 0;
	for(size_t i = 0; i < it_head->second.size(); i++)
	{
		for(size_t j = 0; j < it_tail->second.size(); j++)
		{
			if(it_head->second[i].flank_case != 0 && it_tail->second[j].flank_case != 0  && it_head->second[i].strand != it_tail->second[j].strand)
				continue;
			total_num++;
		}
	}
	return total_num;
}

void JunctionHandler::OutPutJuncSeq(string chrom, string& chromseq, size_t max_seq_thresh, ofstream& out_fs)
{
	map<size_t, vector<Junc_Seq> >::iterator it_head, it_tail;
	for(it_head = m_head_seq.begin(); it_head != m_head_seq.end(); it_head++)
	{
		it_tail = m_tail_seq.find(it_head->first);
		size_t head_bound = it_head->second.size();
		size_t tail_bound = it_tail->second.size();
		if(CheckJuncSeqNum(it_head, it_tail) > max_seq_thresh)  // only output normal path for maxed junction
		{
			head_bound = 1;
			tail_bound = 1;
		}
		for(size_t i = 0; i < head_bound; i++)
		{
			for(size_t j = 0; j < tail_bound; j++)
			{
				if(it_head->second[i].flank_case != 0 && it_tail->second[j].flank_case != 0  && it_head->second[i].strand != it_tail->second[j].strand)
					continue;
				out_fs<<">"<<chrom<<"_"<<it_head->second[i].start<<":";
				for(int jc_index = (int)it_head->second[i].jump_code.size() - 1; jc_index >=0; jc_index-- )
				{
					out_fs << it_head->second[i].jump_code[jc_index].len << (it_head->second[i].jump_code[jc_index].type == M ? 'M':'N');
				}
				for(int jc_index = 0; jc_index < (int)it_tail->second[j].jump_code.size(); jc_index++ )
				{
					out_fs << it_tail->second[j].jump_code[jc_index].len << (it_tail->second[j].jump_code[jc_index].type == M ? 'M':'N');
				}
				out_fs << endl;
				string whole_sequence;
				for(int seq_index = (int)it_head->second[i].sequence.size() - 1; seq_index >=0; seq_index-- )
				{
					whole_sequence.append(chromseq.substr(it_head->second[i].sequence[seq_index].first - 1, it_head->second[i].sequence[seq_index].second - it_head->second[i].sequence[seq_index].first + 1));
				}
				for(int seq_index = 0; seq_index < (int)it_tail->second[j].sequence.size(); seq_index++ )
				{
					whole_sequence.append(chromseq.substr(it_tail->second[j].sequence[seq_index].first - 1, it_tail->second[j].sequence[seq_index].second - it_tail->second[j].sequence[seq_index].first + 1));
				}
				size_t start = 0;
				size_t end= start + NUM_BP_PER_LINE;
				while(true)
				{
					if(end<=whole_sequence.length())
					{
						out_fs<<whole_sequence.substr(start,NUM_BP_PER_LINE)<<endl;
						if(end==whole_sequence.length())
							break;
					}
					else
					{
						out_fs<<whole_sequence.substr(start,whole_sequence.length() - start)<<endl;
						break;
					}
					start=end;
					end=start + NUM_BP_PER_LINE;
				}
			}
		}
	}
}

void JunctionHandler::AddJunSeq(map< size_t, vector<Junc_Seq> >& junc_seq_list, Junc_Seq& my_junc_seq, size_t max_seq_thresh)
{
	map<size_t, vector<Junc_Seq> >::iterator it = junc_seq_list.find(my_junc_seq.junc_id);
	if(it == junc_seq_list.end())
	{
		vector<Junc_Seq> new_junc_seq_vec;
		junc_seq_list.insert(make_pair(my_junc_seq.junc_id, new_junc_seq_vec));
	}
	if(junc_seq_list[my_junc_seq.junc_id].size() < max_seq_thresh)
		junc_seq_list[my_junc_seq.junc_id].push_back(my_junc_seq);
}

void JunctionHandler::GetJuncHead(Junc_Seq& my_junc_seq, vector<JunctionSeed*>& junction_set, size_t curr_junc, size_t min_anchor, size_t max_anchor, size_t max_seq_thresh)
{
	//cout << my_junc_seq.total_len << endl;
	size_t len_left = max_anchor - my_junc_seq.total_len;
	Junc_Seq normal_junc_seq(my_junc_seq);		//normal path
	if(junction_set[curr_junc]->m_start >= len_left)
		normal_junc_seq.start = junction_set[curr_junc]->m_start - len_left + 1;
	else
	{
		normal_junc_seq.start = 1;
		len_left = junction_set[curr_junc]->m_start - normal_junc_seq.start + 1;
	}
	normal_junc_seq.sequence.push_back(make_pair(normal_junc_seq.start, junction_set[curr_junc]->m_start));
	Jump_Code new_jumpcode(len_left, M);
	normal_junc_seq.jump_code.push_back(new_jumpcode);
	normal_junc_seq.total_len += len_left;
	AddJunSeq(m_head_seq, normal_junc_seq, max_seq_thresh);
	for(int i=(int)curr_junc - 1;i>=0; i--)    // more splice
	{
		if(my_junc_seq.flank_case == 0)
		{
			my_junc_seq.flank_case = junction_set[i]->m_flankcase;
			my_junc_seq.strand = junction_set[i]->m_strand;
		}
		if(junction_set[i]->m_flankcase != 0 && junction_set[i]->m_strand != my_junc_seq.strand)
			continue;
		if(junction_set[curr_junc]->m_start <= junction_set[i]->m_end + min_anchor)
			continue;
		else if(junction_set[curr_junc]->m_start >= junction_set[i]->m_end + len_left - 1)
			break;
		else if(m_head_seq[my_junc_seq.junc_id].size() < max_seq_thresh)
		{
			Junc_Seq new_junc_seq(my_junc_seq);
			new_junc_seq.sequence.push_back(make_pair(junction_set[i]->m_end, junction_set[curr_junc]->m_start));
			new_junc_seq.total_len += junction_set[curr_junc]->m_start - junction_set[i]->m_end + 1;
			Jump_Code new_jumpcode1(junction_set[curr_junc]->m_start - junction_set[i]->m_end + 1, M);
			Jump_Code new_jumpcode2(junction_set[i]->m_end - junction_set[i]->m_start - 1 , N);
			new_junc_seq.jump_code.push_back(new_jumpcode1);
			new_junc_seq.jump_code.push_back(new_jumpcode2);
			GetJuncHead(new_junc_seq, junction_set, i, min_anchor, max_anchor, max_seq_thresh);
		}
	}
}

void JunctionHandler::GetJuncTail(Junc_Seq& my_junc_seq, vector<JunctionSeed*>& junction_set, size_t curr_junc, size_t min_anchor, size_t max_anchor, size_t chrom_len, size_t max_seq_thresh)
{
	//cout << my_junc_seq.total_len << endl;
	size_t len_left = max_anchor - my_junc_seq.total_len;
	Junc_Seq normal_junc_seq(my_junc_seq);		//normal path
	size_t seq_end_pos=min(junction_set[curr_junc]->m_end + len_left - 1, chrom_len);
	len_left = seq_end_pos - junction_set[curr_junc]->m_end + 1;
	normal_junc_seq.sequence.push_back(make_pair(junction_set[curr_junc]->m_end, seq_end_pos));
	Jump_Code new_jumpcode(len_left, M);
	normal_junc_seq.jump_code.push_back(new_jumpcode);
	normal_junc_seq.total_len += len_left;
	AddJunSeq(m_tail_seq, normal_junc_seq, max_seq_thresh);
	for(int i = (int)curr_junc+1;i < (int)junction_set.size(); i++)  //more splice
	{
		if(my_junc_seq.flank_case == 0)
		{
			my_junc_seq.flank_case = junction_set[i]->m_flankcase;
			my_junc_seq.strand = junction_set[i]->m_strand;
		}
		if(junction_set[i]->m_flankcase != 0 && junction_set[i]->m_strand != my_junc_seq.strand)
			continue;
		if( junction_set[i]->m_start <= junction_set[curr_junc]->m_end + min_anchor)       // too close
			continue;
		else if(junction_set[i]->m_start >= junction_set[curr_junc]->m_end + len_left - 1)   // too far
			break;
		else if(m_tail_seq[my_junc_seq.junc_id].size() < max_seq_thresh)                   // fit
		{
			Junc_Seq new_junc_seq(my_junc_seq);
			new_junc_seq.sequence.push_back(make_pair(junction_set[curr_junc]->m_end, junction_set[i]->m_start));
			new_junc_seq.total_len += junction_set[i]->m_start - junction_set[curr_junc]->m_end + 1;
			Jump_Code new_jumpcode1(junction_set[i]->m_start - junction_set[curr_junc]->m_end + 1, M);
			Jump_Code new_jumpcode2(junction_set[i]->m_end - junction_set[i]->m_start - 1, N);
			new_junc_seq.jump_code.push_back(new_jumpcode1); 
			new_junc_seq.jump_code.push_back(new_jumpcode2);
			GetJuncTail(new_junc_seq, junction_set, i, min_anchor, max_anchor, chrom_len, max_seq_thresh);
		}
	}	
}

void JunctionHandler::GenerateJunctionSequence(size_t min_anchor, size_t max_anchor, size_t max_seq_thresh, string output_file)
{
	ofstream out_fs(output_file.c_str());
	map<string, vector<JunctionSeed*> >::iterator chrom_iter;
	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		string chromfile = m_chrom_dir + chrom_iter->first;
		chromfile.append(".fa");
		string chromseq;
		readchrom(chromfile.c_str(), chromseq);
		size_t chrom_len = chromseq.length();
		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort_by_end);
		for(size_t i = 0; i< chrom_iter->second.size(); i++)
		{
			Jump_Code new_head_code(chrom_iter->second[i]->m_end - chrom_iter->second[i]->m_start -1, N);
			Junc_Seq new_junc_seq(chrom_iter->second[i]->m_strand, chrom_iter->second[i]->m_flankcase, chrom_iter->second[i]->m_junc_id, new_head_code);
			GetJuncHead(new_junc_seq, chrom_iter->second, i, min_anchor, max_anchor, max_seq_thresh);
		}
		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
		for(size_t i = 0; i< chrom_iter->second.size(); i++)
		{
			Junc_Seq new_junc_seq(chrom_iter->second[i]->m_strand, chrom_iter->second[i]->m_flankcase, chrom_iter->second[i]->m_junc_id, chrom_iter->second[i]->m_end);
			GetJuncTail(new_junc_seq, chrom_iter->second, i, min_anchor, max_anchor, chrom_len, max_seq_thresh);
		}
		OutPutJuncSeq(chrom_iter->first, chromseq, max_seq_thresh,out_fs);
		m_head_seq.clear();
		m_tail_seq.clear();
		chromseq.clear();
	}
	out_fs.close();
}


map<string, vector<JunctionSeed*> >& JunctionHandler::GetJuncSortList()
{
	return m_junc_sort;
}

void JunctionHandler::ConfirmJunction(JunctionHandler& database_handler, size_t range, bool clear_internal_result)
{
	map<string, vector<JunctionSeed*> >& other_junc_sort = database_handler.GetJuncSortList();
	map<string, vector<JunctionSeed*> >::iterator this_chrom_iter, other_chrom_iter;

	string compare_junc = m_junc_files.front(); compare_junc.append(".comp"); ofstream ofs(compare_junc.c_str());
	
	size_t count = 0;
	SortJuncByEnd();
	database_handler.SortJuncByEnd();
	for (this_chrom_iter = m_junc_sort.begin(); this_chrom_iter != m_junc_sort.end(); ++this_chrom_iter)  // search 3 prime 
	{
		// per chrom initialization
		other_chrom_iter = other_junc_sort.find(this_chrom_iter->first);     
		if(other_chrom_iter == other_junc_sort.end())
			continue;
		vector <JunctionSeed* >& this_junc_vec = this_chrom_iter->second;
		vector <JunctionSeed* >& other_junc_vec = other_chrom_iter->second;
		size_t this_index = 0;
		size_t other_index = 0;
		size_t backtrack_index = 0;
		bool has_backtrack_index = false;
		while(this_index < this_junc_vec.size())
		{
			if(other_index >= other_junc_vec.size() || this_junc_vec[this_index]->m_end + range < other_junc_vec[other_index]->m_end)
			{
				this_index ++;
				if(has_backtrack_index)
					other_index = backtrack_index;
				else if(other_index >= other_junc_vec.size())
					break;
				has_backtrack_index = false;
			}
			else if(this_junc_vec[this_index]->m_end > other_junc_vec[other_index]->m_end + range)
			{
				other_index ++;
			}
			else
			{
				if(!has_backtrack_index)
				{
					backtrack_index = other_index;
					has_backtrack_index = true;
				}
				this_junc_vec[this_index]->m_three_prime_known_id.insert(make_pair(other_junc_vec[other_index]->m_junc_id, 3));
				other_index ++;
			}
		}
	}
	SortJunc();
	database_handler.SortJunc();
	for (this_chrom_iter = m_junc_sort.begin(); this_chrom_iter != m_junc_sort.end(); ++this_chrom_iter)  // search 5 prime 
	{
		// per chrom initialization
		other_chrom_iter = other_junc_sort.find(this_chrom_iter->first);     
		if(other_chrom_iter == other_junc_sort.end())
			continue;
		vector <JunctionSeed* >& this_junc_vec = this_chrom_iter->second;
		vector <JunctionSeed* >& other_junc_vec = other_chrom_iter->second;
		size_t this_index = 0;
		size_t other_index = 0;
		size_t backtrack_index = 0;
		bool has_backtrack_index = false;
		while(this_index < this_junc_vec.size() )
		{
			if(other_index >= other_junc_vec.size() || this_junc_vec[this_index]->m_start + range < other_junc_vec[other_index]->m_start)
			{
				this_index ++;
				if(has_backtrack_index)
					other_index = backtrack_index;
				else if(other_index >= other_junc_vec.size())
					break;
				has_backtrack_index = false;
			}
			else if(this_junc_vec[this_index]->m_start > other_junc_vec[other_index]->m_start + range)
			{
				other_index ++;
			}
			else
			{
				if(!has_backtrack_index)
				{
					backtrack_index = other_index;
					has_backtrack_index = true;
				}
				this_junc_vec[this_index]->m_five_prime_known_id.insert(make_pair(other_junc_vec[other_index]->m_junc_id, 5));
				other_index ++;
			}
		}
	}
	for (this_chrom_iter = m_junc_sort.begin(); this_chrom_iter != m_junc_sort.end(); ++this_chrom_iter)  // find confirm pairing
	{
		vector <JunctionSeed* >& this_junc_vec = this_chrom_iter->second;
		for(size_t i = 0 ; i < this_junc_vec.size() ; i ++ )
		{
			if(this_junc_vec[i]->m_five_prime_known_id.size() > 0)
				this_junc_vec[i]->m_five_prime_known = true;
			if(this_junc_vec[i]->m_three_prime_known_id.size() > 0)
				this_junc_vec[i]->m_three_prime_known = true;
			for(map<size_t, size_t>::iterator it = this_junc_vec[i]->m_three_prime_known_id.begin(); it != this_junc_vec[i]->m_three_prime_known_id.end(); ++it)
			{
				if(this_junc_vec[i]->m_five_prime_known_id.find(it->first) != this_junc_vec[i]->m_five_prime_known_id.end())
				{
					this_junc_vec[i]->m_pair_known = true;
					break;
				}
			}
			if(clear_internal_result)
			{
				this_junc_vec[i]->m_five_prime_known_id.clear();
				this_junc_vec[i]->m_three_prime_known_id.clear();
			}

			ofs << this_junc_vec[i]->to_normal_junction(++count) <<'\t' << this_junc_vec[i]->m_left_exon << '\t' << this_junc_vec[i]->m_right_exon <<  "\tcomp:"
				<<this_junc_vec[i]->m_five_prime_known << this_junc_vec[i]->m_three_prime_known << this_junc_vec[i]->m_pair_known << endl;

			//cout << this_junc_vec[i]->m_five_prime_known << this_junc_vec[i]->m_three_prime_known << this_junc_vec[i]->m_pair_known << endl;
		}
	}
}

JuncSeedSp::JuncSeedSp(size_t max_prefix_len, size_t max_suffix_len, size_t start, size_t end, const string& ins_str, SamRec* sam_rec_ptr, bool is_fusion) :
m_max_prefix_len(max_prefix_len), m_max_suffix_len(max_suffix_len), m_start(start), m_end(end), m_ins_str(ins_str), m_sam_rec_ptr(sam_rec_ptr), m_is_fusion(is_fusion)
{
}

void
JuncSeedSp::set(size_t max_prefix_len, size_t max_suffix_len, size_t start, size_t end, const string& ins_str, SamRec* sam_rec_ptr, bool is_fusion)
{
	m_max_prefix_len = max_prefix_len;
	
	m_max_suffix_len = max_suffix_len;
	
	m_start = start;
	
	m_end = end;
	
	m_ins_str = ins_str;
	
	m_sam_rec_ptr = sam_rec_ptr;
	
	m_is_fusion = is_fusion;
}
