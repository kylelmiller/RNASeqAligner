#include "JunctionHandler.h"
//#include "ReadNextTagAlignHandler.h"

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
				return lhs.m_acceptor_end < rhs.m_acceptor_end;
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

JunctionHandler::JunctionHandler(): m_canonical(0), m_semi_canonical(0), m_non_canonical(0), m_fusion(0), m_fusion_canon(0), m_fusion_semicanon(0), m_fusion_noncanon(0), m_delete_len(0), m_small_del(0), m_small_ins(0), m_min_mismatch(0), m_min_anchor(0), m_total_junc(0), m_total_encompass(0)
{
}

bool JunctionHandler::Init(const vector<string>& alignment_files, size_t max_read_width, int insert_len, int delete_len, 
	const string& chrom_dir, size_t min_anchor, size_t min_mismatch, size_t min_junc_anchor, size_t min_fusion_coverage, int do_filter)
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

	m_min_fusion_coverage = min_fusion_coverage;

	m_do_filter = do_filter;

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

	m_fusion = 0;

	m_fusion_canon = 0;
	
	m_fusion_semicanon = 0;
	
	m_fusion_noncanon = 0;
	
	m_small_ins = 0;

	m_min_anchor = 0;

	m_min_mismatch = 0;

	m_min_junc_anchor = 0;

	return true;
}

bool JunctionHandler::ClearJunction()
{
	cout << "clear junction"<<endl;

	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		//cout << chrom_iter->first<<endl;

		JUNC_HASH_COMB::iterator offset_iter;

		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); )
		{
			//cout << offset_iter->second.to_normal_junction(0)<<endl;
			if (offset_iter->second.m_filtered_type != NOT_FILTERED)
			{
				chrom_iter->second.erase(offset_iter++);
			}
			else
			{
				offset_iter->second.clear();
				++offset_iter;
			}
		}
	}

	cout << "cleared junction"<<endl;

	return true;
}

size_t JunctionHandler::ReadJunction(vector<string>& junc_files)
{
	vector<string>::iterator junc_file_iter;

	size_t count = 0;

	for (junc_file_iter = junc_files.begin(); junc_file_iter != junc_files.end(); ++junc_file_iter)
	{
		m_junc_files.push_back(*junc_file_iter);

		cout << *junc_file_iter << endl;

		ifstream ifs(junc_file_iter->c_str());

		if (ifs.is_open())
		{
			//getline(ifs, m_head_line);

			while (!ifs.eof() )
			{
				string line;

				getline(ifs,line);

				if (line == "" || line.find("chr") == string::npos)
				{
					if (line.find("chr") == string::npos)
						cout << line << endl;

					continue;
				}

				char chromname[100], juncname[100], strand, rgb[100], flankchr[10];

				unsigned int hits = 0;

				unsigned short kinds = 0, flankcase = 0;

				size_t juncst, juncend, prefixend, suffixst;

				double rank = 0, lpq = 0, il = 0;

				unsigned short min_mis = 0, max_mis = 0;

				double ave_mis = 0;

				size_t max_prefix_len = 0, max_suffix_len = 0;

				size_t start_block_offset = 0, end_block_offset = 0;

				unsigned int unique_count = 0, multi_count = 0;

				unsigned int paired_count = 0, single_count = 0;

				unsigned int left_paired_count = 0, right_paired_count = 0;

				unsigned int paired_mutiple_count = 0, paired_unique_count = 0;

				unsigned short min_anchor_diff = 0;

				//chr1	14829	14930	JUNC_1	1	-	14829	14930	255,0,0	2	30,20,	0,122,	0.000000	6	CTAC	0.999495	0.999986	1	1	1.000000	1	0	1	0	1	0	1	0	10

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
					//cout << "find same junction"<<endl;
					//cout << junc_hash_comb_iter->second.to_normal_junction(0)<<endl;

					//cout << hits << endl;
					//cout << rank << endl;

					junc_hash_comb_iter->second.m_hits += hits;

					if (junc_hash_comb_iter->second.m_entropy > rank) 
						junc_hash_comb_iter->second.m_entropy = rank;

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

void JunctionHandler::UpdateJuncTable(const string& chrom_name, size_t combined_offset, size_t prefixlen, size_t suffixlen, size_t fusion_prefixlen, size_t fusion_suffixlen, size_t start, size_t end, 
	const string& ins_str, SamRec& samrec, size_t sam_count, bool is_fusion)
{
	//if (((start == 68429995 && end == 29631805) || (start == 29631805 && end == 68429995)) && chrom_name.find("~") != string::npos)
	//	cout << samrec.tostring(0, 0) << endl;

	//string chrom_name;

	//if (is_fusion)
	//	chrom_name = samrec.chrom_name + "~" + samrec.chrom_name2;
	//else
	//	chrom_name = samrec.chrom_name;

	//cout << 11 << endl;

	//cout << chrom_name << endl;

	CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = m_junc_hash.find(chrom_name);

	if (chrom_junc_hash_iter == m_junc_hash.end())
	{
		JUNC_HASH_COMB junc_hash_comb;

		chrom_junc_hash_iter = (m_junc_hash.insert(CHROM_JUNC_HASH_COMB::value_type(chrom_name, junc_hash_comb))).first;							
	}

	//cout << 22 << endl;

	JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

	JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

	if (junc_hash_comb_iter != junc_hash_comb.end())
	{
		//cout << 33 << endl;
		if (samrec.mate_match != "*" && samrec.mate_diff == 0)
			cout <<"mate dist 0"<< endl << samrec.tostring(1, 1) << endl;
		junc_hash_comb_iter->second.inc_hits(prefixlen, suffixlen, fusion_prefixlen, fusion_suffixlen, samrec.tagidx, samrec.mis_match, samrec.strand_t, sam_count, samrec.mate_match, 
			samrec.mate_diff, samrec.left_splice_ways, samrec.right_splice_ways, &samrec, ins_str);
	}
	else
	{
		//cout << 44 << endl;
		++m_total_junc;
		chrom_junc_hash_iter->second.insert(JUNC_HASH_COMB::value_type(combined_offset, JunctionSeed(/*prim, *//*flankseq, */
			prefixlen, suffixlen, fusion_prefixlen, fusion_suffixlen, m_max_read_width, samrec.tagidx, samrec.mis_match, samrec.strand_t, samrec.strand_t2, start, end, samrec.chrom_name, samrec.chrom_name2, 
			sam_count, samrec.mate_match, samrec.mate_diff, samrec.left_splice_ways, samrec.right_splice_ways, is_fusion, &samrec, ins_str)));
	}

	//size_t junc_seed_size = sizeof(unsigned short) * junc_hash_comb.find(combined_offset)->second.m_prefix_count.size();
}


void JunctionHandler::SamRec2Junc(SamRec& samrec, size_t sam_count, hash_map<string, hash_map<size_t, JuncSeedSp> >& cur_read_junc)
{
	vector<pair<size_t, int> >& spliceway_vec = samrec.spliceway_vec;

	vector<pair<size_t, int> >& spliceway_vec2 = samrec.spliceway_vec2;

	const string& readstr = samrec.mapped_seq;

	size_t mappedlen = 0;

	size_t flankidx = 0;

	if (spliceway_vec.size() > 1)
	{
		vector<pair<size_t, int> >::iterator vp_iter;

		//cout << "1"<<endl;
		//cout << "spliceway_vec"<<endl;

		//for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
		//{
		//	cout << vp_iter->first << ':' << vp_iter->second << endl;
		//}

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

			if (ins_str.empty() && m_chrom_dir.empty())
			{
				if (flankidx >= samrec.flankstrings.size())
					cerr << samrec.tostring(0,0)<<endl;

				ins_str = samrec.flankstrings[flankidx++];
			}

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
					//cout << "aaaaa"<<endl;
					//cout << prefixlen << endl;
					//cout << suffixlen << endl;
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
		//cout << "2"<<endl;

		//cout << "spliceway_vec2.size():"<<spliceway_vec2.size()<<endl;

		//for (size_t ii = 0; ii < spliceway_vec.size();++ii)
		//{
		//	cout << spliceway_vec[ii].first <<'\t' << spliceway_vec[ii].second << endl;
		//}

		vector<pair<size_t, int> >::iterator vp_iter;

		size_t flankidx2 = flankidx + 1; //samrec.flankstrings.size() - 1;
		//cout << "spliceway_vec2"<<endl;

		//for (vp_iter =  spliceway_vec2.begin(); vp_iter != spliceway_vec2.end(); ++vp_iter)
		//{
		//	cout << vp_iter->first << ':' << vp_iter->second << endl;
		//}

		//cout << "21"<<endl;

		for (vp_iter =  spliceway_vec2.begin(); vp_iter != spliceway_vec2.end(); ++vp_iter)
		{
			//cout << "22"<<endl;
			size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

			string ins_str;

			//cout << vp_iter->first << '\t' << vp_iter->second << endl;

			if (vp_iter->second < 0)
			{
				//cout << "221"<<endl;
				if (vp_iter == spliceway_vec2.begin() || (vp_iter - 1)->second  < 0)
					prefixlen = 1;
				else
					prefixlen = (vp_iter - 1)->second + 1;

				if (vp_iter == spliceway_vec2.end() - 1 || (vp_iter + 1)->second  < 0)
					suffixlen = 1;
				else
					suffixlen = (vp_iter + 1)->second + 1;

				//cout << "222"<<endl;

				//cout <<"readstr.length():"<< readstr.length() << endl;
				//
				//cout <<"mappedlen:"<< mappedlen << endl;

				//cout <<"-(vp_iter->second)"<< -(vp_iter->second) << endl;
				
				ins_str = readstr.substr(mappedlen, -(vp_iter->second));

				prefixend = vp_iter->first - 1;

				//cout << "223"<<endl;
				suffixst = prefixend;

				//cout << "224"<<endl;
				combined_offset = (prefixend << THIRTY_TWO) + suffixst;
			}
			else if (vp_iter == spliceway_vec2.end() - 1 || (vp_iter + 1)->second < 0)
			{
				//cout << "225"<<endl;
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
			//cout << "23"<<endl;
			if (ins_str.empty() && m_chrom_dir.empty())
			{
				if (flankidx2 >= samrec.flankstrings.size())
					cerr << samrec.tostring(0,0)<<endl;

				//wrong ? 
				ins_str = samrec.flankstrings[flankidx2++];
			}

			//cout << "24"<<endl;
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
					//cout << "bbbbb"<<endl;
					//cout << prefixlen << endl;
					//cout << suffixlen << endl;
					junc_hash_comb.insert(hash_map<size_t, JuncSeedSp>::value_type(combined_offset, JuncSeedSp(prefixlen, suffixlen, prefixend, suffixst, ins_str, &samrec, false)));
				}
			}
				//UpdateJuncTable(combined_offset, prefixlen, suffixlen, prefixend, suffixst, ins_str, samrec, sam_count);
		}
	}

	if (samrec.is_fusion)
	{
		//cout << "3"<<endl;
		//if paired reads and is paired
		
		size_t prefixst, prefixend, suffixst, prefixlen, suffixlen, combined_offset;

		string ins_str;

		if (m_chrom_dir.empty())
		{
			if (flankidx >= samrec.flankstrings.size())
			{
				cerr << flankidx << endl;
				cerr << samrec.flankstrings.size() << endl;
				cerr << "fusion sequence empty:"<<samrec.tostring(0,0)<<endl;
			}

			ins_str = samrec.flankstrings[flankidx];
		}

		//prefixst = samrec.fusion_prefix_st;

		//prefixend = samrec.fusion_prefix_end;

		//suffixst = samrec.fusion_suffix_st;

		//if (prefixend == 98658330 || suffixst == 98658330)
		//	cout << samrec.tostring(0, 0) << endl;

		if (prefixlen >= m_min_anchor && suffixlen >= m_min_anchor && samrec.mis_match < m_min_mismatch)
		{
			string chrom_name;

			if (samrec.need_swap == false)
			{
				chrom_name = samrec.chrom_name + "~" + samrec.chrom_name2;

				prefixend = samrec.fusion_prefix_end;

				suffixst = samrec.fusion_suffix_st;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;

				prefixlen = samrec.mappedlen1;//spliceway_vec.back().second; //samrec.fusion_prefix_len;//spliceway_vec.back().second;//samrec.mappedlen1;

				//if (samrec.strand1 == '-')
				//	prefixlen = spliceway_vec.front().second;

				suffixlen = samrec.mappedlen2;//spliceway_vec2.front().second; //samrec.fusion_suffix_len;//spliceway_vec2.front().second;//samrec.mappedlen2;

				//if (samrec.strand2 == '-')
				//	suffixlen = spliceway_vec2.back().second;
			}
			else
			{
				ins_str = revcomp(ins_str);

				chrom_name = samrec.chrom_name2 + "~" + samrec.chrom_name;

				suffixst = samrec.fusion_prefix_end;

				prefixend = samrec.fusion_suffix_st;

				combined_offset = (prefixend << THIRTY_TWO) + suffixst;//(suffixst << THIRTY_TWO) + prefixend;

				suffixlen = samrec.mappedlen1;//spliceway_vec.back().second; //samrec.fusion_prefix_len;//spliceway_vec.back().second;//samrec.mappedlen1;

				//if (samrec.strand1 == '-')
				//	suffixlen = spliceway_vec.front().second;

				prefixlen = samrec.mappedlen2;//spliceway_vec2.front().second; //samrec.fusion_suffix_len;//spliceway_vec2.front().second;//samrec.mappedlen2;

				//if (samrec.strand2 == '-')
				//	prefixlen = spliceway_vec2.back().second;
			}

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

	//cout << 1<<endl;

	for (sam_rec_iter = samrecs.begin(); sam_rec_iter != samrecs.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->min_anchor < m_min_anchor)
			++small_anchor_count;
	}

	//cout << 2<<endl;

	size_t real_sam_count = sam_count - small_anchor_count;

	//vector<SamRec>::iterator sam_rec_iter;

	hash_map<string, hash_map<size_t, JuncSeedSp> > cur_read_junc;

	for (sam_rec_iter = samrecs.begin(); sam_rec_iter != samrecs.end(); ++sam_rec_iter)
	{
		if ((*sam_rec_iter)->isunmapped)
			continue;

		//cout << (*sam_rec_iter)->tostring(0, 0) << endl;

		if (((*sam_rec_iter)->isspliced || (*sam_rec_iter)->issmallins || (*sam_rec_iter)->is_fusion) && real_sam_count)
			SamRec2Junc(**sam_rec_iter, real_sam_count, cur_read_junc);

		if ((*sam_rec_iter)->encompassed_juncs.size())
		{
			*m_fusion_encompassing_reads_acceptor_ofs<< (*sam_rec_iter)->encompassed_juncs.size()<<endl;

			vector<pair<JunctionSeed*, bool> >::iterator juncseed_iter = (*sam_rec_iter)->encompassed_juncs.begin(), juncseed_end = (*sam_rec_iter)->encompassed_juncs.end();
			
			*m_fusion_encompassing_reads_doner_ofs<< (**sam_rec_iter).tostring(0, 0)<<endl;

			m_fusion_encompassing_reads.push_back((**sam_rec_iter));

			if (sam_count > 1)
				m_fusion_encompassing_reads.back().isunique = false;
			else
				m_fusion_encompassing_reads.back().isunique = true;

			m_fusion_encompassing_reads.back().qual_str.clear();

			m_fusion_encompassing_reads.back().mapped_seq.clear();
			//*m_fusion_encompassing_reads_doner_ofs<< ((juncseed_iter->first))->to_normal_junction(0)<<endl;

			for (; juncseed_iter != juncseed_end; ++juncseed_iter)
			{
				++m_total_encompass;

				if (juncseed_iter->second)
				{
					SamRec& samrec = (**sam_rec_iter);

					//if ( samrec.tag_base_name == "UNC9-SN296_125:8:1:11416:65749/")
					//	cout << "UNC9-SN296_125:8:1:11416:65749/ is being encompmassed doner"<<endl;

					;
					((juncseed_iter->first))->m_fusion_encompassing_reads_doner.push_back(m_fusion_encompassing_reads.size() - 1);
					;

					//(*sam_rec_iter)->tostring(0, 0));

					//if ((((juncseed_iter->first))->m_start == 96968937 || ((juncseed_iter->first))->m_start == 96904172) && (((juncseed_iter->first))->m_end == 96968937 || ((juncseed_iter->first))->m_end == 96904172))
					//{
					//	cout << "96968937_96904172:samrec2junc:" << samrec.tostring(0, 0)<<endl;
					//}

					//((juncseed_iter->first))->left_splice_ways.push_back(SpliceWayTrue(&(samrec.chrom_name), samrec.start, &(samrec.spliceway_vec)));

					++m_fusion_encompassing_reads_doner_count;

					//*m_fusion_encompassing_reads_doner_ofs<< samrec.tostring(0, 0)<<endl;

					//*m_fusion_encompassing_reads_doner_ofs<< ((juncseed_iter->first))->to_normal_junction(0)<<endl;
				}
				else
				{
					SamRec& samrec = (**sam_rec_iter);

					//if ( samrec.tag_base_name == "UNC9-SN296_125:8:1:11416:65749/")
					//	cout << "UNC9-SN296_125:8:1:11416:65749/ is being encompmassed acceptor"<<endl;

					;
					((juncseed_iter->first))->m_fusion_encompassing_reads_acceptor.push_back(m_fusion_encompassing_reads.size() - 1);
					;

					//(*sam_rec_iter)->tostring(0, 0));

					//if ((((juncseed_iter->first))->m_start == 96968937 || ((juncseed_iter->first))->m_start == 96904172) && (((juncseed_iter->first))->m_end == 96968937 || ((juncseed_iter->first))->m_end == 96904172))
					//{
					//	cout << "96968937_96904172:samrec2junc:" << samrec.tostring(0, 0)<<endl;
					//}

					//((juncseed_iter->first))->right_splice_ways.push_back(SpliceWayTrue(&(samrec.chrom_name), samrec.start, &(samrec.spliceway_vec)));

					++m_fusion_encompassing_reads_acceptor_count;

					//*m_fusion_encompassing_reads_acceptor_ofs<< samrec.tostring(0, 0)<<endl;

					//*m_fusion_encompassing_reads_acceptor_ofs<< ((juncseed_iter->first))->to_normal_junction(0)<<endl;
				}

				((juncseed_iter->first))->m_encompass_reads_count++;
			}

			//cout<< "doner_side_spanning_pairs_count:"<< doner_side_spanning_pairs_count<<endl;
			//cout<< "accetpr_side_spanning_pairs_count:"<<accetpr_side_spanning_pairs_count<<endl;
			//cout<< "accetpr_side_spanning_pairs_count:"<<single_spanning_count<<endl;
			//cout<< "spliceway_true_count:"<<spliceway_true_count<<endl;
			//cout<< "m_fusion_encompassing_reads_doner_count:"<<m_fusion_encompassing_reads_doner_count<<endl;
			//cout<< "m_fusion_encompassing_reads_acceptor_count:"<<m_fusion_encompassing_reads_acceptor_count<<endl;

		}
	}

	//cout << 3<<endl;

	hash_map<string, hash_map<size_t, JuncSeedSp> >::iterator chrom_iter;

	for (chrom_iter = cur_read_junc.begin(); chrom_iter != cur_read_junc.end(); ++chrom_iter)
	{
		hash_map<size_t, JuncSeedSp>::iterator offset_iter;

		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{
			UpdateJuncTable(chrom_iter->first, offset_iter->first, offset_iter->second.m_max_prefix_len, offset_iter->second.m_max_suffix_len, 
				offset_iter->second.m_fusion_prefix_len,  offset_iter->second.m_fusion_suffix_len, offset_iter->second.m_start,
				offset_iter->second.m_end, offset_iter->second.m_ins_str, *(offset_iter->second.m_sam_rec_ptr), sam_count, offset_iter->second.m_is_fusion);
				//prefixlen, suffixlen, prefixend, suffixst, ins_str, samrec, sam_count);
		}
	}

	//cout << 4<<endl;
}

/*
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
}*/

map<string, string> JunctionHandler::m_loaded_chromosomes;

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
			string* chromseq;

			if (!m_chrom_dir.empty())
			{
				string chromfile = m_chrom_dir + chm_iter->first;

				chromfile.append(".fa");

				if (m_loaded_chromosomes.find(chromfile) == m_loaded_chromosomes.end())
				{
					readchrom(chromfile.c_str(), m_loaded_chromosomes[chromfile]);
				}

				chromseq = &(m_loaded_chromosomes[chromfile]);

				if (chromseq->empty())
				{
					cout <<"empty chrom: "<<chromfile<<endl;
					exit(1);
				}
			}

			size_t chrom_size = m_chrom_sizes[chm_iter->first];//10000000;//chromseq.size() - 1;

			JUNC_HASH_COMB::iterator iter_conj;

			for (iter_conj = chm_iter->second.begin(); iter_conj != chm_iter->second.end(); ++iter_conj)
			{
				size_t comb_offset = iter_conj->first;

				size_t prefix_end = comb_offset >> THIRTY_TWO;

				size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

				//cout << prefix_end << '\t' << suffix_st << endl;

				iter_conj->second.set_coverage();

				iter_conj->second.set_entropy();

				if (!m_chrom_dir.empty())
				{
					string flankstr;// = chromseq->substr(prefix_end, 2) + chromseq->substr(suffix_st - 3, 2);

					if (suffix_st > 2 && prefix_end < chromseq->length() + 2)
						flankstr = chromseq->substr(prefix_end, 2) + chromseq->substr(suffix_st - 3, 2);
					else
						flankstr = "NNNN";

					for (size_t i = 0; i < flankstr.length(); ++i)
					{
						if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
							flankstr[i] = flankstr[i] + 'A' - 'a';
					}

					iter_conj->second.set_flankstring(flankstr);
				}
				else
					iter_conj->second.set_flankstring();

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
			string* chromseq1;
			string* chromseq2;

			char chr1[1000], chr2[1000];

			sscanf(chm_iter->first.c_str(), "%[^~]~%[^~]", chr1, chr2);

			if (!m_chrom_dir.empty())
			{
				string chromfile1 = m_chrom_dir + chr1;

				chromfile1.append(".fa");

				if (m_loaded_chromosomes.find(chromfile1) == m_loaded_chromosomes.end())
				{
					readchrom(chromfile1.c_str(), m_loaded_chromosomes[chromfile1]);
				}

				chromseq1 = &(m_loaded_chromosomes[chromfile1]);

				if (chromseq1->empty())
				{
					cout <<"empty chrom: "<<chromfile1<<endl;
					exit(1);
				}

				string chromfile2 = m_chrom_dir + chr2;

				chromfile2.append(".fa");

				if (m_loaded_chromosomes.find(chromfile2) == m_loaded_chromosomes.end())
				{
					readchrom(chromfile2.c_str(), m_loaded_chromosomes[chromfile2]);
				}

				chromseq2 = &(m_loaded_chromosomes[chromfile2]);

				if (chromseq2->empty())
				{
					cout <<"empty chrom: "<<chromfile2<<endl;
					exit(1);
				}
			}

			size_t chrom_size1 = m_chrom_sizes[chr1];//chromseq1.size() - 1;

			size_t chrom_size2 = m_chrom_sizes[chr2];//chromseq2.size() - 1;

			JUNC_HASH_COMB::iterator iter_conj;

			for (iter_conj = chm_iter->second.begin(); iter_conj != chm_iter->second.end(); ++iter_conj)
			{
				size_t comb_offset = iter_conj->first;

				size_t prefix_end = comb_offset >> THIRTY_TWO;

				size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

				iter_conj->second.set_coverage();

				iter_conj->second.set_entropy();

				#ifdef DEBUG

				cout << prefix_end << '\t' << suffix_st << endl;
				
				#endif

				if (!m_chrom_dir.empty())
				{
					string flankstr1, flankstr2;

					if (iter_conj->second.m_strand1 == '+' && prefix_end < chromseq1->length() + 2)
						flankstr1 = chromseq1->substr(prefix_end/* + 1*/, 2);
					else if (iter_conj->second.m_strand1 == '-' && prefix_end > 2)
					{
						flankstr1 = chromseq1->substr(prefix_end - 3/* - 2*/, 2);
						flankstr1 = revcomp(flankstr1);
					}
					else
						flankstr1 = "NN";

					if (iter_conj->second.m_strand2 == '+' && suffix_st > 2)
						flankstr2 = chromseq2->substr(suffix_st - 3/* - 2*/, 2);
					else if (iter_conj->second.m_strand2 == '-' && suffix_st < chromseq2->length() + 2)
					{
						flankstr2 = chromseq2->substr(suffix_st/* + 1*/, 2);
						flankstr2 = revcomp(flankstr2);
					}
					else
						flankstr2 = "NN";

					string flankstr = flankstr1 + flankstr2; //chromseq1.substr(prefix_end, 2) + chromseq2.substr(suffix_st - 3, 2);

					for (size_t i = 0; i < flankstr.length(); ++i)
					{
						if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
							flankstr[i] = flankstr[i] + 'A' - 'a';
					}

					iter_conj->second.set_flankstring(flankstr);
				}
				else
					iter_conj->second.set_flankstring();

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

void
JunctionHandler::ReadChromSize(string chrom_size_file)
{
	ifstream ifs(chrom_size_file.c_str());

	if (ifs.is_open())
	{
		vector<bool> dumy;

		while (!ifs.eof() )
		{
			string line;

			getline(ifs, line);

			if (line.empty()/* || line[0] == '@'*/)
				continue;

			char chrom_name[1000];

			size_t chrom_size;

			sscanf(line.c_str(), "%s\t%llu", chrom_name, &chrom_size);

			m_chrom_sizes.insert(make_pair(chrom_name, chrom_size));
		}
	}
	else
	{
		cerr <<"Can't open file: " << chrom_size_file << endl;
		exit(0);
	}

	ifs.close();
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

void JunctionHandler::LoadUnfilteredJuncToSortVec()
{
	m_junc_sort.clear();

	CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter;

	JUNC_HASH_COMB::iterator junc_hash_comb_iter;

	for (chrom_junc_hash_iter = m_junc_hash.begin(); chrom_junc_hash_iter != m_junc_hash.end(); ++chrom_junc_hash_iter)
	{
		for (junc_hash_comb_iter = chrom_junc_hash_iter->second.begin(); junc_hash_comb_iter != chrom_junc_hash_iter->second.end(); ++junc_hash_comb_iter)
		{
			if (junc_hash_comb_iter->second.m_filtered_type == NOT_FILTERED && junc_hash_comb_iter->second.m_is_fusion == false && junc_hash_comb_iter->second.m_start != junc_hash_comb_iter->second.m_end)
				m_junc_sort[chrom_junc_hash_iter->first].push_back(&(junc_hash_comb_iter->second));
		}
	}
}

void JunctionHandler::Load2SimpJunc(hash_map<string, vector<pair<int, int> > >& m_original_junc_locs)
{
	m_original_junc_locs.clear();

	CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter;

	JUNC_HASH_COMB::iterator junc_hash_comb_iter;

	for (chrom_junc_hash_iter = m_junc_hash.begin(); chrom_junc_hash_iter != m_junc_hash.end(); ++chrom_junc_hash_iter)
	{
		for (junc_hash_comb_iter = chrom_junc_hash_iter->second.begin(); junc_hash_comb_iter != chrom_junc_hash_iter->second.end(); ++junc_hash_comb_iter)
		{
			m_original_junc_locs[chrom_junc_hash_iter->first].push_back(make_pair(junc_hash_comb_iter->second.m_start, junc_hash_comb_iter->second.m_end));
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
				if (junc_hash_comb_iter->second.m_filtered_type == NOT_FILTERED)
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

	cout << "m_delete_len:"<< m_delete_len << endl;

	//cout << m_junc_sort.size() << '\t' <<"chroms"<<endl;

	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		//sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
		vector<JunctionSeed*>::iterator junc_sort_iter;

		string chrom = chrom_iter->first;

		vector<JunctionSeed*>& junc_sort_vec = chrom_iter->second;

		//cout << chrom << '\t' << junc_sort_vec.size() << endl;

		for (junc_sort_iter = junc_sort_vec.begin(); junc_sort_iter != junc_sort_vec.end(); ++junc_sort_iter)
		{

			//cout << junc_sort_iter - junc_sort_vec.begin() + 1 << "th"<<endl;

			if ((*junc_sort_iter)->m_hits == 0)
				(*junc_sort_iter)->m_filtered_type = FILTERED_BY_SMALL_ANCHOR;

			//if (((*junc_sort_iter)->m_start == 144952201 && (*junc_sort_iter)->m_end == 149590973) ||
			//		((*junc_sort_iter)->m_end == 144952201 && (*junc_sort_iter)->m_start == 149590973))
			//{
			//	cerr << "write junction"<<endl;
			//	
			//	cerr << (*junc_sort_iter)->to_normal_junction(0) << endl;
			//}

			if ((*junc_sort_iter)->m_is_fusion)
			{
				if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_SMALL_ANCHOR)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_SMALL_ANCHOR"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_SMALL_DELETION)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_SMALL_DELETION"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_LARGE_MULTIPLE_PAIRED)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_LARGE_MULTIPLE_PAIRED"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_LARGE_MIN_ANCHOR_DIFF)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_LARGE_MIN_ANCHOR_DIFF"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_LOW_COVERAGE)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_LOW_COVERAGE"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_UNBALANCED_LEFT_RIGHT_PAIR)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_UNBALANCED_LEFT_RIGHT_PAIR"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_NOPAIRED)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_NOPAIRED"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_INSERTION)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_INSERTION"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_LARGE_MISMATCH)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_LARGE_MISMATCH"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_NONCAN_LEFT_RIGHT_PAIR)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_NONCAN_LEFT_RIGHT_PAIR"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_NONCAN_ENTROPY)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_NONCAN_ENTROPY"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_CAN_ENTROPY)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_CAN_ENTROPY"<<endl;
				}			
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_NONCAN_MULTI)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_NONCAN_MULTI"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_NONCAN_ERROR)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_NONCAN_ERROR"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_ISOLATED_EXON)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_ISOLATED_EXON"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_NO_EXON)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_NO_EXON"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_LARGE_MIN_MISMATCH)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_MIN_MISMATCH"<<endl;
				}
				else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_FUSION_SMALL_ENTROPY)
				{
					++filtered_count;

					ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_FUSION_SMALL_ENTROPY"<<endl;
				}
				else
				{
					++fusion_count;

					ofs_fusion << (*junc_sort_iter)->to_normal_junction(fusion_count)<<endl;
				}
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
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_CAN_ENTROPY)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_CAN_ENTROPY"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_LOW_SUPPORT)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_LOW_SUPPORT"<<endl;
			}
			else if ((*junc_sort_iter)->m_filtered_type == FILTERED_BY_ENTROPY)
			{
				++filtered_count;

				ofs_filtered << (*junc_sort_iter)->to_normal_junction(filtered_count)<< "\tFILTERED_BY_ENTROPY"<<endl;
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

				ofs_deltion << (*junc_sort_iter)->to_normal_junction(del_count)<<endl;
			}
			else
			{
				++normal_count;

				ofs_normal << (*junc_sort_iter)->to_normal_junction(normal_count)<<endl;
			}

			//cout << junc_sort_iter - junc_sort_vec.begin() + 1 << "th finished"<<endl;
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

int JunctionHandler::FindFragmentLengthReverse(vector<vector<int> >& paths, vector<pair<size_t, size_t> >& exons, size_t end_point, hash_set<size_t>& stored_fragment_lengths)
{

	#ifdef DEBUG

	cout << "FindFragmentLengthReverse paths"<<endl;

	cout << "end_point:" << end_point<<endl;

	#endif

	for (size_t j = 0; j < paths.size(); ++j)
	{
		size_t cur_fragment_length = 0;

		vector<int>& cur_path = paths[j];

		#ifdef DEBUG

		for (size_t m = 0; m < cur_path.size(); ++m)
		{
			cout << abs(cur_path[m]) - 1 << "--->";
		}

		cout << endl;

		for (size_t m = 0; m < cur_path.size(); ++m)
		{
			cout << exons[abs(cur_path[m]) - 1].first << '\t' << exons[abs(cur_path[m]) - 1].second << "--->";
		}

		cout << endl;

		#endif

		int k;

		for (k = 0; k < cur_path.size(); ++k)
		{
			if (exons[abs(cur_path[k]) - 1].first <= end_point && exons[abs(cur_path[k]) - 1].second >= end_point)
			{
				cur_fragment_length += exons[abs(cur_path[k]) - 1].second - end_point + 1;

				break;
			}
			else
				cur_fragment_length += exons[abs(cur_path[k]) - 1].second - exons[abs(cur_path[k]) - 1].first + 1;

			if (cur_fragment_length > fragment_length)
				break;
		}

		if (k != cur_path.size() && cur_fragment_length < fragment_length)
		{
			#ifdef DEBUG

			cout << "cur_fragment_length:" <<cur_fragment_length <<  endl;

			#endif

			stored_fragment_lengths.insert(cur_fragment_length);
		}
	}

	return 0; 
}

int JunctionHandler::FindFragmentLengthForward(vector<vector<int> >& paths, vector<pair<size_t, size_t> >& exons, size_t end_point, hash_set<size_t>& stored_fragment_lengths)
{
	#ifdef DEBUG
	cout << "FindFragmentLengthForward paths"<<endl;

	cout << "end_point:" << end_point<<endl;
	#endif
	for (size_t j = 0; j < paths.size(); ++j)
	{
		size_t cur_fragment_length = 0;

		vector<int>& cur_path = paths[j];

		#ifdef DEBUG

		for (size_t m = 0; m < cur_path.size(); ++m)
		{
			cout << abs(cur_path[m]) - 1  << "--->";
		}

		cout << endl;

		for (size_t m = 0; m < cur_path.size(); ++m)
		{
			cout << exons[abs(cur_path[m]) - 1].first << '\t' << exons[abs(cur_path[m]) - 1].second << "--->";
		}

		cout << endl;

		#endif

		int k;

		for (k = 0; k < cur_path.size(); ++k)
		{
			if (exons[abs(cur_path[k]) - 1].first <= end_point && exons[abs(cur_path[k]) - 1].second >= end_point)
			{
				cur_fragment_length += end_point - exons[abs(cur_path[k]) - 1].first + 1;

				break;
			}
			else
				cur_fragment_length += exons[abs(cur_path[k]) - 1].second - exons[abs(cur_path[k]) - 1].first + 1;

			if (cur_fragment_length > fragment_length)
				break;
		}

		if (k != cur_path.size() && cur_fragment_length < fragment_length)
		{
			#ifdef DEBUG

			cout << "cur_fragment_length:" <<cur_fragment_length <<  endl;

			#endif

			stored_fragment_lengths.insert(cur_fragment_length);
		}
	}

	return 0; 
}

double JunctionHandler::PrintExonExpression(vector<vector<int> >& paths, vector<pair<size_t, size_t> >& exons, vector<vector<int> >& exon_express, ofstream* ofs_ptr, double& ks_score)
{
	#ifdef DEBUG
	(*ofs_ptr) << "print exon expression paths"<<endl;
	#endif

	double min_chisquare = 100;

	ks_score = 100000;

	for (size_t j = 0; j < paths.size(); ++j)
	{
		size_t cur_fragment_length = 0;

		vector<int>& cur_path = paths[j];

		#ifdef DEBUG

		for (size_t m = 0; m < cur_path.size(); ++m)
		{
			(*ofs_ptr) << abs(cur_path[m]) - 1  << "--->";
		}

		(*ofs_ptr) << endl;

		for (size_t m = 0; m < cur_path.size(); ++m)
		{
			(*ofs_ptr) << exons[abs(cur_path[m]) - 1].first << '\t' << exons[abs(cur_path[m]) - 1].second << "--->";
		}

		(*ofs_ptr) << endl;

		#endif

		int k;

		vector<int> combined_exon_express;

		for (k = 0; k < cur_path.size(); ++k)
		{
			#ifdef DEBUG

			(*ofs_ptr) << exons[abs(cur_path[k]) - 1].first << '\t' << exons[abs(cur_path[k]) - 1].second<<endl;

			#endif
			for (size_t i = 0; i < exon_express[abs(cur_path[k]) - 1].size(); ++i)
			{
				#ifdef DEBUG
				(*ofs_ptr) << exon_express[abs(cur_path[k]) - 1][i]<<',';
				#endif

				combined_exon_express.push_back(exon_express[abs(cur_path[k]) - 1][i]);
			}

			#ifdef DEBUG
			(*ofs_ptr) << endl;
			#endif
		}

		double mean = 0, chrsquare = 0, max_difference = 0;

		KSTest(combined_exon_express, max_difference, ofs_ptr);

		GenerateMeanAndSD(combined_exon_express, mean, chrsquare);

		if (chrsquare < min_chisquare)
			min_chisquare = chrsquare;

		if (max_difference < ks_score)
			ks_score = max_difference;
	}

	
	return min_chisquare; 
}

void JunctionHandler::GenerateMeanAndSD(vector<int> expressed_regions, double& mean, double& chisqure)
{
	if (expressed_regions.size() > fragment_length)
	{
		expressed_regions.resize(fragment_length);
	}
	else if (expressed_regions.size() < fragment_length)
	{
		expressed_regions.resize(fragment_length, 0);
	}

	mean = 0;

	chisqure = 0;

	double sum = 0;

	for (size_t i = 0; i < expressed_regions.size(); ++i)
	{
		sum += expressed_regions[i];
	}

	mean = sum / (double) (expressed_regions.size());

	for (size_t i = 0; i < expressed_regions.size(); ++i)
	{
		chisqure += ((double)expressed_regions[i] - mean) * ((double)expressed_regions[i] - mean) / mean;
	}

}

void JunctionHandler::KSTest(vector<int> expressed_regions, double& max_distance, ofstream* ofs)
{
	#ifdef DEBUG
	(*ofs) << "KSTest:"<< avearge_fragment_length<<'\t'<<max_distance <<endl;

	for (size_t i = 0; i < expressed_regions.size(); ++i)
		(*ofs) << expressed_regions[i] << ',';

	(*ofs) << endl;
	#endif

	if (expressed_regions.size() > avearge_fragment_length)
		expressed_regions.resize(avearge_fragment_length);

	double sum = 0;

	max_distance = 0;

	vector<double> observed_freq(expressed_regions.size(), 0), expected_freq(expressed_regions.size(), 0), diff(expressed_regions.size(), 0);

	for (size_t i = 0; i < expressed_regions.size(); ++i)
	{
		sum += (double) expressed_regions[i];
	}

	observed_freq[0] = (double) expressed_regions[0] / sum;
	
	for (size_t i = 1; i < expressed_regions.size(); ++i)
	{
		observed_freq[i] = observed_freq[i - 1] + ((double) expressed_regions[i] / sum);
	}

	for (size_t i = 0; i < expressed_regions.size(); ++i)
	{
		expected_freq[i] = (i) / (avearge_fragment_length);
	}

	diff[0] = abs(observed_freq[0] - expected_freq[0]);

	max_distance = diff[0];

	for (size_t i = 1; i < expressed_regions.size(); ++i)
	{
		double A = abs(observed_freq[i] - expected_freq[i]);

		double B = abs(expected_freq[i] - observed_freq[i-1]);

		if (A > B)
			diff[i] = A;
		else
			diff[i] = B;

		if (diff[i] > max_distance)
			max_distance = diff[i];
	}

	//double a0 = sqrt(sum);

	//double c1 = a0 + 0.12 + (0.11/a0);
	//double d15 = 1.138/c1;
	//double d10 = 1.224/c1;
	//double d05 = 1.358/c1;
	//double d025 = 1.480/c1;

	//double t2 = max_distance;

	//max_distance = t2 * c1;
	#ifdef DEBUG
	*ofs << "max_distance:" <<max_distance<<endl;
	#endif

}

void JunctionHandler::WriteFusionJunctionWithEncompassingAlignments(string junction_fusion_file, JunctionHandler* filtered_junction_ptr)
{
	ofstream ofs_fusion(junction_fusion_file.c_str());

	string junction_fusion_file_no_compass = junction_fusion_file; junction_fusion_file_no_compass.append(".no.compass");

	ofstream ofs_fusion_no_compass(junction_fusion_file_no_compass.c_str());


	string junction_fusion_file_no_compass_short_isoformlength = junction_fusion_file; junction_fusion_file_no_compass_short_isoformlength.append(".no.compass.short_isoformlength");

	ofstream ofs_fusion_no_compass_isoformlength(junction_fusion_file_no_compass_short_isoformlength.c_str());


	string junction_fusion_file_no_compass_encompasscount = junction_fusion_file; junction_fusion_file_no_compass_encompasscount.append(".no.compass.encompasscount");

	ofstream ofs_fusion_no_compass_encompasscount(junction_fusion_file_no_compass_encompasscount.c_str());


	ofstream* ofs_fusion_no_compass_ptr;

	map<string, vector<JunctionSeed*> >::iterator chrom_iter;

	size_t fusion_count = 0;

	string fragment_file = junction_fusion_file + ".fusion_fragments.txt";

	ofstream ofs_fragment(fragment_file.c_str());

	vector<size_t> fusion_fragments;

	size_t fusion_fragments_count = 0;

	double total_fragment_length = 0;

	time_t m1, m2;

	time_t p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 = 0, p7 = 0, p8 = 0;

	time_t p9 = 0, p10 = 0, p11 = 0, p12 = 0, p13 = 0, p14 = 0, p15 = 0, p16 = 0;

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
				//stringstream fusion_junc_with_encompass;

				//fusion_junc_with_encompass.clear();

				map<size_t, size_t> doner_fragment_counts;

				size_t max_doner_fragment = 0, max_acceptor_fragment = 0;

				size_t doner_encompass_multiple = 0, doner_encompass_unique = 0;

				size_t acceptor_encompass_multiple = 0, acceptor_encompass_unique = 0;

				vector<size_t> cur_fusion_fragments;

				string chrom_name = chrom_iter->first;

				size_t comb_offset = ((*junc_sort_iter)->m_start << THIRTY_TWO) + (*junc_sort_iter)->m_end;

				hash_map<string, JUNC_HASH_COMB>::iterator junc_chrom_iter = filtered_junction_ptr->m_junc_hash.find(chrom_name);

				if (junc_chrom_iter == filtered_junction_ptr->m_junc_hash.end())
					continue;

				JUNC_HASH_COMB::iterator junc_offset_iter  = junc_chrom_iter->second.find(comb_offset);

				if (junc_offset_iter == junc_chrom_iter->second.end())
					continue;

				if (junc_offset_iter->second.m_hits == 0)
					continue;

				if (junc_offset_iter->second.m_filtered_type != NOT_FILTERED)
					continue;

				JunctionSeed& cur_junc = junc_offset_iter->second;

				m1 = time(NULL);

				(*junc_sort_iter)->generate_fusion_struct(m_fusion_encompassing_reads);

				m2 = time(NULL);

				p1 += m2 - m1;

				set<size_t> mapped_left_exon_id;

				m1 = time(NULL);

				for (size_t i = 0; i < (*junc_sort_iter)->left_paths.size(); ++i)
				{
					vector<int>& cur_path = (*junc_sort_iter)->left_paths[i];

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						mapped_left_exon_id.insert(abs(cur_path[j]) - 1);
					}
				}

				set<size_t> mapped_right_exon_id;

				for (size_t i = 0; i < (*junc_sort_iter)->right_paths.size(); ++i)
				{
					vector<int>& cur_path = (*junc_sort_iter)->right_paths[i];

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						mapped_right_exon_id.insert(abs(cur_path[j]) - 1);
					}
				}

				vector<size_t > left_arrived_exons;

				vector<size_t > right_arrived_exons;

				vector<size_t > left_exons;

				vector<size_t > right_exons;

				vector<vector<int> > left_exon_express;

				vector<vector<int> > right_exon_express;

				set<size_t>::iterator exon_iter;

				for (exon_iter = mapped_left_exon_id.begin(); exon_iter != mapped_left_exon_id.end(); ++exon_iter)
				{
					left_arrived_exons.push_back((*junc_sort_iter)->left_exons[*exon_iter].first);

					left_arrived_exons.push_back((*junc_sort_iter)->left_exons[*exon_iter].second);
				}

				for (exon_iter = mapped_right_exon_id.begin(); exon_iter != mapped_right_exon_id.end(); ++exon_iter)
				{
					right_arrived_exons.push_back((*junc_sort_iter)->right_exons[*exon_iter].first);

					right_arrived_exons.push_back((*junc_sort_iter)->right_exons[*exon_iter].second);
				}

				m2 = time(NULL);

				p2 += m2 - m1;

//#ifdef DEBUG

				ofs_fusion << "all reads before filtering"<<endl;

				ofs_fusion << "encompass reads"<<endl;

//#endif
//
//#ifdef DEBUG

				m1 = time(NULL);

				for (size_t i = 0; i < (*junc_sort_iter)->m_fusion_encompassing_reads_doner.size(); ++i)
				{

					ofs_fusion << m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_doner[i]].tostring(0,0) <<'\t';
					
					ofs_fusion << endl;

					ofs_fusion << m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_acceptor[i]].tostring(0,0) <<'\t';

					ofs_fusion << endl;

				}
//#endif
//
//#ifdef DEBUG
				ofs_fusion << "doner side spanning pairs"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->doner_side_spanning_pairs.size(); ++i)
				{
					if ((*junc_sort_iter)->doner_side_spanning_pairs[i].first.is_fusion)
						ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tostandfusion(1, 1) << endl;
					else
						ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tostring(0,0) << endl;

					//ofs_fusion << endl;

					if ((*junc_sort_iter)->doner_side_spanning_pairs[i].second.is_fusion)
						ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostandfusion(1, 1) << endl;
					else
						ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostring(0,0) <<endl;

					//ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostring(0,0) <<'\t';

					//ofs_fusion << endl;

				}
//#endif
//
//#ifdef DEBUG
				ofs_fusion << "acceptor side spanning pairs"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->accetpr_side_spanning_pairs.size(); ++i)
				{
					if ((*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.is_fusion)
						ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.tostandfusion(1, 1) <<endl;
					else
						ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.tostring(0,0) <<endl;

					//ofs_fusion << endl;

					if ((*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.is_fusion)
						ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.tostandfusion(1, 1) <<endl;
					else
						ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.tostring(0,0) <<endl;

					//ofs_fusion << endl;

				}
//#endif
//
//#ifdef DEBUG
				ofs_fusion << "single spanning"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->single_spanning.size(); ++i)
				{

					ofs_fusion << (*junc_sort_iter)->single_spanning[i].tostandfusion(1, 1) <<endl;
				}

				m2 = time(NULL);

				p3 += m2 - m1;
//#endif

				//ofs_fusion << "init left exon expression"<<endl;

				//for (size_t i = 0; i < (*junc_sort_iter)->left_exons.size(); ++i)
				//{
				//	//vector<int> tmp_express((*junc_sort_iter)->left_exons[i].second - (*junc_sort_iter)->left_exons[i].first + 1, 0);

				//	//ofs_fusion << (*junc_sort_iter)->left_exons[i].second<<'\t' <<(*junc_sort_iter)->left_exons[i].first << '\t'<<(*junc_sort_iter)->left_exons[i].second - (*junc_sort_iter)->left_exons[i].first + 1  <<endl;

				//	//left_exon_express.push_back(tmp_express);

				//	left_exons.push_back((*junc_sort_iter)->left_exons[i].first);

				//	left_exons.push_back((*junc_sort_iter)->left_exons[i].second);
				//}

				//ofs_fusion << "init right exon expression"<<endl;

				//for (size_t i = 0; i < (*junc_sort_iter)->right_exons.size(); ++i)
				//{
				//	//vector<int> tmp_express((*junc_sort_iter)->right_exons[i].second - (*junc_sort_iter)->right_exons[i].first + 1, 0);

				//	//ofs_fusion << (*junc_sort_iter)->right_exons[i].second<<'\t' <<(*junc_sort_iter)->right_exons[i].first << '\t'<<(*junc_sort_iter)->right_exons[i].second - (*junc_sort_iter)->right_exons[i].first + 1  <<endl;

				//	//right_exon_express.push_back(tmp_express);

				//	right_exons.push_back((*junc_sort_iter)->right_exons[i].first);

				//	right_exons.push_back((*junc_sort_iter)->right_exons[i].second);
				//}

				vector<SamRec*> fusion_encompassing_reads_doner, fusion_encompassing_reads_acceptor;

				size_t filtered_count = 0;

				//ofs_fusion << "doner size:" << (*junc_sort_iter)->m_fusion_encompassing_reads_doner.size()<<endl;

				//ofs_fusion << "acceptor size:" << (*junc_sort_iter)->m_fusion_encompassing_reads_acceptor.size()<<endl;

				if ((*junc_sort_iter)->m_fusion_encompassing_reads_doner.size() != (*junc_sort_iter)->m_fusion_encompassing_reads_acceptor.size())
					ofs_fusion << "doner acceptor size not equal"<<endl;


				;
				//wrtie encompass reads spanning reads to after junction

				//ofs_fusion << "encompassing reads"<<endl;

				//for (size_t i = 0; i < fusion_encompassing_reads_doner_filtered.size(); ++i)
				//{
				//	ofs_fusion << fusion_encompassing_reads_doner_filtered[i]->tostring(1, 1)<<endl;

				//	ofs_fusion << fusion_encompassing_reads_acceptor_filtered[i]->tostring(1, 1)<<endl;
				//}

				//ofs_fusion << endl;

				//ofs_fusion << "doner spanning reads"<<endl;

				//for (size_t i = 0; i < doner_side_spanning_pairs_filtered.size(); ++i)
				//{
				//	if (doner_side_spanning_pairs_filtered[i].first->is_fusion)
				//		ofs_fusion << doner_side_spanning_pairs_filtered[i].first->tostandfusion()<<endl;
				//	else
				//		ofs_fusion << doner_side_spanning_pairs_filtered[i].first->tostring(1, 1)<<endl;

				//	if (doner_side_spanning_pairs_filtered[i].second->is_fusion)
				//		ofs_fusion << doner_side_spanning_pairs_filtered[i].second->tostandfusion()<<endl;
				//	else
				//		ofs_fusion << doner_side_spanning_pairs_filtered[i].second->tostring(1, 1)<<endl;
				//}

				//ofs_fusion << endl;

				//ofs_fusion << "acceptor spanning reads"<<endl;

				//for (size_t i = 0; i < accetpr_side_spanning_pairs_filtered.size(); ++i)
				//{
				//	if (accetpr_side_spanning_pairs_filtered[i].first->is_fusion)
				//		ofs_fusion << accetpr_side_spanning_pairs_filtered[i].first->tostandfusion()<<endl;
				//	else
				//		ofs_fusion << accetpr_side_spanning_pairs_filtered[i].first->tostring(1, 1)<<endl;

				//	if (accetpr_side_spanning_pairs_filtered[i].second->is_fusion)
				//		ofs_fusion << accetpr_side_spanning_pairs_filtered[i].second->tostandfusion()<<endl;
				//	else
				//		ofs_fusion << accetpr_side_spanning_pairs_filtered[i].second->tostring(1, 1)<<endl;
				//}

				//ofs_fusion << endl;

				//ofs_fusion << "single spanning reads"<<endl;

				//for (size_t i = 0; i < (*junc_sort_iter)->single_spanning.size(); ++i)
				//{
				//	ofs_fusion << (*junc_sort_iter)->single_spanning[i].tostandfusion()<<endl;
				//}

				//ofs_fusion << endl;


				;
				//ofs_fusion << "filtered not in arrived exons alignments"<<endl;

				m1 = time(NULL);

				ofs_fusion << "fusion before filtering:"<< endl;

				ofs_fusion << cur_junc.to_normal_junction(fusion_count)<<'\t';

				for (size_t i = 0; i < (*junc_sort_iter)->left_paths.size(); ++i)
				{
					size_t struct_length = 0;

					vector<int> cur_path = (*junc_sort_iter)->left_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					ofs_fusion << (*junc_sort_iter)->left_exons[abs(cur_path[0]) - 1].first <<",";

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						ofs_fusion << (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1<<'M';

						//struct_length += (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1;

						if (j + 1 < cur_path.size())
						{
							if (cur_path[j] < 0 && cur_path[j+1] < 0)
							{
								ofs_fusion <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'P';
							}
							else
							{
								ofs_fusion <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'N';
							}
						}
					}

					ofs_fusion << '|';
				}


				ofs_fusion << '\t';

				for (size_t i = 0; i < (*junc_sort_iter)->right_paths.size(); ++i)
				{
					size_t struct_length = 0;

					vector<int> cur_path = (*junc_sort_iter)->right_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					ofs_fusion << (*junc_sort_iter)->right_exons[abs(cur_path[0]) - 1].first <<",";

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						ofs_fusion << (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1<<'M';
						
						//struct_length += (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1;

						if (j + 1 < cur_path.size())
						{
							if (cur_path[j] < 0 && cur_path[j+1] < 0)
							{
								ofs_fusion <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'P';
							}
							else
							{
								ofs_fusion <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'N';
							}
						}
					}

					ofs_fusion << '|';
				}

				ofs_fusion << endl;

				m2 = time(NULL);

				p4 += m2 - m1;

				m1 = time(NULL);

				ofs_fusion << endl << "filter encompassing reads not reached"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->m_fusion_encompassing_reads_doner.size(); ++i)
				{
					if (m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_doner[i]].tag_base_name != m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_acceptor[i]].tag_base_name)
						cout << "warning: tag base name not consistent"<<endl<<m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_doner[i]].tostring(0, 0) 
						<< m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_acceptor[i]].tostring(0, 0)<<endl;
					//cout << "is in region doner"<<endl;

					bool is_in_region = m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_doner[i]].is_in_region(left_arrived_exons);

					if (is_in_region)
					{
						//cout << "is in region acceptor"<<endl;
						is_in_region = m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_acceptor[i]].is_in_region(right_arrived_exons);
					}

					if (is_in_region)
					{
						//cout << "is_in_region"<<endl;

						//ofs_fusion << "is_in_region" << endl;

						fusion_encompassing_reads_doner.push_back(&m_fusion_encompassing_reads[((*junc_sort_iter)->m_fusion_encompassing_reads_doner[i])]);

						fusion_encompassing_reads_acceptor.push_back(&m_fusion_encompassing_reads[((*junc_sort_iter)->m_fusion_encompassing_reads_acceptor[i])]);
					}
					else
					{
						//cout << "is_not_in_region"<<endl;
						//ofs_fusion << "is_not_in_region" << endl;
						//cout <<"read name:"<<  (*junc_sort_iter)->m_fusion_encompassing_reads_doner[i].tag_name<<endl;

						//#ifdef DEBUG

						ofs_fusion << m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_doner[i]].tostring(0, 0) <<endl;

						ofs_fusion << m_fusion_encompassing_reads[(*junc_sort_iter)->m_fusion_encompassing_reads_acceptor[i]].tostring(0, 0) <<endl;

						//#endif

						++filtered_count;
					}

					//cout << "next junc_sort_iter"<<endl;
				}

				m2 = time(NULL);

				p5 += m2 - m1;

				vector<SamRec*> fusion_encompassing_reads_doner_filtered;

				vector<SamRec*> fusion_encompassing_reads_acceptor_filtered;

				vector<pair<SamRec*, SamRec*> > doner_side_spanning_pairs_filtered;

				vector<pair<SamRec*, SamRec*> > accetpr_side_spanning_pairs_filtered;


				//ofs_fusion << "encompass reads"<<endl;

				m1 = time(NULL);

				#ifdef DEBUG

				cout << "estimate encompass reads fragment length "<<endl;

				cout << "fusion_encompassing_reads_doner size:"<< fusion_encompassing_reads_doner.size()<<endl;

				#endif

				ofs_fusion << endl << "filter encompassing reads with too long fragment length"<<endl;

				for (size_t i = 0; i < fusion_encompassing_reads_doner.size(); ++i)
				{
					#ifdef DEBUG

					cout << i <<"th"<<endl;

					cout << fusion_encompassing_reads_doner[i]->tostring(0, 0) << endl;

					cout << fusion_encompassing_reads_acceptor[i]->tostring(0, 0) << endl;

					#endif

					hash_set<size_t> doner_fragments;

					hash_set<size_t> acceptor_fragments;

					vector<size_t> combined_fragments;

					size_t doner_farthest_point;
					
					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						doner_farthest_point = fusion_encompassing_reads_doner[i]->start;

						#ifdef DEBUG

						cout << "doner_farthest_point:"<<doner_farthest_point<<endl;

						#endif

						FindFragmentLengthReverse((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
					}
					else
					{
						doner_farthest_point = fusion_encompassing_reads_doner[i]->end;

						#ifdef DEBUG

						cout << "doner_farthest_point:"<<doner_farthest_point<<endl;

						#endif

						FindFragmentLengthForward((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
					}

					size_t acceptor_farthest_point;

					if ((*junc_sort_iter)->m_strand2 == '-')
					{
						acceptor_farthest_point = fusion_encompassing_reads_acceptor[i]->start;

						#ifdef DEBUG

						cout << "acceptor_farthest_point:"<<acceptor_farthest_point<<endl;

						#endif

						FindFragmentLengthReverse((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
					}
					else
					{
						acceptor_farthest_point = fusion_encompassing_reads_acceptor[i]->end;

						#ifdef DEBUG

						cout << "acceptor_farthest_point:"<<acceptor_farthest_point<<endl;

						#endif

						FindFragmentLengthForward((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
					}

					bool is_fragment_normal = false;

					hash_set<size_t>::iterator doner_fragments_iter;

					hash_set<size_t>::iterator acceptor_fragments_iter;

					for (doner_fragments_iter = doner_fragments.begin(); doner_fragments_iter != doner_fragments.end(); ++doner_fragments_iter)
					{
						for (acceptor_fragments_iter = acceptor_fragments.begin(); acceptor_fragments_iter != acceptor_fragments.end(); ++acceptor_fragments_iter)
						{
							combined_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							#ifdef DEBUG

							fusion_encompassing_reads_doner[i]->fragment_lengths.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							fusion_encompassing_reads_acceptor[i]->fragment_lengths.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							ofs_fusion <<"fragment length:"<< (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl;

							cout << "fragment length:"<< (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl;

							#endif

							if ((*doner_fragments_iter) + (*acceptor_fragments_iter) < fragment_length)
							{
								if (is_fragment_normal == false)
								{
									if (doner_fragment_counts.find((*doner_fragments_iter)) == doner_fragment_counts.end())
										doner_fragment_counts[(*doner_fragments_iter)] = 1;
									else
										++doner_fragment_counts[(*doner_fragments_iter)];
								}

								is_fragment_normal = true;

								ofs_fragment << (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl; //fusion_fragments[i] << endl;
								
								//fusion_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

								total_fragment_length += (*doner_fragments_iter) + (*acceptor_fragments_iter);

								++fusion_fragments_count;

								cur_fusion_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

								if ((*doner_fragments_iter) > max_doner_fragment)
									max_doner_fragment = (*doner_fragments_iter);

								if ((*acceptor_fragments_iter) > max_acceptor_fragment)
									max_acceptor_fragment = (*acceptor_fragments_iter);
							}
						}
					}	

					if (is_fragment_normal)
					{
						fusion_encompassing_reads_doner_filtered.push_back(fusion_encompassing_reads_doner[i]);

						fusion_encompassing_reads_acceptor_filtered.push_back(fusion_encompassing_reads_acceptor[i]);

						if (fusion_encompassing_reads_doner_filtered.back()->isunique)
							++doner_encompass_unique;
						else
							++doner_encompass_multiple;

						if (fusion_encompassing_reads_acceptor_filtered.back()->isunique)
							++acceptor_encompass_unique;
						else
							++acceptor_encompass_multiple;

						//fusion_encompassing_reads_doner_filtered.back()->generate_endpoint_expressison(left_exons, 
						//	doner_farthest_point, left_exon_express);

						//fusion_encompassing_reads_acceptor_filtered.back()->generate_endpoint_expressison(right_exons, 
						//	acceptor_farthest_point, right_exon_express);
					}
					else
					{
						ofs_fusion << fusion_encompassing_reads_doner[i]->tostring(0, 0) << endl;

						ofs_fusion << fusion_encompassing_reads_acceptor[i]->tostring(0, 0) << endl;

						//for (size_t j = 0; j < doner_fragments.size(); ++j)
						//{
						//	for (size_t k = 0; k < acceptor_fragments.size(); ++k)
						//	{
						//		ofs_fusion << "doner_fragments["<<j<<"] + acceptor_fragments["<<k<<"]:" << (*doner_fragments_iter) << "+"<< (*acceptor_fragments_iter) << endl;
						//	}
						//}						

						if (combined_fragments.empty())
							ofs_fusion <<"combined_fragments.empty()"<<endl;
						else
							ofs_fusion <<"fragment length:"<< combined_fragments.back()<< "\tmax fragment length:" << fragment_length <<endl;
					}
					
					//= fusion_encompassing_reads_acceptor[i]->end;
				}

				//ofs_fusion << "doner side spanning pairs"<<endl;

				#ifdef DEBUG

				cout << "estimate doner_side_spanning_pairs fragment length "<<endl;

				cout << "doner_side_spanning_pairs size:" << (*junc_sort_iter)->doner_side_spanning_pairs.size() << endl;

				#endif

				ofs_fusion << endl << "filter doner side reads with too long fragment length"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->doner_side_spanning_pairs.size(); ++i)
				{
					//cout << i <<"th"<<endl;

					//cout << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tag_base_name << endl;

					//cout << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tostring(0, 0) << endl;

					//cout << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostring(0, 0) << endl;

					hash_set<size_t> doner_fragments;

					hash_set<size_t> acceptor_fragments;

					vector<size_t> combined_fragments;

					size_t doner_farthest_point;
					
					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						doner_farthest_point = (*junc_sort_iter)->doner_side_spanning_pairs[i].first.start;

						#ifdef DEBUG

						cout << "doner_farthest_point:"<<doner_farthest_point<<endl;

						#endif

						FindFragmentLengthReverse((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
					}
					else
					{
						doner_farthest_point = (*junc_sort_iter)->doner_side_spanning_pairs[i].first.end;

						#ifdef DEBUG

						cout << "doner_farthest_point:"<<doner_farthest_point<<endl;

						#endif

						FindFragmentLengthForward((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
					}

					size_t acceptor_farthest_point;

					if ((*junc_sort_iter)->m_strand2 == '-')
					{
						if ((*junc_sort_iter)->doner_side_spanning_pairs[i].second.need_swap == false)
							acceptor_farthest_point = (*junc_sort_iter)->doner_side_spanning_pairs[i].second.fusion_suffix_end;
						else
							acceptor_farthest_point = (*junc_sort_iter)->doner_side_spanning_pairs[i].second.fusion_prefix_st;

						#ifdef DEBUG

						cout << "acceptor_farthest_point:"<<acceptor_farthest_point<<endl;

						#endif

						FindFragmentLengthReverse((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
					}
					else
					{
						if ((*junc_sort_iter)->doner_side_spanning_pairs[i].second.need_swap == false)
							acceptor_farthest_point = (*junc_sort_iter)->doner_side_spanning_pairs[i].second.fusion_suffix_end;
						else
							acceptor_farthest_point = (*junc_sort_iter)->doner_side_spanning_pairs[i].second.fusion_prefix_st;

						#ifdef DEBUG

						cout << "acceptor_farthest_point:"<<acceptor_farthest_point<<endl;

						#endif

						FindFragmentLengthForward((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
					}

					bool is_fragment_normal = false;

					hash_set<size_t>::iterator doner_fragments_iter;

					hash_set<size_t>::iterator acceptor_fragments_iter;

					for (doner_fragments_iter = doner_fragments.begin(); doner_fragments_iter != doner_fragments.end(); ++doner_fragments_iter)
					{
						for (acceptor_fragments_iter = acceptor_fragments.begin(); acceptor_fragments_iter != acceptor_fragments.end(); ++acceptor_fragments_iter)
						{
							combined_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							#ifdef DEBUG

							(*junc_sort_iter)->doner_side_spanning_pairs[i].first.fragment_lengths.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							(*junc_sort_iter)->doner_side_spanning_pairs[i].second.fragment_lengths.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							ofs_fusion << "fragment length:" << (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl;

							cout << "fragment length:"<< (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl;

							#endif

							if ((*doner_fragments_iter) + (*acceptor_fragments_iter) < fragment_length)
							{
								if (is_fragment_normal == false)
								{
									if (doner_fragment_counts.find((*doner_fragments_iter)) == doner_fragment_counts.end())
										doner_fragment_counts[(*doner_fragments_iter)] = 1;
									else
										++doner_fragment_counts[(*doner_fragments_iter)];
								}

								is_fragment_normal = true;

								ofs_fragment << (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl; //fusion_fragments[i] << endl;

								//fusion_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

								total_fragment_length += (*doner_fragments_iter) + (*acceptor_fragments_iter);

								++fusion_fragments_count;

								cur_fusion_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

								if ((*doner_fragments_iter) > max_doner_fragment)
									max_doner_fragment = (*doner_fragments_iter);

								if ((*acceptor_fragments_iter) > max_acceptor_fragment)
									max_acceptor_fragment = (*acceptor_fragments_iter);
							}
						}
					}

					if (is_fragment_normal)
					{
						doner_side_spanning_pairs_filtered.push_back(make_pair( &((*junc_sort_iter)->doner_side_spanning_pairs[i].first), &((*junc_sort_iter)->doner_side_spanning_pairs[i].second)));

						//doner_side_spanning_pairs_filtered.back().first->generate_endpoint_expressison(left_exons, 
						//	doner_farthest_point, left_exon_express);

						//doner_side_spanning_pairs_filtered.back().second->generate_endpoint_expressison(right_exons, 
						//	acceptor_farthest_point, right_exon_express);
					}
					else
					{
						if ((*junc_sort_iter)->doner_side_spanning_pairs[i].first.is_fusion)
							ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tostandfusion(1, 1) <<endl;
						else
							ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tostring(0, 0) << endl;

						if ((*junc_sort_iter)->doner_side_spanning_pairs[i].second.is_fusion)
							ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostandfusion(1, 1) <<endl;
						else
							ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostring(0, 0) << endl;

						//for (size_t j = 0; j < doner_fragments.size(); ++j)
						//{
						//	for (size_t k = 0; k < acceptor_fragments.size(); ++k)
						//	{
						//		ofs_fusion << "doner_fragments["<<j<<"] + acceptor_fragments["<<k<<"]:" << (*doner_fragments_iter) << "+"<< (*acceptor_fragments_iter) << endl;
						//	}
						//}	
						
						if (combined_fragments.empty())
							ofs_fusion <<"combined_fragments.empty()"<<endl;
						else
							ofs_fusion <<"fragment length:"<< combined_fragments.back()<< "\tmax fragment length:" << fragment_length <<endl;
					}
				}


				
				//ofs_fusion << "acceptor side spanning pairs"<<endl;

				//cout << "estimate accetpr_side_spanning_pairs fragment length "<<endl;

				ofs_fusion << endl << "filter acceptor side reads with too long fragment length"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->accetpr_side_spanning_pairs.size(); ++i)
				{
					#ifdef DEBUG
					cout << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.tostring(0, 0) << endl;

					cout << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.tostring(0, 0) << endl;
					#endif

					hash_set<size_t> doner_fragments;

					hash_set<size_t> acceptor_fragments;

					vector<size_t> combined_fragments;

					size_t doner_farthest_point;
					
					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						if ((*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.need_swap == false)
							doner_farthest_point = (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fusion_prefix_st;
						else
							doner_farthest_point = (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fusion_suffix_end;

						#ifdef DEBUG

						cout << "doner_farthest_point:"<<doner_farthest_point<<endl;

						#endif

						FindFragmentLengthReverse((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
					}
					else
					{
						if ((*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.need_swap == false)
							doner_farthest_point = (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fusion_prefix_st;
						else
							doner_farthest_point = (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fusion_suffix_end;

						#ifdef DEBUG

						cout << "doner_farthest_point:"<<doner_farthest_point<<endl;

						#endif

						FindFragmentLengthForward((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
					}

					size_t acceptor_farthest_point;

					if ((*junc_sort_iter)->m_strand2 == '-')
					{
						acceptor_farthest_point = (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.start;

						#ifdef DEBUG

						cout << "acceptor_farthest_point:"<<acceptor_farthest_point<<endl;

						#endif

						FindFragmentLengthReverse((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
					}
					else
					{
						acceptor_farthest_point = (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.end;

						#ifdef DEBUG

						cout << "acceptor_farthest_point:"<<acceptor_farthest_point<<endl;

						#endif

						FindFragmentLengthForward((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
					}

					bool is_fragment_normal = false;

					hash_set<size_t>::iterator doner_fragments_iter;

					hash_set<size_t>::iterator acceptor_fragments_iter;

					for (doner_fragments_iter = doner_fragments.begin(); doner_fragments_iter != doner_fragments.end(); ++doner_fragments_iter)
					{
						for (acceptor_fragments_iter = acceptor_fragments.begin(); acceptor_fragments_iter != acceptor_fragments.end(); ++acceptor_fragments_iter)
						{
							combined_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							#ifdef DEBUG

							(*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fragment_lengths.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							(*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.fragment_lengths.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

							ofs_fusion << (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl;

							cout << "fragment length:"<< (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl;

							#endif

							if ((*doner_fragments_iter) + (*acceptor_fragments_iter) < fragment_length)
							{
								if (is_fragment_normal == false)
								{
									if (doner_fragment_counts.find((*doner_fragments_iter)) == doner_fragment_counts.end())
										doner_fragment_counts[(*doner_fragments_iter)] = 1;
									else
										++doner_fragment_counts[(*doner_fragments_iter)];
								}

								is_fragment_normal = true;

								ofs_fragment << (*doner_fragments_iter) + (*acceptor_fragments_iter) << endl; //fusion_fragments[i] << endl;

								//fusion_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

								total_fragment_length += (*doner_fragments_iter) + (*acceptor_fragments_iter);

								++fusion_fragments_count;

								cur_fusion_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));

								if ((*doner_fragments_iter) > max_doner_fragment)
									max_doner_fragment = (*doner_fragments_iter);

								if ((*acceptor_fragments_iter) > max_acceptor_fragment)
									max_acceptor_fragment = (*acceptor_fragments_iter);
							}
						}
					}

					if (is_fragment_normal)
					{
						accetpr_side_spanning_pairs_filtered.push_back(make_pair( &((*junc_sort_iter)->accetpr_side_spanning_pairs[i].first), &((*junc_sort_iter)->accetpr_side_spanning_pairs[i].second)));

						//accetpr_side_spanning_pairs_filtered.back().first->generate_endpoint_expressison(left_exons, 
						//	doner_farthest_point, left_exon_express);

						//accetpr_side_spanning_pairs_filtered.back().second->generate_endpoint_expressison(right_exons, 
						//	acceptor_farthest_point, right_exon_express);
					}
					else
					{
						if ((*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.is_fusion)
							ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.tostandfusion(1, 1) <<endl;
						else
							ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.tostring(0, 0) << endl;

						if ((*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.is_fusion)
							ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.tostandfusion(1, 1) <<endl;
						else
							ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.tostring(0, 0) << endl;

						//for (size_t j = 0; j < doner_fragments.size(); ++j)
						//{
						//	for (size_t k = 0; k < acceptor_fragments.size(); ++k)
						//	{
						//		ofs_fusion << "doner_fragments["<<j<<"] + acceptor_fragments["<<k<<"]:" << (*doner_fragments_iter) << "+"<< (*acceptor_fragments_iter) << endl;
						//	}
						//}	

						if (combined_fragments.empty())
							ofs_fusion <<"combined_fragments.empty()"<<endl;
						else
							ofs_fusion <<"fragment length:"<< combined_fragments.back() << "\tmax fragment length:" << fragment_length <<endl;
					}
				}

				ofs_fusion << "fusion_fragments.size():"<<fusion_fragments_count<< endl;

				ofs_fusion << "total_fragment_length:"<<total_fragment_length<< endl;

				m2 = time(NULL);

				p6 += m2 - m1;

				//ofs_fusion << "single spanning"<<endl;

				//for (size_t i = 0; i < (*junc_sort_iter)->single_spanning.size(); ++i)
				//{
				//	vector<size_t> doner_fragments;

				//	vector<size_t> acceptor_fragments;

				//	vector<size_t> combined_fragments;

				//	size_t doner_farthest_point;
				//	
				//	if ((*junc_sort_iter)->m_strand1 == '+')
				//	{
				//		doner_farthest_point = (*junc_sort_iter)->single_spanning[i].fusion_prefix_st;

				//		FindFragmentLengthReverse((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
				//	}
				//	else
				//	{
				//		doner_farthest_point = (*junc_sort_iter)->single_spanning[i].fusion_prefix_st;

				//		FindFragmentLengthForward((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, doner_farthest_point, doner_fragments);
				//	}

				//	size_t acceptor_farthest_point;

				//	if ((*junc_sort_iter)->m_strand2 == '-')
				//	{
				//		acceptor_farthest_point = (*junc_sort_iter)->single_spanning[i].fusion_suffix_end;

				//		FindFragmentLengthReverse((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
				//	}
				//	else
				//	{
				//		acceptor_farthest_point = (*junc_sort_iter)->single_spanning[i].fusion_suffix_end;

				//		FindFragmentLengthForward((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, acceptor_farthest_point, acceptor_fragments);
				//	}

				//	for (size_t j = 0; j < doner_fragments.size(); ++j)
				//	{
				//		for (size_t k = 0; k < acceptor_fragments.size(); ++k)
				//		{
				//			combined_fragments.push_back((*doner_fragments_iter) + (*acceptor_fragments_iter));
				//		}
				//	}
				//}

				//cout << "encompass alignment filtered"<<endl;

				//ofs_fusion << filtered_count <<" encompass alignment filtered"<<endl;


				//print exon expression

				;

				//////print exon expression

				//double chi_square1, chi_square2;

				//double ks_score1, ks_score2;

				//chi_square1 = PrintExonExpression((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, left_exon_express, &ofs_fusion, ks_score1);

				//chi_square2 = PrintExonExpression((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, right_exon_express, &ofs_fusion, ks_score2);
				//
				//double uniformity1 = alglib::chisquarecdistribution(fragment_length - 1, chi_square1);

				//double uniformity2 = alglib::chisquarecdistribution(fragment_length - 1, chi_square2);

				//cur_junc.m_encompass_reads_count = fusion_encompassing_reads_doner_filtered.size() * 2;

				//cur_junc.m_left_paired_count = doner_side_spanning_pairs_filtered.size();

				//cur_junc.m_right_paired_count = accetpr_side_spanning_pairs_filtered.size();
				//
				//cur_junc.m_hits = cur_junc.m_left_paired_count + cur_junc.m_right_paired_count + cur_junc.m_single_count;

				++fusion_count;

				/////////////////

				#ifdef DEBUG

				ofs_fusion << "encompass reads"<<endl;

				for (size_t i = 0; i < fusion_encompassing_reads_doner.size(); ++i)
				{
					ofs_fusion << fusion_encompassing_reads_doner[i]->tostring(0,0) <<'\t';
					
					ofs_fusion << "ZF:Z:FRAG:";

					for (size_t j = 0; j < fusion_encompassing_reads_doner[i]->fragment_lengths.size(); ++j)
						ofs_fusion << fusion_encompassing_reads_doner[i]->fragment_lengths[j] << ',';

					ofs_fusion << endl;

					ofs_fusion << fusion_encompassing_reads_acceptor[i]->tostring(0,0) <<'\t';

					ofs_fusion << "ZF:Z:FRAG:";

					for (size_t j = 0; j < fusion_encompassing_reads_doner[i]->fragment_lengths.size(); ++j)
						ofs_fusion << fusion_encompassing_reads_doner[i]->fragment_lengths[j] << ',';

					ofs_fusion << endl;
				}

				#endif

				#ifdef DEBUG

				ofs_fusion << "doner side spanning pairs"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->doner_side_spanning_pairs.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.tostring(0,0) <<'\t';

					ofs_fusion << "ZF:Z:FRAG:";

					for (size_t j = 0; j < (*junc_sort_iter)->doner_side_spanning_pairs[i].first.fragment_lengths.size(); ++j)
						ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].first.fragment_lengths[j] << ',';

					ofs_fusion << endl;

					ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.tostring(0,0) <<'\t';

					ofs_fusion << "ZF:Z:FRAG:";

					for (size_t j = 0; j < (*junc_sort_iter)->doner_side_spanning_pairs[i].second.fragment_lengths.size(); ++j)
						ofs_fusion << (*junc_sort_iter)->doner_side_spanning_pairs[i].second.fragment_lengths[j] << ',';

					ofs_fusion << endl;
				}

				#endif

				#ifdef DEBUG

				ofs_fusion << "acceptor side spanning pairs"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->accetpr_side_spanning_pairs.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.tostring(0,0) <<'\t';

					ofs_fusion << "ZF:Z:FRAG:";

					for (size_t j = 0; j < (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fragment_lengths.size(); ++j)
						ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].first.fragment_lengths[j] << ',';

					ofs_fusion << endl;

					ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.tostring(0,0) <<'\t';

					ofs_fusion << "ZF:Z:FRAG:";

					for (size_t j = 0; j < (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.fragment_lengths.size(); ++j)
						ofs_fusion << (*junc_sort_iter)->accetpr_side_spanning_pairs[i].second.fragment_lengths[j] << ',';

					ofs_fusion << endl;
				}

				#endif

				#ifdef DEBUG

				ofs_fusion << "single spanning"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->single_spanning.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->single_spanning[i].tostring(0,0) <<endl;
				}

				#endif

				//reset splice ways

				//cout << "reset splice ways"<<endl;

				m1 = time(NULL);

				(*junc_sort_iter)->reset_splice_ways(fusion_encompassing_reads_doner_filtered,		
													 fusion_encompassing_reads_acceptor_filtered,
													 doner_side_spanning_pairs_filtered,
													 accetpr_side_spanning_pairs_filtered,
													 (*junc_sort_iter)->single_spanning);

				//cout << "left splice ways" << endl;
			
				#ifdef DEBUG

				for (size_t i = 0; i < (*junc_sort_iter)->left_splice_ways.size(); ++i)
				{
					cout << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < (*junc_sort_iter)->left_splice_ways[i].spliceway_vec.size(); ++j)
						cout << (*junc_sort_iter)->left_splice_ways[i].spliceway_vec[j].first<<':' << (*junc_sort_iter)->left_splice_ways[i].spliceway_vec[j].second << '\t';

					cout << endl;
				}

				#endif

				//cout << 
				if ((*junc_sort_iter)->left_splice_ways.size())
				{
					#ifdef DEBUG

					cout << "spliceways2exons_doner"<<endl;

					#endif

					(*junc_sort_iter)->spliceways2exons_doner();

					vector<int> doner_array((*junc_sort_iter)->left_exons.size(), -1);

					#ifdef DEBUG

					cout << "doner_graph"<<endl;

					#endif

					vector<vector<int> > doner_graph((*junc_sort_iter)->left_exons.size(), doner_array);

					#ifdef DEBUG

					cout << "construct_graph"<<endl;

					#endif

					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						(*junc_sort_iter)->construct_graph((*junc_sort_iter)->left_exons, (*junc_sort_iter)->left_splice_ways, doner_graph, '-');

						#ifdef DEBUG

						cout << "DFS"<<endl;

						#endif

						DFS(doner_graph,(*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons.size() - 1);
					}
					else
					{
						(*junc_sort_iter)->construct_graph((*junc_sort_iter)->left_exons, (*junc_sort_iter)->left_splice_ways, doner_graph, '+');

						#ifdef DEBUG

						cout << "DFS"<<endl;

						#endif

						DFS(doner_graph, (*junc_sort_iter)->left_paths, 0);
					}

					#ifdef DEBUG

					cout << "finish" << endl;

					#endif
				}


				#ifdef DEBUG

				cout << "right splice ways" << endl;

				#endif

				for (size_t i = 0; i < (*junc_sort_iter)->right_splice_ways.size(); ++i)
				{
					#ifdef DEBUG

					cout << (*junc_sort_iter)->right_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->right_splice_ways[i].start<<'\t';

					#endif

					#ifdef DEBUG

					for (size_t j = 0; j < (*junc_sort_iter)->right_splice_ways[i].spliceway_vec.size(); ++j)
						cout << (*junc_sort_iter)->right_splice_ways[i].spliceway_vec[j].first<<':' << (*junc_sort_iter)->right_splice_ways[i].spliceway_vec[j].second << '\t';

					cout << endl;

					#endif
				}

				 
				if ((*junc_sort_iter)->right_splice_ways.size())
				{
					//cout << "spliceways2exons_acceptor" << endl;

					#ifdef DEBUG

					cout << "spliceways2exons_doner"<<endl;

					#endif

					(*junc_sort_iter)->spliceways2exons_acceptor();

					vector<int> acceptor_array((*junc_sort_iter)->right_exons.size(), -1);

					#ifdef DEBUG

					cout << "acceptor_graph"<<endl;

					#endif

					vector<vector<int> > acceptor_graph((*junc_sort_iter)->right_exons.size(), acceptor_array);

					//cout << "acceptor_graph" << endl;

					if ((*junc_sort_iter)->m_strand2 == '+')
					{
						//cout << "construct_graph" << endl;

						(*junc_sort_iter)->construct_graph((*junc_sort_iter)->right_exons, (*junc_sort_iter)->right_splice_ways, acceptor_graph, '+');

						#ifdef DEBUG

						cout << "DFS" << endl;

						#endif

						DFS(acceptor_graph, (*junc_sort_iter)->right_paths, 0);

						#ifdef DEBUG

						cout << "finish" << endl;

						#endif
					}
					else
					{
						//cout << "construct_graph" << endl;

						(*junc_sort_iter)->construct_graph((*junc_sort_iter)->right_exons, (*junc_sort_iter)->right_splice_ways, acceptor_graph, '-');

						#ifdef DEBUG

						cout << "DFS" << endl;

						#endif

						DFS(acceptor_graph, (*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons.size() - 1);

						#ifdef DEBUG

						cout << "finish" << endl;

						#endif
					}
				}

				m2 = time(NULL);

				p7 += m2 - m1;

				//ofs_fusion << "init left exon expression"<<endl;

				//for (size_t i = 0; i < (*junc_sort_iter)->left_exons.size(); ++i)
				//{
				//	//vector<int> tmp_express((*junc_sort_iter)->left_exons[i].second - (*junc_sort_iter)->left_exons[i].first + 1, 0);

				//	//ofs_fusion << (*junc_sort_iter)->left_exons[i].second<<'\t' <<(*junc_sort_iter)->left_exons[i].first << '\t'<<(*junc_sort_iter)->left_exons[i].second - (*junc_sort_iter)->left_exons[i].first + 1  <<endl;

				//	//left_exon_express.push_back(tmp_express);

				//	left_exons.push_back((*junc_sort_iter)->left_exons[i].first);

				//	left_exons.push_back((*junc_sort_iter)->left_exons[i].second);
				//}

				//ofs_fusion << "init right exon expression"<<endl;

				//for (size_t i = 0; i < (*junc_sort_iter)->right_exons.size(); ++i)
				//{
				//	//vector<int> tmp_express((*junc_sort_iter)->right_exons[i].second - (*junc_sort_iter)->right_exons[i].first + 1, 0);

				//	//ofs_fusion << (*junc_sort_iter)->right_exons[i].second<<'\t' <<(*junc_sort_iter)->right_exons[i].first << '\t'<<(*junc_sort_iter)->right_exons[i].second - (*junc_sort_iter)->right_exons[i].first + 1  <<endl;

				//	//right_exon_express.push_back(tmp_express);

				//	right_exons.push_back((*junc_sort_iter)->right_exons[i].first);

				//	right_exons.push_back((*junc_sort_iter)->right_exons[i].second);
				//}


				m1 = time(NULL);

				#ifdef DEBUG

				ofs_fusion << "init left exon expression"<<endl;

				#endif

				for (size_t i = 0; i < (*junc_sort_iter)->left_exons.size(); ++i)
				{
					vector<int> tmp_express((*junc_sort_iter)->left_exons[i].second - (*junc_sort_iter)->left_exons[i].first + 1, 0);

					left_exon_express.push_back(tmp_express);

					left_exons.push_back((*junc_sort_iter)->left_exons[i].first);

					left_exons.push_back((*junc_sort_iter)->left_exons[i].second);
				}

				#ifdef DEBUG

				ofs_fusion << "init right exon expression"<<endl;

				#endif

				for (size_t i = 0; i < (*junc_sort_iter)->right_exons.size(); ++i)
				{
					vector<int> tmp_express((*junc_sort_iter)->right_exons[i].second - (*junc_sort_iter)->right_exons[i].first + 1, 0);

					right_exon_express.push_back(tmp_express);

					right_exons.push_back((*junc_sort_iter)->right_exons[i].first);

					right_exons.push_back((*junc_sort_iter)->right_exons[i].second);
				}

				for (size_t i = 0; i < fusion_encompassing_reads_doner_filtered.size(); ++i)
				{
					size_t doner_farthest_point;
					
					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						doner_farthest_point = fusion_encompassing_reads_doner_filtered[i]->start;
					}
					else
					{
						doner_farthest_point = fusion_encompassing_reads_doner_filtered[i]->end;
					}

					size_t acceptor_farthest_point;

					if ((*junc_sort_iter)->m_strand2 == '-')
					{
						acceptor_farthest_point = fusion_encompassing_reads_acceptor_filtered[i]->start;
					}
					else
					{
						acceptor_farthest_point = fusion_encompassing_reads_acceptor_filtered[i]->end;
					}
	
					fusion_encompassing_reads_doner_filtered[i]->generate_endpoint_expressison(left_exons, 
						doner_farthest_point, left_exon_express);

					fusion_encompassing_reads_acceptor_filtered[i]->generate_endpoint_expressison(right_exons, 
						acceptor_farthest_point, right_exon_express);
				}

				#ifdef DEBUG

				ofs_fusion << "doner side spanning pairs"<<endl;

				#endif

				for (size_t i = 0; i < doner_side_spanning_pairs_filtered.size(); ++i)
				{
					size_t doner_farthest_point;
					
					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						doner_farthest_point = doner_side_spanning_pairs_filtered[i].first->start;
					}
					else
					{
						doner_farthest_point = doner_side_spanning_pairs_filtered[i].first->end;
					}

					size_t acceptor_farthest_point;

					if ((*junc_sort_iter)->m_strand2 == '-')
					{
						if (doner_side_spanning_pairs_filtered[i].second->need_swap == false)
							acceptor_farthest_point = doner_side_spanning_pairs_filtered[i].second->fusion_suffix_end;
						else
							acceptor_farthest_point = doner_side_spanning_pairs_filtered[i].second->fusion_prefix_st;
					}
					else
					{
						if (doner_side_spanning_pairs_filtered[i].second->need_swap == false)
							acceptor_farthest_point = doner_side_spanning_pairs_filtered[i].second->fusion_suffix_end;
						else
							acceptor_farthest_point = doner_side_spanning_pairs_filtered[i].second->fusion_prefix_st;
					}

					doner_side_spanning_pairs_filtered[i].first->generate_endpoint_expressison(left_exons, 
						doner_farthest_point, left_exon_express);

					doner_side_spanning_pairs_filtered[i].second->generate_endpoint_expressison(right_exons, 
						acceptor_farthest_point, right_exon_express);
				}

				#ifdef DEBUG

				ofs_fusion << "acceptor side spanning pairs"<<endl;

				#endif

				for (size_t i = 0; i < accetpr_side_spanning_pairs_filtered.size(); ++i)
				{
					size_t doner_farthest_point;
					
					if ((*junc_sort_iter)->m_strand1 == '+')
					{
						if (accetpr_side_spanning_pairs_filtered[i].first->need_swap == false)
							doner_farthest_point = accetpr_side_spanning_pairs_filtered[i].first->fusion_prefix_st;
						else
							doner_farthest_point = accetpr_side_spanning_pairs_filtered[i].first->fusion_suffix_end;
					}
					else
					{
						if (accetpr_side_spanning_pairs_filtered[i].first->need_swap == false)
							doner_farthest_point = accetpr_side_spanning_pairs_filtered[i].first->fusion_prefix_st;
						else
							doner_farthest_point = accetpr_side_spanning_pairs_filtered[i].first->fusion_suffix_end;
					}

					size_t acceptor_farthest_point;

					if ((*junc_sort_iter)->m_strand2 == '-')
					{
						acceptor_farthest_point = accetpr_side_spanning_pairs_filtered[i].second->start;
					}
					else
					{
						acceptor_farthest_point = accetpr_side_spanning_pairs_filtered[i].second->end;
					}

					accetpr_side_spanning_pairs_filtered[i].first->generate_endpoint_expressison(left_exons, 
						doner_farthest_point, left_exon_express);

					accetpr_side_spanning_pairs_filtered[i].second->generate_endpoint_expressison(right_exons, 
						acceptor_farthest_point, right_exon_express);

				}

				//print exon expression

				double chi_square1, chi_square2;

				double ks_score1, ks_score2;

				chi_square1 = PrintExonExpression((*junc_sort_iter)->left_paths, (*junc_sort_iter)->left_exons, left_exon_express, &ofs_fusion, ks_score1);

				chi_square2 = PrintExonExpression((*junc_sort_iter)->right_paths, (*junc_sort_iter)->right_exons, right_exon_express, &ofs_fusion, ks_score2);
				
				m2 = time(NULL);

				p8 += m2 - m1;

				m1 = time(NULL);

				double uniformity1 = 0;//alglib::chisquarecdistribution(fragment_length - 1, chi_square1);

				double uniformity2 = 0;//alglib::chisquarecdistribution(fragment_length - 1, chi_square2);

				cur_junc.m_encompass_reads_count = fusion_encompassing_reads_doner_filtered.size() * 2;

				cur_junc.m_left_paired_count = doner_side_spanning_pairs_filtered.size();

				cur_junc.m_right_paired_count = accetpr_side_spanning_pairs_filtered.size();
				
				cur_junc.m_hits = cur_junc.m_left_paired_count + cur_junc.m_right_paired_count + cur_junc.m_single_count;

				map<size_t, size_t>::iterator doner_fragment_iter;

				size_t total_paired_reads = cur_junc.m_left_paired_count + cur_junc.m_right_paired_count + fusion_encompassing_reads_doner_filtered.size();

				double paired_reads_entropy = 0;

				for (doner_fragment_iter = doner_fragment_counts.begin(); doner_fragment_iter != doner_fragment_counts.end(); ++doner_fragment_iter)
				{
					double pi =  (double) doner_fragment_iter->second / (double)total_paired_reads;

					paired_reads_entropy += pi * log(pi);
				}

				if ( paired_reads_entropy != 0)
					paired_reads_entropy = -paired_reads_entropy;

				//
				double m_pq_score = 0;

				double ppower = pow(0.25, double(max_doner_fragment));

				double pNpower = pow(1.0 - ppower, (double)30000000);

				double qpower = pow(0.25, double(max_acceptor_fragment));

				double pDpower = pow(1.0 - qpower, (double)10000);

				double lpq = 1.0 - (pNpower * pDpower);

				double ppower2 = pow(0.25, double(max_doner_fragment));

				double pNpower2 = pow(1.0 - ppower2, (double)10000 );

				double qpower2 = pow(0.25, double(max_acceptor_fragment));

				double pDpower2 = pow(1.0 - qpower2, (double)30000000);

				double lpq2 = 1.0 - (pNpower2 * pDpower2);

				m_pq_score = 1.0 - (lpq + lpq2) / 2;

				double max_cur_fragment = 0, min_cur_fragment = 1000000, ave_cur_fragment = 0, sum_cur_fragment = 0;

				for (size_t i = 0; i < cur_fusion_fragments.size(); ++i)
				{
					if (max_cur_fragment < cur_fusion_fragments[i])
						max_cur_fragment = cur_fusion_fragments[i];

					if (min_cur_fragment > cur_fusion_fragments[i])
						min_cur_fragment = cur_fusion_fragments[i];

					sum_cur_fragment += cur_fusion_fragments[i];
				}

				if (cur_fusion_fragments.size())
					ave_cur_fragment = sum_cur_fragment / (double)cur_fusion_fragments.size();
				
				//
				#ifdef DEBUG

				cout <<"reset splice way finished"<<endl;
				
				#endif

				//if ((*junc_sort_iter)->m_fusion_encompassing_reads.size())
				//{
				//	(*junc_sort_iter)->m_paired_count += (*junc_sort_iter)->m_fusion_encompassing_reads.size();
				//}

				//

				size_t max_doner_struct = 0, min_doner_struct = -1, max_acceptor_struct = 0, min_acceptor_struct = -1;


				for (size_t i = 0; i < (*junc_sort_iter)->left_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					size_t struct_length = 0;

					vector<int> cur_path = (*junc_sort_iter)->left_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						struct_length += (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1;
					}

					if (struct_length > max_doner_struct)
						max_doner_struct = struct_length;

					if (struct_length < min_doner_struct)
						min_doner_struct = struct_length;
				}


				for (size_t i = 0; i < (*junc_sort_iter)->right_paths.size(); ++i)
				{
					size_t struct_length = 0;

					vector<int> cur_path = (*junc_sort_iter)->right_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						struct_length += (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1;
					}

					if (struct_length > max_acceptor_struct)
						max_acceptor_struct = struct_length;

					if (struct_length < min_acceptor_struct)
						min_acceptor_struct = struct_length;
				}


				if ((max_acceptor_struct < min_isoform_length) || (max_doner_struct < min_isoform_length))
				{
					ofs_fusion_no_compass_ptr = &ofs_fusion_no_compass_isoformlength;
				}
				else if ((*junc_sort_iter)->m_encompass_reads_count < min_encompass_count)
				{
					ofs_fusion_no_compass_ptr = &ofs_fusion_no_compass_encompasscount;
				}
				else
					ofs_fusion_no_compass_ptr = &ofs_fusion_no_compass;

				//#ifdef DEBUG

				ofs_fusion << cur_junc.to_normal_junction(fusion_count)<<'\t';

				//#endif

				(*ofs_fusion_no_compass_ptr) << cur_junc.to_normal_junction(fusion_count)<<'\t';//<<endl;

				/*size_t max_doner_struct = 0, min_doner_struct = -1, max_acceptor_struct = 0, min_acceptor_struct = -1;*/

				for (size_t i = 0; i < (*junc_sort_iter)->left_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					size_t struct_length = 0;

					vector<int> cur_path = (*junc_sort_iter)->left_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					(*ofs_fusion_no_compass_ptr) << (*junc_sort_iter)->left_exons[abs(cur_path[0]) - 1].first <<",";

					//#ifdef DEBUG

					ofs_fusion << (*junc_sort_iter)->left_exons[abs(cur_path[0]) - 1].first <<",";

					//#endif

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						(*ofs_fusion_no_compass_ptr) << (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1<<'M';

						//#ifdef DEBUG

						ofs_fusion << (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1<<'M';

						//#endif

						struct_length += (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].first + 1;

						if (j + 1 < cur_path.size())
						{
							if (cur_path[j] < 0 && cur_path[j+1] < 0)
							{
								(*ofs_fusion_no_compass_ptr) <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'P';

								//#ifdef DEBUG

								ofs_fusion <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'P';

								//#endif
							}
							else
							{
								(*ofs_fusion_no_compass_ptr) <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'N';

								//#ifdef DEBUG

								ofs_fusion <<  (*junc_sort_iter)->left_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->left_exons[abs(cur_path[j]) - 1].second - 1 << 'N';

								//#endif
							}

						}
					}

					if (struct_length > max_doner_struct)
						max_doner_struct = struct_length;

					if (struct_length < min_doner_struct)
						min_doner_struct = struct_length;

					(*ofs_fusion_no_compass_ptr) << '|';

					//#ifdef DEBUG

					ofs_fusion << '|';

					//#endif
				}

				(*ofs_fusion_no_compass_ptr) << '\t';

				//#ifdef DEBUG

				ofs_fusion << '\t';

				//#endif

				for (size_t i = 0; i < (*junc_sort_iter)->right_paths.size(); ++i)
				{
					//ofs_fusion << (*junc_sort_iter)->left_splice_ways[i].chrom_name<<'\t'<< (*junc_sort_iter)->left_splice_ways[i].start<<'\t';

					//for (size_t j = 0; j < (*junc_sort_iter)->right_paths[i].size(); ++j)
					//	ofs_fusion_no_compass << (*junc_sort_iter)->right_exons[(*junc_sort_iter)->right_paths[i][j]].first<<':' << 
					//	(*junc_sort_iter)->right_exons[(*junc_sort_iter)->right_paths[i][j]].second << '=';

					//ofs_fusion_no_compass << '~';

					size_t struct_length = 0;

					vector<int> cur_path = (*junc_sort_iter)->right_paths[i];

					sort(cur_path.begin(), cur_path.end(), less_than);

					(*ofs_fusion_no_compass_ptr) << (*junc_sort_iter)->right_exons[abs(cur_path[0]) - 1].first <<",";

					//#ifdef DEBUG

					ofs_fusion << (*junc_sort_iter)->right_exons[abs(cur_path[0]) - 1].first <<",";

					//#endif

					for (size_t j = 0; j < cur_path.size(); ++j)
					{
						(*ofs_fusion_no_compass_ptr) << (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1<<'M';

						//#ifdef DEBUG

						ofs_fusion << (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1<<'M';

						//#endif
						
						struct_length += (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].first + 1;

						if (j + 1 < cur_path.size())
						{
							if (cur_path[j] < 0 && cur_path[j+1] < 0)
							{
								(*ofs_fusion_no_compass_ptr) <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'P';

								//#ifdef DEBUG

								ofs_fusion <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'P';

								//#endif
							}
							else
							{
								(*ofs_fusion_no_compass_ptr) <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'N';

								//#ifdef DEBUG

								ofs_fusion <<  (*junc_sort_iter)->right_exons[abs(cur_path[j+1]) - 1].first - (*junc_sort_iter)->right_exons[abs(cur_path[j]) - 1].second - 1 << 'N';

								//#endif
							}
						}

					}

					if (struct_length > max_acceptor_struct)
						max_acceptor_struct = struct_length;

					if (struct_length < min_acceptor_struct)
						min_acceptor_struct = struct_length;

					(*ofs_fusion_no_compass_ptr) << '|';

					//#ifdef DEBUG

					ofs_fusion << '|';

					//#endif
				}

				(*ofs_fusion_no_compass_ptr) <<'\t'<< uniformity1 <<'\t' <<uniformity2 <<'\t' <<ks_score1 << '\t' <<ks_score2;

				//#ifdef DEBUG

				ofs_fusion <<'\t'<<uniformity1  <<'\t' <<uniformity2<<'\t' <<ks_score1 << '\t' <<ks_score2;

				//#endif

				(*ofs_fusion_no_compass_ptr) <<'\t'<< max_doner_struct <<'\t' <<min_doner_struct <<'\t' <<max_acceptor_struct << '\t' <<min_acceptor_struct;

				//#ifdef DEBUG

				ofs_fusion <<'\t'<< max_doner_struct <<'\t' <<min_doner_struct <<'\t' <<max_acceptor_struct << '\t' <<min_acceptor_struct;

				(*ofs_fusion_no_compass_ptr) <<'\t'<< paired_reads_entropy<<'\t' << cur_junc.m_ave_mismatch / m_max_read_width <<'\t' << m_pq_score << '\t' << max_doner_fragment<<'\t' << max_acceptor_fragment;

				ofs_fusion <<'\t'<< paired_reads_entropy<<'\t' << cur_junc.m_ave_mismatch / m_max_read_width <<'\t' << m_pq_score << '\t' << max_doner_fragment<<'\t' << max_acceptor_fragment;
				
				//#endif

				(*ofs_fusion_no_compass_ptr) << '\t' <<max_cur_fragment<< '\t' <<  min_cur_fragment << '\t' <<  ave_cur_fragment;

				ofs_fusion << '\t' <<max_cur_fragment<< '\t' <<  min_cur_fragment << '\t' <<  ave_cur_fragment;

				(*ofs_fusion_no_compass_ptr) << '\t' <<doner_encompass_unique << '\t' <<  doner_encompass_multiple << '\t' <<  acceptor_encompass_unique<< '\t' <<  acceptor_encompass_multiple;

				ofs_fusion << '\t' <<doner_encompass_unique << '\t' <<  doner_encompass_multiple << '\t' <<  acceptor_encompass_unique<< '\t' <<  acceptor_encompass_multiple;

				(*ofs_fusion_no_compass_ptr) << endl;

				//#ifdef DEBUG

				ofs_fusion << endl;;

				m2 = time(NULL);

				p9 += m2 - m1;
				//#endif


				//print doner fragment length

				m1 = time(NULL);

				ofs_fusion << "doner fragment length"<<endl;

				for (doner_fragment_iter = doner_fragment_counts.begin(); doner_fragment_iter != doner_fragment_counts.end(); ++doner_fragment_iter)
				{
					ofs_fusion << doner_fragment_iter->first <<'\t' << doner_fragment_iter->second << endl;
				}

				//wrtie encompass reads spanning reads to after junction

				ofs_fusion << "encompassing reads"<<endl;

				for (size_t i = 0; i < fusion_encompassing_reads_doner_filtered.size(); ++i)
				{
					ofs_fusion << fusion_encompassing_reads_doner_filtered[i]->tostring(1, 1)<<endl;

					ofs_fusion << fusion_encompassing_reads_acceptor_filtered[i]->tostring(1, 1)<<endl;
				}

				ofs_fusion << endl;

				ofs_fusion << "doner spanning reads"<<endl;

				for (size_t i = 0; i < doner_side_spanning_pairs_filtered.size(); ++i)
				{
					if (doner_side_spanning_pairs_filtered[i].first->is_fusion)
						ofs_fusion << doner_side_spanning_pairs_filtered[i].first->tostandfusion(1, 1) <<endl;
					else
						ofs_fusion << doner_side_spanning_pairs_filtered[i].first->tostring(1, 1)<<endl;

					if (doner_side_spanning_pairs_filtered[i].second->is_fusion)
						ofs_fusion << doner_side_spanning_pairs_filtered[i].second->tostandfusion(1, 1) <<endl;
					else
						ofs_fusion << doner_side_spanning_pairs_filtered[i].second->tostring(1, 1)<<endl;
				}

				ofs_fusion << endl;

				ofs_fusion << "acceptor spanning reads"<<endl;

				for (size_t i = 0; i < accetpr_side_spanning_pairs_filtered.size(); ++i)
				{
					if (accetpr_side_spanning_pairs_filtered[i].first->is_fusion)
						ofs_fusion << accetpr_side_spanning_pairs_filtered[i].first->tostandfusion(1, 1) <<endl;
					else
						ofs_fusion << accetpr_side_spanning_pairs_filtered[i].first->tostring(1, 1)<<endl;

					if (accetpr_side_spanning_pairs_filtered[i].second->is_fusion)
						ofs_fusion << accetpr_side_spanning_pairs_filtered[i].second->tostandfusion(1, 1) <<endl;
					else
						ofs_fusion << accetpr_side_spanning_pairs_filtered[i].second->tostring(1, 1)<<endl;
				}

				ofs_fusion << endl;

				ofs_fusion << "single spanning reads"<<endl;

				for (size_t i = 0; i < (*junc_sort_iter)->single_spanning.size(); ++i)
				{
					ofs_fusion << (*junc_sort_iter)->single_spanning[i].tostandfusion(1, 1) <<endl;
				}

				ofs_fusion << endl;

				//if ((max_acceptor_struct < min_isoform_length) || (max_doner_struct < min_isoform_length))
				//{
				//	ofs_fusion_no_compass_ptr = &ofs_fusion_no_compass_isoformlength;
				//}
				//else if ((*junc_sort_iter)->m_encompass_reads_count < min_encompass_count)
				//{
				//	ofs_fusion_no_compass_ptr = &ofs_fusion_no_compass_encompasscount;
				//}
				//else
				//	ofs_fusion_no_compass_ptr = &ofs_fusion_no_compass;

				//(*ofs_fusion_no_compass_ptr) << fusion_junc_with_encompass.str();				

				//fusion_junc_with_encompass.clear();

				m2 = time(NULL);

				p10 += m2 - m1;

				#ifdef DEBUG

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

				#endif

				m1 = time(NULL);

				(*junc_sort_iter)->clear_splice_ways();

				m2 = time(NULL);

				p11 += m2 - m1;
			}
		}
	}

	cerr << "p1:"<<p1<<endl;
	cerr << "p2:"<<p2<<endl;
	cerr << "p3:"<<p3<<endl;
	cerr << "p4:"<<p4<<endl;
	cerr << "p5:"<<p5<<endl;
	cerr << "p6:"<<p6<<endl;
	cerr << "p7:"<<p7<<endl;
	cerr << "p8:"<<p8<<endl;
	cerr << "p9:"<<p9<<endl;
	cerr << "p10:"<<p10<<endl;
	cerr << "p11:"<<p11<<endl;
	//for (size_t i = 0; i < fusion_fragments.size(); ++i)
	//{
	//	ofs_fragment << fusion_fragments[i] << endl;
	//}
}



void JunctionHandler::GenerateFusionStruct()
{
	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	JUNC_HASH_COMB::iterator offset_iter;

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		//cout << 1 << endl;

		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{

			//cout << 2 << endl;

			if (offset_iter->second.m_is_fusion)
			{
				#ifdef DEBUG

				cout << "junction "<<endl;

				cout << offset_iter->second.to_normal_junction(0)<<endl;

				cout << "generate splice ways"<<endl;

				#endif

				offset_iter->second.reset_splice_ways(m_fusion_encompassing_reads);

				#ifdef DEBUG

				cout << "left splice ways" << endl;

				#endif

				#ifdef DEBUG

				for (size_t i = 0; i < offset_iter->second.left_splice_ways.size(); ++i)
				{
					cout << offset_iter->second.left_splice_ways[i].chrom_name<<'\t'<< offset_iter->second.left_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < offset_iter->second.left_splice_ways[i].spliceway_vec.size(); ++j)
						cout << offset_iter->second.left_splice_ways[i].spliceway_vec[j].first<<':' << offset_iter->second.left_splice_ways[i].spliceway_vec[j].second << '\t';

					cout << endl;
				}

				#endif

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

				#ifdef DEBUG

				cout << "right splice ways" << endl;

				for (size_t i = 0; i < offset_iter->second.right_splice_ways.size(); ++i)
				{
					cout << offset_iter->second.right_splice_ways[i].chrom_name<<'\t'<< offset_iter->second.right_splice_ways[i].start<<'\t';

					for (size_t j = 0; j < offset_iter->second.right_splice_ways[i].spliceway_vec.size(); ++j)
						cout << offset_iter->second.right_splice_ways[i].spliceway_vec[j].first<<':' << offset_iter->second.right_splice_ways[i].spliceway_vec[j].second << '\t';

					cout << endl;
				}

				#endif
				 
				if (offset_iter->second.right_splice_ways.size())
				{
					//cout << "spliceways2exons_acceptor" << endl;

					offset_iter->second.spliceways2exons_acceptor();

					vector<int> acceptor_array(offset_iter->second.right_exons.size(), -1);

					vector<vector<int> > acceptor_graph(offset_iter->second.right_exons.size(), acceptor_array);

					//cout << "acceptor_graph" << endl;

					if (offset_iter->second.m_strand2 == '+')
					{
						//cout << "construct_graph" << endl;

						offset_iter->second.construct_graph(offset_iter->second.right_exons, offset_iter->second.right_splice_ways, acceptor_graph, '+');

						//cout << "DFS" << endl;

						DFS(acceptor_graph, offset_iter->second.right_paths, 0);
					}
					else
					{
						//cout << "construct_graph" << endl;

						offset_iter->second.construct_graph(offset_iter->second.right_exons, offset_iter->second.right_splice_ways, acceptor_graph, '-');

						//cout << "DFS" << endl;

						DFS(acceptor_graph, offset_iter->second.right_paths, offset_iter->second.right_exons.size() - 1);

						//cout << "finish" << endl;
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
				if (offset_iter->second.m_hits == 0)
					;
				else if (offset_iter->second.m_is_fusion)
				{
					++m_fusion;

					if (offset_iter->second.m_flankcase >= 5)
						++m_fusion_canon;
					else if (offset_iter->second.m_flankcase >= 1)
						++m_fusion_semicanon;
					else if (offset_iter->second.m_flankcase == 0)
						++m_fusion_noncanon;
				}
				else if (offset_iter->second.m_start == offset_iter->second.m_end)
					++m_small_ins;
				else if (offset_iter->second.m_end - offset_iter->second.m_start < (size_t) m_delete_len )
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

	m_fusion = 0;

	m_fusion_canon = 0;
	
	m_fusion_semicanon = 0;
	
	m_fusion_noncanon = 0;

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

	ofs << "non_canonical_junction\t"<< m_non_canonical << endl;

	ofs << "fusion\t"<< m_fusion << endl;

	ofs << "fusion_canon\t"<< m_fusion_canon << endl;

	ofs << "fusion_semi_canon\t"<< m_fusion_semicanon << endl;

	ofs << "fusion_non_canon\t"<< m_fusion_noncanon << endl;
}

void JunctionHandler::MarkFiltered(bool paired, JunctionHandler& junc_db)
{
	CHROM_JUNC_HASH_COMB::iterator chrom_iter;

	JUNC_HASH_COMB::iterator offset_iter;

	hash_map<string, hash_map<size_t, int> > doner_not_filtered;

	hash_map<string, hash_map<size_t, int> > acceptor_not_filtered;

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{
			if (offset_iter->second.m_is_fusion)
			{
				//if ((offset_iter->second.m_start == 144952201 && offset_iter->second.m_end == 149590973) ||
				//	(offset_iter->second.m_end == 144952201 && offset_iter->second.m_start == 149590973))
				//	cerr << offset_iter->second.to_normal_junction(0) << endl;

				if ((offset_iter->second.m_max_prefix_len < m_min_junc_anchor || offset_iter->second.m_max_suffix_len < m_min_junc_anchor)
					/* &&(offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) */)
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_SMALL_ANCHOR;
				//else if (offset_iter->second.m_end - offset_iter->second.m_start <= (size_t) m_delete_len && 
				//	offset_iter->second.m_end - offset_iter->second.m_start > 0/* && offset_iter->second.m_flankcase < 5*/)
				//	offset_iter->second.m_filtered_type = FILTERED_BY_SMALL_DELETION;
				//else if (offset_iter->second.m_min_anchor_difference > max_anchor_diff && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
				//	offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MIN_ANCHOR_DIFF;
				//			else if (paired && offset_iter->second.m_paired_count && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) )
				//			{
				///*				if ((double)offset_iter->second.m_paired_unique_count / (double) offset_iter->second.m_paired_count < 0.6)
				//					offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MULTIPLE_PAIRED;			
				//				else*/ //if (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)
				//					offset_iter->second.m_filtered_type = FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR;
				//			}
				//else if (offset_iter->second.m_paired_count == 0 && paired)
				//{
				//	offset_iter->second.m_filtered_type = FILTERED_BY_NOPAIRED;
				//}
				else if (offset_iter->second.m_ave_mismatch > 1 &&  (offset_iter->second.m_entropy < 2)/*(offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)*/)
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_LARGE_MISMATCH;
				//else if (!isindel && paired && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_left_paired_count < 2 || offset_iter->second.m_right_paired_count < 2))
				//{
				//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_LEFT_RIGHT_PAIR;
				//}
				else if (/*!isindel &&*/ offset_iter->second.m_flankcase > 4 && (offset_iter->second.m_entropy < 2.5) && ((offset_iter->second.m_paired_count == 0 && paired))
					/* && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)*/)
				{
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_CAN_ENTROPY;
				}
				else if (/*!isindel && */offset_iter->second.m_flankcase <= 4 && (offset_iter->second.m_entropy < 2.5) && ((offset_iter->second.m_paired_count == 0 && paired))
					)
				{
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_NONCAN_ENTROPY;
				}
				else if (/*!isindel && */offset_iter->second.m_hits <= m_min_fusion_coverage)
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_LOW_COVERAGE;
				else if (/*!isindel && */offset_iter->second.m_min_mismatch > 1 && offset_iter->second.m_entropy < 2)
				{
					//if((offset_iter->second.m_start == 46414792 || offset_iter->second.m_start == 64085601) && (offset_iter->second.m_end == 46414792 || offset_iter->second.m_end == 64085601))
					//	cout << 46414792 <<"\t" <<64085601 << "\tmismatch:" <<offset_iter->second.m_min_mismatch  <<"\tentropy:" <<offset_iter->second.m_entropy <<endl  ;
						
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_LARGE_MIN_MISMATCH;
				}
				else if (/*!isindel && */offset_iter->second.m_entropy < 0.001)
					offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_SMALL_ENTROPY;
				//else if ( (offset_iter->second.m_chrom.find("chrM") != string::npos) ||  (offset_iter->second.m_chrom2.find("chrM") != string::npos) || 
				//	      (offset_iter->second.m_chrom.find("chrUn") != string::npos) ||  (offset_iter->second.m_chrom2.find("chrUn") != string::npos))
				//	offset_iter->second.m_filtered_type = FILTERED_BY_FUSION_SMALL_ENTROPY;

				//else if (/*!isindel &&*/ offset_iter->second.m_flankcase <= 4 && (offset_iter->second.m_entropy < 2.5) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
				//{
				//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ENTROPY;
				//}
				//else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_multi_count > 0) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
				//{
				//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_MULTI;
				//}
				//else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_ave_mismatch > 1) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
				//{
				//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ERROR;
				//}

				continue;
			}
				

			if (m_do_filter & 1024)
				continue;

			bool isindel = false;

			if (offset_iter->second.m_start == offset_iter->second.m_end || offset_iter->second.m_end - offset_iter->second.m_start <= (size_t) m_delete_len)
				isindel = true;

			bool is_annotated = false;

			CHROM_JUNC_HASH_COMB::iterator cjhc_iter = junc_db.m_junc_hash.find(chrom_iter->first);

			if (cjhc_iter != junc_db.m_junc_hash.end())
			{
				JUNC_HASH_COMB::iterator jhc_iter = cjhc_iter->second.find(offset_iter->first);

				if (jhc_iter != cjhc_iter->second.end())
					is_annotated = true;
			}

/*			if (offset_iter->second.m_end == offset_iter->second.m_start)
				offset_iter->second.m_filtered_type = FILTERED_BY_INSERTION;
			else */
			if (is_annotated)
				;
			else if ((offset_iter->second.m_max_prefix_len < m_min_junc_anchor || offset_iter->second.m_max_suffix_len < m_min_junc_anchor)/* && 
				(offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) */)
				offset_iter->second.m_filtered_type = FILTERED_BY_SMALL_ANCHOR;
			//else if (offset_iter->second.m_end - offset_iter->second.m_start <= (size_t) m_delete_len && 
			//	offset_iter->second.m_end - offset_iter->second.m_start > 0/* && offset_iter->second.m_flankcase < 5*/)
			//	offset_iter->second.m_filtered_type = FILTERED_BY_SMALL_DELETION;
			//else if (offset_iter->second.m_min_anchor_difference > max_anchor_diff && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
			//	offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MIN_ANCHOR_DIFF;
//			else if (paired && offset_iter->second.m_paired_count && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) )
//			{
///*				if ((double)offset_iter->second.m_paired_unique_count / (double) offset_iter->second.m_paired_count < 0.6)
//					offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MULTIPLE_PAIRED;			
//				else*/ //if (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)
//					offset_iter->second.m_filtered_type = FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR;
//			}
			//else if (offset_iter->second.m_paired_count == 0 && paired)
			//{
			//	offset_iter->second.m_filtered_type = FILTERED_BY_NOPAIRED;
			//}
			else if (offset_iter->second.m_ave_mismatch > 2 &&  (offset_iter->second.m_entropy < 2.5)/*(offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)*/)
				offset_iter->second.m_filtered_type = FILTERED_BY_LARGE_MISMATCH;
			else if (offset_iter->second.m_entropy <= min_entropy)
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_ENTROPY;
			}
			//else if (!isindel && paired && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_left_paired_count < 2 || offset_iter->second.m_right_paired_count < 2))
			//{
			//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_LEFT_RIGHT_PAIR;
			//}
			/*else if (!isindel && offset_iter->second.m_flankcase > 4 && (offset_iter->second.m_entropy < 2.5) && 
				(( ((m_do_filter & 256) == 0 || m_max_read_width >= 75) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)) ||
				(((m_do_filter & 256)) && (offset_iter->second.m_left_paired_count == 0 && offset_iter->second.m_right_paired_count == 0)))	&& 
				paired)*/
			else if (paired && !isindel && offset_iter->second.m_flankcase > 4 && (offset_iter->second.m_entropy < 2.5) && 
				offset_iter->second.m_left_paired_count == 0 && offset_iter->second.m_right_paired_count == 0)
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_CAN_ENTROPY;
			}
			else if((m_do_filter & 2048) && offset_iter->second.m_hits < low_support_threshold && (offset_iter->second.m_max_prefix_len < 15 || offset_iter->second.m_max_suffix_len < 15 || offset_iter->second.m_ave_mismatch > 0))
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_LOW_SUPPORT;
			}
			else if (!isindel && offset_iter->second.m_flankcase <= 4 && (offset_iter->second.m_entropy < 2.5))
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ENTROPY;
			}
			else if (!isindel && offset_iter->second.m_flankcase <= 4 && (offset_iter->second.m_entropy < 2.5) && 
				((((m_do_filter & 256) == 0 || m_max_read_width >= 75 ) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0)) ||
				(((m_do_filter & 256)) && (offset_iter->second.m_left_paired_count == 0 && offset_iter->second.m_right_paired_count == 0)))	&& 
				/*(offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0) && */
				paired)
			{
				offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ENTROPY;
			}
			//else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_multi_count > 0) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
			//{
			//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_MULTI;
			//}
			//else if (!isindel && offset_iter->second.m_flankcase < 4 && (offset_iter->second.m_ave_mismatch > 1) && (offset_iter->second.m_left_paired_count == 0 || offset_iter->second.m_right_paired_count == 0))
			//{
			//	offset_iter->second.m_filtered_type = FILTERED_BY_NONCAN_ERROR;
			//}

			if (!isindel && offset_iter->second.m_filtered_type == NOT_FILTERED)
			{
				doner_not_filtered[offset_iter->second.m_chrom][offset_iter->second.m_start] = 0;

				acceptor_not_filtered[offset_iter->second.m_chrom][offset_iter->second.m_end] = 0;
			}

			//offset_iter->second.m_five_prime_known_id.clear();

			//offset_iter->second.m_mapped_idx.clear();

			//offset_iter->second.m_prefix_count.clear();

			//offset_iter->second.m_three_prime_known_id.clear();
		}
	}

	for (chrom_iter = m_junc_hash.begin(); chrom_iter != m_junc_hash.end(); ++chrom_iter)
	{
		for (offset_iter = chrom_iter->second.begin(); offset_iter != chrom_iter->second.end(); ++offset_iter)
		{
			if (offset_iter->second.m_is_fusion)
				continue;

			if (offset_iter->second.m_filtered_type != NOT_FILTERED)
			{
				//if ((offset_iter->second.m_max_prefix_len < m_min_junc_anchor || offset_iter->second.m_max_suffix_len < m_min_junc_anchor) || offset_iter->second.m_entropy == 0)
				//	continue;

				bool doner_confirmed = false;

				bool acceptor_confirmed = false;

				{
					hash_map<string, hash_map<size_t, int> >::iterator c_iter = doner_not_filtered.find(offset_iter->second.m_chrom);

					if (c_iter != doner_not_filtered.end())
					{
						hash_map<size_t, int>::iterator off_iter = c_iter->second.find(offset_iter->second.m_start);

						if (off_iter != c_iter->second.end())
						{
							doner_confirmed = true;
							//offset_iter->second.m_filtered_type = NOT_FILTERED;
						}
						else
							continue;
					}
					else
						continue;
				}

				{
					hash_map<string, hash_map<size_t, int> >::iterator c_iter = acceptor_not_filtered.find(offset_iter->second.m_chrom);

					if (c_iter != acceptor_not_filtered.end())
					{
						hash_map<size_t, int>::iterator off_iter = c_iter->second.find(offset_iter->second.m_end);

						if (off_iter != c_iter->second.end())
						{
							acceptor_confirmed = true;
							//offset_iter->second.m_filtered_type = NOT_FILTERED;
						}
						else
							continue;
					}
					else
							continue;
				}

				offset_iter->second.m_filtered_type = NOT_FILTERED;
			}
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

		sam_rec->filter_type = junc_hash_comb_iter->second.m_filtered_type;

		//sam_rec->corresponding_juncs.push_back(&junc_hash_comb_iter->second);

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

void JunctionHandler::OutPutJuncSeq(size_t junc_id, string chrom, string& chromseq, size_t max_seq_thresh, ofstream& out_fs)
{
	map<size_t, vector<Junc_Seq> >::iterator it_head = m_head_seq.find(junc_id);
	map<size_t, vector<Junc_Seq> >::iterator it_tail = m_tail_seq.find(junc_id);
	if(it_head != m_head_seq.end() && it_tail != m_tail_seq.end())
	{
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
		it_head->second.clear();
		it_tail->second.clear();
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
	map<size_t, vector<Junc_Seq> >::iterator it = m_head_seq.find(my_junc_seq.junc_id);
	if(it != m_head_seq.end() && m_head_seq[my_junc_seq.junc_id].size() >= max_seq_thresh)
		return;
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
	map<size_t, vector<Junc_Seq> >::iterator it = m_tail_seq.find(my_junc_seq.junc_id);
	if(it != m_tail_seq.end() && m_tail_seq[my_junc_seq.junc_id].size() >= max_seq_thresh)
		return;
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
	map<string, vector<JunctionSeed*> > m_junc_sort2(m_junc_sort);
	map<string, vector<JunctionSeed*> >::iterator chrom_iter, chrom_iter2;
	for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
	{
		string chromfile = m_chrom_dir + chrom_iter->first;
		chromfile.append(".fa");
		string chromseq;
		readchrom(chromfile.c_str(), chromseq);
		size_t chrom_len = chromseq.length();
		chrom_iter2 = m_junc_sort2.find(chrom_iter->first);
		sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
		sort(chrom_iter2->second.begin(), chrom_iter2->second.end(), comp_junc_sort_by_end);
		for(size_t i = 0; i< chrom_iter->second.size(); i++)
		{
			Junc_Seq new_junc_seq(chrom_iter->second[i]->m_strand, chrom_iter->second[i]->m_flankcase, chrom_iter->second[i]->m_junc_id, chrom_iter->second[i]->m_end);
			GetJuncTail(new_junc_seq, chrom_iter->second, i, min_anchor, max_anchor, chrom_len, max_seq_thresh);
			
			Jump_Code new_head_code(chrom_iter2->second[i]->m_end - chrom_iter2->second[i]->m_start -1, N);
			Junc_Seq new_junc_seq2(chrom_iter2->second[i]->m_strand, chrom_iter2->second[i]->m_flankcase, chrom_iter2->second[i]->m_junc_id, new_head_code);
			GetJuncHead(new_junc_seq2, chrom_iter2->second, i, min_anchor, max_anchor, max_seq_thresh);
			
			OutPutJuncSeq(chrom_iter->second[i]->m_junc_id, chrom_iter->first, chromseq, max_seq_thresh,out_fs);
			OutPutJuncSeq(chrom_iter2->second[i]->m_junc_id, chrom_iter->first, chromseq, max_seq_thresh,out_fs);
		}
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
	if (m_is_fusion)
	{
		if (m_sam_rec_ptr->need_swap)
		{
			m_fusion_suffix_len = m_sam_rec_ptr->fusion_prefix_len;

			m_fusion_prefix_len = m_sam_rec_ptr->fusion_suffix_len;
		}
		else
		{
			m_fusion_prefix_len = m_sam_rec_ptr->fusion_prefix_len;

			m_fusion_suffix_len = m_sam_rec_ptr->fusion_suffix_len;
		}
	}
	else
	{
		m_fusion_prefix_len = 0;

		m_fusion_suffix_len = 0;
	}
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

	if (m_is_fusion)
	{
		m_fusion_prefix_len = sam_rec_ptr->fusion_prefix_len;

		m_fusion_suffix_len = sam_rec_ptr->fusion_suffix_len;
	}
	else
	{
		m_fusion_prefix_len = 0;

		m_fusion_suffix_len = 0;
	}
}
