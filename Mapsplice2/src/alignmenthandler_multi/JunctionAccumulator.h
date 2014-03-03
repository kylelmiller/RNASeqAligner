#ifndef JUNCTIONACCUMULATOR_H
#define JUNCTIONACCUMULATOR_H

#include "sharedlib.h"
#include "SamRec.h"

struct JuncAccum{
	JuncAccum(/*int pm, *//*const string fs, */size_t loc, size_t suffix_len, size_t rw, size_t tagidx, unsigned short mis, size_t strand, string insert = "") : 
/*prim(pm), *//*flankstr(fs),*/ p(rw - 1, 0), positive_count(0), negative_count(0), entropy(0), coverage(0), flankcase(0), il_score(0)
	{
		++p[loc-1];
		max_prefix_len = loc;
		max_suffix_len = suffix_len;

		max_mismatch = mis;

		min_mismatch = mis;

		sum_mismatch = mis;

		m[tagidx] = 1;

		if (insert == "")
			;
		else
			ins[insert] = 1;

		if (strand & IS_REVERSE)
			++negative_count;
		else
			++positive_count;

	}
	bool inc_hits(size_t idx, size_t suffix_len, size_t tagidx, unsigned short mis, size_t strand, string insert = "")
	{
		if (tagidx != -1 && m.find(tagidx) != m.end())
			return false;

		if (insert == "")
			;
		else if (ins.find(insert) == ins.end())
			ins[insert] = 1;
		else
			++ins[insert];			

		++p[idx-1];

		m[tagidx] = 1;

		if (max_prefix_len < idx)
			max_prefix_len = idx;

		if (max_suffix_len < suffix_len)
			max_suffix_len = suffix_len;

		if (mis > max_mismatch)
			max_mismatch = mis;
		
		if (mis < min_mismatch)
			min_mismatch = mis;

		sum_mismatch += mis;

		if (strand & IS_REVERSE)
			++negative_count;
		else
			++positive_count;

		return true;
	}

	//char get_strand() const
	//{
	//	return positive_count > negative_count ? '+' : '-';
	//}

	void set_coverage()
	{
		coverage = positive_count + negative_count;
	}

	void set_entropy()
	{
		for (size_t i = 0; i < p.size(); ++i)
		{
			if (p[i] > 0)
			{
				double pi = p[i] / (double)coverage;
				entropy += pi * log(pi);
			}
		}

		if (entropy != 0)
			entropy = -entropy;
	}

	void set_flankstring(const string& flankstring)
	{
		//char strand = '+';

		flankstr = flankstring;

		if (flankstr == "ATAC")
		{
			flankcase = 1;
			strand = '+';
		}
		else if (flankstr == "CTAC")
		{
			flankcase = 6;
			strand = '-';
		}
		else if (flankstr == "CTGC")
		{
			flankcase = 3;
			strand = '-';
		}
		else if (flankstr == "GCAG")
		{
			flankcase = 4;
			strand = '+';
		}
		else if (flankstr == "GTAG")
		{
			flankcase = 5;
			strand = '+';
		}
		else if (flankstr == "GTAT")
		{
			flankcase = 2;
			strand = '-';
		}
	}

	void set_pq_score(size_t intron_len, size_t chrom_size)
	{
		double ppower = pow(0.25, double(max_prefix_len));

		double pNpower = pow(1.0 - ppower, (double)chrom_size);

		double qpower = pow(0.25, double(max_suffix_len));

		double pDpower = pow(1.0 - qpower, (double)intron_len);

		double lpq = 1.0 - (pNpower * pDpower);

		double ppower2 = pow(0.25, double(max_prefix_len));

		double pNpower2 = pow(1.0 - ppower2, (double)intron_len );

		double qpower2 = pow(0.25, double(max_suffix_len));

		double pDpower2 = pow(1.0 - qpower2, (double)chrom_size);

		double lpq2 = 1.0 - (pNpower2 * pDpower2);

		lpqave = 1.0 - (lpq + lpq2) / 2;
	}

	void set_il_score(size_t junc_st, size_t junc_end)
	{
		il_score = 1.0 - (((double) (junc_end - junc_st - 1 + 1)) / (double (200000 - 1 + 2)));
	}

	void set_ave_mis()
	{
		ave_mismatch = (double)sum_mismatch / (double)coverage;
	}

	unsigned short flankcase;
	string flankstr;
	size_t max_prefix_len;
	size_t max_suffix_len;
	vector<unsigned short> p;

	double lpqave;

	double il_score;

	char strand;

	double entropy;

	unsigned int positive_count, negative_count;

	unsigned int coverage;

	unsigned short max_mismatch;

	unsigned short min_mismatch;

	unsigned short sum_mismatch;

	double ave_mismatch;

	map<size_t, int> m;

	map<string, int> ins;

};

typedef hash_map<size_t, JuncAccum> JUNC_ACCUM_HASH_COMB;

typedef hash_map<string, JUNC_ACCUM_HASH_COMB> CHROM_JUNC_ACCUM_HASH_COMB;

struct JuncSort {

	JuncSort(size_t st, size_t ed, JuncAccum* jptr) : start(st), end(ed), junc_accum_ptr(jptr) 
	{
	}

	char* to_normal_junction(size_t junc_id)
	{
		char buf[5000];

		sprintf(buf, "%llu\t%llu\tJUNC_%llu\t%u\t%c\t%llu\t%llu\t255,0,0\t2\t%llu,%llu,\t0,%llu,\t%lf\t%hu\t%s\t%lf\t%lf\t%hu\t%hu\t%lf", start, end, junc_id, junc_accum_ptr->coverage, junc_accum_ptr->strand, start, end, 
			junc_accum_ptr->max_prefix_len, junc_accum_ptr->max_suffix_len, end + junc_accum_ptr->max_suffix_len - start + 1, junc_accum_ptr->entropy, junc_accum_ptr->flankcase, junc_accum_ptr->flankstr.c_str(),
			junc_accum_ptr->il_score, junc_accum_ptr->lpqave, junc_accum_ptr->min_mismatch, junc_accum_ptr->max_mismatch, junc_accum_ptr->ave_mismatch);

	}

	char* to_normal_junction_bed(size_t junc_id)
	{
	}

	char* to_insert_junction(size_t junc_id)
	{
		char buf[5000];

		string ins_str;

		map<string, int>::iterator ins_iter;
		for (ins_iter = junc_accum_ptr->ins.begin(); ins_iter != junc_accum_ptr->ins.end(); ++ins_iter)
		{
			ins_str.append(ins_iter->first);
			ins_str.append("-");
			char intbuf[10];
			ins_str.append(itoa(ins_iter->second, intbuf, 10));
			ins_str.append(",");
		}

		sprintf(buf, "%llu\t%llu\tJUNC_%llu\t%u\t%c\t%llu\t%llu\t255,0,0\t2\t%llu,%llu,\t0,%llu,\t%lf\t%hu\t%s\t%lf\t%lf\t%hu\t%hu\t%lf\t%s", start, end, junc_id, junc_accum_ptr->coverage, junc_accum_ptr->strand, start, end, 
			junc_accum_ptr->max_prefix_len, junc_accum_ptr->max_suffix_len, end + junc_accum_ptr->max_suffix_len - start + 1, junc_accum_ptr->entropy, junc_accum_ptr->flankcase, junc_accum_ptr->flankstr.c_str(),
			junc_accum_ptr->il_score, junc_accum_ptr->lpqave, junc_accum_ptr->min_mismatch, junc_accum_ptr->max_mismatch, junc_accum_ptr->ave_mismatch, ins_str.c_str());
	}
	
	size_t start;
	size_t end;
	JuncAccum* junc_accum_ptr;
};

bool comp_junc_sort(const JuncSort& lhs, const JuncSort& rhs)
{
	if (lhs.start == rhs.start)
		return lhs.end < rhs.end;
	else
		return lhs.start < rhs.start;
}

class JunctionAccumulator {

private:

	vector<string> m_alignment_files;

	string m_junction_file;

	string m_junction_ins_file;

	string m_head_line;

	size_t m_max_read_width;

	string m_chrom_dir;

	CHROM_JUNC_ACCUM_HASH_COMB m_junc_accum_hash;

	map<string, vector<JuncSort> > m_junc_sort;

public:
	JunctionAccumulator(const vector<string>& alignment_files, string junction_file, string junction_ins_file, size_t max_read_width, const string& chrom_dir) : m_alignment_files(alignment_files), 
		m_junction_file(junction_file), m_junction_ins_file(junction_ins_file), m_max_read_width(max_read_width), m_chrom_dir(chrom_dir)
	{
	}

	bool Init(const vector<string>& alignment_files, string junction_file, string junction_ins_file, size_t max_read_width, const string& chrom_dir)
	{
		Clear();

		m_alignment_files = alignment_files;

		m_junction_file = junction_file;

		m_junction_ins_file = junction_ins_file;

		m_max_read_width = max_read_width;

		m_chrom_dir = chrom_dir;
	}

	bool Clear()
	{
		m_alignment_files.clear();

		m_head_line.clear();

		m_junction_file.clear();

		m_junc_accum_hash.clear();

		m_max_read_width = 0;

		m_chrom_dir.clear();

		m_junction_ins_file.clear();

		m_junc_sort.clear();

		return true;
	}

	void UpdateJuncTable(size_t combined_offset, size_t prefixlen, size_t suffixlen, const string& ins_str, const SamRec& samrec)
	{
		CHROM_JUNC_ACCUM_HASH_COMB::iterator chrom_junc_hash_iter = m_junc_accum_hash.find(samrec.chrom_name);

		if (chrom_junc_hash_iter == m_junc_accum_hash.end())
		{
			JUNC_ACCUM_HASH_COMB junc_hash_comb;

			chrom_junc_hash_iter = (m_junc_accum_hash.insert(CHROM_JUNC_ACCUM_HASH_COMB::value_type(samrec.chrom_name, junc_hash_comb))).first;							
		}

		JUNC_ACCUM_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

		JUNC_ACCUM_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

		if (junc_hash_comb_iter != junc_hash_comb.end())
		{
			junc_hash_comb_iter->second.inc_hits(prefixlen, suffixlen, samrec.tagidx, samrec.mis_match, samrec.strand_t, ins_str);
		}
		else
		{
			chrom_junc_hash_iter->second.insert(JUNC_ACCUM_HASH_COMB::value_type(combined_offset, JuncAccum(/*prim, *//*flankseq, */
				prefixlen, suffixlen, m_max_read_width, samrec.tagidx, samrec.mis_match, samrec.strand_t, ins_str)));
		}
	}

	void SamRec2Junc(SamRec& samrec)
	{
		vector<pair<size_t, int> >& spliceway_vec = samrec.spliceway_vec;

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

				UpdateJuncTable(combined_offset, prefixlen, suffixlen, ins_str, samrec);
			}
		}
	}

	void ReadSam()
	{
		vector<string>::iterator alignment_file_iter;

		for (alignment_file_iter = m_alignment_files.begin(); alignment_file_iter != m_alignment_files.end(); ++alignment_file_iter)
		{
			ifstream ifs(alignment_file_iter->c_str());

			if (ifs.is_open())
			{
				while (!ifs.eof() )
				{
					string line;

					getline(ifs,line);

					if (line.empty() || line[0] == '@')
						continue;

					SamRec samrec(line);

					if (samrec.isspliced)
						SamRec2Junc(samrec);
				}

				ifs.close();
			}
			else 
				cerr << "Can't open file: "<<*alignment_file_iter<<endl;
		}
	}

	void LoadFlankString()
	{
		CHROM_JUNC_ACCUM_HASH_COMB::iterator chm_iter;

		for (chm_iter = m_junc_accum_hash.begin(); chm_iter != m_junc_accum_hash.end(); ++chm_iter)
		{
			if (m_junc_sort.find(chm_iter->first) == m_junc_sort.end())
				m_junc_sort[chm_iter->first];

			vector<JuncSort>& junc_sort_vec = m_junc_sort.find(chm_iter->first)->second;

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

			JUNC_ACCUM_HASH_COMB::iterator iter_conj;

			for (iter_conj = chm_iter->second.begin(); iter_conj != chm_iter->second.end(); ++iter_conj)
			{
				size_t comb_offset = iter_conj->first;

				size_t prefix_end = comb_offset >> THIRTY_TWO;

				size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

				iter_conj->second.set_coverage();

				iter_conj->second.set_entropy();

				string flankstr = chromseq.substr(prefix_end, 2) + chromseq.substr(suffix_st - 3, 2);

				for (size_t i = 0; i < flankstr.length(); ++i)
				{
					if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
						flankstr[i] = flankstr[i] + 'A' - 'a';
				}

				iter_conj->second.set_flankstring(flankstr);

				//ofs <<chm_iter->first <<'\t'<< prefix_end << '\t' << suffix_st<< '\t'<<juncidstr<<juncid<<'\t';

				//ofs << hits << '\t'<<strand <<'\t'<<prefix_end << '\t' << suffix_st<< "\t255,0,0\t2\t"<<iter_conj->second.max_prefix_len 
				//	<< ','<< iter_conj->second.max_suffix_len << ",\t0,"<<suffix_st + iter_conj->second.max_suffix_len - prefix_end + 1<<",\t";

				size_t intron_len = suffix_st - prefix_end - 1;

				iter_conj->second.set_pq_score(intron_len, chrom_size);

				iter_conj->second.set_il_score(prefix_end, suffix_st);

				iter_conj->second.set_ave_mis();

				junc_sort_vec.push_back(JuncSort(prefix_end, suffix_st, &iter_conj->second));

				//if (suffix_st == prefix_end)
				//{
				//}
				//else
				//{
				//}

				//ofs << rank << '\t'<< flankcase <<'\t'<<flankstr<< '\t'<< lpqave<<'\t'<<iter_conj->second.min_mismatch<<'\t' <<iter_conj->second.max_mismatch<<'\t'
				//	<< (double)iter_conj->second.sum_mismatch / (double)hits << '\t';

				//map<string, int>::const_iterator m_iter;

				//for (m_iter = iter_conj->second.ins.begin(); m_iter != iter_conj->second.ins.end(); ++m_iter)
				//{
				//	ofs << m_iter->first <<'-'<<m_iter->second<<';';
				//}

				//ofs<<endl;

				//++juncid;
			}
		}
	}

	void SortJunc()
	{
		map<string, vector<JuncSort> >::iterator chrom_iter;

		for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
		{
			sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
			//vector<JuncSort>::iterator junc_sort_iter;

			//vector<JuncSort>& junc_sort_vec = chrom_iter->second;

			//for (junc_sort_iter = junc_sort_vec.begin(); junc_sort_iter != junc_sort_vec.end(); ++junc_sort_iter)
			//{
			//	sort(
			//}
		}

	}

	void WriteJunction()
	{
		ofstream ofs_normal(m_junction_file.c_str());

		ofstream ofs_insert(m_junction_ins_file.c_str());

		map<string, vector<JuncSort> >::iterator chrom_iter;

		size_t normal_count = 0, ins_count = 0;

		for (chrom_iter = m_junc_sort.begin(); chrom_iter != m_junc_sort.end(); ++chrom_iter)
		{
			//sort(chrom_iter->second.begin(), chrom_iter->second.end(), comp_junc_sort);
			vector<JuncSort>::iterator junc_sort_iter;

			string chrom = chrom_iter->first;

			vector<JuncSort>& junc_sort_vec = chrom_iter->second;

			for (junc_sort_iter = junc_sort_vec.begin(); junc_sort_iter != junc_sort_vec.end(); ++junc_sort_iter)
			{
				if (junc_sort_iter->start == junc_sort_iter->end)
				{
					++ins_count;

					ofs_insert << chrom<< '\t' << junc_sort_iter->to_insert_junction(ins_count)<<endl;

					//ofs_insert<<chrom<<'\t'<<junc_sort_iter->start<<'\t'<<junc_sort_iter->end<<'\t'<<"JUNC_"<<ins_count<<'\t'<<junc_sort_iter->junc_accum_ptr->coverage
					//	<<'\t'<<junc_sort_iter->junc_accum_ptr->strand<<'\t'<<junc_sort_iter->start<<'\t'<<junc_sort_iter->end<<"\t255,0,0\t"
					//	<<2<<'\t'<<junc_sort_iter->junc_accum_ptr->max_prefix_len<<','<<junc_sort_iter->junc_accum_ptr->max_suffix_len
					//	<< ",\t0,"<<junc_sort_iter->end + junc_sort_iter->junc_accum_ptr->max_suffix_len - junc_sort_iter->start + 1<<",\t"
					//	<< junc_sort_iter->junc_accum_ptr->entropy << '\t'<< junc_sort_iter->junc_accum_ptr->flankcase <<'\t'
					//	<<junc_sort_iter->junc_accum_ptr->flankstr <<'\t'<<junc_sort_iter->junc_accum_ptr->il_score<<'\t'<<junc_sort_iter->junc_accum_ptr->lpqave
					//	<<'\t' << junc_sort_iter->junc_accum_ptr->min_mismatch<<'\t'<<junc_sort_iter->junc_accum_ptr->max_mismatch<<'\t'
					//	<<(double)junc_sort_iter->junc_accum_ptr->sum_mismatch/(double)junc_sort_iter->junc_accum_ptr->coverage<<endl;


					//ofs_insert<< chrom << '\t' << junc_sort_iter->start << 
				}
				else
				{
					++normal_count;

					ofs_normal << chrom<< '\t' << junc_sort_iter->to_normal_junction(normal_count)<<endl;
				}
			}
		}
	}
};

#endif