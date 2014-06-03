#ifndef __SAM_2_JUNC_H_INCLUDED__
#define __SAM_2_JUNC_H_INCLUDED__

/*
 *    fusionsam2junc_filteranchor_newfmt.cpp
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Kai Wang
 *
 *    Authors: Kai Wang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


//#define VS

#include <iostream>
#include <vector>

#include <string>

#include <ext/hash_map> //g++ only

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
#include <dirent.h>
#include <iomanip>
#include <map>
#include <queue>
#include <list>
#include <ext/hash_map> //g++ only

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;

namespace __gnu_cxx
{
	template<typename Traits, typename Allocator>
	struct hash<std::basic_string<char, Traits, Allocator> >
	{
		size_t operator()(const std::basic_string<char, Traits, Allocator>& __s) const
		{
			return __stl_hash_string(__s.c_str());
		}
	};
}

using namespace std;


#define IS_PAIRED 0x0001
#define IS_PAIRED_MAPPED 0x0002
#define IS_UNMAPPED 0x0004
#define MATE_UNMAPPED 0x0008
#define IS_REVERSE 0x0010
#define IS_MATE_REVERSE 0x0020
#define IS_FIRST_END 0x040
#define IS_SECOND_END 0x0080
#define IS_PRIMARY 0x0100
#define IS_FAILED_QUAL_CHECK 0x0200
#define IS_PCR_DUP 0x0400

const size_t THIRTY_TWO = 32;
const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;

struct JuncInfo{
	JuncInfo(/*int pm, *//*const string fs, */size_t loc, size_t suffix_len, size_t rw, size_t tagidx, unsigned short mis, size_t strand, string insert = "") : /*prim(pm), *//*flankstr(fs),*/ p(rw - 1, 0), positive_count(0), negative_count(0)
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

	char get_strand() const
	{
		return positive_count > negative_count ? '+' : '-';
	}

	//int prim;
	//string flankstr;
	size_t max_prefix_len;
	size_t max_suffix_len;
	vector<unsigned short> p;

	unsigned short positive_count, negative_count;

	unsigned short max_mismatch;

	unsigned short min_mismatch;

	unsigned short sum_mismatch;

	map<size_t, int> m;

	map<string, int> ins;
};

typedef hash_map<size_t, JuncInfo> JUNC_SUFFIX;
typedef JUNC_SUFFIX::iterator JUNC_SUFFIX_ITER;
typedef JUNC_SUFFIX::const_iterator JUNC_SUFFIX_CITER;

typedef hash_map<size_t, JUNC_SUFFIX > JUNC_HASH;
typedef JUNC_HASH::iterator JUNC_HASH_ITER;
typedef JUNC_HASH::const_iterator JUNC_HASH_CITER;

typedef map<string, JUNC_HASH > CONJ_HASH_MAP;

typedef hash_map<size_t, JuncInfo> JUNC_HASH_COMB;

typedef hash_map<string, JUNC_HASH_COMB> CHROM_JUNC_HASH_COMB;

struct JuncForSort{
	int juncst;
	int juncend;
	int hits;
	int kinds;
	string blocks;
	string blocksoffset;
	string rank;
	string lpq;
	unsigned short min_mismatch;
	unsigned short max_mismatch;
	double ave_mismatch;

	string flankstr;

	size_t flankcase;

	char strand;

	string ins_str;

	JuncForSort(const int& jst, const int& jend, const int& hts, char sd, const int& kds, const string& blks, const string& blksoft,
		const string& rk, /*const string& fs, const string& fc, */const string& l, unsigned short min_mis, unsigned short max_mis, double ave_mis) :
	juncst(jst), juncend(jend), hits(hts), strand(sd), kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), /*flankstr(fs), flankcase(fc),*/
		lpq(l), min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis) {}

	JuncForSort(const int& jst, const int& jend, const int& hts, char sd, const int& kds, const string& blks, const string& blksoft,
		const string& rk, const string& fs, size_t fc, const string& l, unsigned short min_mis, unsigned short max_mis, double ave_mis, const string& ins) :
	juncst(jst), juncend(jend), hits(hts), strand(sd), kinds(kds), blocks(blks), blocksoffset(blksoft), rank(rk), flankstr(fs), flankcase(fc),
		lpq(l), min_mismatch(min_mis), max_mismatch(max_mis), ave_mismatch(ave_mis), ins_str(ins) {}

};

bool compjunc(const JuncForSort& lhs, const JuncForSort& rhs)
{
	if (lhs.juncst == rhs.juncst)
		return lhs.juncend < rhs.juncend;
	else
		return lhs.juncst < rhs.juncst;
}

void readchrom(const char* filename, string& longseq)
{
	size_t size;

	ifstream longfile(filename);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	if (longfile.is_open())
	{
		string skipline;
		getline(longfile,skipline);

		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);

			if (line.empty())
				continue;
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else
	{
		cout << "Unable to open file";

		cerr << "can not open chromosome file:" <<filename << endl;
		exit(1);
	}
}

void SortJuncComb(const char* juncfile, string chrom_dir, size_t min_intron, size_t max_intron)
{
	map<string, vector<JuncForSort> > m_p;

    ifstream ifs(juncfile);

	if (ifs.is_open())
	{
		string skipline;
		getline(ifs,skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			char chromname[100], juncname[100], strand, rgb[100], blocks[1000], blocksoffset[1000], rank[100], lpq[100], flankchr[100], ins_chr[100000];
			int juncst, juncend, prefixend, suffixst, kinds, hits, flankcase;

			if (line.length() > 100000)
			{
				cout << line << endl;

				continue;
			}

			unsigned short min_mis, max_mis;

			double ave_mis;

			//chr21	33637598	33637711	JUNC_1	1	+	33637598	33637711	255,0,0	2	48,52,	0	5	GTAG	1	0	0	0
			int read_count = sscanf(line.c_str(), "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%hu\t%hu\t%lf\t%s", chromname,  &prefixend, &suffixst,
				juncname, &hits, &strand, &juncst, &juncend, rgb, &kinds, blocks, blocksoffset, rank, &flankcase, flankchr, lpq, &min_mis, &max_mis, &ave_mis, ins_chr);



			string chromstr = chromname;

			string ins_str;
			if (read_count == 20)
				ins_str = ins_chr;

			m_p[chromstr].push_back(JuncForSort(juncst, juncend, hits, strand, kinds, blocks, blocksoffset,	rank, flankchr, flankcase, lpq, min_mis, max_mis, ave_mis, ins_str));
		}

		ifs.close();
	}
	map<string, vector<JuncForSort> >::iterator mpit;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		sort(mpit->second.begin(), mpit->second.end(), compjunc);
	}

	ofstream ofs(juncfile);

	string insertfile = juncfile;insertfile.append("ins");

	ofstream ins_ofs(insertfile.c_str());

	string juncidstr = "JUNC_";
	size_t juncid = 1;

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	for (mpit = m_p.begin(); mpit != m_p.end(); ++mpit)
	{
		vector<JuncForSort>::iterator vit;
		for (vit = mpit->second.begin(); vit != mpit->second.end(); ++vit)
		{

			double il = 1.0 - (((double) (vit->juncend - vit->juncst - min_intron + 1)) / (double (max_intron - min_intron + 2)));

			if (vit->juncst != vit->juncend)
				ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<vit->strand<<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< vit->flankcase <<'\t'<<vit->flankstr <<'\t'<<il<<'\t'<<vit->lpq
				<<'\t' << vit->min_mismatch<<'\t'<<vit->max_mismatch<<'\t'<<vit->ave_mismatch<<endl;
			else
				ins_ofs<<mpit->first<<'\t'<<vit->juncst<<'\t'<<vit->juncend<<'\t'<<juncidstr<<juncid<<'\t'<<vit->hits<<'\t'<<vit->strand<<'\t'<<vit->juncst<< '\t'<<vit->juncend<<"\t255,0,0\t"
				<<vit->kinds<<'\t'<<vit->blocks<<'\t'<< vit->blocksoffset<<'\t'<< vit->rank << '\t'<< vit->flankcase <<'\t'<<vit->flankstr <<'\t'<<il<<'\t'<<vit->lpq
				<<'\t' << vit->min_mismatch<<'\t'<<vit->max_mismatch<<'\t'<<vit->ave_mismatch<<'\t' <<vit->ins_str<< endl;

			++juncid;
		}
	}

	cout << juncid-1 <<" junctions"<<endl;
}

void
WriteJuncComb(const CHROM_JUNC_HASH_COMB& conj_hash_map, ofstream& ofs, size_t m_read_width, double m_max_rank, string chrom_dir, size_t max_intron)
{
	string juncidstr = "JUNC_";
	size_t juncid = 1;

	CHROM_JUNC_HASH_COMB::const_iterator chm_iter;

	for (chm_iter = conj_hash_map.begin(); chm_iter != conj_hash_map.end(); ++chm_iter)
	{
		string chromfile = chrom_dir + chm_iter->first;

		chromfile.append(".fa");

		string chromseq;

		//cout << chromfile << endl;

		readchrom(chromfile.c_str(), chromseq);

		if (chromseq.empty())
		{
			cout <<"empty chrom: "<<chromfile<<endl;

			cerr << "can not open chromosome file:" << chromfile<<endl;

			exit(1);
		}

		size_t chrom_size = chromseq.size() - 1;

		JUNC_HASH_COMB::const_iterator iter_conj;

		for (iter_conj = chm_iter->second.begin(); iter_conj != chm_iter->second.end(); ++iter_conj)
		{
			size_t comb_offset = iter_conj->first;

			size_t prefix_end = comb_offset >> THIRTY_TWO;

			size_t suffix_st = comb_offset & LOWER_THIRTY_TWO_MASK;

			int kinds = 0;
			int hits = 0;
			for (size_t i = 0; i < iter_conj->second.p.size(); ++i)
			{
				if (iter_conj->second.p[i] > 0)
				{
					++kinds;
					hits += iter_conj->second.p[i];
				}
			}

			double rank = 0;

			for (size_t i = 0; i < iter_conj->second.p.size(); ++i)
			{
				if (iter_conj->second.p[i] > 0)
				{
					double pi = (double)iter_conj->second.p[i] / (double)hits;
					rank += pi * log(pi);
				}
			}

			if (rank!=0)
				rank = -rank;

			if (rank < m_max_rank)
			{
				//iter_conj_pre++;
				continue;
			}

			if (chromseq.length() < prefix_end + 2 || chromseq.length() < suffix_st)
			{
				cerr << "chromosome length less than junction start or end position" << endl;

				cerr << "chromosome:"<< chromfile << endl;

				cerr << "chromsoome length: "<< chromseq.length();

				cerr << "prefix_end:" << prefix_end<< endl;

				cerr << "suffix_st:" << suffix_st<< endl;

				exit(1);
			}

			string flankstr = chromseq.substr(prefix_end, 2) + chromseq.substr(suffix_st - 3, 2);

			for (size_t i = 0; i < flankstr.length(); ++i)
			{
				if (flankstr[i] >= 'a' && flankstr[i] <= 'z' )
					flankstr[i] = flankstr[i] + 'A' - 'a';
			}

			int flankcase = 0;

			char strand = '+';

			if (flankstr == "ATAC")
				flankcase = 1;
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
				 flankcase = 4;
			else if (flankstr == "GTAG")
				 flankcase = 5;
			else if (flankstr == "GTAT")
			{
				 flankcase = 2;
				 strand = '-';
			}

			ofs <<chm_iter->first <<'\t'<< prefix_end << '\t' << suffix_st<< '\t'<<juncidstr<<juncid<<'\t';

			ofs << hits << '\t'<<strand <<'\t'<<prefix_end << '\t' << suffix_st<< "\t255,0,0\t2\t"<<iter_conj->second.max_prefix_len
				<< ','<< iter_conj->second.max_suffix_len << ",\t0,"<<suffix_st + iter_conj->second.max_suffix_len - prefix_end + 1<<",\t";

			size_t intron_len = suffix_st - prefix_end - 1;

			double ppower = pow(0.25, double(iter_conj->second.max_prefix_len));

			double pNpower = pow(1.0 - ppower, (double)chrom_size);

			double qpower = pow(0.25, double(iter_conj->second.max_suffix_len));

			double pDpower = pow(1.0 - qpower, (double)intron_len);

			double lpq = 1.0 - (pNpower * pDpower);

			double ppower2 = pow(0.25, double(iter_conj->second.max_prefix_len));

			double pNpower2 = pow(1.0 - ppower2, (double)intron_len );

			double qpower2 = pow(0.25, double(iter_conj->second.max_suffix_len));

			double pDpower2 = pow(1.0 - qpower2, (double)chrom_size);

			double lpq2 = 1.0 - (pNpower2 * pDpower2);

			double lpqave = 1.0 - (lpq + lpq2) / 2;

			ofs << rank << '\t'<< flankcase <<'\t'<<flankstr<< '\t'<< lpqave<<'\t'<<iter_conj->second.min_mismatch<<'\t' <<iter_conj->second.max_mismatch<<'\t'
				<< (double)iter_conj->second.sum_mismatch / (double)hits << '\t';

			map<string, int>::const_iterator m_iter;

			for (m_iter = iter_conj->second.ins.begin(); m_iter != iter_conj->second.ins.end(); ++m_iter)
			{
				ofs << m_iter->first <<'-'<<m_iter->second<<';';
			}

			ofs<<endl;

			++juncid;
		}
	}
}

void
Covert2JuncComb(const char* junc_filename, vector<string>& m_mapreads_files, int m_read_width, string chrom_dir, int m_max_rank, int min_intron, int max_intron, int min_anchor)
{
	ofstream ofs(junc_filename);

	string headline = "track name=junctions description=\"Mapsplice junctions\"";

	ofs << headline << endl;

	//sprintf(wig_filename, "%scoverage.wig.txt", m_path.c_str());
	//ofstream ofs_wig(wig_filename);

	//string headline2 = "track type=bedGraph name=\"Mapsplice - read coverage\"";

	//ofs_wig << headline2 << endl;

	CHROM_JUNC_HASH_COMB conj_hash_map;

	hash_map<string, JUNC_SUFFIX> conj_insert_map;

	for (size_t mi = 0; mi < m_mapreads_files.size(); ++mi)
	{
		ifstream ifs(m_mapreads_files[mi].c_str());

		//VEC_COVER_BLOCK v_coverage_block;

		if (ifs.is_open())
		{
			while (!ifs.eof() )
			{
				string line;
				getline(ifs,line);
				if (line == "" || line[0] == '@')
					continue;

				//cout <<line<<endl;

				char chromname[1000], readname[1000], /*flankseq[10], */chromseq[1000], qualseq[1000], spliceway[2000];
				char strand = '+';
				//size_t /*prim, */
				int prefixst, /*prefixlen,suffixend,*/ strand_t, incorrect/*, suffixlen*//*, spliceoutlen*//*, score*/, mate_offest;

				int mate_diff;

				//unsigned short
				int mis_match;

				char mate_match;
				                                       //%dM\t%dN\t%dM
				int firstTabPos = line.find("\t");
				line = line.substr(firstTabPos - 1);
				//cout << "line: " << endl << line << endl;
				sscanf(line.c_str(), "%s\t%d\t%s\t%d\t%d\t%s\t%c\t%d\t%d\t%s\t%s\tNM:i:%d",
					readname, &strand_t, chromname, &prefixst, &incorrect, spliceway/*, &prefixlen, &spliceoutlen, &suffixlen,*/ , &mate_match, &mate_offest, &mate_diff, chromseq, qualseq, &mis_match);
				//cout << "readname = " << readname << endl;
				if (spliceway[0] == '*')
					continue;
				string cigarString = spliceway;

				string tagnamestr = readname;
				size_t last_idx = tagnamestr.find_first_of("~");
				size_t tagidx = -1;

				if (last_idx != string::npos)
				{
					string tagidxstr = tagnamestr.substr(0, last_idx);
					tagidx = atoi(tagidxstr.c_str()) - 1;
				}

				//cout << line<<endl;

				vector<pair<size_t, int> > spliceway_vec;

				string splicewaystr = spliceway;

				size_t index = 0;

				string flag_str = " ";

				bool exceed = false;

				if ((splicewaystr.find("S") != string::npos))
				{
					//clipe tail
					if (splicewaystr[splicewaystr.length() - 1] == 'S')
					{
						size_t m_idx = splicewaystr.find_last_of("M");

						splicewaystr = splicewaystr.substr(0, m_idx + 1);
					}

					size_t s_index = splicewaystr.find("S");

					if (s_index != string::npos)
					{
						splicewaystr = splicewaystr.substr(s_index + 1, splicewaystr.length() - s_index - 1);
					}

					//ori_splice_way = splice_way;
				}


				for (size_t i = 0; i < splicewaystr.length(); ++i)
				{
					if (splicewaystr[i] == 'D')
						splicewaystr[i] = 'N';
				}

				while (true)
				{
					if (index >= splicewaystr.length())
						break;

					int maplen;

					char flag;

					sscanf(splicewaystr.c_str() + index, "%d%c", &maplen, &flag);

					if (flag_str[0] == ' ')
					{
						if (flag == 'I')
							spliceway_vec.push_back(make_pair(prefixst, -maplen));
						else if (flag == 'M')
							spliceway_vec.push_back(make_pair(prefixst, maplen));
						else if (flag == 'N')
						{
							cout<<"start with N?"<<endl;
							spliceway_vec.push_back(make_pair(prefixst + maplen, 0));
						}
					}
					else if (flag_str[0] == 'M')
					{
						if (flag == 'I')
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second, -maplen));
						else if (flag == 'M')
						{
							cout << "continue Ms?"<<endl;
							spliceway_vec.back().second += maplen;
						}
						else if (flag == 'N')
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));
					}
					else if (flag_str[0] == 'N')
					{
						if (flag == 'I')
							spliceway_vec.back().second = -maplen;
						else if (flag == 'M')
							spliceway_vec.back().second = maplen;
						else if (flag == 'N')
						{
							cout << "continue Ns?"<<endl;
							spliceway_vec.back().first += maplen;
						}
					}
					else if (flag_str[0] == 'I')
					{
						if (flag == 'I')
						{
							cout << "continue Is?"<<endl;
							spliceway_vec.back().second += -maplen;
						}
						else if (flag == 'M')
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));
						else if (flag == 'N')
							spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					}

					flag_str[0] = flag;

					index = splicewaystr.find(flag_str, index) + 1;
				}

				string readstr = chromseq;

				size_t mappedlen = 0;

				if (spliceway_vec.size() > 1 && !exceed)
				{
					vector<pair<size_t, int> >::iterator vp_iter;
					for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
					{
						size_t prefixend, suffixst, prefixlen, suffixlen, combined_offset;

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

							if (ins_str != "A")
							{
								int temp = 0;
							}

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
						CHROM_JUNC_HASH_COMB::iterator chrom_junc_hash_iter = conj_hash_map.find(chromname);

						if (chrom_junc_hash_iter == conj_hash_map.end())
						{
							JUNC_HASH_COMB junc_hash_comb;

							chrom_junc_hash_iter = (conj_hash_map.insert(CHROM_JUNC_HASH_COMB::value_type(chromname, junc_hash_comb))).first;
						}

						JUNC_HASH_COMB& junc_hash_comb = chrom_junc_hash_iter->second;

						JUNC_HASH_COMB::iterator junc_hash_comb_iter = junc_hash_comb.find(combined_offset);

						if (junc_hash_comb_iter != junc_hash_comb.end())
						{
							junc_hash_comb_iter->second.inc_hits(prefixlen, suffixlen, tagidx, mis_match, strand_t, ins_str);
						}
						else
						{
							chrom_junc_hash_iter->second.insert(JUNC_HASH_COMB::value_type(combined_offset, JuncInfo(/*prim, *//*flankseq, */prefixlen, suffixlen, m_read_width, tagidx, mis_match, strand_t, ins_str)));
						}
					}
				}
			}
			ifs.close();
		}
		else cout << "Unable to open file";
	}

	cout <<"write junc"<<endl;
	WriteJuncComb(conj_hash_map, ofs, m_read_width, m_max_rank, chrom_dir, max_intron);

	ofs.close();

	cout <<"sort junc"<<endl;

	SortJuncComb(junc_filename, chrom_dir, min_intron, max_intron);
}

#endif
