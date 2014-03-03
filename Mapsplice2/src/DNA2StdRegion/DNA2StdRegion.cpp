#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <iterator>
#include <map>

using namespace std;

//#define VS

#ifdef VS
#include <hash_map> //vc only
#include <hash_set>
//#define _CRT_SECURE_NO_WARNINGS 
//#pragma warning(disable:_CRT_SECURE_NO_WARNINGS)
 
#else
#include <ext/hash_map> //g++ only
#include <ext/hash_set>

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_set;

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

#endif

#ifdef VS
using namespace stdext;
#endif

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

struct fusion_db_rec {
	size_t start;
	size_t end;
	size_t strand1;
	size_t strand2;
	string chr1;
	string chr2;
	string line;
	string type;

	string gene1;
	string gene2;

	fusion_db_rec(size_t st, size_t ed, size_t std1, size_t std2, const string& c1, const string& c2, const string& tp, const string& g1, const string& g2, string& ln) 
		: start(st), end(ed), strand1(std1), strand2(std2), chr1(c1), chr2(c2), type(tp), gene1(g1), gene2(g2), line(ln) {}

	string tostring(size_t ext) 
	{
		char buf[10000];

		char std1, std2;

		if (strand1 == 0)
			std1 = '+';
		else
			std1 = '-';

		if (strand2 == 0)
			std2 = '-';
		else
			std2 = '+';

		size_t st1, end1, st2, end2;

		if (start > ext)
		{
			st1 = start - ext;
		}
		else
		{
			st1 = 1;
		}

		end1 = start + ext;

		if (end > ext)
		{
			st2 = end - ext;
		}
		else
		{
			st2 = 1;
		}

		end2 = end + ext;

		sprintf(buf, "%s\t%c\t%llu\t%llu\t%s\t%c\t%llu\t%llu\n%s\t%c\t%llu\t%llu\t%s\t%c\t%llu\t%llu", chr1.c_str(), std1, st1, end1, chr2.c_str(), std2, st2, end2,
			chr2.c_str(), std2, st2, end2, chr1.c_str(), std1, st1, end1);

		return buf;
	}
};

inline char
complement(int i) {
	static const int b2c_size = 20;
	static const char b2c[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
	};
	static const char b2cl[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
	};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

inline string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}


size_t boundary = 1000;

bool comp_mate_dist(int lhs, int rhs)
{
	if (abs(abs(lhs) - abs(rhs)) < 1000)
		return true;

	return false;
}

string
basename2(string filename) {

	//cout << "bef: "<<filename<<endl;
	const string s(filename.substr(0, filename.find_last_of(".")));
	size_t final_slash = s.find_last_of("/");

	if (final_slash == string::npos)
		final_slash = s.find_last_of("\\");
	if (final_slash != string::npos)
	{
		//cout << "aft 1: "<<s.substr(final_slash + 1)<<endl;
		return s.substr(final_slash + 1);
	}
	else
	{
		//cout << "aft 2: "<<s<<endl;
		return s;
	}
}


void 
read_db(const char* fusion_junction_file, const char* std_region, map<string, vector<fusion_db_rec> > & fusion_junc_map, map<string, string> fusion_db_set, size_t ext)
{
	ifstream ifs(fusion_junction_file);

	ofstream ofs(std_region);

	string std_region_skip_type = std_region; std_region_skip_type.append(".skiptype");

	ofstream ofs_skip_type(std_region_skip_type.c_str());

	if (ifs.is_open())
	{
		string line;

		getline(ifs,line);

		while (getline(ifs,line))
		{
			if (line.empty())
				continue;

			char gene1[1000], gene2[1000], skip[1000], chr1[1000], chr2[1000], type[1000];

			size_t start, end;

			size_t strand1, strand2;

			for (size_t i = 0; i < line.size(); ++i)
			{
				if (line[i] == ' ')
					line[i] = '_';
			}

			//individual	num	chr1	str1	pos1	chr2	str2	pos2	class
			sscanf(line.c_str(), "%s\t%s\t%s\t%llu\t%llu\t%s\t%llu\t%llu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
				skip, skip, chr1, &strand1, &start, chr2, &strand2, &end, type, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, gene1, skip, gene2);

			string chrstr1(chr1), chrstr2(chr2);

			if (chrstr1 == "23")
			{
				chrstr1 = "X";
			}

			if (chrstr2 == "23")
			{
				chrstr2 = "X";
			}

			if (chrstr1 == "24")
			{
				chrstr1 = "Y";
			}

			if (chrstr2 == "24")
			{
				chrstr2 = "Y";
			}

			if (chrstr1 == "25")
			{
				chrstr1 = "M";
			}

			if (chrstr2 == "25")
			{
				chrstr2 = "M";
			}

			string comb;

			if (chrstr1 > chrstr2 || (chrstr1 == chrstr2 && start > end ) )
			{
				comb =  "chr" + chrstr2 + "_" + "chr" + chrstr1;
				swap(start, end);
				swap(strand1, strand2);
				swap(chrstr1, chrstr2);
			}
			else
			{
				comb =  "chr" + chrstr1 + "_" + "chr" + chrstr2;
			}

			chrstr1 = "chr" + chrstr1;

			chrstr2 = "chr" + chrstr2;

			//comb = chrstr1 + "_" + chrstr2;

			char fusion_key[10000];

			sprintf(fusion_key, "%s_%llu_%llu", comb.c_str(), start, end);

			//fusion_junc_map[comb].push_back(fusion_db_rec(start, end, 0, 0, chrstr1, chrstr2, type, line));

			//fusion_junc_map[comb].push_back(fusion_db_rec(start, end, 0, 1, chrstr1, chrstr2, type, line));

			//fusion_junc_map[comb].push_back(fusion_db_rec(start, end, 1, 0, chrstr1, chrstr2, type, line));

			//fusion_junc_map[comb].push_back(fusion_db_rec(start, end, 1, 1, chrstr1, chrstr2, type, line));

			fusion_junc_map[comb].push_back(fusion_db_rec(start, end, strand1, strand2, chrstr1, chrstr2, type, gene1, gene2, line));

			if (fusion_db_set.find(fusion_key) != fusion_db_set.end())
			{
				cout << fusion_db_set.find(fusion_key)->second << endl;

				cout << line << endl;
			}
			else
				fusion_db_set.insert(make_pair(fusion_key, line));
		}

		ifs.close();

	}
	else
	{
		cout << "can't open file "<< fusion_junction_file<<endl; exit(1);
	}

	map<string, vector<fusion_db_rec> >::iterator chr_iter;

	vector<fusion_db_rec>::iterator loc_iter;

	for (chr_iter = fusion_junc_map.begin(); chr_iter != fusion_junc_map.end(); ++chr_iter)
	{
		for (loc_iter = chr_iter->second.begin(); loc_iter != chr_iter->second.end(); ++loc_iter)
		{
			if (loc_iter->type == "long_range" || loc_iter->type == "inter_chr" || loc_iter->gene1 != loc_iter->gene2)
			{
				ofs << loc_iter->tostring(ext) << endl;
			}
			else
			{
				ofs_skip_type << loc_iter->tostring(ext) << '\t' << loc_iter->type << endl;
			}
		}
	}
}


void splitcomma(string coma_str, vector<string>& splitted)
{
	size_t pre_index = 0, index = 0;

	index = coma_str.find(",", pre_index);

	while (index != string::npos)
	{
		string tobe_fixed_hmer_str = coma_str.substr(pre_index, index - pre_index);

		splitted.push_back(tobe_fixed_hmer_str);

		pre_index = index + 1;

		index = coma_str.find(",", pre_index);
	}

}

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "DNA_fusion std_region ext" <<endl;
		exit(0);
	}

	const char* fusion_db = argv[1];

	const char* std_region_file = argv[2];

	size_t ext = atoi(argv[3]);

	cout << fusion_db <<endl;

	map<string, vector<fusion_db_rec> >  fusion_db_map;

	map<string, string> fusion_db_set;

	read_db(fusion_db, std_region_file, fusion_db_map, fusion_db_set, ext);
}