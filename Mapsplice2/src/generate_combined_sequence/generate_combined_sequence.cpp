#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string.h>
#include <iterator>
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

struct PairedRegions
{
	string chrom1;

	char strand1;

	size_t start1, end1;

	hash_set<size_t> mapped_reads_index;

	string chrom2;

	char strand2;

	size_t start2, end2;

	//vector<size_t> mapped_reads_index2;	

	PairedRegions(const string& chr1, char sd1, size_t st1, size_t en1, const string& chr2,  char sd2, size_t st2, size_t en2) : chrom1(chr1), strand1(sd1), start1(st1), end1(en1), chrom2(chr2), strand2(sd2), start2(st2), end2(en2)
	{
	}

	string tostring()
	{
		char output_chr[5000];

		sprintf(output_chr, "%s\t%c\t%llu\t%llu\t%s\t%c\t%llu\t%llu", chrom1.c_str(), strand1, start1, end1, chrom2.c_str(), strand2, start2, end2);

		return output_chr;		
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

typedef vector<PairedRegions> PAIRED_REGIONS_VEC;
typedef vector<PairedRegions>::iterator PAIRED_REGIONS_ITERATOR;
typedef PairedRegions* PAIRED_REGIONS_PTR;
typedef vector<PAIRED_REGIONS_PTR> PAIRED_REGIONS_PTR_VEC;
typedef vector<PAIRED_REGIONS_PTR>::iterator PAIRED_REGIONS_PTR_VEC_ITERATOR;
typedef hash_map<string, PAIRED_REGIONS_PTR_VEC > HASH_CHROM_PAIRED_REGIONS_PTR_VEC;
typedef hash_map<string, PAIRED_REGIONS_PTR_VEC >::iterator HASH_CHROM_PAIRED_REGIONS_PTR_VEC_ITERATOR;

size_t max_dist = 0;

int ReadRegions(const char* region_file, PAIRED_REGIONS_VEC& paired_regions, size_t ext)
{
	ifstream ifs(region_file);

	int count = 0;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line.empty())
				continue;

			char chrom1[1000];

			char strand1;
			
			size_t start1, end1;

			char chrom2[1000];

			char strand2;

			size_t start2, end2;

			//chr10	-	46959544	46969544	chr10	-	88971413	88981413
			sscanf(line.c_str(), "%s\t%c\t%llu\t%llu\t%s\t%c\t%llu\t%llu", chrom1, &strand1, &start1, &end1, chrom2, &strand2, &start2, &end2);

			string chrom1str = chrom1;

			string chrom2str = chrom2;

			if (start1 > ext)
				start1 -= ext;
			else
				start1 = 1;

			end1 += ext;

			if (start2 > ext)
				start2 -= ext;
			else
				start2 = 1;

			end2 += ext;

			paired_regions.push_back(PairedRegions(chrom1, strand1, start1, end1, chrom2, strand2, start2, end2));

			if (max_dist <  (end1 - start1))
				max_dist = end1 - start1;

			if (max_dist <  (end2 - start2))
				max_dist = end2 - start2;

			++count;
		}
	}
	else
		cout <<"can't open file: " << region_file << endl;

	return count;
}

bool comp1(const PairedRegions* lhs, const PairedRegions* rhs)
{
	if (lhs->start1 == rhs->start1)
		return lhs->end1 < rhs->end1 ;

	return lhs->start1 < rhs->start1;
}

bool comp2(const PairedRegions* lhs, const PairedRegions* rhs)
{
	if (lhs->start2 == rhs->start2)
		return lhs->end2 < rhs->end2;

	return lhs->start2 < rhs->start2;
}

void SortRegion(PAIRED_REGIONS_VEC& paired_regions, HASH_CHROM_PAIRED_REGIONS_PTR_VEC& paired_regions_sorted1, HASH_CHROM_PAIRED_REGIONS_PTR_VEC& paired_regions_sorted2)
{
	PAIRED_REGIONS_ITERATOR paired_regions_iter;

	for (paired_regions_iter = paired_regions.begin(); paired_regions_iter != paired_regions.end(); ++paired_regions_iter)
	{
		paired_regions_sorted1[paired_regions_iter->chrom1].push_back(&(*paired_regions_iter));

		paired_regions_sorted2[paired_regions_iter->chrom2].push_back(&(*paired_regions_iter));
	}

	HASH_CHROM_PAIRED_REGIONS_PTR_VEC_ITERATOR chrom_paired_regions_iter;
	for (chrom_paired_regions_iter = paired_regions_sorted1.begin(); chrom_paired_regions_iter != paired_regions_sorted1.end(); ++chrom_paired_regions_iter)
		sort(chrom_paired_regions_iter->second.begin(), chrom_paired_regions_iter->second.end(), comp1);

	for (chrom_paired_regions_iter = paired_regions_sorted2.begin(); chrom_paired_regions_iter != paired_regions_sorted2.end(); ++chrom_paired_regions_iter)
		sort(chrom_paired_regions_iter->second.begin(), chrom_paired_regions_iter->second.end(), comp2);
}

void
readchrom(const char* filename, string& longseq)
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

			if (line.length() == 0)
				continue;
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";
}

int CombineSequence(PAIRED_REGIONS_VEC& paired_regions, string out_path, string chrom_dir)
{
	int count = 0;

	string* chromseq1;
	string* chromseq2;
	string prevchrom1 = "";
	string prevchrom2 = "";

	PAIRED_REGIONS_VEC::iterator paired_regions_iter;

	int region_count = 0;

	string sub_dir = out_path + "/splitchromosome";

	string mkdir_cmd1 = "mkdir ";

	mkdir_cmd1.append(out_path);

	system(mkdir_cmd1.c_str());

	string mkdir_cmd2 = "mkdir ";

	mkdir_cmd2.append(sub_dir);

	system(mkdir_cmd2.c_str());

	string sub_dir2 = out_path + "/splitindex";

	string mkdir_cmd3 = "mkdir ";

	mkdir_cmd3.append(sub_dir2);

	system(mkdir_cmd3.c_str());

	ofstream combined_ofs;

	string combined_sequence_str = out_path + "combined_sequence.fa";

	combined_ofs.open(combined_sequence_str.c_str());

	map<string, string> chrom_map;

	for (paired_regions_iter = paired_regions.begin(); paired_regions_iter != paired_regions.end(); ++paired_regions_iter)
	{
		++region_count;

		cout << paired_regions_iter->tostring() << endl;

		string chrom1 = paired_regions_iter->chrom1;

		string chrom2 = paired_regions_iter->chrom2;

		size_t start1 = paired_regions_iter->start1;

		size_t start2 = paired_regions_iter->start2;

		size_t end1 = paired_regions_iter->end1;

		size_t end2 = paired_regions_iter->end2;

		char strand1 = paired_regions_iter->strand1;

		char strand2 = paired_regions_iter->strand2;

		//cout << "read chromosome" << endl;

		if (chrom_map.find(chrom1) == chrom_map.end())
		{
			string chrom_file = chrom_dir;

			chrom_file.append(chrom1);

			chrom_file.append(".fa");

			chromseq1 = &(chrom_map[chrom1]);

			readchrom(chrom_file.c_str(), *chromseq1);

			prevchrom1 = chrom1;
		}
		else
			chromseq1 = &(chrom_map[chrom1]);


		if (chrom_map.find(chrom2) == chrom_map.end())
		{
			string chrom_file = chrom_dir;

			chrom_file.append(chrom2);

			chrom_file.append(".fa");

			chromseq2 = &(chrom_map[chrom2]);

			readchrom(chrom_file.c_str(), *chromseq2);

			prevchrom2 = chrom2;
		}
		else
			chromseq2 = &(chrom_map[chrom2]);

		if (end1 > (*chromseq1).length())
			end1 = (*chromseq1).length();

		if (end2 > (*chromseq2).length())
			end2 = (*chromseq2).length(); 

		char path_chr[1000];

		cout << "generated chromosome name" << endl;

		sprintf(path_chr, "%s=%c=%llu=%llu=%s=%c=%llu=%llu", chrom1.c_str(), strand1, start1, end1, chrom2.c_str(), strand2, start2, end2);

		string cur_path = sub_dir;

		string cur_chrom_file = cur_path + "/";

		cur_chrom_file.append(path_chr); cur_chrom_file.append(".fa");

		ofstream cur_chrom_file_ofs(cur_chrom_file.c_str());

		cur_chrom_file_ofs << '>' <<path_chr<<endl;

		combined_ofs <<'>' <<path_chr<<endl;

		string comb_chrom_seq, chrom1_seq, chrom2_seq;
		
		//cout << "chromosome size"<<endl;
		//cout << "chromseq1.size():" <<chromseq1.size() <<endl;
		//cout << "start1:" <<start1 << endl;
		//cout << "end1 - start1:" <<end1 - start1 << endl;
		//cout << "chromseq2.size():" <<chromseq2.size() <<endl;
		//cout << "start2:" <<start2 << endl;
		//cout << "end2- start2:" <<end2- start2 << endl;

		if (strand1 == '+')
			chrom1_seq = (*chromseq1).substr(start1, end1 - start1);
		else
			chrom1_seq = revcomp((*chromseq1).substr(start1, end1 - start1));
		
		if (strand2 == '+')
			chrom2_seq = (*chromseq2).substr(start2, end2 - start2);
		else
			chrom2_seq = revcomp((*chromseq2).substr(start2, end2 - start2));

		comb_chrom_seq = chrom1_seq + chrom2_seq;

		cout << "write file"<<endl;
		for (size_t i = 0; i < comb_chrom_seq.length(); i = i + 50)
		{
			cur_chrom_file_ofs << comb_chrom_seq.substr(i, 50)<<endl;
			combined_ofs << comb_chrom_seq.substr(i, 50)<<endl;
		}
	}

	return count;
}

int main(int argc, char* argv[])
{
	if (argc < 5)
	{
		cout << "region_file output_path chromosome_dir extension"<<endl;

		exit(1);
	}

	const char* region_file = argv[1];

	string out_path = argv[2];

	string chrom_dir = argv[3];

	size_t ext = atoi(argv[4]);

	vector<PairedRegions> paired_regions;

	cout << "read regions" << endl;

	ReadRegions(region_file, paired_regions, ext);

	cout << max_dist<<endl;

	cout << "combine sequence"<<endl;

	CombineSequence(paired_regions, out_path, chrom_dir);

}