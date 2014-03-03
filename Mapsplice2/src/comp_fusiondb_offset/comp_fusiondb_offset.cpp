#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <stdlib.h>
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

struct transcript;

struct transcript {
	//chr14	23067146	23073161	GENE.369246	1	+	23067146	23073161	0,0,0	3	92,89,867	0,3425,5148
	string chromosome;
	size_t start, end;
	string trans_name;
	string m_line;

	char strand;

	vector<size_t> block_sizes;

	vector<size_t> offsets;

	transcript(string& line)
	{
		m_line = line;

		char chr[1000], trans[1000], block_chr[10000], offset_chr[10000];

		char skip[1000];

		size_t kind;

		sscanf(line.c_str(), "%s\t%llu\t%llu\t%s\t%s\t%c\t%s\t%s\t%s\t%llu\t%s\t%s", chr, &start, &end, trans, skip, &strand, skip, skip, skip, &kind, block_chr, offset_chr);

		chromosome = chr;

		trans_name = trans;

		string block_str = block_chr;

		string offset_str = offset_chr;

		size_t pre_index = 0, index = 0;

		while(true)
		{
			index = block_str.find(",", pre_index);

			if (index == string::npos)
			{
				//if (pre_index != 

				string block = block_str.substr(pre_index, block_str.length() - pre_index);

				block_sizes.push_back(atoi(block.c_str()));

				break;
			}

			string block = block_str.substr(pre_index, index - pre_index);

			block_sizes.push_back(atoi(block.c_str()));

			pre_index = index + 1;
		}

		pre_index = 0; index = 0;

		while(true)
		{
			index = offset_str.find(",", pre_index);

			if (index == string::npos)
			{
				//if (pre_index != 

				string block = offset_str.substr(pre_index, offset_str.length() - pre_index);

				offsets.push_back(atoi(block.c_str()));

				break;
			}

			string block = offset_str.substr(pre_index, index - pre_index);

			offsets.push_back(atoi(block.c_str()));

			pre_index = index + 1;
		}
	}
};

struct fusion_db_rec {
	size_t start;
	size_t end;
	string line;
	bool matched;
	fusion_db_rec(size_t st, size_t ed, string& ln) : start(st), end(ed), line(ln) { matched = false; }
};

struct SamRec {
	string tag_name;
	unsigned short strand_t;
	string chrom_name;
	string qual_str;
	string alters;
	size_t start;
	size_t end;
	unsigned short confid;
	string splice_way;
	string mapped_seq;
	unsigned short mis_match;

	vector<pair<size_t, int> > spliceway_vec;

	char paired;

	size_t mate_offset;

	int mate_diff;

	bool wrong_format;

	size_t intron_size;

	size_t mappedlen;

	string line;

	SamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, 
		unsigned short mismatch, const string& alt, const string& qualstr, const string& ln) : 
	tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
		mis_match(mismatch), qual_str(qualstr), alters(alt), paired('*'), mate_offset(0), mate_diff(0), 
	intron_size(0), mappedlen(0), line(ln)
	{
		size_t index = 0;

		wrong_format = false;

		string flag_str = " ";

		while (true)
		{
			if (index >= splice_way.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way.c_str() + index, "%d%c", &maplen, &flag);

			if (flag_str[0] == ' ')
			{
				if (flag == 'I')
					spliceway_vec.push_back(make_pair(start, -maplen));
				else if (flag == 'M')
				{
					spliceway_vec.push_back(make_pair(start, maplen));
					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec.push_back(make_pair(start + maplen, 0));
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

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str[0] == 'N')
			{
				if (flag == 'I')
					spliceway_vec.back().second = -maplen;
				else if (flag == 'M')
				{
					spliceway_vec.back().second = maplen;
					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec.back().first += maplen;

					intron_size += maplen;
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
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first, maplen));
					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));

					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;
	}

	string tostring()
	{
		char sam_rec_char[5000];

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s"/*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
			chrom_name.c_str(), start, confid, splice_way.c_str(), paired, mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, alters.c_str()/*, best, filter_score*/);

		return sam_rec_char;
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


bool comp(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.mis_match < rhs.mis_match;
}

struct PairedSamRec {

	PairedSamRec(int dist, unsigned short tm, size_t is, size_t maplen, SamRec* sam_rec1, SamRec* sam_rec2) : mate_dist(dist), intron_size(is), mappedlen(maplen), total_mismatch(tm) 
	{
		paired_sam_rec = make_pair(sam_rec1, sam_rec2);
	}

	int mate_dist;
	size_t intron_size;
	size_t mappedlen;
	unsigned short total_mismatch;
	pair<SamRec*, SamRec*> paired_sam_rec;
};


size_t boundary = 1000;

bool comp_mate_dist(int lhs, int rhs)
{
	if (abs(abs(lhs) - abs(rhs)) < 1000)
		return true;

	return false;
}

bool comp_dist(const PairedSamRec& lhs, const PairedSamRec& rhs)
{
	if (lhs.mappedlen == rhs.mappedlen)
	{
		if (comp_mate_dist(lhs.mate_dist, rhs.mate_dist))
		{
			if (lhs.total_mismatch == rhs.total_mismatch)
				return lhs.intron_size < rhs.intron_size;
			else
				return lhs.total_mismatch < rhs.total_mismatch;
		}
		else
			return abs(lhs.mate_dist) < abs(rhs.mate_dist);
	}
	else
		return lhs.mappedlen > rhs.mappedlen;
}

string generate_trans_align(map<int, vector<SamRec> >& mapped_reads, hash_map<string, transcript>& trans_map, hash_map<string, string>& read_trans_map, char* chromsome_dir)
{
	string generated_align;

	string base_tag_name = mapped_reads[1].front().tag_name;

	base_tag_name = base_tag_name.substr(0, base_tag_name.length() - 2);

	size_t start1 = mapped_reads[1].front().start;

	size_t end1 = mapped_reads[1].front().end;

	size_t start2 = mapped_reads[2].front().start;

	size_t end2 = mapped_reads[2].front().end;
	
	size_t start = start1 < start2 ? start1 : start2;

	size_t end = end1 > end2 ? end1 : end2;

	if (read_trans_map.find(base_tag_name) == read_trans_map.end())
		cout << "trans not found: " << base_tag_name << endl;

	string cur_trans_name = (read_trans_map.find(base_tag_name))->second;

	if (trans_map.find(cur_trans_name) == trans_map.end())
		cout << "trans not found: " << cur_trans_name << endl;

	transcript& cur_trans = (trans_map.find(cur_trans_name))->second;

	size_t i, j, m;

	ofstream ofs("test.region");

	for (m = 0; m < cur_trans.block_sizes.size(); ++m)
	{
		ofs << cur_trans.start +  cur_trans.offsets[m] << '\t' << cur_trans.start +  cur_trans.offsets[m] + cur_trans.block_sizes[m] - 1 << endl;

	}

	ofs.close();

	for (i = 0; i < cur_trans.block_sizes.size(); ++i)
	{
		if (start >= cur_trans.start +  cur_trans.offsets[i] + 1 && start <= cur_trans.start +  cur_trans.offsets[i] + cur_trans.block_sizes[i])
			break;
	}

	if (i == cur_trans.block_sizes.size())
	{
		cout << "not start block found: "   << endl;
		cout << base_tag_name << endl;
		cout << cur_trans_name << endl;
		cout << start << '\t' << end << endl;
		cout << mapped_reads[1].front().line<<endl;
		cout << mapped_reads[2].front().line<<endl;
		cout << cur_trans.m_line << endl;

		for (m = 0; m < cur_trans.block_sizes.size(); ++m)
		{
			cout << cur_trans.start +  cur_trans.offsets[m] << '\t' << cur_trans.start +  cur_trans.offsets[m] + cur_trans.block_sizes[m] - 1 << endl;

		}

		

		return "";
	}

	for (j = 0; j < cur_trans.block_sizes.size(); ++j)
	{
		if (end >= cur_trans.start +  cur_trans.offsets[j] + 1 && end <= cur_trans.start +  cur_trans.offsets[j] + cur_trans.block_sizes[j])
			break;
	}

	if (j == cur_trans.block_sizes.size())
	{
		cout << "not end block found: "   << endl;
		cout << base_tag_name << endl;
		cout << cur_trans_name << endl;
		cout << cur_trans.m_line << endl;
		cout << start << '\t' << end << endl;
		cout << mapped_reads[1].front().line<<endl;
		cout << mapped_reads[2].front().line<<endl;

		for (m = 0; m < cur_trans.block_sizes.size(); ++m)
		{
			cout << cur_trans.start +  cur_trans.offsets[m] << '\t' << cur_trans.start +  cur_trans.offsets[m] + cur_trans.block_sizes[m] - 1 << endl;

		}

		

		return "";
	}

	generated_align.append(base_tag_name);

	generated_align.append("\t");

	if (cur_trans.strand == '+')
		generated_align.append("0\t");
	else
		generated_align.append("16\t");

	generated_align.append(cur_trans.chromosome);

	generated_align.append("\t");

	char start_chr[1000];

	sprintf(start_chr, "%llu", start);

	generated_align.append(start_chr);

	generated_align.append("\t");

	generated_align.append("255\t");
	//seq.13/1        16      chr16   30080866        255     75M     *       0       0       TGGAAGGCACCTTGCTGAAGCCCAACATGGTCACCCCAGGCCATGCTTGCACTCAGAAGTTTTCTCATGAGGAGA     
	//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    NM:i:0  MD:Z:75

	for (size_t k = i; k <= j; ++k)
	{
		char jump_code[1000];

		if (i == j)
		{
			sprintf(jump_code, "%lluM", end - start + 1);

			generated_align.append(jump_code);
		}
		else if (k == i )
		{
			sprintf(jump_code, "%lluM%lluN", cur_trans.start +  cur_trans.offsets[k] + cur_trans.block_sizes[k] - start + 1, 
				cur_trans.start +  cur_trans.offsets[k+1] - (cur_trans.start +  cur_trans.offsets[k] + cur_trans.block_sizes[k] - 1) - 1);

			generated_align.append(jump_code);
		}
		else if (k == j)
		{
			sprintf(jump_code, "%lluM", end - (cur_trans.start +  cur_trans.offsets[k] + 1) + 1);

			generated_align.append(jump_code);
		}
		else
		{
			sprintf(jump_code, "%lluM%lluN", cur_trans.block_sizes[k], 
				cur_trans.start +  cur_trans.offsets[k+1] - (cur_trans.start +  cur_trans.offsets[k] + cur_trans.block_sizes[k]));

			generated_align.append(jump_code);
		}
	}

	generated_align.append("\t*\t0\t0\t");

	generated_align.append(mapped_reads[1].front().mapped_seq);

	generated_align.append(mapped_reads[2].front().mapped_seq);

	generated_align.append("\t");

	//string Is = string(mapped_reads[1].front().mapped_seq.length() + mapped_reads[2].front().mapped_seq.length(), 'I');

	generated_align.append(string(mapped_reads[1].front().mapped_seq.length() + mapped_reads[2].front().mapped_seq.length(), 'I'));

	generated_align.append("\t");

	generated_align.append("NM:i:");

	char mis[1000];

	sprintf(mis, "%llu", mapped_reads[1].front().mis_match + mapped_reads[2].front().mis_match);

	generated_align.append(mis);

	return generated_align;
}

void
GenerateAlignmentsFromTranscript(const char* infile, const char* out_align, hash_map<string, transcript>& trans_map, hash_map<string, string>& read_trans_map, const char* chromsome_dir)
{
	ifstream ifs(infile);

	ofstream ofs(out_align);

	string prev_tagname = "";

	//chrom pairid offset line
	map<int, vector<SamRec> > mapped_reads;

	//map<int, vector<pair<char, string> > > paired_reads;

	//vector<PairedSamRec> paired_reads_ptr;

	//mapped_reads.

	size_t count = 0, paired_count = 0, unpaired_count = 0;

	bool paired = false;

	bool end1_mapped = false, end2_mapped = false;

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "")
				continue;

			char tagname[1000], chrom[100], mapped[100], seq[1000], qual_str[1000], alters[1000];

			unsigned short strand, something, mismatch;

			size_t offset;

			string alterstr = "";

			string qualstr = "I";

			char mate_match;

			size_t mate_offest, mate_diff;

			//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
			//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

			size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%c\t%llu\t%llu\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand, chrom, &offset, &something, mapped, &mate_match, &mate_offest, &mate_diff, seq, qual_str, &mismatch, alters);

			if (read_count == 10)
				alterstr = alters;

			string tagnamestr = tagname;

			int pairid = atoi(tagnamestr.substr(tagnamestr.length() - 1, 1).c_str());

			size_t st_idx = 0;//tagnamestr.find("~");			

			tagnamestr = tagnamestr.substr(0, tagnamestr.length() - 2);			

			++count;

			if (prev_tagname.empty() || prev_tagname == tagnamestr)
			{
				mapped_reads[pairid].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, alterstr, qual_str, line));

				if (pairid == 1)
					end1_mapped = true;

				if (pairid == 2)
					end2_mapped = true;

				if (prev_tagname.empty())
					prev_tagname = tagnamestr;

				//tag_count++;
			}
			else
			{
				//hash_map<string, int>::iterator msi_iter;

				map<int, vector<SamRec> >::iterator pair_iter1 = mapped_reads.begin();

				if (end1_mapped && end2_mapped)
				{
					for (pair_iter1 = mapped_reads.begin(); pair_iter1 != mapped_reads.end(); ++pair_iter1)
					{
						vector<SamRec>::iterator sam_iter;

						if (pair_iter1->second.size() > 1)
						{
							cout << "size large" << endl;
							cout << pair_iter1->second.front().line<<endl;
						}

						//for (sam_iter = pair_iter1->second.begin(); sam_iter != pair_iter1->second.end(); ++sam_iter)
						//{
						//	if (sam_iter->strand_t && 16)
						//		paired_reads_ofs << '>' << sam_iter->tag_name<<endl << revcomp(sam_iter->mapped_seq)<<endl;
						//	else
						//		paired_reads_ofs << '>' << sam_iter->tag_name<<endl << sam_iter->mapped_seq<<endl;

						//	paired_reads_sam_ofs << sam_iter->line << endl;
						//}
					}

					string generated_align = generate_trans_align(mapped_reads, trans_map, read_trans_map, "");


					if (!generated_align.empty())
						ofs << generated_align << endl;
				}
				else
				{
					for (pair_iter1 = mapped_reads.begin(); pair_iter1 != mapped_reads.end(); ++pair_iter1)
					{
						cout << pair_iter1->second.front().line << endl;

						vector<SamRec>::iterator sam_iter;

						if (pair_iter1->second.size() > 1)
						{
							cout << "size large" << endl;
							cout << pair_iter1->second.front().line<<endl;
						}

						//for (sam_iter = pair_iter1->second.begin(); sam_iter != pair_iter1->second.end(); ++sam_iter)
						//{
						//	if (sam_iter->strand_t && 16)
						//		single_reads_ofs << '>' << sam_iter->tag_name<<endl << revcomp(sam_iter->mapped_seq)<<endl;
						//	else
						//		single_reads_ofs << '>' << sam_iter->tag_name<<endl << sam_iter->mapped_seq<<endl;

						//	single_reads_sam_ofs << sam_iter->line << endl;
						//}
					}
				}

				//paired_reads.clear();

				//paired_reads_ptr.clear();

				mapped_reads.clear();

				end1_mapped = false, end2_mapped = false;

				mapped_reads[pairid].push_back(SamRec(tagname, strand, chrom, offset, something, mapped, seq, mismatch, alterstr, qual_str, line));

				if (pairid == 1)
					end1_mapped = true;

				if (pairid == 2)
					end2_mapped = true;

				paired = false;				

				prev_tagname = tagnamestr;

			}		
		}
		ifs.close();
	}
	else
	{
		cout << "can't open file "<< infile<<endl; exit(1);
	}

	map<int, vector<SamRec> >::iterator pair_iter1 = mapped_reads.begin();

	if (end1_mapped && end2_mapped)
	{
		for (pair_iter1 = mapped_reads.begin(); pair_iter1 != mapped_reads.end(); ++pair_iter1)
		{
			vector<SamRec>::iterator sam_iter;

			if (pair_iter1->second.size() > 1)
			{
				cout << "size large" << endl;
				cout << pair_iter1->second.front().line<<endl;
			}

			

			//for (sam_iter = pair_iter1->second.begin(); sam_iter != pair_iter1->second.end(); ++sam_iter)
			//{
			//	if (sam_iter->strand_t && 16)
			//		paired_reads_ofs << '>' << sam_iter->tag_name<<endl << revcomp(sam_iter->mapped_seq)<<endl;
			//	else
			//		paired_reads_ofs << '>' << sam_iter->tag_name<<endl << sam_iter->mapped_seq<<endl;

			//	paired_reads_sam_ofs << sam_iter->line << endl;
			//}
		}

		string generated_align = generate_trans_align(mapped_reads, trans_map, read_trans_map, "");

		if (!generated_align.empty())
			ofs << generated_align << endl;
	}
	else
	{
		for (pair_iter1 = mapped_reads.begin(); pair_iter1 != mapped_reads.end(); ++pair_iter1)
		{
			cout << pair_iter1->second.front().line << endl;

			vector<SamRec>::iterator sam_iter;

			if (pair_iter1->second.size() > 1)
			{
				cout << "size large" << endl;
				cout << pair_iter1->second.front().line<<endl;
			}

			//for (sam_iter = pair_iter1->second.begin(); sam_iter != pair_iter1->second.end(); ++sam_iter)
			//{
			//	if (sam_iter->strand_t && 16)
			//		single_reads_ofs << '>' << sam_iter->tag_name<<endl << revcomp(sam_iter->mapped_seq)<<endl;
			//	else
			//		single_reads_ofs << '>' << sam_iter->tag_name<<endl << sam_iter->mapped_seq<<endl;

			//	single_reads_sam_ofs << sam_iter->line << endl;
			//}
		}
	}
}



void 
read_transcript(const char* transcript_file, hash_map<string, transcript>& trans_map)
{
	ifstream ifs(transcript_file);

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "")
				continue;

			transcript ts(line);

			trans_map.insert(make_pair(ts.trans_name, ts));
		}
		ifs.close();
	}
	else
	{
		cout << "can't open file "<< transcript_file<<endl; exit(1);
	}
}

void 
read_read2transcript(const char* read_transcritp_file, hash_map<string, string>& read_trans_map)
{
	ifstream ifs(read_transcritp_file);

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "")
				continue;

			char read[1000], trans[1000];

			sscanf(line.c_str(), "%s\t%s", read, trans);

			read_trans_map[read] = trans;
		}
		ifs.close();
	}
	else
	{
		cout << "can't open file "<< read_transcritp_file<<endl; exit(1);
	}
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
read_junction(const char* fusion_junction_file, hash_map<string, hash_map<size_t, hash_map<size_t, pair<string, bool> > > > & fusion_junc_map)
{
	ifstream ifs(fusion_junction_file);

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "" || line.find("track") != string::npos)
				continue;

			char read[1000];

			size_t start, end;

			//chr11chr17	125778183	56787220	JUNC_4	18	++	255,0,0	2	60,63	0,-68990963	2.73634	5	GTAG	CTAAAGTCCAAGAGATGGCTCTCCCTATGATGCTGGCACATCCGTTCGACTAGTGATAGTGGATGGTATTGCTTT	0	0	0	32	36	3

			sscanf(line.c_str(), "%s\t%llu\t%llu", read, &start, &end);

			fusion_junc_map[read][start][end] = make_pair(line, false);
		}
		ifs.close();
	}
	else
	{
		cout << "can't open file "<< fusion_junction_file<<endl; exit(1);
	}
}

void 
read_db(const char* fusion_junction_file, hash_map<string, vector<fusion_db_rec> > & fusion_junc_map)
{
	ifstream ifs(fusion_junction_file);

	if (ifs.is_open())
	{
		string line;

		getline(ifs,line);

		while (getline(ifs,line))
		{
			char gene1[1000], gene2[1000], skip[1000], chr1[1000], chr2[1000];

			size_t start, end;

			//ISL1	ISL1	LUSC-60-2695	2138	5	0	50619521	5	0	50620655	inversion	1134	93	1	1	1	
			//BR-MEX-050:1,25_samples:0	50619185	50619521	337	72	50620335	50620655	321	71	IGR: 59Kb before ISL1(+)	IGR: 58Kb before ISL1(+)	-	0	0	3.00E-02	3.00E-02	3	0	2	1	1.64E+00	1.54E+00	6.84E-01	6.36E+01	1	6.36E+01	1		


			//individual	num	chr1	str1	pos1	chr2	str2	pos2	class	span	tumreads	normreads	normpanelreads	normpanelsamps	normpaneldetails	
			//min1	max1	range1	stdev1	min2	max2	range2	stdev2	gene1	site1	gene2	site2	fusion	fmapqzT1	fmapqzN1	fmapqzT2	fmapqzN2	nuwpT1	nuwpN1	nuwpT2	nuwpN2	zstdev1	zstdev2	quality	score	somatic	somatic_score	BPtry


			//individual	num	chr1	str1	pos1	chr2	str2	pos2	class	span	tumreads	normreads	normpanelreads	normpanelsamps	normpaneldetails	
			//min1	max1	range1	stdev1	min2	max2	range2	stdev2	gene1	site1	gene2	site2	fusion	fmapqzT1	fmapqzN1	fmapqzT2	fmapqzN2	nuwpT1	
			//nuwpN1	nuwpT2	nuwpN2	zstdev1	zstdev2	quality	score	somatic	somatic_score	BPtry

//LUSC-66-2756	4252	9	0	25220040	9	0	25309748	inversion	89708	443	2	0	0	26_samples:0	25219774	25220040	267	49	25309388	25309748	361	55	TUSC1	IGR: 456Kb before TUSC1(-)	TUSC1	IGR: 367Kb before TUSC1(-)	-	0	0	2.00E-02	1.20E-01	29	2	6	3	-4.15E-01	9.78E-02	9.07E-01	4.02E+02	1	4.02E+02	1

			for (size_t i = 0; i < line.size(); ++i)
			{
				if (line[i] == ' ')
					line[i] = '_';
			}

			sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%llu\t%s\t%s\t%llu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
				skip, skip, chr1, skip, &start, chr2, skip, &end, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, skip, gene1, skip, gene2);

			string chrstr1(chr1), chrstr2(chr2);

			string comb;

			comb =  "chr" + chrstr1 + "~" + "chr" + chrstr2;

			fusion_junc_map[comb].push_back(fusion_db_rec(start, end, line));

			//if (chrstr1 > chrstr2 || (chrstr1 == chrstr2 && start > end ) )
			//{
			comb =  "chr" + chrstr2 + "~" + "chr" + chrstr1;

			swap(start, end);

			fusion_junc_map[comb].push_back(fusion_db_rec(start, end, line));
			//}
			//else
			//{
				
			//}

			
		}

		ifs.close();

	}
	else
	{
		cout << "can't open file "<< fusion_junction_file<<endl; exit(1);
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

void 
read_junction_comp2db(const char* fusion_junction_file, const char* comp_fusion, hash_map<string, vector<fusion_db_rec> > & fusion_junc_map, size_t range, string sampleid, string dataname)
{
	ifstream ifs(fusion_junction_file);
	
	ofstream ofs(comp_fusion);
//	inversion
//long_range
//inter_chr
//deletion
//tandem_dup

	string stats = comp_fusion; stats.append(".stats"); ofstream stats_ofs(stats.c_str());

	//string inversion = comp_fusion; inversion.append("inversion"); ofstream inversion_ofs(inversion.c_str());

	//string long_range = comp_fusion; long_range.append("long_range"); ofstream long_range_ofs(long_range.c_str());

	//string inter_chr = comp_fusion; inter_chr.append("inter_chr"); ofstream inter_chr_ofs(inter_chr.c_str());

	//string deletion = comp_fusion; deletion.append("deletion"); ofstream deletion_ofs(deletion.c_str());

	//string tandem_dup = comp_fusion; tandem_dup.append("tandem_dup"); ofstream tandem_dup_ofs(tandem_dup.c_str());

	string mpsonly = comp_fusion; mpsonly.append(".mpsonly"); ofstream mpsonly_ofs(mpsonly.c_str());

	string dbonly = comp_fusion; dbonly.append(".drangeonly"); ofstream dbonly_ofs(dbonly.c_str());

	size_t inversion_t = 0, long_range_t = 0, inter_chr_t = 0, deletion_t = 0, tandem_dup_t = 0;

	size_t shared;

	size_t mps_count = 0;

	map<string, int> matchedlines;

	if (ifs.is_open())
	{
		string line;
		while (getline(ifs,line))
		{
			if (line == "" || line.find("track") != string::npos)
				continue;

			++mps_count;

			char skip[1000], gene1[10000], gene2[10000], chrname[1000];

			size_t start, end;

			//chr11chr17	125778183	56787220	JUNC_4	18	++	255,0,0	2	60,63	0,-68990963	2.73634	5	GTAG	CTAAAGTCCAAGAGATGGCTCTCCCTATGATGCTGGCACATCCGTTCGACTAGTGATAGTGGATGGTATTGCTTT	0	0	0	32	36	3


			//chrX_chrX	153608727	153608637	JUNC_55394	1	++	255,0,0	2	115,17,	0,206,	0.000000	6	GTAG	0	0	
			//0.000000	0	17	98	1	0	0	153608612	153608654	153608637,17M|	153608614,114M|

			sscanf(line.c_str(), "%s\t%llu\t%llu", chrname, &start, &end);

			//string chrstr = chrname;

			//string chr1 = chrstr.substr(0, chrstr.find("~"));

			//string chr2 = chrstr.substr(chrstr.find("~") + 1, chrstr.length() - 1);
			//if (start > 
			//from_fusion	overlapping	SAMD11,NOC2L,	SAMD11,NOC2L,	chr1_chr1	879733	879783	JUNC_23683	0	+-	255,0,0	2	
			//35,113,	4391856,360197048,	0.000000	0		0	0	0.000000	35	0	62	2	0	2879672:879734=~	879672:879784=~

			if (fusion_junc_map.find(chrname) != fusion_junc_map.end())
			{
				/*cout << chrname <<endl;*/

				vector<fusion_db_rec>& cur_vec = (fusion_junc_map.find(chrname))->second;

				bool matched = false;

				for (size_t i = 0; i < cur_vec.size(); ++i)
				{
					string fusion_db_line = cur_vec[i].line;

					if (((cur_vec[i].start > start && cur_vec[i].start - start < range) ||
						(cur_vec[i].start <= start && start - cur_vec[i].start < range)) &&
						((cur_vec[i].end > end && cur_vec[i].end - end < range) ||
						(cur_vec[i].end <= end && end - cur_vec[i].end < range)))
					{

						if (fusion_db_line.find("inversion") != string::npos)
						{
							ofs << fusion_db_line << "\t----\t----\t";

							ofs << line << endl;

							++inversion_t;
						}

						if (fusion_db_line.find("long_range") != string::npos)
						{
							ofs << fusion_db_line << "\t----\t----\t";

							ofs << line << endl;
							
							++long_range_t;
						}

						if (fusion_db_line.find("inter_chr") != string::npos)
						{
							ofs << fusion_db_line << "\t----\t----\t";

							ofs << line << endl;
							
							++inter_chr_t;
						}

						if (fusion_db_line.find("deletion") != string::npos)
						{
							ofs << fusion_db_line << "\t----\t----\t";

							ofs << line << endl;
							
							++deletion_t;
						}

						if (fusion_db_line.find("tandem_dup") != string::npos)
						{
							ofs << fusion_db_line << "\t----\t----\t";

							ofs << line << endl;
							
							++tandem_dup_t;
						}

						matchedlines[fusion_db_line] = 0;

						cur_vec[i].matched = true;

						matched = true;

						break;
					}
				}

				if (matched == false)
					mpsonly_ofs << line << endl;

			}
			else
				mpsonly_ofs << line << endl;
		}

		ifs.close();
	}
	else
	{
		cout << "can't open file "<< fusion_junction_file<<endl; exit(1);
	}

	hash_map<string, vector<fusion_db_rec> >::iterator chr_iter;

	vector<fusion_db_rec>::iterator fd_iter;

	size_t notmatch = 0;

	size_t dbcout = 0;

	size_t dbmatched = 0;

	size_t dbmatchedlongintershared = 0;

	map<string, int> matched_dbs;

	for (chr_iter = fusion_junc_map.begin(); chr_iter != fusion_junc_map.end(); ++chr_iter)
	{
		for (fd_iter = chr_iter->second.begin(); fd_iter != chr_iter->second.end(); ++fd_iter)
		{
			if (fd_iter->matched == false && matchedlines.find(fd_iter->line) == matchedlines.end())
				dbonly_ofs << fd_iter->line << endl;
			else
			{
				if (matched_dbs.find(fd_iter->line) == matched_dbs.end())
					matched_dbs[fd_iter->line] = 1;
				else
					++matched_dbs[fd_iter->line];

				//++dbmatched;

				
			}

			++dbcout;
		}
	}

	dbmatched = matched_dbs.size();

	map<string, int>::iterator dbs_iter;

	for (dbs_iter = matched_dbs.begin(); dbs_iter != matched_dbs.end(); ++dbs_iter)
	{
		if (dbs_iter->first.find("long_range") != string::npos || dbs_iter->first.find("inter_chr") != string::npos)
			++dbmatchedlongintershared;
	}

	shared = inversion_t + long_range_t + inter_chr_t + deletion_t + tandem_dup_t;

	size_t longintershared = long_range_t + inter_chr_t;

	stats_ofs << sampleid <<'\t'<< dataname << '\t' <<dbcout / 2<< '\t' << dbmatched << '\t' << dbmatchedlongintershared << '\t' << mps_count - shared + dbmatched << '\t' << mps_count<<endl;

	cout << sampleid <<'\t'<< dataname << endl;
	cout << "inversion:\t"<<inversion_t<<endl;
	cout << "long_range:\t"<<long_range_t<<endl;
	cout << "inter_chr:\t"<<inter_chr_t<<endl;
	cout << "deletion:\t"<<deletion_t<<endl;
	cout << "tandem_dup:\t"<<tandem_dup_t<<endl;
}



int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "infile unspliced_reads unique_spliced_reads multiple_spliced_reads" <<endl;
		exit(0);
	}
	
	const char* fusion_file = argv[1];

	const char* fusion_db = argv[2];

	const char* fusion_matched = argv[3];

	size_t range = atoi(argv[4]);

	cout << fusion_db << '\t' << fusion_file << '\t' <<fusion_matched<< endl;

	hash_map<string, vector<fusion_db_rec> >  fusion_db_map;

	read_db(fusion_db, fusion_db_map);

	read_junction_comp2db(fusion_file, fusion_matched, fusion_db_map, range, argv[5], argv[6]);
}