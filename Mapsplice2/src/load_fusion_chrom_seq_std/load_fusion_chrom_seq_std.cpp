#include <algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <iterator>
#include <string.h>

using namespace std;

bool filter = false;

bool fulllength = false;

void
readchrom(const char* filename, string& longseq)
{
	long size;  

	ifstream longfile(filename);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	if (longfile.is_open())
	{
		string skipline;
		getline(longfile,skipline);
		int count = 0;
		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);
			longseq.append(line);

			if (line.length() == 0)
				continue;

			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);

		}

		longfile.close();
	}
	else cout << "Unable to open file";
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

void load_fusion_chrom_seq(string chrom_dir, string fusion_junc, string output_fusion_junc, size_t length, size_t anchor_diff, size_t min_anchor, double minentropy)
{
	map<string, string> loaded_chromos;

	ifstream ifs(fusion_junc.c_str());

	ofstream ofs(output_fusion_junc.c_str());

	string output_fusion_junc_filtered = output_fusion_junc; output_fusion_junc_filtered.append(".filtered");

	ofstream ofs_filtered(output_fusion_junc_filtered.c_str());

	ofstream ofs_1_fa((output_fusion_junc + ".1.fa").c_str());

	ofstream ofs_2_fa((output_fusion_junc + ".2.fa").c_str());

	size_t count = 0;

	if (ifs.is_open())
	{
		int count = 0;
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);

			if (line.length() == 0)
				continue;
			
			//chr10~chr10	100903995	100481585	--	1	4	6	4	4	6	6	1.03972	1.03972	42	31	6	
			//100481444,142M|	100903995,162M87948N15M|	100481444,142M|	100903995,162M87948N15M|	doner_exact_matched	acceptor_exact_matched

			char chromname[1000], strand1, strand2;

			size_t start, end;

			//size_t number_of_dataset_supported, sum_of_support_reads;

			//size_t sum_of_encompass_reads, max_support;

			//size_t min_support, max_encompass;

			//size_t min_encompass;
			
			//double max_entropy, min_entropy;

			size_t max_doner, max_acceptor;

			size_t min_anchor_diff;

			//Chromosome	Fusion_end_1_location	Fusion_end_2_location	Fusion_Strand	number_of_dataset_supported	sum_of_support_reads	
			//sum_of_encompass_reads	max_support	min_support	max_encompass	min_encompass	max_entropy	min_entropy	max_doner	max_acceptor
			//chr6~chr6	29856519	31323583	+-	4	23	1550	7	4	488	236	1.54983	1.33218	44	41

			;

			//chr10~chr10	51276994	51402778	JUNC_8	11	-+	255,0,0	2	46,35,6327,14330,	0,125820,	1.540306	5	GTAG	0	2	0.272727	
			
			//24	22	3	11	0	9	6	3	0	9	2	34	51277040	51402813	
			
			//51276994,111M3294N59M2452N85M229N96M|51276994,111M5805P85M229N96M|51276994,111M6119P96M|	51402778,121M14027P181M|	0	0	0.24	0.157172	351	207	302	302

			char skip[1000];

			//sprintf(buf, "%s~%s\t%llu\t%llu\tJUNC_%llu\t%u\t%c%c\t255,0,0\t2\t%llu,%llu,%llu,%llu,\t%llu,%llu,\t%lf\t%hu\t%s\t%hu\t%hu\t%lf\t%hu\t%hu\t%hu\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u",
			//	m_chrom2.c_str(), m_chrom.c_str(), 
			//	m_end, m_start, junc_id, m_hits, strand1, strand2, m_max_suffix_len, m_max_prefix_len, m_max_fusion_suffix_len, m_max_fusion_prefix_len, m_start_blockoffset, m_end_blockoffset, m_entropy, 
			//	6, revcomp(m_flankstring).c_str(), m_min_mismatch, m_max_mismatch, m_ave_mismatch, m_max_min_suffix, m_max_min_prefix, m_min_anchor_difference, 
			//	m_unique_count, m_multi_count, 
			//	m_paired_count, m_left_paired_count, m_right_paired_count, m_paired_mutiple_count, m_paired_unique_count, 
			//	m_single_count, 
			//	m_encompass_reads_count, fusion_doner_start, fusion_acceptor_end);

			//sscanf(line.c_str(), "%s\t%llu\t%llu\t%c%c\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%lf\t%lf\t%llu\t%llu\t%llu", chromname, &start, &end, &strand1, &strand2, 
			//	&number_of_dataset_supported, &sum_of_support_reads, &sum_of_encompass_reads, &max_support, &min_support, &max_encompass, &min_encompass, 
			//	&max_entropy, &min_entropy, &max_doner, &max_acceptor, &min_anchor_diff);

			double entropy;

			double min_mis, max_mis, ave_mis;

			char block[1000];

			//new added columns

			size_t /*max_min_suffix, max_min_prefix, min_anchor_difference,*/ unique_count;

			//"%llu\t";

			//, &unique_count;

			//21	22	23	24	25	26	27	28	29	30
			//multi_count	paired_count	left_paired_count	right_paired_count	paired_multi_count	paired_unique_count	single_count	encompass_count	doner_start	acceptor_end

			size_t multi_count, paired_count, left_paired_count, right_paired_count, paired_multi_count, paired_unique_count, single_count, encompass_count, doner_start, acceptor_end;

			//"%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t"

			//&multi_count, &paired_count, &left_paired_count, &right_paired_count, &paired_multi_count, &paired_unique_count, &single_count, &encompass_count, &doner_start, &acceptor_end;

			//31	32	33	34	35	36	37	38	39	40
			//doner_iosforms	acceptor_isoforms	notused	notused	notused	notused	minimal_doner_isoform_length	maximal_doner_isoform_length	minimal_acceptor_isoform_length	maximal_acceptor_isoform_length

			char doner_iosforms[50000], acceptor_isoforms[50000];

			double uni_score1, uni_score2, ks_score1, ks_score2;

			size_t minimal_doner_isoform_length, maximal_doner_isoform_length, minimal_acceptor_isoform_length, maximal_acceptor_isoform_length;

			//"%s\t%s\t%lf\t%lf\t%lf\t%lf\t%llu\t%llu\t%llu\t%llu\t"

			//doner_iosforms, acceptor_isoforms, &uni_score1, &uni_score2, &ks_score1, &ks_score2, &minimal_doner_isoform_length, &maximal_doner_isoform_length, &minimal_acceptor_isoform_length, &maximal_acceptor_isoform_length

			//41	42	43	44	45	46	47	48	49	50
			//fusion_pair_entropy	ave_mis_per_bp	anchor_score	max_doner	max_acceptor	max_fragment	min_fragment	ave_fragment	doner_match_to_normal	acceptor_match_to_normal

			double fusion_pair_entropy, ave_mis_per_bp,	anchor_score;

			size_t max_doner_frag, max_acceptor_frag, max_fragment, min_fragment;

			double ave_fragment;

			size_t doner_encompass_unique, doner_encompass_multiple, acceptor_encompass_unique, acceptor_encompass_multiple;

			//"%lf\t%lf\t%lf\t%llu\t%llu\t%llu\t%llu\t%lf\t%llu\t%llu\t%llu\t%llu"

			//&fusion_pair_entropy, &ave_mis_per_bp, &anchor_score, &max_doner, &max_acceptor, &max_fragment, &min_fragment, &ave_fragment, &doner_encompass_unique, &doner_encompass_multiple, &acceptor_encompass_unique, &acceptor_encompass_multiple;

			//"\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%s\t%s\t%lf\t%lf\t%lf\t%lf\t%llu\t%llu\t%llu\t%llu\t%lf\t%lf\t%lf\t%llu\t%llu\t%llu\t%llu\t%lf\t%llu\t%llu\t%llu\t%llu"

			//, &unique_count, &multi_count, &paired_count, &left_paired_count, &right_paired_count, &paired_multi_count, &paired_unique_count, &single_count, &encompass_count, &doner_start, &acceptor_end, doner_iosforms, acceptor_isoforms, &uni_score1, &uni_score2, &ks_score1, &ks_score2, &minimal_doner_isoform_length, &maximal_doner_isoform_length, &minimal_acceptor_isoform_length, &maximal_acceptor_isoform_length, &fusion_pair_entropy, &ave_mis_per_bp, &anchor_score, &max_doner, &max_acceptor, &max_fragment, &min_fragment, &ave_fragment, &doner_encompass_unique, &doner_encompass_multiple, &acceptor_encompass_unique, &acceptor_encompass_multiple;
			
			sscanf(line.c_str(), "%s\t%llu\t%llu\t%s\t%s\t%c%c\t%s\t%s\t%s\t%s\t%lf\t%s\t%s\t%lf\t%lf\t%lf\t%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%s\t%s\t%lf\t%lf\t%lf\t%lf\t%llu\t%llu\t%llu\t%llu\t%lf\t%lf\t%lf\t%llu\t%llu\t%llu\t%llu\t%lf\t%llu\t%llu\t%llu\t%llu", 
				chromname, &start, &end, skip, skip, &strand1, &strand2, skip, skip, /*&max_doner, &max_acceptor, skip,*/ block, skip,
				&entropy, skip, skip, &min_mis, &max_mis, &ave_mis, skip, skip, &min_anchor_diff, &unique_count, &multi_count, &paired_count, &left_paired_count, &right_paired_count, &paired_multi_count, 
				&paired_unique_count, &single_count, &encompass_count, &doner_start, &acceptor_end, doner_iosforms, acceptor_isoforms, &uni_score1, &uni_score2, &ks_score1, 
				&ks_score2, &minimal_doner_isoform_length, &maximal_doner_isoform_length, &minimal_acceptor_isoform_length, &maximal_acceptor_isoform_length, &fusion_pair_entropy, 
				&ave_mis_per_bp, &anchor_score, &max_doner_frag, &max_acceptor_frag, &max_fragment, &min_fragment, &ave_fragment, &doner_encompass_unique, &doner_encompass_multiple, &acceptor_encompass_unique, &acceptor_encompass_multiple);

			size_t max_doner_extend, max_acceptor_extend;
			
			sscanf(block, "%llu,%llu,%llu,%llu", &max_doner, &max_acceptor, &max_doner_extend, &max_acceptor_extend);

			string chrnamestr = chromname;

			string chr1, chr2;

			char chr1chr[100], chr2chr[100];

			sscanf(chrnamestr.c_str(), "%[^~]~%[^~]", chr1chr, chr2chr);

			//size_t idx = chrnamestr.find("chr", 4);

			//idx = idx;

			chr1 = chr1chr; //chrnamestr.substr(0, idx - 1);

			chr2 = chr2chr;//.substr(idx, chrnamestr.length() - idx);

			if (loaded_chromos.find(chr1) == loaded_chromos.end())
			{
				string chromfile = chrom_dir + "/" + chr1 + ".fa";

				readchrom(chromfile.c_str(), loaded_chromos[chr1]);
			}

			if (loaded_chromos.find(chr2) == loaded_chromos.end())
			{
				string chromfile = chrom_dir + "/" + chr2 + ".fa";

				readchrom(chromfile.c_str(), loaded_chromos[chr2]);
			}

			string& chr1seq = loaded_chromos[chr1];

			string& chr2seq = loaded_chromos[chr2];

			string chr1substr, chr2substr;

			size_t length1 = length, length2 = length;

			if (length1 > max_doner && fulllength == false)
				length1 = max_doner;

			if (length2 > max_acceptor && fulllength == false)
				length2 = max_acceptor;

			if (filter)
			{
				if ((length1 < min_anchor && max_doner_extend < min_anchor) || (fulllength && max_doner < min_anchor && max_doner_extend < min_anchor) )
				{
					ofs_filtered << line << "\tdoner short anchor:" << max_doner << endl;

					continue;
				}

				if ((length2 < min_anchor && max_acceptor_extend < min_anchor) || (fulllength && max_acceptor < min_anchor && max_acceptor_extend < min_anchor) )
				{
					ofs_filtered << line << "\tacceptor short anchor:" << max_acceptor << endl;

					continue;
				}

				if (min_anchor_diff > anchor_diff)
				{
					ofs_filtered << line << "\tlarge anchor diff:" << min_anchor_diff << endl;

					continue;
				}

				if (minentropy > entropy && minentropy > fusion_pair_entropy && ((fusion_pair_entropy + entropy) < (2 * minentropy) ))
				{
					ofs_filtered << line << "\tmin entropy:" << entropy << '\t' << fusion_pair_entropy<< endl;

					continue;
				}
				

				if (chr1.find("M") != string::npos || chr2.find("M") != string::npos ||
					chr1.find("Un") != string::npos || chr2.find("Un") != string::npos ||
					chr1.find("random") != string::npos || chr2.find("random") != string::npos)
				{
					ofs_filtered << line << "\tfusion in chromosome M Un or random" << endl;

					continue;
				}
			}

			++count;

			if (strand1 == '-')
			{
				if (start + length1 > chr1seq.length())
				{
					//cout << line << endl;

					ofs_filtered << line << "\tdoner anchor exceed chromosome boundary:" << start + length1 << "\tVS\t"<<chr1seq.length() << endl;

					continue;
				}

				chr1substr = chr1seq.substr(start - 1, length1);

				chr1substr = revcomp(chr1substr);
			}
			else
			{
				if (start < length1 )
				{
					//cout << line << endl;

					ofs_filtered << line << "\tdoner anchor exceed chromosome boundary:" << start - length1 << "\tVS\t"<<chr1seq.length() << endl;

					continue;
				}

				chr1substr = chr1seq.substr(start - length1, length1);
			}

			if (strand2 == '-')
			{
				if (end < length2 )
				{
					//cout << line << endl;

					ofs_filtered << line << "\tacceptor anchor exceed chromosome boundary:" << end - length2 << "\tVS\t"<<chr2seq.length() << endl;

					continue;
				}

				chr2substr = chr2seq.substr(end - length2, length2);

				chr2substr = revcomp(chr2substr);
			}
			else
			{
				if (end + length2 > chr2seq.length())
				{
					//cout << line << endl;

					ofs_filtered << line << "\tacceptor anchor exceed chromosome boundary:" << end + length2 << "\tVS\t"<<chr2seq.length() << endl;

					continue;
				}

				chr2substr = chr2seq.substr(end - 1, length2);

			}

			ofs << line << '\t' << chr1substr << '\t' << revcomp(chr2substr) << endl;

			ofs_1_fa << ">" << count << ":"<<chr1<<":"<<start<<":"<<strand1<< ":" << chr2<<":"<<end<<":"<<strand2<< "/1"<<endl;

			ofs_1_fa << chr1substr <<endl;

			ofs_2_fa << ">" << count << ":"<<chr1<<":"<<start<<":"<<strand1<< ":" << chr2<<":"<<end<<":"<<strand2<< "/2"<< endl;

			ofs_2_fa << revcomp(chr2substr) <<endl;
		}

		ifs.close();
	}
	else cout << "Unable to open file :"<<fusion_junc<<endl;

}

int 
main(int argc, const char** argv)
{
	if (argc < 4)
	{
		cout << "chr start length"<<endl;
		return 0;
	}

	if (argc > 6)
		filter = true;

	if (argc > 9)
		fulllength = true;

	string longseq;

	size_t length = atoi(argv[2]);

	size_t anchor_diff = atoi(argv[3]);

	size_t min_anchor = atoi(argv[4]);

	double minentropy;

	char* stop;

	minentropy = strtod(argv[8], &stop);

	load_fusion_chrom_seq(argv[1], argv[5], argv[6], length, anchor_diff, min_anchor, minentropy);

	return 0;
}