#include <algorithm>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <iterator>
#include <string.h>

using namespace std;


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

void load_fusion_chrom_seq(string input_junc, string output_fusion_junc, string output_normal_junc, string output_circular_RNAs)
{
	ifstream ifs(input_junc.c_str());

	ofstream ofs(output_fusion_junc.c_str());

	ofstream ofs_normal(output_normal_junc.c_str());

	ofstream ofs_circular_RNAs(output_circular_RNAs.c_str());

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

			++count;
			//chr10~chr10	100903995	100481585	--	1	4	6	4	4	6	6	1.03972	1.03972	42	31	6	
			//100481444,142M|	100903995,162M87948N15M|	100481444,142M|	100903995,162M87948N15M|	doner_exact_matched	acceptor_exact_matched


			//from_fusion	fusion	-,+	MTOR,	RBP7,	chr1~chr1	11288725	10071334	
			//--	1	5	34	5	5	34	34	1.05492	1.05492	34	34	14	10071167,168M|	
			//11288725,226M|	10071167,168M|	11288725,226M|	doner_exact_matched	not_matched	ATTCGAGTCTGTGATGGGGCCATCCGGGAA	tttgtattagagaaggaagaggttcacctg

			char chromname[1000], strand1, strand2;

			char from[1000], type[1000], genestrand[100];

			char gene1[1000], gene2[1000];

			size_t start, end;

			//cout << line << endl;
			//from_fusion	intergenic		-,	DDX11L1,	chr1~chr1	10507	12010	-+	1	18	28	18	18	28	28	1.95053	1.95053	43	39	2	
			//10507,107M1435P160M|10507,107M1997P111M499N17M|	12010,202M|	10507,107M1435P160M|10507,107M1997P111M499N17M|	12010,202M|	
			//not_matched	acceptor_exact_matched	actctgaaggcggagcacagttctcctcag	ACAGGCCAGCAGTTGCTGGAAGTCAGACAC	0

			//from	type	gene_strand	gene1	gene2	Chromosome	Fusion_end_1_location	Fusion_end_2_location	Fusion_Strand	
			//number_of_dataset_supported	sum_of_support_reads	sum_of_encompass_reads	max_support	min_support	max_encompass	
			//min_encompass	max_entropy	min_entropy	max_doner	max_acceptor	min_anchor_difference

			size_t number_of_dataset_supported, sum_of_support_reads, sum_of_encompass_reads, max_support, min_support, max_encompass;

			size_t min_encompass, max_doner, max_acceptor, min_anchor_difference;

			double max_entropy, min_entropy;

			char skip[1000];
			//from_fusion	fusion	-,+	DST,	TRIM38,	chr6~chr6	56707853	25963277	JUNC_1405	4	-+	255,0,0	2	45,34,301,151,	0,30744622,	1.386294	
			//6	GTAG	0	1	0.250000	16	19	12	4	0	4	1	3	0	4	0	24	56707898	25963311	25963277,234M|	56707853,175M75P50M|	
			//0	0	0.337778	0.275278	234	234	225	225	doner_exact_matched	acceptor_exact_matched	GAGGATGTCCTGGAGAGGTACAAAG	AATCTTCTGAGATTTGGTGTAAAGT	1

			sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%s\t%s\t%c%c", from, type, gene1, gene2, chromname, &start, &end, skip, skip, &strand1, &strand2);

			//ofs_normal << line.substr(juncidx, line.length() - juncidx);

			//if (start == 104310906 && end == 104312004)
			//	int i = 0;

			string typestr = type;

			string genestrandstr = genestrand;

			string gene1str = gene1, gene2str = gene2;

			size_t firstdot = line.find(",	");

			size_t seconddot = line.find(",	", firstdot + 1);

			size_t chridx = seconddot + 2;

			string prefix = line.substr(0, chridx);

			string suffix = line.substr(chridx);

			string swapline = suffix + "\t" + prefix;

			if (typestr == "normal" && strand1 != strand2)
			{
				char chr1[1000];

				sscanf(chromname, "%[^~]", chr1);

				ofs_normal << /*from << '\t' << type<< '\t' <<gene1<< '\t' << gene2<< '\t' << */chr1<< '\t' <<start << '\t' <<end << '\t';// <<"+" << '\t';

				size_t juncidx = line.find("JUNC");

				ofs_normal << line.substr(juncidx, line.length() - juncidx) <<'\t' << prefix << endl;

				continue;
			}
			else if (typestr == "overlapping"/* && strand1 == strand2*/)
			{
				vector<string> donergenes;

				vector<string> acceptorgenes;

				size_t startidx = 0;
				
				size_t comma = 0;

				while (true)
				{
					comma = gene1str.find(',', startidx + 1);

					string cur_gene = gene1str.substr(startidx, comma - startidx);

					donergenes.push_back(cur_gene);

					startidx = comma + 1;

					if (startidx >= gene1str.length())
						break;
				}

				startidx = 0;
				
				comma = 0;

				while (true)
				{
					comma = gene2str.find(',', startidx + 1);

					string cur_gene = gene2str.substr(startidx, comma - startidx);

					acceptorgenes.push_back(cur_gene);

					startidx = comma + 1;

					if (startidx >= gene2str.length())
						break;
				}

				bool normal = false;

				for (int i = 0; i < donergenes.size(); ++i)
				{
					for (int j = 0; j < acceptorgenes.size(); ++j)
					{
						if (donergenes[i] == acceptorgenes[j])
						{
							normal = true;

							break;
						}
					}
				}

				if (normal)
				{
					if ( strand1 != strand2)
					{
						char chr1[1000];

						sscanf(chromname, "%[^~]", chr1);

						ofs_normal << /*from << '\t' << type<< '\t' <<gene1<< '\t' << gene2<< '\t' << */chr1<< '\t' <<start << '\t' <<end << '\t';// <<"+" << '\t';

						size_t juncidx = line.find("JUNC");

						ofs_normal << line.substr(juncidx, line.length() - juncidx) <<'\t' << prefix << endl;

						continue;
					}
					else if (strand1 == '+' && start > end)
					{
						ofs_circular_RNAs << swapline <<endl; 
						//ofs << swapline <<endl;
						continue;
					}
					else if (strand1 == '-' && end > start)
					{
						ofs_circular_RNAs << swapline <<endl; 
						//ofs << swapline <<endl;
						continue;
					}
					else if (strand1 == '-' && start > end)
					{
						char chr1[1000];

						sscanf(chromname, "%[^~]", chr1);

						ofs_normal << /*from << '\t' << type<< '\t' <<gene1<< '\t' << gene2<< '\t' << */chr1<< '\t'  <<start << '\t'<<end << '\t';// <<"+" << '\t';

						size_t juncidx = line.find("JUNC");

						ofs_normal << line.substr(juncidx, line.length() - juncidx) <<'\t' << prefix << endl;

						continue;
					}
					else if (strand1 == '+' && end > start)
					{
						char chr1[1000];

						sscanf(chromname, "%[^~]", chr1);

						ofs_normal << /*from << '\t' << type<< '\t' <<gene1<< '\t' << gene2<< '\t' << */chr1<< '\t' <<start  << '\t' << end << '\t';// <<"-" << '\t';

						size_t juncidx = line.find("JUNC");

						ofs_normal << line.substr(juncidx, line.length() - juncidx) <<'\t' << prefix << endl;

						continue;
					}
				}
				else
				{
					ofs << swapline <<endl;
						continue;
				}
			}
			else if (typestr != "normal" || gene1str.empty() || gene2str.empty() || gene1str != gene2str || strand1 != strand2)
			{
				ofs << swapline <<endl;
				continue;
			}
			else if (strand1 == '+' && start > end)
			{
				ofs_circular_RNAs << swapline <<endl;
				continue;
			}
			else if (strand1 == '-' && end > start)
			{
				ofs_circular_RNAs << swapline <<endl;
				continue;
			}
			else if (strand1 == '-' && start > end)
			{
				char chr1[1000];

				sscanf(chromname, "%[^~]", chr1);

				ofs_normal << /*from << '\t' << type<< '\t' <<gene1<< '\t' << gene2<< '\t' << */chr1<< '\t' <<start << '\t' <<end << '\t';// <<"+" << '\t';

				size_t juncidx = line.find("JUNC");

				ofs_normal << line.substr(juncidx, line.length() - juncidx) << '\t' << prefix <<endl;

				//ofs_normal << number_of_dataset_supported << '\t' << sum_of_support_reads<< '\t' <<max_support<< '\t' << min_support<< '\t' << max_entropy<< '\t' <<min_entropy << '\t' <<max_doner << '\t'; 

				//ofs_normal << max_acceptor << '\t' << min_anchor_difference<< endl;

				continue;
			}
			else if (strand1 == '+' && end > start)
			{
				char chr1[1000];

				sscanf(chromname, "%[^~]", chr1);

				ofs_normal << /*from << '\t' << type<< '\t' <<gene1<< '\t' << gene2<< '\t' << */chr1<< '\t' <<start  << '\t' << end << '\t';// <<"-" << '\t';

				size_t juncidx = line.find("JUNC");

				ofs_normal << line.substr(juncidx, line.length() - juncidx) << '\t' << prefix <<endl;

				//ofs_normal << number_of_dataset_supported << '\t' << sum_of_support_reads<< '\t' <<max_support<< '\t' << min_support<< '\t' << max_entropy<< '\t' <<min_entropy << '\t' <<max_doner << '\t'; 

				//ofs_normal << max_acceptor << '\t' << min_anchor_difference<< endl;

				continue;
			}
		}

		ifs.close();
	}
	else cout << "Unable to open file :"<<input_junc<<endl;

}

int 
main(int argc, const char** argv)
{
	if (argc < 3)
	{
		cout << "chr start length"<<endl;
		return 0;
	}

	load_fusion_chrom_seq(argv[1], argv[2], argv[3], argv[4]);

	return 0;
}