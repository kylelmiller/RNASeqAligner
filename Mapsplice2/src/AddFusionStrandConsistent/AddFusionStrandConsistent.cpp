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

void load_fusion_chrom_seq(string input_junc, string output_fusion_junc)
{
	ifstream ifs(input_junc.c_str());

	string match = output_fusion_junc + ".match"; ofstream ofs_match(match.c_str());

	string notmatch = output_fusion_junc + ".notmatch"; ofstream ofs_notmatch(notmatch.c_str());

	ofstream ofs(output_fusion_junc.c_str());

	size_t count = 0;

	if (ifs.is_open())
	{
		//string headline;

		//getline(ifs, headline);

		//ofs << headline << endl;

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

			char skip[1000];


			//from_fusion	normal		POGK,	POGK,	chr1~chr1	166823269	166823352	JUNC_779	3	--	255,0,0	2	34,32,161,210,	0,116,	1.098612	
			//5	GTAG	0	0	0	21	16	8	3	0	3	1	2	0	3	0	2	166823303	166823320	166823269,62M48P50M|	166823144,84M53P72M|	
			//0	0	0.506667	0.528889	112	112	156	156	not_matched	not_matched	GTAACTAGGACAGATGAAAAATTAG	GGGAGCTTAAAAGATTTTACAAGAC


			if (line.find("JUNC") != string::npos)
				sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%s\t%s\t%c%c", from, type, genestrand, gene1, gene2, chromname, &start, &end, skip, skip, &strand1, &strand2);
			else
				sscanf(line.c_str(), "%s\t%s\t%s\t%s\t%s\t%s\t%llu\t%llu\t%c%c", from, type, genestrand, gene1, gene2, chromname, &start, &end, &strand1, &strand2);

			string typestr = type;

			string genestrandstr = genestrand;

			if (typestr != "fusion" || genestrandstr.empty())
			{
				ofs_notmatch << line << '\t' << 0<<endl;
				//ofs_match << line << '\t' << 0<<endl;

				ofs << line << '\t' << 0<<endl;

				continue;
			}
			string doner_strands, acceptor_strands;

			size_t index = genestrandstr.find(",");

			doner_strands = genestrandstr.substr(0, index);

			acceptor_strands = genestrandstr.substr(index + 1, genestrandstr.length() - index - 1);

			if (doner_strands.length() > 1 || acceptor_strands.length() > 1)
			{
				bool donermatch = false, acceptormatch = false;

				for (size_t i = 0; i < doner_strands.length(); ++i)
				{
					if (doner_strands[i] == strand1)
						donermatch = true;
				}

				for (size_t i = 0; i < acceptor_strands.length(); ++i)
				{
					if (acceptor_strands[i] == strand2)
						acceptormatch = true;
				}

				if (donermatch && acceptormatch)
				{
					ofs_match << line << '\t' << 1<<endl;

					ofs << line << '\t' << 1<<endl;
				}
				else
				{
					ofs_notmatch << line << '\t' << 0<<endl;

					ofs << line << '\t' << 0<<endl;
				}
			}
			else
			{
				if (doner_strands[0] != strand1 || acceptor_strands[0] != strand2)
				{
					ofs_notmatch << line << '\t' << 0<<endl;

					ofs << line << '\t' << 0<<endl;
				}
				else
				{
					ofs_match << line << '\t' << 1<<endl;

					ofs << line << '\t' << 1<<endl;
				}
			}
			//bool doner_plus = false, doner_minus = false, acceptor_plus = false, acceptor_minus = false;

			//for (size_t i = 0; i < doner_strands.length(); ++i)
			//{
			//	if (doner_strands[i] == '+')
			//		doner_plus = true;
			//	else if (doner_strands[i] == '-')
			//		doner_minus = true;
			//}

			//for (size_t i = 0; i < acceptor_strands.length(); ++i)
			//{
			//	if (acceptor_strands[i] == '+')
			//		acceptor_plus = true;
			//	else if (acceptor_strands[i] == '-')
			//		acceptor_minus = true;
			//}

			//bool same_gene_strand = false, diff_gene_strand = false, same_fusion_strand = false, diff_fusion_strand = false;

			//if ((doner_plus && acceptor_plus) || (doner_minus && acceptor_minus))
			//{
			//	same_gene_strand = true;
			//}

			//if ((doner_plus == strand2))
			//{
			//	same_fusion_strand = true;
			//}
			//else
			//{
			//	diff_fusion_strand = true;
			//}
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

	load_fusion_chrom_seq(argv[1], argv[2]);

	return 0;
}