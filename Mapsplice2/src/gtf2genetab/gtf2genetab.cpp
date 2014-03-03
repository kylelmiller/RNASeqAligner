#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string.h>
#include <sstream>
#include <stdlib.h>
using namespace std;

//#define _DEBUG

void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false)
{
	stringstream my_ss(line);
	string item;
	while (getline(my_ss, item, delim))
	{
#ifdef _DEBUG
cerr << "parsed_item: " << item << endl;
#endif		
		if (ignore_empty_item && item.empty())
			continue;    // use to skip empty item
		parsed_item.push_back(item);
	}
}

void
readchrom(const char* filename, string& longseq)
{
	cout << " read chrom: " << filename << endl;
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
	else cout << "Unable to open file";

	cout <<"chrom size:"<< longseq.size() << endl;
}


char
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

string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

struct exon {

	int start, end;

	string chr;

	string gene;

	string transcript;

	int strand;

	exon() {}

	exon(int _start, int _end, string _chr, string _gene, string _transcript, int _strand) : start(_start), end(_end), chr(_chr), gene(_gene), transcript(_transcript), strand(_strand)
	{
	}

	string tostring(int id)
	{
		char buf[10000];

		sprintf(buf, "%s\t%d\t%s\t%s\tE%d\t%d\t%d", chr.c_str(), strand, gene.c_str(), transcript.c_str(), id, start, end);

		return buf;
	}
};

bool comp_exon(const exon& lhs, const exon& rhs)
{
	if (lhs.chr == rhs.chr)
	{
		if (lhs.start == rhs.start)
		{
			return lhs.end < rhs.end;
		}
		else
		{
			return lhs.start < rhs.start;
		}
	}
	else
	{
		return lhs.chr < rhs.chr;
	}
}

struct transcript {

	string transcript_id;

	char transcript_strand;

	map<int, exon> exons;

	vector<pair<int, int> > exon_region;

	transcript (string& _transcript_id, char _transcript_strand) : transcript_id(_transcript_id), transcript_strand(_transcript_strand)
	{
	}
};

struct gene {

	string chrname;

	string geneid;

	string protein_type;

	char gene_strand;

	map<string, transcript> transcripts;


	gene (string _chrname, string _geneid, string _protein_type, char _gene_strand) : chrname(_chrname), geneid(_geneid), protein_type(_protein_type), gene_strand(_gene_strand)
	{
	}
};


void read_gtf_file(string gtf_file, string genetabout)
{
	ifstream ifs_gtf(gtf_file.c_str());

	ofstream ofs(genetabout.c_str());

	vector<exon> exons;

	if (ifs_gtf.is_open())
	{
		while(!ifs_gtf.eof())
		{
			string line;

			getline(ifs_gtf, line);

			if (line.empty() || line[0] == '#')
				continue;

			vector<string> parsed_item;

			split(line, '\t', parsed_item);

			string chr = parsed_item[0];

			string annotation_source = parsed_item[1];

			string feature_type = parsed_item[2];

			int start = atoi(parsed_item[3].c_str());

			int end = atoi(parsed_item[4].c_str());

			char strand = parsed_item[6][0];

			string gene_id, transcript_id, gene_type, gene_status, gene_name, transcript_type, transcript_status, transcript_name;

			vector<string> parsed_optional_item;
			
			split(parsed_item[8], ' ', parsed_optional_item);
			
			
			for(size_t i = 0; i < parsed_optional_item.size(); i++)
			{
				if(parsed_optional_item[i] == "gene_id")
					gene_id = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);
				else if(parsed_optional_item[i] == "transcript_id")
					transcript_id = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);	
				else if(parsed_optional_item[i] == "gene_type" || parsed_optional_item[i] == "gene_biotype")
					gene_type = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);		
				else if(parsed_optional_item[i] == "gene_status" || parsed_optional_item[i] == "gene_biotype")
					gene_status = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);				
				else if(parsed_optional_item[i] == "gene_name")
					gene_name = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);	
				else if(parsed_optional_item[i] == "transcript_type")
					transcript_type = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);	
				else if(parsed_optional_item[i] == "transcript_status")
					transcript_status = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);				
				else if(parsed_optional_item[i] == "transcript_name")
					transcript_name = parsed_optional_item[i + 1].substr(1, parsed_optional_item[i + 1].length() - 3);		
			}

			if(feature_type != "exon")
				continue;
			if ((gene_type != "protein_coding" && 
				gene_type != "processed_transcript" && 
				gene_type != "IG_C_gene" && 
				gene_type != "IG_D_gene" && 
				gene_type != "IG_J_gene" && 
				gene_type != "IG_V_gene"))
				continue;

			int strand_int = (strand == '+' ? 1 : 0);
			exons.push_back(exon(start, end, chr, gene_name, transcript_name, strand_int));
		}
	}
	else
	{
		cout <<"can't open:" << gtf_file<<endl;
	}

	cout << "sort"<<endl;

	sort(exons.begin(), exons.end(), comp_exon);

	vector<exon>::iterator exon_iter;

	int id = 0;

	cout << "write"<<endl;

	for (exon_iter = exons.begin(); exon_iter != exons.end(); ++exon_iter)
	{
		ofs << exon_iter->tostring(++id) << endl;
	}
}

int main(int argc, char** argv)
{
	cout << "gtf genetab_out"<<endl;

	read_gtf_file(argv[1], argv[2]);

}





