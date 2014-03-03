#include <string.h>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
using namespace std;

#pragma warning(disable:4996)

bool chrexist = false;

enum Junc_Type{normal, fusion, intergenic,overlapping};
string type_str[4]={"normal", "fusion", "intergenic", "overlapping"};
enum Junc_Source{from_normal, from_fusion};
string source_str[2]={"from_normal", "from_fusion"};

class Exon
{
public:
	int start;
	int end;

	Exon()
	{}

	~Exon()
	{}
};

class Gene
{
public:
	string gene_name;
	string chrom;
	string strand;
	int start;
	int end;
	vector<Exon> exons;

	Gene()
	{
	}

	~Gene()
	{}
};

class Fusion_Junction
{
public:
	string chrom1;
	string chrom2;
	int pos1;
	int pos2;
	int support;
	int original_support;
	string strand1;
	string strand2;
	string flank_seq;
	vector<Gene> start_support_gene;
	vector<Gene> end_support_gene;
	string line;
	Junc_Source source;

	Fusion_Junction()
	{
	}

	~Fusion_Junction()
	{}
};

inline bool sort_junc_start(const Fusion_Junction& junc1, const Fusion_Junction& junc2)
{
	if(junc1.chrom1 != junc2.chrom1)
		return junc1.chrom1 < junc2.chrom1;
	else 
		return junc1.pos1 <junc2.pos1;
}

inline bool sort_junc_end(const Fusion_Junction& junc1, const Fusion_Junction& junc2)
{
	if(junc1.chrom2 != junc2.chrom2)
		return junc1.chrom2 < junc2.chrom2;
	else 
		return junc1.pos2 <junc2.pos2;
}

inline bool sort_gene(const Gene& gene1, const Gene& gene2)
{
	if(gene1.chrom != gene2.chrom)
		return gene1.chrom < gene2.chrom;
	else if(gene1.start != gene2.start)
		return gene1.start < gene2.start;
	else
		return gene1.end < gene2.end;
}

int compare_gene_and_pos(Gene& gene, string chrom, int pos, int range)
{
	if(gene.chrom < chrom)
		return -1;
	else if(gene.chrom > chrom)
		return 1;
	else
	{
		if(gene.end + range < pos)
			return -1;
		else if(gene.start - range > pos)
			return 1;
		else 
			return 0;
	}
}


class Search_Fusion_Gene
{
public:
	vector<Fusion_Junction> junc_set;
	vector<Gene> gene_set;

	void load_junc(char* junc_file, bool header)
	{
		ifstream junc_fs(junc_file);
		string line;
		if(header)
			getline(junc_fs, line);  //skip header
		while(getline(junc_fs, line))
		{
			char chrom_check[1000];
			sscanf(line.c_str(), "%s", chrom_check);
			Fusion_Junction new_junc;
			new_junc.line=line;
			char* pos =strchr(chrom_check,'~');
			if(pos) //fusion
			{
				char chrom_tmp1[1000], chrom_tmp2[1000], strand_tmp1[1000], strand_tmp2[1000];
				sscanf(line.c_str(), "%[^~]%*[~]%s\t%d\t%d\t%*s%d%1s%1s", chrom_tmp1, chrom_tmp2, &(new_junc.pos1), &(new_junc.pos2), &(new_junc.original_support), strand_tmp1, strand_tmp2);
				string chrstr = chrom_tmp1;

                            if (chrstr.find("chr") != string::npos)
					chrexist = true;

				new_junc.chrom1=chrom_tmp1;
				new_junc.chrom2=chrom_tmp2;
				new_junc.strand1=strand_tmp1;
				new_junc.strand2=strand_tmp2;
				new_junc.source=from_fusion;
			}
			else//normal junc
			{
				char chrom_tmp1[1000], strand_tmp[1000];
				sscanf(line.c_str(), "%s\t%d\t%d\t%*s%d%1s", chrom_tmp1, &(new_junc.pos1), &(new_junc.pos2), &(new_junc.original_support), strand_tmp);
				new_junc.chrom1=chrom_tmp1;

				string chrstr = chrom_tmp1;

                            if (chrstr.find("chr") != string::npos)
					chrexist = true;

				new_junc.chrom2=chrom_tmp1;
				new_junc.strand1=strand_tmp;
				new_junc.strand2=strand_tmp;
				new_junc.source=from_normal;
			}
			junc_set.push_back(new_junc);
		}
		junc_fs.close();
		cout<<"load junc complete"<<endl;
	}

	void load_gene_gff(char* gene_file)
	{
		ifstream gene_fs(gene_file);
		string line;
		while(getline(gene_fs, line))
		{
			if(line.empty() || line.substr(0,2)=="##")
				continue;
			char chrom_tmp[1000], type_tmp[1000], name_tmp[1000];
			int start, end;
			sscanf(line.c_str(), "%s\t%*s\t%s\t%d\t%d\t%*s\t%*s\t%*s\t%*3s%[^;]", chrom_tmp, type_tmp, &start, &end, name_tmp);
			string type=type_tmp;
			if(type=="match")
			{
				Gene new_gene;
				new_gene.chrom=chrom_tmp;
				new_gene.gene_name=name_tmp;
				new_gene.start=start;
				new_gene.end=end;
				gene_set.push_back(new_gene);
			}
		}
		gene_fs.close();
	}

	void load_gene_tab(char* gene_file)
	{
		ifstream gene_fs(gene_file);
		string line;
		string curGene="";
		string curChrom;
		string curStrand;
		int GeneStart=-1;
		int GeneEnd=-1;
		while(getline(gene_fs, line))
		{
			if(line.empty())
				continue;
			char chrom_tmp[1000], name_tmp[1000];
			int start, end, strand_flag;
			sscanf(line.c_str(), "%s\t%d\t%s\t%*s\t%*s\t%d\t%d", chrom_tmp, &strand_flag, name_tmp, &start, &end);
			string newGene=name_tmp;
			string newChrom="";

			string chromtmpstr = chrom_tmp;

			if (chrexist && chromtmpstr.find("chr") == string::npos)
				newChrom = "chr";

			string newStrand = (strand_flag == 1 ? "+" : "-");
			newChrom.append(chrom_tmp);
			if(newGene==curGene)
			{
				GeneStart=min(GeneStart, start);
				GeneEnd=max(GeneEnd, end);
			}
			else
			{
				if(curGene!="")
				{
					Gene new_gene;
					new_gene.chrom=curChrom;
					new_gene.strand = curStrand;
					new_gene.gene_name=curGene;
					new_gene.start=GeneStart;
					new_gene.end=GeneEnd;
					gene_set.push_back(new_gene);
				}
				curChrom=newChrom;
				curStrand = newStrand;
				curGene=newGene;
				GeneStart=start;
				GeneEnd=end;
			}
		}
		Gene new_gene;
		new_gene.chrom=curChrom;
		new_gene.strand = curStrand;
		new_gene.gene_name=curGene;
		new_gene.start=GeneStart;
		new_gene.end=GeneEnd;
		gene_set.push_back(new_gene);
		gene_fs.close();
		cout<<"load gene complete"<<endl;
	}

	void search(int range)
	{
		sort(gene_set.begin(), gene_set.end(), sort_gene);
		sort(junc_set.begin(),junc_set.end(),sort_junc_start); //search support for junc start
		size_t check_point = 0;
		bool check_point_set;
		for(size_t i=0;i<junc_set.size();i++)
		{
			check_point_set = false;
			for(size_t j=check_point;j<gene_set.size();j++)
			{
				int relation=compare_gene_and_pos(gene_set[j], junc_set[i].chrom1, junc_set[i].pos1, range);
				if(relation == -1)
					continue;
				if(!check_point_set)
				{
					check_point=j;
					check_point_set=true;
				}
				if(relation ==0)
				{
					junc_set[i].start_support_gene.push_back(gene_set[j]);
					continue;
				}
				if(relation == 1)
					break;
			}
		}

		sort(junc_set.begin(),junc_set.end(),sort_junc_end); //search support for junc start
		check_point = 0;
		for(size_t i=0;i<junc_set.size();i++)
		{
			check_point_set = false;
			for(size_t j=check_point;j<gene_set.size();j++)
			{
				int relation=compare_gene_and_pos(gene_set[j], junc_set[i].chrom2, junc_set[i].pos2, range);
				if(relation == -1)
					continue;
				if(!check_point_set)
				{
					check_point=j;
					check_point_set=true;
				}
				if(relation ==0)
				{
					junc_set[i].end_support_gene.push_back(gene_set[j]);
					continue;
				}
				if(relation == 1)
					break;
			}
		}
	}

	bool check_overlapping_gene(vector<Gene> &gene_set1, vector<Gene> &gene_set2)
	{
		for(size_t i=0;i<gene_set1.size();i++)
		{
			for(size_t j=0;j<gene_set2.size();j++)
			{
				if(gene_set1[i].gene_name==gene_set2[j].gene_name)
					return true;
			}
		}
		return false;
	}

	void output(char* outfile)
	{
		ofstream out_fs(outfile);
		Junc_Type junc_type;
		for(size_t i=0;i<junc_set.size();i++)
		{
			out_fs<<source_str[junc_set[i].source]<<"\t";
			if(junc_set[i].start_support_gene.size()==0 || junc_set[i].end_support_gene.size()==0)
				junc_type=intergenic;
			else if(junc_set[i].start_support_gene.size()==1 && junc_set[i].end_support_gene.size()==1 && junc_set[i].start_support_gene[0].gene_name == junc_set[i].end_support_gene[0].gene_name)
				junc_type=normal;
			else if(check_overlapping_gene(junc_set[i].start_support_gene, junc_set[i].end_support_gene))
				junc_type=overlapping;
			else
				junc_type=fusion;
			out_fs<<type_str[junc_type]<<"\t";
			if(junc_type==fusion)
			{
				for(size_t j=0;j<junc_set[i].start_support_gene.size();j++)
				{
					out_fs<<junc_set[i].start_support_gene[j].strand;
				}
				out_fs << ",";
				for(size_t j=0;j<junc_set[i].end_support_gene.size();j++)
				{
					out_fs<<junc_set[i].end_support_gene[j].strand;
				}
			}
			out_fs << "\t";
			if(junc_set[i].start_support_gene.size()==0)
			{
				out_fs<<"-"<<",";
			}
			else
			{
				for(size_t j=0;j<junc_set[i].start_support_gene.size();j++)
				{
					out_fs<<junc_set[i].start_support_gene[j].gene_name<<",";
				}
			}
			out_fs<<"\t";
			if(junc_set[i].end_support_gene.size()==0)
			{
				out_fs<<"-"<<",";
			}
			else
			{
				for(size_t j=0;j<junc_set[i].end_support_gene.size();j++)
				{
					out_fs<<junc_set[i].end_support_gene[j].gene_name<<",";
				}
			}
			out_fs<<"\t"<<junc_set[i].line<<endl;
		}
	}
};


void print_usage()
{
	fprintf(stderr,"-g:          gene annotation file name\n");
	fprintf(stderr,"-f:          fusion junction file name\n");
	fprintf(stderr,"-f_header:   skip header of fusion junction file\n");
	fprintf(stderr,"-n:          normal junction file name\n");
	fprintf(stderr,"-n_header:   skip header of normal junction file\n");
	fprintf(stderr,"-o:          output file name\n");
	fprintf(stderr,"-r:          comparison range, default [0]\n");
	exit(1);
}

int main(int argc, char* argv[])
{
	if(argc==1)
		print_usage();
	char* gene_file=NULL;
	char* fusion_junc_file=NULL;
	char* normal_junc_file=NULL;
	char* output_file=NULL;
	bool fusion_junc_header=false;
	bool normal_junc_header=false;
	int range=0;
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"-g")==0)
			gene_file=argv[++i];
		if(strcmp(argv[i],"-f")==0)
			fusion_junc_file=argv[++i];
		if(strcmp(argv[i],"-n")==0)
			normal_junc_file=argv[++i];
		if(strcmp(argv[i],"-o")==0)
			output_file=argv[++i];
		if(strcmp(argv[i],"-r")==0)
			range=atoi(argv[++i]);
		if(strcmp(argv[i],"-f_header")==0)
			fusion_junc_header=true;
		if(strcmp(argv[i],"-n_header")==0)
			normal_junc_header=true;
	}	
	Search_Fusion_Gene sfg;

	if(fusion_junc_file!=NULL)
		sfg.load_junc(fusion_junc_file, fusion_junc_header);
	if(normal_junc_file!=NULL)
		sfg.load_junc(normal_junc_file, normal_junc_header);
	sfg.load_gene_tab(gene_file);

	sfg.search(range);
	sfg.output(output_file);
	exit(0);


	/*Search_Fusion_Gene sfg;
	sfg.load_gene_tab("D://test//human_hg19_transcript.tab");
	sfg.load_junc("D://test//matched_with_dlx_all_normal_filtered.matched_junc.add_chrom_seq.filtered_repeats", false);
	sfg.search(0);
	sfg.output("D://test//test.txt");*/




}