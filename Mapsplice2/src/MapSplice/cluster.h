#ifndef ClUSTER_H
#define ClUSTER_H

#include <string.h>
using namespace std;

class Cluster
{
public:
	Cluster(string _chrom1, string _strand1, int _start1, int _end1, string _chrom2, string _strand2, int _start2, int _end2, int _support, bool _same_strand)
	{
		chrom1 = _chrom1;
		strand1 = _strand1;
		start1 = _start1;
		end1 = _end1;
		chrom2 = _chrom2;
		strand2 = _strand2;
		start2 = _start2;
		end2 = _end2;
		support = _support;
		same_strand = _same_strand;
	}

	Cluster(string value_line)
	{
		char chrom1_tmp[1000], chrom2_tmp[1000];
		int strand1_tmp, strand2_tmp;
		sscanf(value_line.c_str(), "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d",chrom1_tmp, &strand1_tmp, &start1, &end1, chrom2_tmp, &strand2_tmp, &start2, &end2, &support);
		strand1 = (strand1_tmp == 0 ? "+" : "-");
		strand2 = (strand2_tmp == 0 ? "+" : "-");	
		chrom1 = chrom1_tmp;
		chrom2 = chrom2_tmp;
		same_strand = (strand1 == strand2);
	}
	
	~Cluster()
	{
	}

	string chrom1;
	string strand1;
	int start1;
	int end1;
	string chrom2;
	string strand2;
	int start2;
	int end2;
	int support;
	bool same_strand;
};

#endif

