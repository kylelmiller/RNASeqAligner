#include <iostream>
#include <vector>
#include "SamRec.h"
#include <string>

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
//#include <dirent.h>
#include <iomanip>
#include <map>
#include <queue>
#include <list>
#include <iterator>

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>

using namespace std;

size_t convert_to_abs_offset(vector<pair<size_t, int> >& spliceway_vec, size_t st, size_t end, char strand)
{
	if (strand == '-')
	{
		size_t length = end - st;

		for (size_t i = 0; i < spliceway_vec.size(); ++i)
		{
			size_t temp = length - (spliceway_vec[i].first + spliceway_vec[i].second - 1) + 1;

			spliceway_vec[i].first = temp + st;
		}

		reverse(spliceway_vec.begin(), spliceway_vec.end());
	}
	else
	{
		for (size_t i = 0; i < spliceway_vec.size(); ++i)
			spliceway_vec[i].first += st;
	}

	return spliceway_vec.size();
}

int main(int argc, char* argv[])
{
	ifstream ifs(argv[1]);

	ofstream ofs(argv[2]);

	ofstream ofs_fusion(argv[3]);

	//size_t combined_length = atoi(argv[5]);

	SamRec samrec;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "")
				continue;

			if (line.find("ZF:Z:FUS") != string::npos)
				continue;

			samrec.Clear();

			samrec.Set(line, 3);

			if (samrec.issmallins)
			{
				cout << line << endl;
				continue;
			}

			char chrom1[1000], chrom2[1000];

			char strand1, strand2;

			size_t start1, end1, start2, end2;

			sscanf(samrec.chrom_name.c_str(), "%[^=]=%c=%llu=%llu=%[^=]=%c=%llu=%llu", chrom1, &strand1, &start1, &end1, chrom2, &strand2, &start2, &end2);

			string chrom1str = chrom1, chrom2str = chrom2;

			size_t length1 = end1 - start1;

			size_t length2 = end2 - start2;

			size_t end_pos = samrec.end;

			size_t prefixst = samrec.start;

			vector<pair<size_t, int> > spliceway_vec_p1, spliceway_vec_p2;
			
			vector<pair<size_t, int> >& spliceway_vec = samrec.spliceway_vec;

			if (end_pos <= length1)
			{
				spliceway_vec_p1 = spliceway_vec;
			}
			else if (prefixst > length1)
			{
				spliceway_vec_p2 = spliceway_vec;

				for (size_t i = 0; i < spliceway_vec_p2.size(); ++i)
				{
					spliceway_vec_p2[i].first = spliceway_vec_p2[i].first - length1;
				}
			}
			else if (prefixst <= length1 && end_pos > length1)
			{
				size_t i;

				for (i = 0; i < spliceway_vec.size(); ++i)
				{
					if (spliceway_vec[i].first + spliceway_vec[i].second - 1 <= length1)
						spliceway_vec_p1.push_back(make_pair(spliceway_vec[i].first, spliceway_vec[i].second));
					else
						break;
				}

				for (; i < spliceway_vec.size(); ++i)
				{
					if (spliceway_vec[i].first <= length1 && spliceway_vec[i].first + spliceway_vec[i].second - 1 > length1)
					{
						spliceway_vec_p1.push_back(make_pair(spliceway_vec[i].first, length1 - spliceway_vec[i].first + 1));

						spliceway_vec_p2.push_back(make_pair(1, spliceway_vec[i].second - (length1 - spliceway_vec[i].first + 1)));
					}
					else
						break;
				}

				for (; i < spliceway_vec.size(); ++i)
				{
					if (spliceway_vec[i].first > length1)
						spliceway_vec_p2.push_back(make_pair(spliceway_vec[i].first - length1, spliceway_vec[i].second));
					else
						break;
				}
			}

			if (spliceway_vec_p1.size())
				convert_to_abs_offset(spliceway_vec_p1, start1, end1, strand1);

			if (spliceway_vec_p2.size())
				convert_to_abs_offset(spliceway_vec_p2, start2, end2, strand2);

			size_t strand_t = samrec.strand_t;

			SamRec samrec1(samrec), samrec2(samrec);

			if (spliceway_vec_p1.size() && !spliceway_vec_p2.size())
			{
				if (strand1 == '-')
				{
					if (strand_t == 0)
						strand_t = 16;
					else
						strand_t = 0;
				}

				samrec1.strand_t = strand_t;

				samrec1.spliceway_vec = spliceway_vec_p1;

				samrec1.start = spliceway_vec_p1.front().first;

				samrec1.chrom_name = chrom1;

				vector<pair<size_t, int> >::iterator spliceway_iter;

				string newjunccode;

				for (spliceway_iter = spliceway_vec_p1.begin(); spliceway_iter != spliceway_vec_p1.end() - 1; ++spliceway_iter)
				{
					char buf[1000];

					sprintf(buf, "%dM%lluN", spliceway_iter->second, (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second);

					newjunccode.append(buf);
				}

				char buf[1000];

				sprintf(buf, "%dM", spliceway_iter->second);

				newjunccode.append(buf);

				samrec1.ori_splice_way = newjunccode;

				if (strand1 == '-')
				{
					samrec1.mapped_seq = revcomp(samrec1.mapped_seq);

					reverse(samrec1.qual_str.begin(), samrec1.qual_str.end());
				}

				ofs<<samrec1.tostring(0, 0)<<endl;
			}
			else if (!spliceway_vec_p1.size() && spliceway_vec_p2.size())
			{
				if (strand2 == '-')
				{
					if (strand_t == 0)
						strand_t = 16;
					else
						strand_t = 0;
				}

				samrec2.strand_t = strand_t;

				samrec2.spliceway_vec = spliceway_vec_p2;

				samrec2.start = spliceway_vec_p2.front().first;

				samrec2.chrom_name = chrom2;

				vector<pair<size_t, int> >::iterator spliceway_iter;

				string newjunccode;

				for (spliceway_iter = spliceway_vec_p2.begin(); spliceway_iter != spliceway_vec_p2.end() - 1; ++spliceway_iter)
				{
					char buf[1000];

					sprintf(buf, "%dM%lluN", spliceway_iter->second, (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second);

					newjunccode.append(buf);
				}
				
				char buf[1000];

				sprintf(buf, "%dM", spliceway_iter->second);

				newjunccode.append(buf);

				samrec2.ori_splice_way = newjunccode;

				if (strand2 == '-')
				{
					samrec2.mapped_seq = revcomp(samrec2.mapped_seq);

					reverse(samrec2.qual_str.begin(), samrec2.qual_str.end());
				}

				ofs<< samrec2.tostring(0, 0)<<endl;
			}
			else if (spliceway_vec_p1.size() && spliceway_vec_p2.size())
			{
				//1
				size_t strand_t1 = strand_t, strand_t2 = strand_t;

				char std1, std2;

				if (strand1 == '-')
				{
					if ((strand_t1 & IS_REVERSE) == 0)
						strand_t1 = 16;
					else
						strand_t1 = 0;
				}

				if (strand2 == '-')
				{
					if ((strand_t2 & IS_REVERSE) == 0)
						strand_t2 = 16;
					else
						strand_t2 = 0;
				}

				if (strand_t1 & IS_REVERSE)
				{
					std1 = '-';
				}
				else
				{
					std1 = '+';
				}

				if (strand_t2 & IS_REVERSE)
				{
					std2 = '-';
				}
				else
				{
					std2 = '+';
				}

				string mapped_seq = samrec1.mapped_seq;

				string mapped_seq1, mapped_seq2;

				string qual_str = samrec1.qual_str;

				string qual_str1, qual_str2;

				{
					samrec1.strand_t = strand_t1;

					samrec1.spliceway_vec = spliceway_vec_p1;

					samrec1.start = spliceway_vec_p1.front().first;

					samrec1.chrom_name = chrom1;

					vector<pair<size_t, int> >::iterator spliceway_iter;

					string newjunccode;

					int mappedlen = 0;

					for (spliceway_iter = spliceway_vec_p1.begin(); spliceway_iter != spliceway_vec_p1.end() - 1; ++spliceway_iter)
					{
						char buf[1000];

						sprintf(buf, "%dM%lluN", spliceway_iter->second, (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second);

						mappedlen += spliceway_iter->second;

						newjunccode.append(buf);
					}

					char buf[1000];

					sprintf(buf, "%dM", spliceway_iter->second);

					mappedlen += spliceway_iter->second;

					newjunccode.append(buf);

					char clip[1000];

					sprintf(clip, "%lluS", samrec1.mapped_seq.length() - mappedlen);

					if (strand_t & IS_REVERSE)
						newjunccode = clip + newjunccode;
					else
						newjunccode.append(clip);

					samrec1.ori_splice_way = newjunccode;

					mapped_seq1 = mapped_seq.substr(0, mappedlen);

					qual_str1 = qual_str.substr(0, mappedlen);

					if (strand1 == '-')
					{						
						mapped_seq1 = revcomp(mapped_seq1);

						reverse(qual_str1.begin(), qual_str1.end());
					}
				}

				{
					samrec2.strand_t = strand_t2;

					samrec2.spliceway_vec = spliceway_vec_p2;

					samrec2.start = spliceway_vec_p2.front().first;

					samrec2.chrom_name = chrom2;

					vector<pair<size_t, int> >::iterator spliceway_iter;

					string newjunccode;

					int mappedlen = 0;

					for (spliceway_iter = spliceway_vec_p2.begin(); spliceway_iter != spliceway_vec_p2.end() - 1; ++spliceway_iter)
					{
						char buf[1000];

						sprintf(buf, "%dM%lluN", spliceway_iter->second, (spliceway_iter + 1)->first - spliceway_iter->first - spliceway_iter->second);

						mappedlen += spliceway_iter->second;

						newjunccode.append(buf);
					}

					char buf[1000];

					sprintf(buf, "%dM", spliceway_iter->second);

					mappedlen += spliceway_iter->second;

					newjunccode.append(buf);

					char clip[1000];

					sprintf(clip, "%lluS", samrec1.mapped_seq.length() - mappedlen);

					if (strand_t & IS_REVERSE)
						newjunccode.append(clip);
					else
						newjunccode = clip + newjunccode;

					samrec2.ori_splice_way = newjunccode;

					mapped_seq2 = mapped_seq.substr(mapped_seq.length() - mappedlen, mappedlen);

					qual_str2 = qual_str.substr(mapped_seq.length() - mappedlen, mappedlen);

					if (strand2 == '-')
					{						
						mapped_seq2 = revcomp(mapped_seq2);

						reverse(qual_str2.begin(), qual_str2.end());
					}					
				}				

				if (strand_t & IS_REVERSE)
				{
					samrec1.qual_str = qual_str2 + qual_str1;

					samrec1.mapped_seq = mapped_seq2 + mapped_seq1;

					samrec2.qual_str = qual_str2 + qual_str1;

					samrec2.mapped_seq = mapped_seq2 + mapped_seq1;

					ofs_fusion<< samrec2.tostring(0, 0)<<"\tZF:Z:FUS"<<endl;
					ofs<< samrec2.tostring(0, 0)<<"\tZF:Z:FUS"<<endl;

					ofs_fusion<<samrec1.tostring(0, 0)<< "\tZF:Z:FUS"<< endl;
					ofs<<samrec1.tostring(0, 0)<<"\tZF:Z:FUS"<<endl;						
				}
				else
				{
					samrec1.qual_str = qual_str1 + qual_str2;

					samrec1.mapped_seq = mapped_seq1 + mapped_seq2;

					samrec2.qual_str = qual_str1 + qual_str2;

					samrec2.mapped_seq = mapped_seq1 + mapped_seq2;

					ofs_fusion<<samrec1.tostring(0, 0)<< "\tZF:Z:FUS"<< endl;
					ofs<<samrec1.tostring(0, 0)<<"\tZF:Z:FUS"<<endl;

					ofs_fusion<< samrec2.tostring(0, 0)<<"\tZF:Z:FUS"<<endl;
					ofs<< samrec2.tostring(0, 0)<<"\tZF:Z:FUS"<<endl;
				}
			}
		}

		ifs.close();
	}
}