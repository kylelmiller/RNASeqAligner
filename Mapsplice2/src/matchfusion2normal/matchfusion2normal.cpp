#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
using namespace std;

bool paircomp ( const pair<int, int>& lhs, const pair<int, int>& rhs)
{
	return lhs.first < rhs.first;
}

int
readjunc(const char* true_juncfile, vector<string>& junc_vec, map<string, map<size_t, size_t> >& start_map, map<string, map<size_t, size_t> >& end_map)
{
	ifstream ifs(true_juncfile);

	if (ifs.is_open())
	{
		size_t count = 0;
		while (!ifs.eof() )
		{
			string line;
			getline(ifs, line);

			if (line == "")
				continue;

			char chrom[100];
			size_t start, end;

			sscanf(line.c_str(), "%s\t%lu\t%lu", chrom, &start, &end);

			junc_vec.push_back(line);

			start_map[chrom][start] = junc_vec.size();

			end_map[chrom][end] = junc_vec.size();			
		}

		ifs.close();
	}
	return 0;
}

size_t abs_diff(size_t lhs, size_t rhs)
{
	if (lhs <= rhs)
		return rhs - lhs;
	else 
		return lhs - rhs;
}

int
MatchStartEnd(const char* verify_junc, const char* out_root, vector<string>& junc_vec, map<string, map<size_t, size_t> >& start_map, map<string, map<size_t, size_t> >& end_map, size_t range)
{
	ifstream ifs(verify_junc);

	//string out_matched = out_root; out_matched.append(".matched"); ofstream ofs_matched(out_matched.c_str());

	string out_matched_junc = out_root; out_matched_junc.append(".matched_junc"); ofstream ofs_matched_junc(out_matched_junc.c_str());

	//string out_rangematched = out_root; out_rangematched.append(".range_matched"); ofstream ofs_rangematched(out_rangematched.c_str());

	//string out_rangematched_junc = out_root; out_rangematched_junc.append(".range_matched_junc"); ofstream ofs_rangematched_junc(out_rangematched_junc.c_str());

	//string out_sidematched = out_root; out_sidematched.append(".sidematched"); ofstream ofs_sidematched(out_sidematched.c_str());

	//string out_sidematched_junc = out_root; out_sidematched_junc.append(".sidematched_junc"); ofstream ofs_sidematched_junc(out_sidematched_junc.c_str());
 //
	//string out_rangesidematched = out_root; out_rangesidematched.append(".range_sidematched"); ofstream ofs_rangesidematched(out_rangesidematched.c_str());

	//string out_rangesidematched_junc = out_root; out_rangesidematched_junc.append(".range_sidematched_junc"); ofstream ofs_rangesidematched_junc(out_rangesidematched_junc.c_str());

	//string out_notmatched = out_root; out_notmatched.append(".notmatched"); ofstream ofs_notmatched(out_notmatched.c_str());

	int count_match = 0, count_sidematch = 0, count_rangematch = 0, count_rangesidematch = 0, count_notmatch = 0, count_total = 0;

	if (ifs.is_open())
	{
		//string skipline;

		//getline(ifs, skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs, line);

			if (line == "" || line.find("JUNC") == string::npos)
				continue;

			char chrom1[100], chrom2[100];
			size_t start, end;

			char juncname[100], strand[10], rgb[100], blocks[1000], blocksoffset[1000], rank[100],/* il[100], */flankchr[100], synseq[1000];
			int /*juncst, juncend,*/ kinds, hits, flankcase;

//			double lpq;

			unsigned short min_mis, max_mis;

			double ave_mis;

			//chr1_chr1	569482	38077348	JUNC_8733	94	++	255,0,0	2	26,58	0,37507866	2.03725	6	CTAC	
			//TTTCGCTCTAAGATTAAAAATGCCCTAGCCCACTTCTTACCACAAGGCACACCTACACCCCTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTG	1	3	1.85106

			sscanf(line.c_str(), "%[^~]~%s\t%llu\t%llu\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%hu\t%hu\t%lf", chrom1, chrom2, &start, &end, 
				juncname, &hits, strand, rgb, &kinds, blocks, blocksoffset, rank, &flankcase, flankchr, synseq, &min_mis, &max_mis, &ave_mis);

			size_t start_junc_id = -1, end_junc_id = -1, cur_start_junc_start_id = -1, cur_start_junc_end_id = -1, cur_end_junc_start_id = -1, cur_end_junc_end_id = -1;

			size_t cur_start_junc_start = 0, cur_start_junc_end = 0, cur_end_junc_start = 0, cur_end_junc_end = 0;

			if (true)
			{
				map<size_t, size_t>::const_iterator cur_start_junc_start_RcIter, cur_start_junc_end_RcIter, cur_end_junc_start_RcIter, cur_end_junc_end_RcIter;

				//find start
				if (start_map.find(chrom1) != start_map.end())
				{
					//cur start junc start
					cur_start_junc_start_RcIter = start_map[chrom1].lower_bound(start);

					if (cur_start_junc_start_RcIter !=  start_map[chrom1].end())
					{
						if (cur_start_junc_start_RcIter->first == start)
						{
							cur_start_junc_start = 1;

							cur_start_junc_start_id = cur_start_junc_start_RcIter->second;
						}
						else
						{
							if (cur_start_junc_start_RcIter != start_map[chrom1].begin())
								--cur_start_junc_start_RcIter;

							while (start + range >= cur_start_junc_start_RcIter->first && cur_start_junc_start_RcIter != start_map[chrom1].end())
							{
								if (abs_diff(cur_start_junc_start_RcIter->first, start) <= range)
								{
									cur_start_junc_start = 2;

									cur_start_junc_start_id = cur_start_junc_start_RcIter->second;

									break;
								}

								++cur_start_junc_start_RcIter;
							}
						}
					}

					//cur start junc end
					cur_start_junc_end_RcIter = end_map[chrom1].lower_bound(start);

					if (cur_start_junc_end_RcIter != end_map[chrom1].end())
					{
						if (cur_start_junc_end_RcIter->first == start)
						{
							cur_start_junc_end = 1;

							cur_start_junc_end_id = cur_start_junc_end_RcIter->second;
						}
						else
						{
							if (cur_start_junc_end_RcIter != end_map[chrom1].begin())
								--cur_start_junc_end_RcIter;

							while (start + range >= cur_start_junc_end_RcIter->first && cur_start_junc_end_RcIter != end_map[chrom1].end())
							{
								if (abs_diff(cur_start_junc_end_RcIter->first, start) <= range)
								{
									cur_start_junc_end = 2;

									cur_start_junc_end_id = cur_start_junc_end_RcIter->second;

									break;
								}

								++cur_start_junc_end_RcIter;
							}
						}
					}

				}

				//find end

				if (end_map.find(chrom2) != end_map.end())
				{
					//cur end junc start
					cur_end_junc_start_RcIter = start_map[chrom2].lower_bound(end);

					if (cur_end_junc_start_RcIter != start_map[chrom2].end())
					{
						if (cur_end_junc_start_RcIter->first == end)
						{
							cur_end_junc_start = 1;

							cur_end_junc_start_id = cur_end_junc_start_RcIter->second;
						}
						else
						{
							if (cur_end_junc_start_RcIter != start_map[chrom2].begin())
								--cur_end_junc_start_RcIter;

							while (end + range >= cur_end_junc_start_RcIter->first)
							{
								if (abs_diff(cur_end_junc_start_RcIter->first, end) <= range && cur_end_junc_start_RcIter != start_map[chrom2].end())
								{
									cur_end_junc_start = 2;

									cur_end_junc_start_id = cur_end_junc_start_RcIter->second;

									break;
								}

								++cur_end_junc_start_RcIter;
							}
						}
					}

					//cur end junc end
					cur_end_junc_end_RcIter = end_map[chrom2].lower_bound(end);

					if (cur_end_junc_end_RcIter != end_map[chrom2].end())
					{
						if (cur_end_junc_end_RcIter->first == end)
						{
							cur_end_junc_end = 1;

							cur_end_junc_end_id = cur_end_junc_end_RcIter->second;
						}
						else
						{
							if (cur_end_junc_end_RcIter != end_map[chrom2].begin())
								--cur_end_junc_end_RcIter;

							while (end + range >= cur_end_junc_end_RcIter->first)
							{
								if (abs_diff(cur_end_junc_end_RcIter->first, end) <= range && cur_end_junc_end_RcIter != end_map[chrom2].end())
								{
									cur_end_junc_end = 2;

									cur_end_junc_end_id = cur_end_junc_end_RcIter->second;

									break;
								}

								++cur_end_junc_end_RcIter;
							}
						}
					}
				}

				count_total++;

				ofs_matched_junc << line;

				if (cur_start_junc_start == 1 || cur_start_junc_end == 1)
				{
					ofs_matched_junc << "\tdoner_exact_matched";
				}
				else if (cur_start_junc_start >= 1 || cur_start_junc_end >= 1)
				{
					ofs_matched_junc << "\tdoner_range_matched";
				}
				else
					ofs_matched_junc << "\tnot_matched";

				if (cur_end_junc_start == 1 || cur_end_junc_end == 1)
				{
					ofs_matched_junc << "\tacceptor_exact_matched";
				}
				else if (cur_end_junc_start >= 1 || cur_end_junc_end >= 1)
				{
					ofs_matched_junc << "\tacceptor_range_matched";
				}
				else
					ofs_matched_junc << "\tnot_matched";

				ofs_matched_junc << endl;

				if ((cur_start_junc_start == 1 || cur_start_junc_end == 1) && (cur_end_junc_start == 1 || cur_end_junc_end == 1) )
				{
					//start_junc_id = cur_start_junc_start == 1 ? cur_start_junc_start_id : cur_start_junc_end_id;

					//end_junc_id = cur_end_junc_start == 1 ? cur_end_junc_start_id : cur_end_junc_end_id;

					//ofs_matched << line << endl <<junc_vec[start_junc_id - 1]<<endl<< junc_vec[end_junc_id - 1] << endl;

					count_match++;
				}
				else if ((cur_start_junc_start >= 1 || cur_start_junc_end >= 1) && (cur_end_junc_start >= 1 || cur_end_junc_end >= 1) )
				{
					//start_junc_id = cur_start_junc_start >= 1 ? cur_start_junc_start_id : cur_start_junc_end_id;

					//end_junc_id = cur_end_junc_start >= 1 ? cur_end_junc_start_id : cur_end_junc_end_id;

					//ofs_rangematched << line << endl <<junc_vec[start_junc_id - 1]<<endl<< junc_vec[end_junc_id - 1] << endl;

					//ofs_matched_junc << line << endl;

					count_rangematch++;
				}
				else if ((cur_start_junc_start == 1 || cur_start_junc_end == 1) || (cur_end_junc_start == 1 || cur_end_junc_end == 1) )
				{
					if (cur_start_junc_start == 1 || cur_start_junc_end == 1)
					{
						//start_junc_id = cur_start_junc_start == 1 ? cur_start_junc_start_id : cur_start_junc_end_id;

						//ofs_sidematched << line << endl <<junc_vec[start_junc_id - 1]<<endl;

						//ofs_sidematched_junc << line << endl;
					}
					else
					{
						//end_junc_id = cur_end_junc_start == 1 ? cur_end_junc_start_id : cur_end_junc_end_id;

						//ofs_sidematched << line << endl <<junc_vec[end_junc_id - 1]<<endl;

						//ofs_sidematched_junc << line << endl;
					}

					count_sidematch++;
				}
				else if ((cur_start_junc_start >= 1 || cur_start_junc_end >= 1) || (cur_end_junc_start >= 1 || cur_end_junc_end >= 1) )
				{
					if (cur_start_junc_start >= 1 || cur_start_junc_end >= 1)
					{
						//start_junc_id = cur_start_junc_start >= 1 ? cur_start_junc_start_id : cur_start_junc_end_id;
						//
						//ofs_rangesidematched << line << endl << junc_vec[start_junc_id - 1] << endl;

						//ofs_rangesidematched_junc << line << endl;
					}
					else
					{
						//end_junc_id = cur_end_junc_start >= 1 ? cur_end_junc_start_id : cur_end_junc_end_id;
						//
						//ofs_rangesidematched << line << endl << junc_vec[end_junc_id - 1] << endl;

						//ofs_rangesidematched_junc << line << endl;
					}

					count_rangesidematch++;
				}
				else
				{
					//ofs_notmatched << line << endl; 
					count_notmatch++;
				}
			}
		}

	}

	cout << "count_match:\t"<<count_match << endl;
	cout << "count_rangematch:\t"<<count_rangematch<< endl;
	cout << "count_sidematch:\t"<<count_sidematch << endl;
	cout << "count_rangesidematch:\t"<<count_rangesidematch << endl;
	cout << "count_notmatch:\t"<<count_notmatch << endl;
	cout << "count_total:\t"<<count_total << endl;

	return 0;
}

int
main(int argc, const char** argv)
{
	if (argc < 4)
	{
		cout << " need three files\n";
		return 0;
	}

	const char* true_juncfile = argv[1];
	const char* verify_junc = argv[2];
	const char* out_root = argv[3];
	int range = atoi(argv[4]);

	vector<string> junc_vec;
	
	map<string, map<size_t, size_t> > start_map;
	
	map<string, map<size_t, size_t> > end_map;

	readjunc(true_juncfile, junc_vec, start_map, end_map);

	MatchStartEnd(verify_junc, out_root, junc_vec, start_map, end_map, range);
}
