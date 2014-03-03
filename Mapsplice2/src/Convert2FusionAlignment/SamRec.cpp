#include "SamRec.h"
//
SamRec::SamRec(const string& line, int min_ins) : matched_id(-1),  wrong_format(false), tagidx(-1), filter_score(0), best(0), 
	isunique(false), isexonic(false), isspliced(false), ave_intron_len(0), ave_junc_mis(0), canon_count(0), 
	noncanon_count(0), canon_rate(0), issmallins(false), issmalldel(false), iscanonical(false), 
	issemicanonical(false), isnoncanoical(false), isunmapped(false), isclipped(false), cur_line(line),
	/*is_insert(false),*/ mappedlen(0), intron_size(0), min_anchor(-1), max_anchor(0), pair_rate(0), m_contig_len(0), 
	fusion_prefix_len(0), fusion_suffix_len(0), is_swapped(false), mate_offset(0), mate_offset2(0), mate_diff(0), 
	mate_diff2(0), mate_match("*"), mate_match2("*"), paired_type(SINGLE), start(0), start2(0), strand_t(0), strand_t2(0),
	forward_count(0), reverse_count(0), fusion_prefix_paired(false), fusion_suffix_paired(false)
{
	if (line.find("seq.10001997/") != string::npos)
	{
		cout << line << endl;
	}

	char tagname[1000], chrom[100], mapped[100], seq[1000], qual_chr[1000], alters_chr[1000];

	char chrom2[100], mapped2[100];

	char mate_match_chr[1000];

	string alterstr = "";

	char strand[100];

	sscanf(line.c_str(), "%s\t%s\t%s", tagname, chrom, strand);

	if (string(strand) == "+" || string(strand) == "-")
		is_fusion = true;
	else
		is_fusion = false;

	if (is_fusion)
	{
		//14~R_7/2	chr1	+	11184637	54M	chr11	+	61646784	21M	GGCGCAGATCTTCATGGCCTGTTAGAAGGAAAACACACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGGATGAA	
		//283348985863448157653149274159715381432428925324118868345815774332128375427~CTGCATCACACGCTCATCCTGGCGCAGATCTTCATGGCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTG
		//CCTGAGAGATGGCCAGGATGAAGGCGGCCAGGGCACTGGGCACCCAGCCAGGACCCAGGAGGTAGATAAGGAGCC	NM:i:2

		size_t read_count = sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tagname, chrom, &strand1, &start, mapped, 
			chrom2, &strand2, &start2, mapped2, seq, qual_chr, &mis_match);

		sscanf(tagname, "%llu", &tagidx);

		mapped_seq = seq;

		qual_str = qual_chr; 

		syn_seq = qual_chr;

		size_t last_idx = syn_seq.find_last_of("~");

		qual_str = qual_str.substr(0, last_idx);

		syn_seq = syn_seq.substr(last_idx + 1, syn_seq.length() - last_idx - 1);

		tag_name = tagname;

		if (strand1 == '+')
			strand_t = 0;
		else
			strand_t = 16;

		if (strand2 == '+')
			strand_t2 = 0;
		else
			strand_t2 = 16;

		if (string(chrom) < string(chrom2) || (string(chrom) == string(chrom2) && start <= start2))
		{
			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;
			
			chrom_name = chrom;
			splice_way = mapped;
			ori_splice_way = mapped;

			splice_way2 = mapped2;
			ori_splice_way2 = mapped2;
			chrom_name2 = chrom2;
		}
		else
		{
			if (strand1 == '+')
				strand1 = '-';
			else
				strand1 = '+';

			if (strand2 == '+')
				strand2 = '-';
			else
				strand2 = '+';

			swap(strand1, strand2);

			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;

			
			swap(start, start2);
			
			chrom_name = chrom2;//chrom;
			splice_way = mapped2;//mapped;
			ori_splice_way = mapped2;//mapped;

			splice_way2 = mapped;//mapped2;
			ori_splice_way2 = mapped;//mapped2;
			chrom_name2 = chrom;//chrom2;

			is_swapped = true;
		}

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}
		}	

		size_t index = 0;

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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		mappedlen1 = mappedlen;

		fusion_prefix_len = mappedlen1;

		//segment 2	

		if (mapped2[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way2.find("I") != string::npos))
			issmallins = true;
		if ((splice_way2.find("D") != string::npos))
			issmalldel = true;
		if (splice_way2.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way2.length(); ++i)
			if (splice_way2[i] == 'D')
				splice_way2[i] = 'N';

		if ((splice_way2.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way2[splice_way2.length() - 1] == 'S')
			{
				size_t m_idx = splice_way2.find_last_of("M");

				splice_way = splice_way2.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way2.find("S");

			if (s_index != string::npos)
			{
				splice_way2.substr(s_index + 1, splice_way2.length() - s_index - 1);
			}
		}

		size_t index2 = 0;

		string flag_str2 = " ";

		while (true)
		{
			if (index2 >= splice_way2.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way2.c_str() + index2, "%d%c", &maplen, &flag);

			if (flag_str2[0] == ' ')
			{
				if (flag == 'I')
					spliceway_vec2.push_back(make_pair(start2, -maplen));
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(start2, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec2.push_back(make_pair(start2 + maplen, 0));
				}
			}
			else if (flag_str2[0] == 'M')
			{
				if (flag == 'I')
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second, -maplen));
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec2.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'N')
			{
				if (flag == 'I')
					spliceway_vec2.back().second = -maplen;
				else if (flag == 'M')
				{
					spliceway_vec2.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec2.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec2.back().second += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str2[0] = flag;

			index2 = splice_way2.find(flag_str2, index2) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end2 = spliceway_vec2.back().first + spliceway_vec2.back().second - 1;

		mappedlen2 = mappedlen - mappedlen1;

		fusion_suffix_len = mappedlen2;
		//fusion_prefix_st, fusion_prefix_end, fusion_suffix_st, fusion_suffix_end;

		if (strand1 == '+')
		{
			fusion_prefix_st = start;

			fusion_prefix_end = end;
		}
		else
		{
			fusion_prefix_st = end;

			fusion_prefix_end = start;
		}

		if (strand2 == '+')
		{
			fusion_suffix_st = start2;
			
			fusion_suffix_end = end2;
		}
		else
		{
			fusion_suffix_st = end2;
			
			fusion_suffix_end = start2;
		}
	}
	else
	{
		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       
		//TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
		//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		//seq.1000000/2   0       chr22   51213810        255     75M     *       0       0       TCTACTCAGTCCCAGCATGTGCCCCCTAGGGAAAGGCAGCAGGCAACAGTGGAGGGAGGCCTTGGGAGCTTGTCT     
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    NM:i:2  NM:i:2  MD:Z:16G3C55

		//seq.46/1	0	chr9	123606192	255	75M	*	0	0	GACCCAGTTAACCTTTTAGGTGGTAGGCCAGATGACAGAGGGCGGCATCTGACCCGGAAGGCCCTCTGCAACCAG	
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	MD:Z:A23GA5G45G1

		//JOHNLENNON_0006:1:1:18984:1271#0/2	435	chr7	52225290	255	50S19M	=	18546378	33678863	
		//CTCGGAGCTGATGCTTGGCGGAACCAACACACTAGTGCTGNACAACACGTGTCAGGACACGCTGCTGGC	eeeeLfffffhhhehWffffhhfghgghffhhhRh\`b^^Ba`^`bbbd_fbJdd`^`WZNX[[\^``_	NM:i:4	ZF:Z:FUS_18546378_52225308(--)

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand_t, chrom, &start, &confid, 
			mapped, mate_match_chr, &mate_offset, &mate_diff, seq, qual_chr, &mis_match, alters_chr);

		strand_t = strand_t & IS_REVERSE;

		if (mis_match)
		{
			size_t mis_idx = line.find("MD:Z:");

			if (mis_idx != string::npos)
			{
				mis_info_str = line.substr(mis_idx + 5);

				//bool was_mis = false;

				size_t current_offset = 0;

				size_t current_num = 0;

				for (size_t i = 0; i < mis_info_str.length(); ++i)
				{
					if (mis_info_str[i] >= '0' && mis_info_str[i] <= '9')
					{
						current_num = current_num * 10 + mis_info_str[i] - '0';

						if (i + 1 < mis_info_str.length() && (mis_info_str[i + 1] < '0' || mis_info_str[i + 1] > '9'))
							current_offset += current_num;
					}
					else
					{
						current_num = 0;

						mis_info.push_back(make_pair(current_offset++, mis_info_str[i]));
					}
				}
				//mis_info;
			}
		}

		if (line.find("ZF:Z:FUS") != string::npos)
			is_fusion_newfmt = true;

		sscanf(tagname, "%llu", &tagidx);

		tag_name = tagname;
		chrom_name = chrom;
		splice_way = mapped;
		ori_splice_way = mapped;
		mapped_seq = seq;
		qual_str = qual_chr; 
		mate_match = mate_match_chr;

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		size_t index = 0;

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		if (isexonic)
		{
			filter_score = 5;

			junc_anchor_len = 500;

			ave_junc_mis = 0;

			pair_rate = 5;

			junc_hits = 0;
		}

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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		if (mis_info.size() && mis_info.back().first > mappedlen)
		{
			cout << "mismatch position large than mapped length?"<<endl;
		}
	}
}

SamRec::SamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& alt, const string& qualstr, 
	const string& matematch, size_t mateoffest, int matediff, string line) : tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
	mis_match(mismatch), matched_id(-1), mate_match(matematch), mate_offset(mateoffest), mate_diff(matediff), wrong_format(false), tagidx(tidx), filter_score(0), best(0), alters(alt), qual_str(qualstr), 
	isunique(false), isexonic(false), isspliced(false), ave_intron_len(100000000), ave_junc_mis(100000000), canon_count(0), noncanon_count(0), canon_rate(0), issmallins(false), issmalldel(false), cur_line(line),
	/*is_insert(false),*/ mappedlen(0), mappedlen1(0), mappedlen2(0), intron_size(0)
{
	sscanf(tag_name.c_str(), "%llu", &tagidx);

	size_t index = 0;

	if ((splice_way.find("I") != string::npos))
		issmallins = true;
	else if ((splice_way.find("D") != string::npos))
		issmalldel = true;
	else if (splice_way.find("N") != string::npos)
		isspliced = true;
	else
		isexonic = true;		

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
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + spliceway_vec.back().second + maplen, 0));
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
				spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
		}

		flag_str[0] = flag;

		index = splice_way.find(flag_str, index) + 1;
	}

	end = spliceway_vec.back().first + spliceway_vec.back().second - 1;
}


void SamRec::Set(const string& line, int min_ins) 
{
	Clear();

	cur_line = line;

	//if (line.find("seq.10001997/") != string::npos)
	//{
	//	cout << line << endl;
	//}

	char tagname[1000], chrom[100], mapped[100], seq[1000], qual_chr[1000], alters_chr[1000];

	char chrom2[100], mapped2[100];

	char mate_match_chr[1000];

	string alterstr = "";

	char strand[100];

	sscanf(line.c_str(), "%s\t%s\t%s", tagname, chrom, strand);

	if (string(strand) == "+" || string(strand) == "-")
		is_fusion = true;
	else
		is_fusion = false;

	if (line.find("ZF:Z:FUS") != string::npos)
		is_fusion_newfmt = true;
	else
		is_fusion_newfmt = false;

	if (is_fusion)
	{
		//14~R_7/2	chr1	+	11184637	54M	chr11	+	61646784	21M	GGCGCAGATCTTCATGGCCTGTTAGAAGGAAAACACACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGGATGAA	
		//283348985863448157653149274159715381432428925324118868345815774332128375427~CTGCATCACACGCTCATCCTGGCGCAGATCTTCATGGCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTG
		//CCTGAGAGATGGCCAGGATGAAGGCGGCCAGGGCACTGGGCACCCAGCCAGGACCCAGGAGGTAGATAAGGAGCC	NM:i:2

		size_t read_count = sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tagname, chrom, &strand1, &start, mapped, 
			chrom2, &strand2, &start2, mapped2, seq, qual_chr, &mis_match);

		sscanf(tagname, "%llu", &tagidx);

		mapped_seq = seq;

		qual_str = qual_chr; 

		syn_seq = qual_chr;

		size_t last_idx = syn_seq.find_last_of("~");

		qual_str = qual_str.substr(0, last_idx);

		syn_seq = syn_seq.substr(last_idx + 1, syn_seq.length() - last_idx - 1);

		tag_name = tagname;

		tag_base_name = tag_name.substr(0, tag_name.length() - 1);

		if (tag_name[tag_name.length() - 1] == '1' || tag_name[tag_name.length() - 1] == 'a')
			end_id = 1;
		else if (tag_name[tag_name.length() - 1] == '2' || tag_name[tag_name.length() - 1] == 'b')
			end_id = 2;

		if (strand1 == '+')
			strand_t = 0;
		else
			strand_t = 16;

		if (strand2 == '+')
			strand_t2 = 0;
		else
			strand_t2 = 16;

		if (string(chrom) < string(chrom2) || (string(chrom) == string(chrom2) && start <= start2))
		{
			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;
			
			chrom_name = chrom;
			splice_way = mapped;
			ori_splice_way = mapped;

			splice_way2 = mapped2;
			ori_splice_way2 = mapped2;
			chrom_name2 = chrom2;
		}
		else
		{
			if (strand1 == '+')
				strand1 = '-';
			else
				strand1 = '+';

			if (strand2 == '+')
				strand2 = '-';
			else
				strand2 = '+';

			swap(strand1, strand2);

			if (strand1 == '+')
				strand_t = 0;
			else
				strand_t = 16;

			if (strand2 == '+')
				strand_t2 = 0;
			else
				strand_t2 = 16;
			
			swap(start, start2);
			
			chrom_name = chrom2;//chrom;
			splice_way = mapped2;//mapped;
			ori_splice_way = mapped2;//mapped;

			splice_way2 = mapped;//mapped2;
			ori_splice_way2 = mapped;//mapped2;
			chrom_name2 = chrom;//chrom2;

			is_swapped = true;
		}

		//strand_t |= IS_PRIMARY;

		//strand_t2 |= IS_PRIMARY;

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}

			//ori_splice_way = splice_way;
		}	

		size_t index = 0;

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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		mappedlen1 = mappedlen;

		fusion_prefix_len = mappedlen1;

		//segment 2	

		if (mapped2[0] == '*')
		{
			isunmapped = true;
			return;
		}

		if ((splice_way2.find("I") != string::npos))
			issmallins = true;
		if ((splice_way2.find("D") != string::npos))
			issmalldel = true;
		if (splice_way2.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		for (size_t i = 0; i < splice_way2.length(); ++i)
			if (splice_way2[i] == 'D')
				splice_way2[i] = 'N';

		if ((splice_way2.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way2[splice_way2.length() - 1] == 'S')
			{
				size_t m_idx = splice_way2.find_last_of("M");

				splice_way2 = splice_way2.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way2.find("S");

			if (s_index != string::npos)
			{
				splice_way2 = splice_way2.substr(s_index + 1, splice_way2.length() - s_index - 1);
			}

			//ori_splice_way2 = splice_way2;
		}

		size_t index2 = 0;

		string flag_str2 = " ";

		while (true)
		{
			if (index2 >= splice_way2.length())
				break;

			int maplen;

			char flag;

			sscanf(splice_way2.c_str() + index2, "%d%c", &maplen, &flag);

			if (flag_str2[0] == ' ')
			{
				if (flag == 'I')
					spliceway_vec2.push_back(make_pair(start2, -maplen));
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(start2, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout<<"start with N?"<<endl;
					spliceway_vec2.push_back(make_pair(start2 + maplen, 0));
				}
			}
			else if (flag_str2[0] == 'M')
			{
				if (flag == 'I')
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second, -maplen));
				else if (flag == 'M')
				{
					cout << "continue Ms?"<<endl;
					spliceway_vec2.back().second += maplen;

					mappedlen += maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + spliceway_vec2.back().second + maplen, 0));

					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'N')
			{
				if (flag == 'I')
					spliceway_vec2.back().second = -maplen;
				else if (flag == 'M')
				{
					spliceway_vec2.back().second = maplen;

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					cout << "continue Ns?"<<endl;
					spliceway_vec2.back().first += maplen;
					intron_size += maplen;
				}
			}
			else if (flag_str2[0] == 'I')
			{
				if (flag == 'I')
				{
					cout << "continue Is?"<<endl;
					spliceway_vec2.back().second += -maplen;
				}
				else if (flag == 'M')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first, maplen));

					mappedlen += maplen;

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec2.push_back(make_pair(spliceway_vec2.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str2[0] = flag;

			index2 = splice_way2.find(flag_str2, index2) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end2 = spliceway_vec2.back().first + spliceway_vec2.back().second - 1;

		mappedlen2 = mappedlen - mappedlen1;

		fusion_suffix_len = mappedlen2;
		//fusion_prefix_st, fusion_prefix_end, fusion_suffix_st, fusion_suffix_end;

		if (strand1 == '+')
		{
			fusion_prefix_st = start;

			fusion_prefix_end = end;
		}
		else
		{
			fusion_prefix_st = end;

			fusion_prefix_end = start;
		}

		if (strand2 == '+')
		{
			fusion_suffix_st = start2;
			
			fusion_suffix_end = end2;
		}
		else
		{
			fusion_suffix_st = end2;
			
			fusion_suffix_end = start2;
		}
	}
	else
	{
		//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       
		//TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
		//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

		//seq.1000000/2   0       chr22   51213810        255     75M     *       0       0       TCTACTCAGTCCCAGCATGTGCCCCCTAGGGAAAGGCAGCAGGCAACAGTGGAGGGAGGCCTTGGGAGCTTGTCT     
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII    NM:i:2  NM:i:2  MD:Z:16G3C55

		//seq.46/1	0	chr9	123606192	255	75M	*	0	0	GACCCAGTTAACCTTTTAGGTGGTAGGCCAGATGACAGAGGGCGGCATCTGACCCGGAAGGCCCTCTGCAACCAG	
		//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	MD:Z:A23GA5G45G1

		size_t read_count = sscanf(line.c_str(), "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\t%s", tagname, &strand_t, chrom, &start, &confid, 
			mapped, mate_match_chr, &mate_offset, &mate_diff, seq, qual_chr, &mis_match, alters_chr);

		strand_t = strand_t & IS_REVERSE;

		//strand_t |= IS_PRIMARY;

		if (mis_match)
		{
			size_t mis_idx = line.find("MD:Z:");

			if (mis_idx != string::npos)
			{
				mis_info_str = line.substr(mis_idx + 5);

				//bool was_mis = false;

				size_t current_offset = 0;

				size_t current_num = 0;

				for (size_t i = 0; i < mis_info_str.length(); ++i)
				{
					if (mis_info_str[i] >= '0' && mis_info_str[i] <= '9')
					{
						current_num = current_num * 10 + mis_info_str[i] - '0';

						if (i + 1 < mis_info_str.length() && (mis_info_str[i + 1] < '0' || mis_info_str[i + 1] > '9'))
							current_offset += current_num;
					}
					else
					{
						current_num = 0;

						mis_info.push_back(make_pair(current_offset++, mis_info_str[i]));
					}
				}
				//mis_info;
			}
		}

		if (line.find("ZF:Z:FUS") != string::npos)
			is_fusion_newfmt = true;

		sscanf(tagname, "%llu", &tagidx);

		tag_name = tagname;

		tag_base_name = tag_name.substr(0, tag_name.length() - 1);

		if (tag_name[tag_name.length() - 1] == '1' || tag_name[tag_name.length() - 1] == 'a')
			end_id = 1;
		else if (tag_name[tag_name.length() - 1] == '2' || tag_name[tag_name.length() - 1] == 'b')
			end_id = 2;

		chrom_name = chrom;
		splice_way = mapped;
		ori_splice_way = mapped;
		mapped_seq = seq;
		qual_str = qual_chr; 
		mate_match = mate_match_chr;

		if (read_count == 13)
			alters = alters_chr;

		if (mapped[0] == '*')
		{
			isunmapped = true;
			return;
		}

		size_t index = 0;

		for (size_t i = 0; i < splice_way.length(); ++i)
			if (splice_way[i] == 'D')
				splice_way[i] = 'N';

		if ((splice_way.find("S") != string::npos))
		{
			isclipped = true;

			//clipe tail
			if (splice_way[splice_way.length() - 1] == 'S')
			{
				size_t m_idx = splice_way.find_last_of("M");

				splice_way = splice_way.substr(0, m_idx + 1);
			}

			size_t s_index = splice_way.find("S");

			if (s_index != string::npos)
			{
				splice_way = splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
			}
		}

		if ((splice_way.find("I") != string::npos))
			issmallins = true;
		if ((splice_way.find("D") != string::npos))
			issmalldel = true;
		if (splice_way.find("N") != string::npos)
			isspliced = true;

		if (!issmallins && !issmalldel && !isspliced)
			isexonic = true;	

		if (isexonic)
		{
			filter_score = 5;

			junc_anchor_len = 500;

			ave_junc_mis = 0;

			junc_hits = 0;

			pair_rate = 5;
		}

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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
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

					if (maplen < min_anchor)
						min_anchor = maplen;

					if (maplen > max_anchor)
						max_anchor = maplen;
				}
				else if (flag == 'N')
				{
					spliceway_vec.push_back(make_pair(spliceway_vec.back().first + maplen, 0));
					intron_size += maplen;
				}
			}

			flag_str[0] = flag;

			index = splice_way.find(flag_str, index) + 1;

			if (flag == 'N' && maplen <= min_ins)
			{
				//splice_way[index-1] = 'D';
				issmalldel = true;
			}
		}

		end = spliceway_vec.back().first + spliceway_vec.back().second - 1;

		if (mis_info.size() && mis_info.back().first > mappedlen)
		{
			cout << "mismatch position large than mapped length?"<<endl;
		}
	}
}

void 
SamRec::Clear()
{
	spliceway_vec.clear();

	spliceway_vec2.clear();

	corresponding_juncs.clear();

	encompassed_juncs.clear();

	left_splice_ways.clear();
	
	right_splice_ways.clear();

	mis_info.clear();

	filter_type.clear();

	matched_id = -1;
	
	wrong_format = false;
	
	tagidx = -1;
	
	filter_score = 0;
	
	best = 0;

	isunique = false;
	
	isexonic = false;
	
	isspliced = false;
	
	ave_intron_len = 0;
	
	ave_junc_mis = 0;
	
	canon_count = 0;

	noncanon_count = 0;
	
	canon_rate = 0;
	
	issmallins = false;
	
	issmalldel = false;
	
	iscanonical = false;

	issemicanonical = false;
	
	isnoncanoical = false; 
	
	isunmapped = false;
	
	isclipped = false; 
	/*is_insert(false),*/ 
	
	mappedlen = 0;
	
	intron_size = 0;
	
	min_anchor = -1;
	
	max_anchor = 0;
	
	pair_rate = 0;
	
	m_contig_len = 0;

	fusion_prefix_len = 0;
	
	fusion_suffix_len = 0;
	
	is_swapped = false;
	
	mate_offset = 0;
	
	mate_offset2 = 0;
	
	mate_diff = 0;

	mate_diff2 = 0;

	//m = NOT_FILTERED;
	
	mate_match = "*";
	
	mate_match2 = "*";

	paired_type = SINGLE;

	chrom_name.clear();
	splice_way.clear();
	ori_splice_way.clear();

	splice_way2.clear();
	ori_splice_way2.clear();
	chrom_name2.clear();

	start = 0;

	start2 = 0;

	strand_t = 0;
	
	strand_t2 = 0;

	forward_count = 0;
	
	reverse_count = 0;

	is_fusion_newfmt = false;

	is_fusion = false;

	fusion_prefix_paired = false;
	
	fusion_suffix_paired = false;
}

string 
SamRec::tostring(size_t total, size_t idx) const
{
	//if (!cur_line.empty())
	//	return cur_line;

	char sam_rec_char[5000];

	if (is_fusion)
		sprintf(sam_rec_char, "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tag_name.c_str(), chrom_name.c_str(), strand1, start, ori_splice_way.c_str(), 
			chrom_name2.c_str(), strand2, start2, ori_splice_way2.c_str(), mapped_seq.c_str(), qual_str.c_str(), mis_match);
	else
	{
		if (forward_count || reverse_count)
		{

			if (forward_count >= reverse_count)
				sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu\tXS:A:+"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
				chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);
			else 
				sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu\tXS:A:-"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
				chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);
		}
		else
			sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
				chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);
	}

	return sam_rec_char;
}

string
SamRec::tostandfusion() const 
{
	//R_112:1	115	11	61646784	255	38S37M	1	11184653	0	
	//GCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGTATGAAGGCGGCCAGGGCACTG	181166788356889826893782431399631578517864278199618228323647256272457757777	ZF:Z:FUS_11184689_446432191(++)	ZT:Z:Seed

	//R_112:1	179	1	11184653	255	38M37S	11	61646784	0	
	//GCCTTTTAGAAGGAAAACAAACTCATGTCCGTTGCTGCCTGAGAGATGGCCAGTATGAAGGCGGCCAGGGCACTG	181166788356889826893782431399631578517864278199618228323647256272457757777	ZF:Z:FUS_11184689_446432191(++)	ZT:Z:Seed

	char sam_rec_char[5000];

	if (is_fusion)
	{
		char append_chr1[1000], append_chr2[1000];

		sprintf(append_chr1, "%lluS", mappedlen2);

		sprintf(append_chr2, "%lluS", mappedlen1);

		string sp1, sp2;

		if (is_fusion_newfmt)
		{
			sp1 = ori_splice_way;

			sp2 = ori_splice_way2;
		}
		else
		{
			if (is_swapped)
			{
				sp1 = append_chr1; sp1.append(splice_way);

				sp2 = splice_way2; sp2.append(append_chr2);
			}
			else
			{
				sp1 = splice_way; sp1.append(append_chr1);

				sp2 = append_chr2; sp2.append(splice_way2);
			}
		}

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tZF:Z:FUS_%llu_%llu(%c%c)\n%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tZF:Z:FUS_%llu_%llu(%c%c)", 
			tag_name.c_str(), strand_t, chrom_name.c_str(), start, sp1.c_str(), 
			mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, fusion_prefix_end, fusion_suffix_st, strand1, strand2,
			tag_name.c_str(), strand_t2, chrom_name2.c_str(), start2, sp2.c_str(), 
			mate_match2.c_str(), mate_offset2, mate_diff2, mapped_seq.c_str(), qual_str.c_str(), mis_match, fusion_prefix_end, fusion_suffix_st, strand1, strand2);
	}
	else
		return "";

	return sam_rec_char;
}


void SamRec::genearte_fusion_str()
{
	/*char sam_rec_char[5000];

	if (is_fusion)
	{
		char append_chr1[1000], append_chr2[1000];

		sprintf(append_chr1, "%lluS", mappedlen2);

		sprintf(append_chr2, "%lluS", mappedlen1);

		string sp1, sp2;

		if (is_swapped)
		{
			sp1 = append_chr1; sp1.append(ori_splice_way);

			sp2 = ori_splice_way2; sp2.append(append_chr2);
		}
		else
		{
			sp1 = ori_splice_way; sp1.append(append_chr1);

			sp2 = append_chr2; sp2.append(ori_splice_way2);
		}

		sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tZF:Z:FUS_%llu_%llu(%c%c)\n%s\t%hu\t%s\t%llu\t255\t%s\t%s\t%llu\t%ld\t%s\t%s\tNM:i:%hu\tZF:Z:FUS_%llu_%llu(%c%c)", 
			tag_name.c_str(), strand_t, chrom_name.c_str(), start, sp1.c_str(), 
			mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, fusion_prefix_end, fusion_suffix_st, strand1, strand2,
			tag_name.c_str(), strand_t2, chrom_name2.c_str(), start2, sp2.c_str(), 
			mate_match2.c_str(), mate_offset2, mate_diff2, mapped_seq.c_str(), qual_str.c_str(), mis_match, fusion_prefix_end, fusion_suffix_st, strand1, strand2);
	}
	else
		return "";

	return sam_rec_char;*/
}

void SamRec::generate_normal_str()
{
	//char sam_rec_char[5000];

	//if (is_fusion)
	//	sprintf(sam_rec_char, "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tag_name.c_str(), chrom_name.c_str(), strand1, start, ori_splice_way.c_str(), 
	//		chrom_name2.c_str(), strand2, start2, ori_splice_way2.c_str(), mapped_seq.c_str(), qual_str.c_str(), mis_match);
	//else
	//	sprintf(sam_rec_char, "%s\t%hu\t%s\t%llu\t%hu\t%s\t%s\t%llu\t%d\t%s\t%s\tNM:i:%hu\tIH:i:%llu\tHI:i:%llu"/*\t%s",*//*\t%d\t%lf"*/, tag_name.c_str(), strand_t, 
	//		chrom_name.c_str(), start, confid, ori_splice_way.c_str(), mate_match.c_str(), mate_offset, mate_diff, mapped_seq.c_str(), qual_str.c_str(), mis_match, total, idx/*, clippled_way.c_str()*//*, best, filter_score*/);

	//return sam_rec_char;
}

void SamRec::set_unmapped()
{
	mate_match = "*";

	mate_offset = 0;

	mate_diff = 0;

	mate_match2 = "*";

	mate_offset2 = 0;

	mate_diff2 = 0;
}

bool SamRec::clip_by_small_anchor(bool add_S/*JunctionSeed* filtered_junc*/)
{
	const string& readstr = mapped_seq;

	if (max_anchor >= 50)
	{
		vector<pair<size_t, int> >::iterator vp_iter, max_vp_iter;

		int max_frame_len = 0;
	
		for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
		{
			if (vp_iter->second > max_frame_len)
			{
				max_frame_len = vp_iter->second;
				max_vp_iter = vp_iter;
			}
		}

		pair<size_t, int> best = *max_vp_iter;

		if (max_frame_len >= 50)
		{
			size_t prev_map_len = 0;

			for (vp_iter =  spliceway_vec.begin(); vp_iter != max_vp_iter; ++vp_iter)
			{
				prev_map_len += abs(vp_iter->second);
			}

			if (!add_S)
			{
				qual_str = qual_str.substr(prev_map_len, max_vp_iter->second);

				mapped_seq = mapped_seq.substr(prev_map_len, max_vp_iter->second);
			}

			mappedlen = max_vp_iter->second;

			start = max_vp_iter->first;

			end = max_vp_iter->first + max_vp_iter->second - 1;

			min_anchor = max_vp_iter->second;

			max_anchor = max_vp_iter->second;

			spliceway_vec.clear();

			spliceway_vec.push_back(best);

			char clipped_jmp[100];

			if (true)
			{
				size_t tail_clip = mapped_seq.length() - prev_map_len - best.second;

				if (add_S && tail_clip > mapped_seq.length() )
					cerr <<" tail clip length"<<endl;

				if (!add_S)
					sprintf(clipped_jmp, "%dM", best.second);
				else if (add_S && prev_map_len > 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%lluS%dM%lluS", prev_map_len, best.second, tail_clip);
				else if (add_S && prev_map_len == 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%dM%lluS", best.second, tail_clip);
				else if (add_S && (prev_map_len > 0 && tail_clip == 0))
					sprintf(clipped_jmp, "%lluS%dM", prev_map_len, best.second);
				else 
				{
					cerr <<"clip not correct"<<endl;
					cerr << this->tostring(0, 0)<<endl;
				}
			}

			ori_splice_way = clipped_jmp;
			 
			//unsigned short mis_match;

			canon_rate = 0;

			ave_intron_len = 500;

			filter_score = 5;

			ave_junc_mis = 0;

			pair_rate = 5;

			isexonic = true;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = true;

			canon_count = 0;			
			noncanon_count = 0;

			return true;
		}
	}

	return false;
}

bool SamRec::clip_fusion_alignment(bool add_S/*JunctionSeed**/)
{
	const string& readstr = mapped_seq;

	size_t fusion_prefix_len, fusion_suffix_len;

	if (fusion_prefix_len > 50 || fusion_suffix_len > 50)
	{
		if (fusion_prefix_len >= fusion_suffix_len)
		{
			if (!add_S)
			{
				qual_str = qual_str.substr(0, fusion_prefix_len);

				mapped_seq = mapped_seq.substr(0, fusion_prefix_len);
			}

			mappedlen = fusion_prefix_len;

			start = spliceway_vec.front().first;

			end = spliceway_vec.back().first + spliceway_vec.back().second;

			char clipped_jmp[100];

			sprintf(clipped_jmp, "%lluS", fusion_suffix_len);

			isexonic = false;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = false;
			is_fusion = false;
			is_fusion_newfmt = false;

			if ((ori_splice_way.find("I") != string::npos))
				issmallins = true;
			if ((ori_splice_way.find("D") != string::npos))
				issmalldel = true;
			if (ori_splice_way.find("N") != string::npos)
				isspliced = true;

			if (add_S)
				ori_splice_way.append(clipped_jmp);

			//unsigned short mis_match;

			canon_rate = 0;

			ave_intron_len = 500;

			filter_score = 5;

			ave_junc_mis = 0;

			pair_rate = 5;

			canon_count = 0;

			noncanon_count = 0;

			ori_splice_way2.clear();

			spliceway_vec2.clear();

			splice_way2.clear();

			return true;
		}
		else
		{
			spliceway_vec = spliceway_vec2;

			ori_splice_way = ori_splice_way2;

			splice_way = splice_way2;

			strand1 = strand2;

			strand_t = strand_t2;

			mappedlen1 = mappedlen2;

			if (!add_S)
			{
				qual_str = qual_str.substr(fusion_prefix_len + 1, fusion_suffix_len);

				mapped_seq = mapped_seq.substr(fusion_prefix_len + 1, fusion_suffix_len);
			}

			mappedlen = fusion_suffix_len;

			start = spliceway_vec.front().first;

			end = spliceway_vec.back().first + spliceway_vec.back().second;

			char clipped_jmp[100];

			sprintf(clipped_jmp, "%lluS", fusion_prefix_len);

			isexonic = false;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = false;
			is_fusion = false;
			is_fusion_newfmt = false;

			if ((ori_splice_way.find("I") != string::npos))
				issmallins = true;
			if ((ori_splice_way.find("D") != string::npos))
				issmalldel = true;
			if (ori_splice_way.find("N") != string::npos)
				isspliced = true;

			if (add_S)
				ori_splice_way = clipped_jmp + ori_splice_way;

			//unsigned short mis_match;

			canon_rate = 0;

			ave_intron_len = 500;

			filter_score = 5;

			ave_junc_mis = 0;

			pair_rate = 5;

			canon_count = 0;

			noncanon_count = 0;

			ori_splice_way2.clear();

			spliceway_vec2.clear();

			splice_way2.clear();

			return true;
		}
	}

	return false;
}

bool SamRec::clip_by_end_mismatch(bool add_S, size_t max_insert_len/*JunctionSeed* filtered_junc*/)
{
	vector<pair<size_t, char> >::iterator mis_info_iter;

	for (mis_info_iter = mis_info.begin(); mis_info_iter != mis_info.end(); ++mis_info_iter)
	{
		if (mis_info_iter->first >= max_insert_len)
			break;
	}

	string clipped_str;

	if (mis_info_iter != mis_info.begin())
	{
		spliceway_vec.front().first += (mis_info_iter - 1)->first + 1;

		start += (mis_info_iter - 1)->first + 1;

		spliceway_vec.front().second -= (mis_info_iter - 1)->second;

		mis_match = mis_info_iter - mis_info.begin();

		char clipped_chr[1000];

		sprintf(clipped_chr, "%lluI", (mis_info_iter - 1)->first + 1);

		clipped_str.append(clipped_chr);
		//spliceway_vec.pu
	}

	const string& readstr = mapped_seq;

	if (max_anchor >= 50)
	{
		vector<pair<size_t, int> >::iterator vp_iter, max_vp_iter;

		int max_frame_len = 0;
	
		for (vp_iter =  spliceway_vec.begin(); vp_iter != spliceway_vec.end(); ++vp_iter)
		{
			if (vp_iter->second > max_frame_len)
			{
				max_frame_len = vp_iter->second;
				max_vp_iter = vp_iter;
			}
		}

		pair<size_t, int> best = *max_vp_iter;

		if (max_frame_len >= 50)
		{
			size_t prev_map_len = 0;

			for (vp_iter =  spliceway_vec.begin(); vp_iter != max_vp_iter; ++vp_iter)
			{
				prev_map_len += abs(vp_iter->second);
			}

			if (!add_S)
			{
				qual_str = qual_str.substr(prev_map_len, max_vp_iter->second);

				mapped_seq = mapped_seq.substr(prev_map_len, max_vp_iter->second);
			}

			mappedlen = max_vp_iter->second;

			start = max_vp_iter->first;

			end = max_vp_iter->first + max_vp_iter->second - 1;

			min_anchor = max_vp_iter->second;

			max_anchor = max_vp_iter->second;

			spliceway_vec.clear();

			spliceway_vec.push_back(best);

			char clipped_jmp[100];

			if (true)
			{
				size_t tail_clip = mapped_seq.length() - prev_map_len - best.second;

				if (add_S && tail_clip > mapped_seq.length() )
					cerr <<" tail clip length"<<endl;

				if (!add_S)
					sprintf(clipped_jmp, "%dM", best.second);
				else if (add_S && prev_map_len > 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%lluS%dM%lluS", prev_map_len, best.second, tail_clip);
				else if (add_S && prev_map_len == 0 && tail_clip > 0)
					sprintf(clipped_jmp, "%dM%lluS", best.second, tail_clip);
				else if (add_S && (prev_map_len > 0 && tail_clip == 0))
					sprintf(clipped_jmp, "%lluS%dM", prev_map_len, best.second);
				else 
				{
					cerr <<"clip not correct"<<endl;
					cerr << this->tostring(0, 0)<<endl;
				}
			}

			ori_splice_way = clipped_jmp;
			 
			//unsigned short mis_match;

			canon_rate = 0;

			ave_intron_len = 0;

			filter_score = 0;

			ave_junc_mis = 0;

			isexonic = true;
			isspliced = false;
			issmallins = false;
			issmalldel = false;
			iscanonical = false;
			issemicanonical = false;
			isnoncanoical = false;
			isunmapped = false;
			isclipped = true;

			canon_count = 0;			
			noncanon_count = 0;

			return true;
		}
	}

	return false;
}

void 
SamRec::setfusionbit()
{
	strand_t |= IS_PAIRED;

	strand_t |= IS_PAIRED_MAPPED;

	if (strand_t2 & IS_REVERSE)
		strand_t |= IS_MATE_REVERSE;

	strand_t |= IS_FIRST_END;

	//if (mate_match == "*")
	//{
	//	if (chrom_name2 == chrom_name)
	//		mate_match = "=";
	//	else
	//		mate_match = chrom_name2;

	//	mate_offset = start2;

	//	mate_diff = fusion_prefix_st - fusion_suffix_end;
	//}

	strand_t2 |= IS_PAIRED;

	strand_t2 |= IS_PAIRED_MAPPED;

	if (strand_t & IS_REVERSE)
		strand_t2 |= IS_MATE_REVERSE;

	strand_t2 |= IS_SECOND_END;

	//if (mate_match2 == "*")
	//{
	//	if (chrom_name2 == chrom_name)
	//		mate_match2 = "=";
	//	else
	//		mate_match2 = chrom_name2;

	//	mate_offset2 = start;

	//	mate_diff2 = fusion_suffix_end - fusion_prefix_st;
	//}

}

bool SamRec::modify_jumpcode_by_filtered_junc()
{
	vector<JunctionSeed*>::iterator junc_iter;

	bool filtered = false;

	if (issmallins)
		filtered = true;

	for (junc_iter = corresponding_juncs.begin(); !filtered && junc_iter != corresponding_juncs.end(); ++junc_iter)
	{
		//if ((*junc_iter)->m_filtered_type != NOT_FILTERED)
		//	filtered = true;

		if ((*junc_iter)->m_filtered_type == FILTERED_BY_SMALL_ANCHOR)
		{

		}

		//if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_SMALL_ANCHOR)
		//{
		//	if ((*junc_iter)->m_max_prefix_len < m_min_junc_anchor)
		//		;
		//}

		//if ((*sam_rec_ptr_iter)->filter_type.empty())
		//	filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);
		//else if ((*sam_rec_ptr_iter)->filter_type.size() == 1)
		//{
		//	if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == NOT_FILTERED)
		//		filtered_sam_rec_ptr.push_back(*sam_rec_ptr_iter);
		//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_SMALL_ANCHOR)
		//		;//trim
		//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_SMALL_DELETION)
		//		;//convert
		//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_LARGE_MIN_ANCHOR_DIFF)
		//		;//filter
		//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_LARGE_MULTIPLE_PAIRED)
		//		;//filter
		//	else if ((*(*sam_rec_ptr_iter)->filter_type.begin()) == FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR)
		//		;
		//}
		//else
		//{
		//}
	}

	return true;
}

bool SamRec::Swap()
{
	swap(splice_way2, splice_way);

	swap(ori_splice_way2, ori_splice_way);

	swap(chrom_name2, chrom_name);

	swap(splice_way2, splice_way);

	swap(start2, start);

	swap(end2, end);

	swap(strand_t2, strand_t);

	swap(mappedlen1, mappedlen2);

	swap(fusion_prefix_len, fusion_suffix_len);

	swap(spliceway_vec2, spliceway_vec);

	swap(mappedlen1, mappedlen2);

	swap(strand1, strand2);

	is_swapped = true;

	return true;
}

void SamRec::fixfusionseq(SamRec& samrec2)
{
	string combseq;

	string newqualseq;

	if (ori_splice_way[ori_splice_way.length() - 1] == 'S')
	{
		string fusionseq1 = mapped_seq.substr(0, mappedlen), fusionseq2 = mapped_seq.substr(mappedlen, mapped_seq.length() - mappedlen);
		
		string qualseq1 = qual_str.substr(0, mappedlen), qualseq2 = qual_str.substr(mappedlen, qual_str.length() - mappedlen);

		if (strand_t & IS_REVERSE)
		{
			fusionseq1 = revcomp(fusionseq1);

			reverse(qualseq1.begin(), qualseq1.end());
		}

		if (samrec2.strand_t & IS_REVERSE)
		{
			fusionseq2 = revcomp(fusionseq2);

			reverse(qualseq2.begin(), qualseq2.end());
		}

		combseq = fusionseq1 + fusionseq2;

		newqualseq = qualseq1 + qualseq2;
	}
	else
	{
		string fusionseq2 = mapped_seq.substr(0, mapped_seq.length() - mappedlen), fusionseq1 = mapped_seq.substr(mapped_seq.length() - mappedlen, mappedlen);

		string qualseq2 = qual_str.substr(0, qual_str.length() - mappedlen), qualseq1 = qual_str.substr(qual_str.length() - mappedlen, mappedlen);

		if (!(strand_t & IS_REVERSE))
		{
			fusionseq1 = revcomp(fusionseq1);

			reverse(qualseq1.begin(), qualseq1.end());
		}

		if (!(samrec2.strand_t & IS_REVERSE))
		{
			fusionseq2 = revcomp(fusionseq2);

			reverse(qualseq2.begin(), qualseq2.end());
		}

		combseq = fusionseq2 + fusionseq1;

		newqualseq = qualseq2 + qualseq1;
	}

	mapped_seq = combseq;

	samrec2.mapped_seq = combseq;

	qual_str = newqualseq;

	samrec2.qual_str = newqualseq;
}