#include "FusionSamRec.h"
//
FusionSamRec::FusionSamRec(const string& line, int min_ins) : matched_id(-1),  wrong_format(false), tagidx(-1), filter_score(0), best(0), 
	isunique(false), isexonic(false), isspliced(false), ave_intron_len(100000000), ave_junc_mis(100000), canon_count(0), 
	noncanon_count(0), canon_rate(0), issmallins(false), issmalldel(false), iscanonical(false), 
	issemicanonical(false), isnoncanoical(false), isunmapped(false), isclipped(false), cur_line(line),
	/*is_insert(false),*/ mappedlen(0), intron_size(0), min_anchor(-1), max_anchor(0)
{
	char tagname[1000], chrom[100], mapped[100], seq[1000], qual_chr[1000], alters_chr[1000];

	char chrom2[100], mapped2[100];

	string alterstr = "";

	//TRAN00000027662:59:252  16      chr14   56062712        255     8M157375N92M    *       0       0       TATTATTTTCCGCTTTCCCTGGGCTTACAGAGAATCCTTGCCCTTCTTGTACTGTGTCACTTTATGGGGTTGGTGCTTGCCACACTTCTTACAGAAAGTC 
	//#######$#%'*,*++*,/121222012122233455766655666555666677989:;<<=>>>>>==>??>>>>>>>>>>>>>>>>>>>=>=>>>>> NM:i:5  8:G>T,24:C>T,33:T>A,63:G>A,81:T>A

	//2661~87323~seq.43662/1	chr17	-	40970259	28M	chr17	-	40971609	19M1164N28M	CGAGAAGGTCCAGGCTGAGGCTGAAAGCCCAGGAGGAAGAGACTAACTCAGGTGAGGAGCCATTTATTGAAACTC	
	//IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII~AGCTGGAAGACGTGGAAAAGAACCGCAAGATAGTGGCAGAAAATCTCGAGAAGGTCCAGGCTGAGGCTGAGAGACCAGGAG
	//GAAGAGACTAACTCAGGAGAGGAGCCATTTATTGAAACTCCTCGCCAGGATGGTGTCTCTCGCAGAT	NM:i:3

	//char strand1, strand2;

	size_t read_count = sscanf(line.c_str(), "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tagname, chrom, &strand1, &start, mapped, 
		chrom2, &strand2, &start2, mapped2, seq, qual_chr, &mis_match);

	sscanf(tagname, "%llu", &tagidx);

	tag_name = tagname;
	chrom_name = chrom;
	splice_way = mapped;
	ori_splice_way = mapped;
	
	mapped_seq = seq;
	qual_str = qual_chr; 

	splice_way2 = mapped2;
	ori_splice_way2 = mapped2;
	chrom_name2 = chrom2;


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
			splice_way.substr(s_index + 1, splice_way.length() - s_index - 1);
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

	//segment 2	

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
}

FusionSamRec::FusionSamRec(const string& tname, unsigned short strand, const string& cname, size_t st, unsigned short conf, const string& spliceway, const string& mapseq, unsigned short mismatch, size_t tidx, const string& alt, const string& qualstr, 
	char matematch, size_t mateoffest, int matediff, string line) : tag_name(tname), strand_t(strand), chrom_name(cname), start(st), confid(conf), splice_way(spliceway), mapped_seq(mapseq), 
	mis_match(mismatch), matched_id(-1), mate_match(matematch), mate_offset(mateoffest), mate_diff(matediff), wrong_format(false), tagidx(tidx), filter_score(0), best(0), alters(alt), qual_str(qualstr), 
	isunique(false), isexonic(false), isspliced(false), ave_intron_len(100000000), ave_junc_mis(100000000), canon_count(0), noncanon_count(0), canon_rate(0), issmallins(false), issmalldel(false), cur_line(line),
	/*is_insert(false),*/ mappedlen(0), intron_size(0)
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

string 
FusionSamRec::tostring(size_t total, size_t idx) const
{
	//if (!cur_line.empty())
	//	return cur_line;

	char sam_rec_char[5000];

	//cout << tag_name.c_str()<<'\t' << strand_t<<'\t' <<  chrom_name.c_str()<<'\t' <<  start<<'\t' <<  confid<<'\t' <<  splice_way.c_str()<<"\t*\t0\t0\t" <<  mapped_seq.c_str()<<'\t' <<  Is.substr(0, mapped_seq.length()).c_str()<<'\t' <<  mis_match << endl;

	// IH:i:2  HI:i:2
	sprintf(sam_rec_char, "%s\t%s\t%c\t%llu\t%s\t%s\t%c\t%llu\t%s\t%s\t%s\tNM:i:%hu", tag_name.c_str(), chrom_name.c_str(), strand1, start, ori_splice_way.c_str(), 
		chrom_name2.c_str(), strand2, start2, ori_splice_way2.c_str(), mapped_seq.c_str(), qual_str.c_str(), mis_match);

	return sam_rec_char;
}


void FusionSamRec::set_unmapped()
{
	mate_match = '*';

	mate_offset = 0;

	mate_diff = 0;
}

bool FusionSamRec::clip_by_small_anchor(bool add_S/*JunctionSeed* filtered_junc*/)
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

			ave_intron_len = 0;

			filter_score = 0;

			ave_junc_mis = 10000;

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

bool FusionSamRec::modify_jumpcode_by_filtered_junc()
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
	}

	return true;
}
