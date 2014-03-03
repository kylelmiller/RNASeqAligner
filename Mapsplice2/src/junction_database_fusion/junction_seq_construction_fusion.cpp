/*    
 *    junction_seq_construction.cpp		
 *    MapSplice
 *
 *    Copyright (C) 2010 University of Kentucky and
 *                       Zeng Zheng
 *
 *    Authors: Zeng Zheng
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm> 
#include <map>
#include "junction_data_struct.h"
#include "junction_file.h"
#include "junction_sequence.h"
#include "sequence_file.h"

const int NUM_BP_PER_LINE=70;

inline bool sort_junc_by_start(const Junction* junc1, const Junction* junc2)
{
	if(junc1->start != junc2->start)
		return junc1->start < junc2->start;
	else
		return junc1->end < junc2->end;
}

inline bool sort_junc_by_end(const Junction* junc1, const Junction* junc2)
{
	if(junc1->end != junc2->end)
		return junc1->end < junc2->end;	
	else
		return junc1->start < junc2->start;
}

inline bool sort_junc_seq_by_name(const Junc_Seq& junc1, const Junc_Seq& junc2)
{
	return junc1.junc_name < junc2.junc_name;
}

inline char complement(int i) 
{
	static const int b2c_size = 20;
	static const char b2c[] = {'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'};
	static const char b2cl[] = {'t','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

inline string revcomp(const string& s) 
{
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

class Junction_Seq_Construction
{
public:
	int min_anchor;
	int max_anchor;
	int max_seq_thresh_each;
	int max_seq_thresh_total;
	int total_seq_num_for_site;
	string normal_jump_code;
	Seq_File sf;
	ofstream output_fs;

	vector<Junction> all_junction;
	map<string, vector<Junction*> > junc_sort_start;
	map<string, vector<Junction*> > junc_sort_end;

	size_t junc_count;
	size_t junc_seq_count;
	int curr_chrom_len;
	map <int, vector<Junc_Seq> > junction_head_set;
	map <int, vector<Junc_Seq> > junction_tail_set;

	Junction_Seq_Construction(int min_anchor_len, int max_anchor_len, int max_sequence_threshold_each, int max_sequence_threshold_total)
	{
		min_anchor=min_anchor_len;
		max_anchor=max_anchor_len;
		junc_count=0;
	    junc_seq_count=0;
		max_seq_thresh_each = max_sequence_threshold_each;
		max_seq_thresh_total = max_sequence_threshold_total;
		normal_jump_code.append(int_to_str(max_anchor_len) + "M" );
	}

	void load_junc_sort_table()
	{
		string curChrom = "";
		for(size_t i = 0;i < all_junction.size(); i++)
		{
			if(curChrom != all_junction[i].chrom)
			{
				junc_sort_start[all_junction[i].chrom].push_back( &(all_junction[i]) );
				junc_sort_end[all_junction[i].chrom].push_back( &(all_junction[i]) );
			}
		}
	}

	int check_link(int pos1, int pos2, int min_dis, int max_dis)
	{
		if(pos2-pos1 <= min_dis)
			return -1;
		else if(pos2-pos1 >= max_dis -1)
			return 1;
		else
			return 0;
	}

	string int_to_str(int numerical)
	{
		char c[100];
		sprintf(c,"%d",numerical);
		string str(c);
		return str;
	}

	bool push_back_seq(Fusion_Type type, Junc_Seq& my_junc_seq)
	{
		bool pushed = false;
		int head_size = 0;
		int tail_size = 0;
		if(junction_head_set.find(my_junc_seq.junc_name) != junction_head_set.end())
			head_size = (int)(junction_head_set[my_junc_seq.junc_name].size());
		if(junction_tail_set.find(my_junc_seq.junc_name) != junction_tail_set.end())
			tail_size = (int)(junction_tail_set[my_junc_seq.junc_name].size());
		if(type == START)
		{
			if(head_size < max_seq_thresh_each)
			{
				junction_head_set[my_junc_seq.junc_name].push_back(my_junc_seq);
				pushed = true;
			}
		}
		else
		{
			if(tail_size < max_seq_thresh_each)
			{	
				junction_tail_set[my_junc_seq.junc_name].push_back(my_junc_seq);
				pushed = true;
			}
		}
		return pushed;
	}

	void generate_junc_head_forward(vector<Junction*>& junction_set, Junc_Seq& my_junc_seq, size_t curr_junc)
	{
		int len_left=max_anchor - (int)my_junc_seq.sequence.length();
		Junc_Seq normal_junc_seq(my_junc_seq);  		                                        // normal path
		normal_junc_seq.start= max(junction_set[curr_junc]->start - len_left + 1, 1);
		string start_sequence=sf.get_sequence(normal_junc_seq.start, junction_set[curr_junc]->start);
		normal_junc_seq.sequence.insert(0,start_sequence);
		normal_junc_seq.jump_code.insert(0,int_to_str((int)start_sequence.length())+"M");
		if(!push_back_seq(START, normal_junc_seq))
			return;
		for(int i = (int)curr_junc - 1; i>= 0; i--)                                     // jump path
		{
			if(junction_set[i]->junc_type==FUSION)
				continue;
			int link_status=check_link(junction_set[i]->end, junction_set[curr_junc]->start, min_anchor, len_left);
			if(link_status == -1)
				continue;
			else if(link_status == 1)
				break;
			else
			{
				Junc_Seq new_junc_seq(my_junc_seq);
				string inter_sequence=sf.get_sequence(junction_set[i]->end, junction_set[curr_junc]->start);
				new_junc_seq.sequence.insert(0,inter_sequence);
				new_junc_seq.jump_code.insert(0, int_to_str((int)inter_sequence.length()) + "M" ); 
				new_junc_seq.jump_code.insert(0, int_to_str(junction_set[i]->end - junction_set[i]->start - 1) + "N" );
				generate_junc_head_forward(junction_set, new_junc_seq, i);
			}
		}	
	}

	void generate_junc_head_backward(vector<Junction*>& junction_set, Junc_Seq& my_junc_seq, size_t curr_junc)
	{
		int len_left=max_anchor - (int)my_junc_seq.sequence.length();
		Junc_Seq normal_junc_seq(my_junc_seq);  		//normal path
		int seq_end_pos=min(normal_junc_seq.end + len_left - 1, curr_chrom_len);
		string start_sequence=sf.get_sequence(normal_junc_seq.end, seq_end_pos);
		normal_junc_seq.sequence.insert(0, revcomp(start_sequence));
		normal_junc_seq.jump_code.insert(0, int_to_str((int)start_sequence.length())+"M");
		normal_junc_seq.end=seq_end_pos;
		if(!push_back_seq(START, normal_junc_seq))
			return;
		for(int i= (int)curr_junc + 1;i < (int)(junction_set.size());i++)
		{
			if(junction_set[i]->junc_type==FUSION)
				continue;
			int link_status=check_link(my_junc_seq.end, junction_set[i]->start, min_anchor, len_left);
			if(link_status == -1)
				continue;
			else if(link_status == 1)
				break;
			else
			{
				Junc_Seq new_junc_seq(my_junc_seq);
				string inter_sequence=sf.get_sequence(new_junc_seq.end , junction_set[i]->start);
				new_junc_seq.sequence.insert(0, revcomp(inter_sequence));
				new_junc_seq.jump_code.insert(0, int_to_str((int)inter_sequence.length()) + "M" ); 
				new_junc_seq.jump_code.insert(0, int_to_str(junction_set[i]->end - junction_set[i]->start - 1) + "N" );
				new_junc_seq.end+=(int)inter_sequence.length()+junction_set[i]->end - junction_set[i]->start - 1;
				generate_junc_head_backward(junction_set, new_junc_seq, i);
			}
		}	
	}

	void generate_junc_tail_forward(vector<Junction*>& junction_set, Junc_Seq& my_junc_seq, size_t curr_junc)
	{
		int len_left=max_anchor - (int)my_junc_seq.sequence.length();
		Junc_Seq normal_junc_seq(my_junc_seq);  		//normal path
		int seq_end_pos=min(junction_set[curr_junc]->end + len_left - 1, curr_chrom_len);
		string end_sequence=sf.get_sequence(junction_set[curr_junc]->end, seq_end_pos);
		normal_junc_seq.sequence.append(end_sequence);
		normal_junc_seq.jump_code.append( int_to_str((int)end_sequence.length())+"M");
		normal_junc_seq.end=seq_end_pos;
		if(!push_back_seq(END, normal_junc_seq))
			return;
		for(int i= (int)curr_junc + 1 ; i < (int)(junction_set.size());i++)
		{
			if(junction_set[i]->junc_type==FUSION)
				continue;
			int link_status=check_link(junction_set[curr_junc]->end, junction_set[i]->start, min_anchor, len_left);
			if(link_status == -1)       // too close
				continue;
			else if(link_status == 1)   // too far
				break;
			else                        // fit
			{
				Junc_Seq new_junc_seq(my_junc_seq);
				string inter_sequence=sf.get_sequence(junction_set[curr_junc]->end, junction_set[i]->start);
				new_junc_seq.sequence.append(inter_sequence);
				new_junc_seq.jump_code.append( int_to_str((int)inter_sequence.length()) + "M" ); 
				new_junc_seq.jump_code.append(int_to_str(junction_set[i]->end - junction_set[i]->start - 1) + "N" );
				generate_junc_tail_forward(junction_set, new_junc_seq, i);
			}
		}	
	}

	void generate_junc_tail_backward(vector<Junction*>& junction_set, Junc_Seq& my_junc_seq, size_t curr_junc)
	{
		int len_left=max_anchor - (int)my_junc_seq.sequence.length();
		Junc_Seq normal_junc_seq(my_junc_seq);  		//normal path
		int seq_start_pos =  max(normal_junc_seq.start - len_left + 1, 1);
		string end_sequence=sf.get_sequence(seq_start_pos, normal_junc_seq.start);
		normal_junc_seq.sequence.append(revcomp(end_sequence));
		normal_junc_seq.jump_code.append(int_to_str((int)end_sequence.length())+"M");
		normal_junc_seq.start=seq_start_pos;
		if(!push_back_seq(END, normal_junc_seq))
			return;
		for(int i = (int)curr_junc - 1; i >= 0; i--)
		{
			if(junction_set[i]->junc_type==FUSION)
				continue;
			int link_status=check_link(junction_set[i]->end, my_junc_seq.start, min_anchor, len_left);
			if(link_status == -1)       // too close
				continue;
			else if(link_status == 1)   // too far
				break;
			else                        // fit
			{
				Junc_Seq new_junc_seq(my_junc_seq);
				string inter_sequence=sf.get_sequence(junction_set[i]->end, new_junc_seq.start);
				new_junc_seq.sequence.append(revcomp(inter_sequence));
				new_junc_seq.jump_code.append(int_to_str((int)inter_sequence.length()) + "M" ); 
				new_junc_seq.jump_code.append(int_to_str(junction_set[i]->end - junction_set[i]->start - 1) + "N" );
				new_junc_seq.start-=(int)inter_sequence.length()+junction_set[i]->end - junction_set[i]->start - 1;
				generate_junc_tail_backward(junction_set, new_junc_seq, i);
			}
		}	
	}

	void output_junc_seq()
	{
		map <int, vector<Junc_Seq> >::iterator it_head, it_tail;
		for(it_head = junction_head_set.begin(); it_head != junction_head_set.end(); it_head ++)
		{
			it_tail = junction_tail_set.find(it_head->first);
			int head_bound = (int)(it_head->second.size());
			int tail_bound = (int)(it_tail->second.size());
			bool exceed_max_total = false;
			if(head_bound * tail_bound > max_seq_thresh_total)
			{
				head_bound = 1;
				tail_bound = 1;
				exceed_max_total = true;
			}
			for(int i = 0; i < head_bound; i ++)
			{
				for(int j = 0; j < tail_bound; j ++)
				{
					int start1 = (it_head->second[i].strand == "+" ? it_head->second[i].start : it_head->second[i].end);
					int start2 = (it_tail->second[j].strand == "+" ? it_tail->second[j].start : it_tail->second[j].end);
					output_fs<<">"<<it_head->second[i].chrom<<"_"<<it_head->second[i].strand<<"_"<<start1<<":"<<it_head->second[i].jump_code<<";"<<it_tail->second[j].chrom<<"_"<<it_tail->second[j].strand<<"_"<<start2<<":"<<it_tail->second[j].jump_code<<endl;	
					string whole_sequence=it_head->second[i].sequence + it_tail->second[j].sequence;
					size_t start = 0;
					size_t end= start + NUM_BP_PER_LINE;
					while(true)
					{
						if(end<=whole_sequence.length())
						{
							output_fs<<whole_sequence.substr(start,NUM_BP_PER_LINE)<<endl;
							if(end==whole_sequence.length())
								break;
						}
						else
						{
							output_fs<<whole_sequence.substr(start,whole_sequence.length() - start)<<endl;
							break;
						}
						start=end;
						end=start + NUM_BP_PER_LINE;
					}
				}
			}
		}
	}

	void construct_junc_seq(char* junc_file_normal, bool normal_header, char* junc_file_fusion, bool fusion_header, char* refseq_path, char* output_file)
	{
		output_fs.open(output_file);
		if( !output_fs ) 
		{
			fprintf(stderr,"error: write junction sequence file error\n");exit(1);
		} 		
		cout << "loading normal junction file " << junc_file_normal <<endl;
		Junction_File junc_normal(junc_file_normal, normal_header);
		int normal_junc_end_id = junc_normal.load_all_normal_junctions(all_junction); //load all junctions
		cout << normal_junc_end_id - 1 << " normal junction loaded" << endl;
		junc_normal.close();
		cout << "loading fusion junction file " << junc_file_fusion <<endl;
		Junction_File junc_fusion(junc_file_fusion, fusion_header, normal_junc_end_id);
		int fusion_junc_end_id = junc_fusion.load_all_fusion_junctions(all_junction); //load all junctions
		cout << (fusion_junc_end_id - normal_junc_end_id) << " fusion junction loaded" << endl;
		junc_fusion.close();
		load_junc_sort_table();
		for(map<string, vector<Junction*> >::iterator it = junc_sort_start.begin(); it != junc_sort_start.end(); it++)
		{
			sf.init(refseq_path + it->first + ".fa");
			sf.load_next_chrom_seq();
			cout<< it->first << " sequence load complete"<<endl;
			curr_chrom_len=sf.get_seq_len();
			sort(junc_sort_start[it->first].begin(), junc_sort_start[it->first].end(), sort_junc_by_start);
			sort(junc_sort_end[it->first].begin(), junc_sort_end[it->first].end(), sort_junc_by_end);
			for(size_t i = 0; i < junc_sort_start[it->first].size(); i++)
			{
				if(junc_sort_start[it->first][i]->junc_type != FUSION)
					continue;
				if(junc_sort_start[it->first][i]->fusion_type == START && junc_sort_start[it->first][i]->strand == "-")
				{
					Junc_Seq head_junc_seq;   
					head_junc_seq.junc_name = junc_sort_start[it->first][i]->junc_name;
					head_junc_seq.strand = junc_sort_start[it->first][i]->strand;
					head_junc_seq.chrom = junc_sort_start[it->first][i]->chrom;
					head_junc_seq.start = junc_sort_start[it->first][i]->start;
					head_junc_seq.end = junc_sort_start[it->first][i]->start;
					generate_junc_head_backward(junc_sort_start[it->first], head_junc_seq, i);
				}
				if(junc_sort_start[it->first][i]->fusion_type == END && junc_sort_start[it->first][i]->strand == "+")
				{
					Junc_Seq tail_junc_seq;   
					tail_junc_seq.junc_name = junc_sort_start[it->first][i]->junc_name;
					tail_junc_seq.strand = junc_sort_start[it->first][i]->strand;
					tail_junc_seq.chrom = junc_sort_start[it->first][i]->chrom;
					tail_junc_seq.start =junc_sort_start[it->first][i]->end;
					generate_junc_tail_forward(junc_sort_start[it->first], tail_junc_seq, i);
				}
			}
			for(size_t i = 0; i < junc_sort_end[it->first].size(); i++)
			{
				if(junc_sort_end[it->first][i]->junc_type != FUSION)
					continue;
				if(junc_sort_end[it->first][i]->fusion_type == START && junc_sort_end[it->first][i]->strand == "+")
				{
					Junc_Seq head_junc_seq;   
					head_junc_seq.junc_name = junc_sort_end[it->first][i]->junc_name;
					head_junc_seq.strand = junc_sort_end[it->first][i]->strand;
					head_junc_seq.chrom = junc_sort_end[it->first][i]->chrom;
					head_junc_seq.end = junc_sort_end[it->first][i]->start;
					generate_junc_head_forward(junc_sort_end[it->first], head_junc_seq, i);
				}
				if(junc_sort_end[it->first][i]->fusion_type == END && junc_sort_end[it->first][i]->strand == "-")
				{
					Junc_Seq tail_junc_seq;   
					tail_junc_seq.junc_name = junc_sort_end[it->first][i]->junc_name;
					tail_junc_seq.strand = junc_sort_end[it->first][i]->strand;
					tail_junc_seq.chrom = junc_sort_end[it->first][i]->chrom;
					tail_junc_seq.start =junc_sort_end[it->first][i]->end;
					tail_junc_seq.end = junc_sort_end[it->first][i]->end;
					generate_junc_tail_backward(junc_sort_end[it->first], tail_junc_seq, i);
				}
			}
		}
		sf.close();
		cout << "output synthetic sequence" << endl;
		output_junc_seq();
		output_fs.close();
	}
};

void print_usage()
{
	fprintf(stderr,"junc_db_fusion [min_anchor] [max_anchor] [max_sequence_threshold_for_each_end] [max_sequence_threshold_total] [normal_junction_file] [normal_junction_file_header 0/1] [fusion_junction_file] [fusion_junction_file_header 0/1] [sequence_path] [output_file]\n");
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc ==1 )
		print_usage();
	else if (argc < 11)
	{
		fprintf(stderr,"error: too few arguments\n");
		print_usage();
		exit(1);
	} 
	else if (argc > 11)
	{
		fprintf(stderr,"error: too many arguments\n");
		print_usage();
		exit(1);
	}
	bool normal_header = (atoi(argv[6]) != 0);
	bool fusion_header = (atoi(argv[8]) != 0);
	Junction_Seq_Construction my_construct(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
	my_construct.construct_junc_seq(argv[5], normal_header, argv[7], fusion_header, argv[9],argv[10]);
}


