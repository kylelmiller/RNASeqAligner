#include "SamRec.h"
#include "ReadNextTagAlignHandler.h"

size_t convert_point_coordinate(size_t cur_point, size_t x, size_t y, size_t x_prime, size_t y_prime, char strand)
{
	size_t cur_point_prime;

	if (strand == '+')
	{
		cur_point_prime = cur_point - x + x_prime;
	}
	else if (strand == '-')
	{
		cur_point_prime = y - cur_point + x_prime;
	}
	else
	{
		cerr << "strand?"<<endl;

		exit(0);
	}

	return cur_point_prime;
}


bool MergeStdFusion(vector<SamRec>::iterator sam_rec_ptr1, vector<SamRec>::iterator sam_rec_ptr2)
{
	//fusion 1
	sam_rec_ptr1->splice_way2 = sam_rec_ptr2->splice_way;

	sam_rec_ptr1->ori_splice_way2 = sam_rec_ptr2->ori_splice_way;

	sam_rec_ptr1->chrom_name2 = sam_rec_ptr2->chrom_name;

	sam_rec_ptr1->start2 = sam_rec_ptr2->start;

	sam_rec_ptr1->end2 = sam_rec_ptr2->end;

	sam_rec_ptr1->strand_t2 = sam_rec_ptr2->strand_t;

	sam_rec_ptr1->mappedlen1 = sam_rec_ptr1->mappedlen;

	sam_rec_ptr1->mappedlen2 = sam_rec_ptr2->mappedlen;

	sam_rec_ptr1->mappedlen = sam_rec_ptr1->mappedlen1 + sam_rec_ptr1->mappedlen2;

	sam_rec_ptr1->fusion_prefix_len = sam_rec_ptr1->mappedlen1;

	sam_rec_ptr1->fusion_suffix_len = sam_rec_ptr1->mappedlen2;

	//cout << "sam_rec_ptr1->fusion_prefix_len:"<<sam_rec_ptr1->fusion_prefix_len<<endl;

	//cout << "sam_rec_ptr1->fusion_suffix_len:"<<sam_rec_ptr1->fusion_suffix_len<<endl;

	sam_rec_ptr1->intron_size = sam_rec_ptr1->intron_size + sam_rec_ptr2->intron_size;

	sam_rec_ptr1->spliceway_vec2 = sam_rec_ptr2->spliceway_vec;

	sam_rec_ptr1->is_fusion = true;

	if (sam_rec_ptr1->strand_t & IS_REVERSE)
	{
		sam_rec_ptr1->strand1 = '-';
	}
	else
	{
		sam_rec_ptr1->strand1 = '+';
	}

	if (sam_rec_ptr1->strand_t2 & IS_REVERSE)
	{
		sam_rec_ptr1->strand2 = '-';
	}
	else
	{
		sam_rec_ptr1->strand2 = '+';
	}

	if (sam_rec_ptr1->strand1 == '+')
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->start;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->end;
	}
	else
	{
		sam_rec_ptr1->fusion_prefix_st = sam_rec_ptr1->end;

		sam_rec_ptr1->fusion_prefix_end = sam_rec_ptr1->start;
	}

	if (sam_rec_ptr1->strand2 == '+')
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->start2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->end2;
	}
	else
	{
		sam_rec_ptr1->fusion_suffix_st = sam_rec_ptr1->end2;

		sam_rec_ptr1->fusion_suffix_end = sam_rec_ptr1->start2;
	}

	if (sam_rec_ptr1->min_anchor > sam_rec_ptr2->min_anchor)
	{
		sam_rec_ptr1->min_anchor = sam_rec_ptr2->min_anchor;
	}

	if (sam_rec_ptr1->max_anchor < sam_rec_ptr2->max_anchor)
	{
		sam_rec_ptr1->max_anchor = sam_rec_ptr2->max_anchor;
	}

	if (sam_rec_ptr1->chrom_name < sam_rec_ptr1->chrom_name2 || (sam_rec_ptr1->chrom_name == sam_rec_ptr1->chrom_name2 && sam_rec_ptr1->start <= sam_rec_ptr1->start2))
	{
	}
	else
	{
		sam_rec_ptr1->need_swap = true;

		//sam_rec_ptr1->Swap();
	}

	return true;
}

bool comp_sam_rec(const SamRec& lhs, const SamRec& rhs)
{
	return lhs.zfzfusstr < rhs.zfzfusstr;

}

bool comp_sam_rec_by_idx(const SamRec& lhs, const SamRec& rhs)
{
	if (lhs.fusion_cur_idx == rhs.fusion_cur_idx)
		return lhs.fusion_cur_seg_idx < rhs.fusion_cur_seg_idx;
	else
		return lhs.fusion_cur_idx < rhs.fusion_cur_idx;
}

void fusiostd2oneline(pair<vector<SamRec>, vector<SamRec> >& stored_tag_aligns_pe, pair<vector<SamRec*>, vector<SamRec*> >& merged_tag_aligns_pe)
{
	vector<SamRec>::iterator align_iter;

	sort(stored_tag_aligns_pe.first.begin(), stored_tag_aligns_pe.first.end(), comp_sam_rec_by_idx);

	sort(stored_tag_aligns_pe.second.begin(), stored_tag_aligns_pe.second.end(), comp_sam_rec_by_idx);

	for (align_iter = stored_tag_aligns_pe.first.begin(); align_iter != stored_tag_aligns_pe.first.end(); ++align_iter)
	{
		if (align_iter->is_fusion_newfmt)
		{
			if (align_iter + 1 == stored_tag_aligns_pe.first.end() || (align_iter + 1)->is_fusion_newfmt == false)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			if ((align_iter + 1)->fusion_cur_seg_idx != align_iter->fusion_cur_seg_idx + 1)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			MergeStdFusion(align_iter, align_iter + 1);

			merged_tag_aligns_pe.first.push_back(&(*align_iter));

			++align_iter;
		}
		else
		{
			merged_tag_aligns_pe.first.push_back(&(*align_iter));
		}
	}

	for (align_iter = stored_tag_aligns_pe.second.begin(); align_iter != stored_tag_aligns_pe.second.end(); ++align_iter)
	{
		if (align_iter->is_fusion_newfmt)
		{
			if (align_iter + 1 == stored_tag_aligns_pe.second.end() || (align_iter + 1)->is_fusion_newfmt == false)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			if ((align_iter + 1)->fusion_cur_seg_idx != align_iter->fusion_cur_seg_idx + 1)
			{
				cout << "new format fusion alignment not together:"<<endl <<align_iter->tostring(0, 0 ) << endl;

				continue;
			}

			MergeStdFusion(align_iter, align_iter + 1);

			merged_tag_aligns_pe.second.push_back(&(*align_iter));

			++align_iter;
		}
		else
		{
			merged_tag_aligns_pe.second.push_back(&(*align_iter));
		}
	}
}

string reverse_jump_code(vector<pair<size_t, int> >& spliceway_vec)
{
	string jump_code;

	for (int i = spliceway_vec.size() - 1; i >= 0; --i)
	{
		char jump_code_chr[100];

		if (spliceway_vec[i].second > 0)
		{
			if ((i != spliceway_vec.size() - 1) && spliceway_vec[i + 1].second > 0)
			{
				int intron = spliceway_vec[i + 1].first - spliceway_vec[i].first - spliceway_vec[i].second;

				sprintf(jump_code_chr, "%dN%dM", intron, spliceway_vec[i].second);
			}
			else
			{
				sprintf(jump_code_chr, "%dM", spliceway_vec[i].second);
			}
		}
		else
		{
			int insertlen = -spliceway_vec[i].second;

			sprintf(jump_code_chr, "%dI", insertlen);
		}

		jump_code.append(jump_code_chr);
	}

	return jump_code;
}

void set_mate_pair_distance(SamRec& samrec1, SamRec& samrec2)
{
	int st1 = samrec1.start, end1 = samrec1.end, st2 = samrec2.start, end2 = samrec2.end;

	int mate_dist1, mate_dist2;

	mate_dist1 = ((st1 - end2));

	mate_dist2 = ((end1 - st2));

	int outter_dist;

	if (abs(mate_dist1) < abs(mate_dist2))
	{
		if (mate_dist2 > 0)
			outter_dist = mate_dist2 + 1;
		else
			outter_dist = mate_dist2 - 1;
	}
	else
	{
		if (mate_dist1 > 0)
			outter_dist = mate_dist1 + 1;
		else
			outter_dist = mate_dist1 - 1;
	}

	samrec1.mate_offset = samrec2.start;

	samrec1.mate_diff = -outter_dist;

	samrec1.mate_match = "=";

	if (samrec2.strand_t & IS_REVERSE)
		samrec1.strand_t |= IS_MATE_REVERSE;

	samrec1.strand_t |= IS_PAIRED_MAPPED;

	samrec1.strand_t |= IS_PAIRED;

	samrec1.strand_t |= IS_FIRST_END;
	//

	samrec2.mate_offset = samrec1.start;

	samrec2.mate_diff = outter_dist;

	samrec2.mate_match = "=";

	if (samrec1.strand_t & IS_REVERSE)
		samrec2.strand_t |= IS_MATE_REVERSE;

	samrec2.strand_t |= IS_PAIRED_MAPPED;

	samrec2.strand_t |= IS_PAIRED;

	samrec2.strand_t |= IS_FIRST_END;
}

void extract_and_convert_sam(const char* fusion_alignment_file, char* ordered_fusion_alignment_file, char* fastq_file)
{
	ofstream ofs(ordered_fusion_alignment_file);

	string fastq_1 = fastq_file; fastq_1.append(".1"); ofstream fastq_1_ofs(fastq_1.c_str());

	string fastq_2 = fastq_file; fastq_2.append(".2"); ofstream fastq_2_ofs(fastq_2.c_str());

	vector<pair<string, size_t> >::iterator chriter;

	ReadNextTagAlignHandler<SamRec> readnextalignmenthandler(fusion_alignment_file, 5);

	pair<vector<SamRec>, vector<SamRec> > stored_tag_aligns_pe;

	size_t cur_read_id;

	while(readnextalignmenthandler.ReadNextTagAlignPE(stored_tag_aligns_pe, cur_read_id))
	{
		pair<vector<SamRec*>, vector<SamRec*> > merged_tag_aligns_pe;

		fusiostd2oneline(stored_tag_aligns_pe, merged_tag_aligns_pe);

		vector<SamRec*>::iterator samrec_iter;
		/*
		if (merged_tag_aligns_pe.first.empty() && merged_tag_aligns_pe.second.size())
		{
		cout << "the other end of reads is not in sam file:" <<endl;

		cout << merged_tag_aligns_pe.second.front()->tostring(0, 0)<< endl;
		}

		if (merged_tag_aligns_pe.second.empty() && merged_tag_aligns_pe.first.size())
		{
		cout << "the other end of reads is not in sam file:" <<endl;

		cout << merged_tag_aligns_pe.first.front()->tostring(0, 0)<< endl;
		}
		*/
		for (samrec_iter = merged_tag_aligns_pe.first.begin(); samrec_iter != merged_tag_aligns_pe.first.end(); ++samrec_iter)
		{
			if (samrec_iter == merged_tag_aligns_pe.first.begin())
				fastq_1_ofs << (*samrec_iter)->tofastq()<< endl;

			if ((*samrec_iter)->is_fusion)
				ofs << (*samrec_iter)->tostandfusion()<< endl;
			else
				ofs << (*samrec_iter)->tostring(merged_tag_aligns_pe.first.size(), samrec_iter - merged_tag_aligns_pe.first.begin() + 1)<< endl;
		}

		for (samrec_iter = merged_tag_aligns_pe.second.begin(); samrec_iter != merged_tag_aligns_pe.second.end(); ++samrec_iter)
		{
			if (samrec_iter == merged_tag_aligns_pe.second.begin())
				fastq_2_ofs << (*samrec_iter)->tofastq()<< endl;

			if ((*samrec_iter)->is_fusion)
				ofs << (*samrec_iter)->tostandfusion()<< endl;
			else
				ofs << (*samrec_iter)->tostring(merged_tag_aligns_pe.second.size(), samrec_iter - merged_tag_aligns_pe.second.begin() + 1)<< endl;
		}

	}
}

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "Executable <Input_Sam_File> <Output_Sorted_Sam_File> <Output_Fastq_File> <temp_directory>" <<endl;
		exit(0);
	}

	char* fusion_alignment_file = argv[1];

	char* converted_fusion_alignment_file = argv[2];

	char* fastq_file = argv[3];
	string output = fastq_file;

	string tmpdir = ".";
	//char *tmpdir_char;
	string tmpdir_sort;


	if (argc >= 5)
		tmpdir = argv[4];
	tmpdir.append("/tmp_MPSSam2fq");
	//tmpdir_char = tmpdir.c_str();
	tmpdir_sort = tmpdir + "/tmp_sort";

	if (!tmpdir.empty())
	{
		string mktmpdircmd = "mkdir -p " + tmpdir;
		system(mktmpdircmd.c_str());
	}

	char* alignment_file = argv[1];
	string orig_file = alignment_file; 
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////check the orig_file single-end or pair-ends/////////////////
	/////////////////////////////////////////////////////////////////////////////////////        
	printf("checking input SAM file single or pair end...\n");
	char single_pair_end = 'n'; 
	char *file1 = argv[1];
	FILE *fp;
	fp = fopen(file1, "r");
	char ch[500];
	char read_name[100], flag[3], rname[10],pos[20],mapq[3],cigar[50],rnext[10];
	char pnext[20], tlen[4],seq[200],qual[200],others[500];
	char NM[20],IH[20],HI[20],YI[20],YH[20],XS[20],XF[20],ZF[50];
	while(!feof(fp))
	{   

		strcpy(ZF," ");
		fgets(ch,sizeof(ch),fp);
		if (ch[0]=='@')
			continue;
		if(feof(fp))
			break;
		sscanf(ch,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
			read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,NM,IH,HI,YI,YH,XS,XF,ZF);    

		if(ZF[0] != 'Z')
		{
			if((flag[strlen(flag)-1]=='1')||
				(flag[strlen(flag)-1]=='3')||
				(flag[strlen(flag)-1]=='5')||
				(flag[strlen(flag)-1]=='7')||
				(flag[strlen(flag)-1]=='9'))
			{printf("pair end file\n");single_pair_end = 'p';break;}

			else if((flag[strlen(flag)-1]=='2')||
				(flag[strlen(flag)-1]=='4')||
				(flag[strlen(flag)-1]=='6')||
				(flag[strlen(flag)-1]=='8')||
				(flag[strlen(flag)-1]=='0'))
			{ printf("single end file\n");single_pair_end = 's';break;}  

			else
			{printf("not sure single or pair end file!\n");  
			exit(0);}
		} 
	}
	fclose(fp);
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	if (single_pair_end == 'p') ///////////////////////////pair end
	{ 
		printf("input sam is pair end file\n");
		//printf("add123 begin\n");

		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		///////////////////////           add 123        ////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////


		FILE *fp1, *fp_add123;
		int flag_value,flag_bit[8];
		int tmp,len,sam_name_size=500;
		int i,j,k;

		char* ori_file = argv[1];
		fp1 = fopen(ori_file, "r");

		string add123_full_path = tmpdir + "/ori_add123.sam";
		fp_add123 = fopen(add123_full_path.c_str(),"w");
		//fp_add123 = fopen("ori_add123.sam", "w");

		while(!feof(fp1))
		{ 
			fgets(ch,sizeof(ch),fp1);
			if (ch[0]=='@')
				continue;
			if(feof(fp1))
				break;
			sscanf(ch,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
				read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,NM,IH,HI,YI,YH,XS,XF,ZF); 

			if (flag[2] <= '9' && flag[2] >= '0')
				flag_value = (flag[0]-'0')*100 + (flag[1]-'0')*10+(flag[2]-'0');
			else 
				flag_value = (flag[0]-'0')*10 + (flag[1]-'0');

			tmp = flag_value;

			for (k = 0; k<=7; k++)
			{  
				flag_bit[k] = tmp%2;
				tmp = (tmp - flag_bit[k])/2; 
			}      
			len = strlen(read_name);
			if (flag_bit[0] == 1 && flag_bit[6] == 1 && flag_bit[7] == 0)
			{             
				for (i=sam_name_size; i>=len+2; i--)
					ch[i]=ch[i-2];
				ch[len] = '/'; 
				ch[len+1] = '1';
				fputs(ch,fp_add123);  
			} 
			else if (flag_bit[0] == 1 && flag_bit[6] == 0 && flag_bit[7] == 1)
			{             
				for (i=sam_name_size; i>=len+2; i--)
					ch[i]=ch[i-2];

				ch[len] = '/'; 
				ch[len+1] = '2';
				fputs(ch,fp_add123);     
			} 
			else 
			{               
				for (i=sam_name_size; i>=len+2; i--)
					ch[i]=ch[i-2];

				ch[len] = '/'; 
				ch[len+1] = '3';        
				fputs(ch,fp_add123);              
			} 
		}
		fclose(fp1);
		fclose(fp_add123);


		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////     sort       ///////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////

		//printf("add123 end\n");
		string sort_tmp_add123 = "mkdir -p "+ tmpdir + "/sort_tmp_add123";

		system(sort_tmp_add123.c_str());

		//printf("sort begin\n");  
		string add123_sort_full_path = tmpdir + "/ori_add123_sorted.sam";

		string add123_sort_cmd = "LC_ALL=C sort -k1,1 -T ./" + tmpdir + "/sort_tmp_add123" + " -o ./" + add123_sort_full_path + " ./" + add123_full_path; 
		system(add123_sort_cmd.c_str()); 



		//   system("LC_ALL=C sort -k1,1 -T ./sort_tmp_add123 -o ori_add123_sorted.sam ori_add123.sam"); 
		//printf("sort end\n");

		//   system("rm ori_add123.sam");
		//   system("rm -r sort_tmp_add123");

		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////       delete          ///////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////



		//	string delete_cmd = "./delete ori_add123_sorted_delete";// + orig_file;
		//printf("delete begin\n");
		//	system(delete_cmd.c_str());
		//	printf("delete end\n");



		FILE *fp_sort, *fp_final;
		int num=0;
		int end;
		char ch_ref[500]="";

		char read_name_ref[20][100];
		char read_name_1[100];
		char read_name_2[100];
		char seq_ref[20][200];

		fp_sort = fopen(add123_sort_full_path.c_str(), "r");
		string add123_sort_delete_full_path = tmpdir + "/ori_add123_sorted_delete";
		fp_final = fopen(add123_sort_delete_full_path.c_str(), "w");

		while(!feof(fp_sort))
		{
			fgets(ch,sizeof(ch),fp_sort);
			if(feof(fp_sort))
				break;
			sscanf(ch,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,others);

			end = strlen(read_name)-1;
			if(read_name[end]=='1'||read_name[end]=='2')
			{   
				num = num +1;
				strcpy(read_name_ref[num%20],read_name);
				strcpy(seq_ref[num%20],seq);
				fputs(ch,fp_final);
			}
			else if(ch[end]=='3')
			{
				for (j=0; j<20; j++)
				{  
					if(strcmp(seq_ref[j], seq) == 0) 
					{
						sscanf(ch,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",read_name,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,others);
						strcpy(read_name_1, read_name); read_name_1[end]='1';
						strcpy(read_name_2, read_name); read_name_2[end]='2';

						if (strcmp(read_name_ref[j],read_name_1) == 0)		 
						{ 
							ch[end]='1';
							fputs(ch,fp_final);
							break;
						}

						else if(strcmp(read_name_ref[j],read_name_2) == 0)
						{
							ch[end]='2';
							fputs(ch,fp_final);
							break;
						}

					}
				}
			}
			else
			{}

		}

		fclose(fp_sort);
		fclose(fp_final);

		//printf("delete end\n");

		//system("rm ori_add123_sorted.sam");

		string sorted_sam = fusion_alignment_file;

		if (!tmpdir_sort.empty())
		{
			string mkdircmd = "mkdir -p " + tmpdir_sort;

			system(mkdircmd.c_str());

			sorted_sam = tmpdir_sort + "/tmp.sorted.sam";

			string sortcmd = "LC_ALL=C sort -k1,1 -T " + tmpdir_sort + " -o " + sorted_sam + " " + add123_sort_delete_full_path;

			system(sortcmd.c_str());
		}

		extract_and_convert_sam(sorted_sam.c_str(), converted_fusion_alignment_file, fastq_file);

		if (!tmpdir_sort.empty())
		{
			string rmcmd = "rm " + sorted_sam;

			system(rmcmd.c_str());
		}

		//  system("rm ori_add123_sorted_delete");
	}

	else if(single_pair_end == 's')///////////////////////////////////// single end
	{    
		printf("input SAM file is single end file\n");
		string sorted_sam = fusion_alignment_file;

		if (!tmpdir_sort.empty())
		{
			string mkdircmd = "mkdir -p " + tmpdir_sort;

			system(mkdircmd.c_str());

			sorted_sam = tmpdir_sort + "/tmp.sorted.sam";

			string sortcmd = "LC_ALL=C sort -k1,1 -T " + tmpdir_sort + " -o " + sorted_sam + " " + fusion_alignment_file;

			system(sortcmd.c_str());
		}

		extract_and_convert_sam(sorted_sam.c_str(), converted_fusion_alignment_file, fastq_file);

		if (!tmpdir_sort.empty())
		{
			string rmcmd = "rm " + sorted_sam;

			system(rmcmd.c_str());
		}
		string cat_cmd = "cat " + output + ".1" + " " + output + ".2"+ " > " + output;
		system(cat_cmd.c_str());
		string delete_1_cmd = "rm " + output + ".1";
		string delete_2_cmd = "rm " + output + ".2";
		system(delete_1_cmd.c_str());
		system(delete_2_cmd.c_str());
	}
	else
	{
		printf("not sure single end or pair end\n");
		exit(0);
	}

	string rm_tmp_cmd = "rm -r " + tmpdir;
	system(rm_tmp_cmd.c_str());

	printf("conversion finished\n");
}