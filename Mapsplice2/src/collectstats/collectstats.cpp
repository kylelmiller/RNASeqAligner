#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
//#include <hash_map>
#include <iomanip>
#include <cmath>
#include <iterator>
using namespace std;

void read_stats_sep_by_coulmn(string stats_file, map<string, size_t>& stored_stats)
{
	ifstream ifs(stats_file.c_str());

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;

			getline(ifs,line);

			if (line == "" || line.find("\t") == string::npos)
				continue;

			char name[1000];

			size_t count;

			sscanf(line.c_str(), "%s\t%llu", name, &count);

			cout << stats_file << "\t" << name << "\t" << count << endl;

			if (stored_stats.find(name) != stored_stats.end())
				cerr << "name already exists:" <<name;
			else
				stored_stats[name] = count;
		}
		ifs.close();
	}
	else cout << "Unable to open file:" << stats_file << endl;
}

size_t read_by_line(string infile)
{
	ifstream ifs(infile.c_str());

	size_t count = 0;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;

			getline(ifs,line);

			if (line == "" || line[0] == '@')
				continue;

			++count;
		}

		ifs.close();
	}
	else 
		cout << "Unable to open file:" << infile << endl;

	return count;
}

int main(int argc, const char** argv)
{
	if (argc < 4)
	{
		cout << "normal_stats fusion_stats unmapped_stats is_paired is_fusion"<<endl;

		exit(0);
	}

	const char* output_stats = argv[1];

	const char* normal_stats = argv[2];

	const char* fusion_stats = argv[3];

	const char* unmapped_sttats = argv[4];

	int is_paired = atoi(argv[5]);

	int is_fusion = atoi(argv[6]);

	int do_repeat = atoi(argv[7]);

	int do_annot = atoi(argv[8]);

	const char* candidate_fusion_file = argv[9];

	const char* well_annotated_fusion_file = argv[10];

	const char* not_well_annotated_fusion_file = argv[11];

	const char* circular_fusion_file = argv[12];

	map<string, size_t> normal_stats_map;

	map<string, size_t> fusion_stats_map;

	map<string, size_t> unmapped_stats_map;

	read_stats_sep_by_coulmn(normal_stats, normal_stats_map);

	read_stats_sep_by_coulmn(fusion_stats, fusion_stats_map);

	read_stats_sep_by_coulmn(unmapped_sttats, unmapped_stats_map);

	size_t paired = 0, single = 0, fusion_paired = 0, unmapped = 0;

	size_t multiple = 0, unique = 0;

	size_t unspliced = 0, spliced = 0, deletion = 0, insertion = 0, fusion_reads = 0, clipped = 0;

	size_t normal_canonical = 0, normal_semi_canonical = 0, normal_non_canonical = 0;

	size_t fusion_canonical = 0, fusion_semi_canonical = 0, fusion_non_canonical = 0;

	size_t fusion_multiple = 0, fusion_unique = 0;

	size_t normal_junc_canonical = 0, normal_junc_semi_canonical = 0, normal_junc_non_canonical = 0;

	size_t small_del_junction = 0, small_ins_junction = 0;

	size_t fusion_junc_canonical = 0, fusion_junc_semi_canonical = 0, fusion_junc_non_canonical = 0;

	size_t candidate_fusion = 0, well_annotated_fusion = 0, notwell_annotated_fusion = 0, circular_RNAs = 0;

	if (is_paired && is_fusion)
	{
		paired = normal_stats_map["paired"] + normal_stats_map["both_unspliced_paired_multiple"] + normal_stats_map["both_unspliced_paired_unique"] + fusion_stats_map["paired"];

		single = fusion_stats_map["single"];

		fusion_paired = fusion_stats_map["fusion_pair"];

		unmapped = unmapped_stats_map["unmapped_read"];

		multiple = normal_stats_map["both_unspliced_paired_multiple"] + normal_stats_map["paired_multiple"] + fusion_stats_map["paired_multiple"] + fusion_stats_map["single_multiple"] + fusion_stats_map["fusion_paired_multiple"];

		unique = normal_stats_map["both_unspliced_paired_unique"] + normal_stats_map["paired_unique"] + fusion_stats_map["paired_unique"] + fusion_stats_map["single_unique"] + fusion_stats_map["fusion_paired_unique"];

		unspliced = normal_stats_map["both_unspliced_paired_multiple"] + normal_stats_map["both_unspliced_paired_unique"] + normal_stats_map["unspliced_paired"] + fusion_stats_map["unspliced"];

		spliced = normal_stats_map["spliced_paired"] + fusion_stats_map["spliced"];

		insertion = normal_stats_map["insertion_paired"] + fusion_stats_map["insertion"];

		deletion = normal_stats_map["deletion_paired"] + fusion_stats_map["deletion"];

		fusion_reads = fusion_stats_map["fusion_reads"];

		clipped = normal_stats_map["clipped_paired"] + fusion_stats_map["clipped"];

		normal_canonical = normal_stats_map["canoical_paired"] + fusion_stats_map["canoical"];

		normal_semi_canonical = normal_stats_map["semi_canonical_paired"] + fusion_stats_map["semi_canonical"];

		normal_non_canonical = normal_stats_map["non_canonical_paired"] + fusion_stats_map["non_canonical"];

		fusion_canonical = fusion_stats_map["fusion_canonical"];

		fusion_semi_canonical = fusion_stats_map["fusion_semi_canonical"];

		fusion_non_canonical = fusion_stats_map["fusion_non_canonical"];

		fusion_multiple = fusion_stats_map["fusion_multiple"];

		fusion_unique = fusion_stats_map["fusion_unique"];

		normal_junc_canonical = normal_stats_map["canonical_junction"];

		normal_junc_semi_canonical = normal_stats_map["semi_canonical_junction"];

		normal_junc_non_canonical = normal_stats_map["non_canonical_junction"];

		small_del_junction = normal_stats_map["small_deletion_junction"];

		small_ins_junction = normal_stats_map["small_insertion_junction"];

		fusion_junc_canonical = fusion_stats_map["fusion_canon"];

		fusion_junc_semi_canonical = fusion_stats_map["fusion_semi_canon"];

		fusion_junc_non_canonical = fusion_stats_map["fusion_non_canon"];
	}
	else if (is_paired && is_fusion == 0)
	{
		paired = normal_stats_map["paired"] + normal_stats_map["both_unspliced_paired_multiple"] + normal_stats_map["both_unspliced_paired_unique"];

		single = normal_stats_map["single"] + normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["both_unspliced_single_unique"];

		fusion_paired = normal_stats_map["fusion_pair"] + normal_stats_map["both_unspliced_fusion_paired_multiple"] + normal_stats_map["both_unspliced_fusion_paired_unique"];

		unmapped = unmapped_stats_map["unmapped_read"];

		multiple = normal_stats_map["both_unspliced_multiple"] + normal_stats_map["multiple"];

		unique = normal_stats_map["both_unspliced_unique"] + normal_stats_map["unique"];

		unspliced = normal_stats_map["both_unspliced_unique"] + normal_stats_map["both_unspliced_multiple"] + normal_stats_map["unspliced"];

		spliced = normal_stats_map["spliced"];

		insertion = normal_stats_map["insertion"];

		deletion = normal_stats_map["deletion"];

		//fusion_reads = fusion_stats_map["fusion_reads"];

		clipped = normal_stats_map["clipped"];

		normal_canonical = normal_stats_map["canoical"];

		normal_semi_canonical = normal_stats_map["semi_canonical"];

		normal_non_canonical = normal_stats_map["non_canonical"];

		//fusion_canonical = fusion_stats_map["fusion_canonical"];

		//fusion_semi_canonical = fusion_stats_map["fusion_semi_canonical"];

		//fusion_non_canonical = fusion_stats_map["fusion_non_canonical"];

		//fusion_multiple = fusion_stats_map["fusion_multiple"];

		//fusion_unique = fusion_stats_map["fusion_unique"];

		normal_junc_canonical = normal_stats_map["canonical_junction"];

		normal_junc_semi_canonical = normal_stats_map["semi_canonical_junction"];

		normal_junc_non_canonical = normal_stats_map["non_canonical_junction"];

		small_del_junction = normal_stats_map["small_deletion_junction"];

		small_ins_junction = normal_stats_map["small_insertion_junction"];

		//fusion_junc_canonical = fusion_stats_map["fusion_canon"];

		//fusion_junc_semi_canonical = fusion_stats_map["fusion_semi_canon"];

		//fusion_junc_non_canonical = fusion_stats_map["fusion_non_canon"];
	}
	else if (is_paired == 0 && is_fusion)
	{
		//paired = normal_stats_map["paired"] + normal_stats_map["both_unspliced_paired_multiple"] + normal_stats_map["both_unspliced_paired_unique"] + fusion_stats_map["paired"];

		single = normal_stats_map["single"] + normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["both_unspliced_single_unique"] + fusion_stats_map["single"];

		//fusion_paired = fusion_stats_map["fusion_pair"];

		//unmapped = unmapped_stats_map["unmapped_read"];

		multiple = normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["single_multiple"] + fusion_stats_map["multiple"];

		unique = normal_stats_map["both_unspliced_single_unique"] + normal_stats_map["single_unique"] + fusion_stats_map["unique"];

		unspliced = normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["both_unspliced_single_unique"] + normal_stats_map["spliced"];

		spliced = normal_stats_map["spliced"] + fusion_stats_map["spliced"];

		insertion = normal_stats_map["insertion"] + fusion_stats_map["insertion"];

		deletion = normal_stats_map["deletion"] + fusion_stats_map["deletion"];

		fusion_reads = fusion_stats_map["fusion_reads"];

		clipped = normal_stats_map["clipped"] + fusion_stats_map["clipped"];

		normal_canonical = normal_stats_map["canoical"];

		normal_semi_canonical = normal_stats_map["semi_canonical"];

		normal_non_canonical = normal_stats_map["non_canonical"];

		fusion_canonical = fusion_stats_map["fusion_canonical"];

		fusion_semi_canonical = fusion_stats_map["fusion_semi_canonical"];

		fusion_non_canonical = fusion_stats_map["fusion_non_canonical"];

		fusion_multiple = fusion_stats_map["fusion_multiple"];

		fusion_unique = fusion_stats_map["fusion_unique"];

		normal_junc_canonical = normal_stats_map["canonical_junction"];

		normal_junc_semi_canonical = normal_stats_map["semi_canonical_junction"];

		normal_junc_non_canonical = normal_stats_map["non_canonical_junction"];

		small_del_junction = normal_stats_map["small_deletion_junction"];

		small_ins_junction = normal_stats_map["small_insertion_junction"];

		fusion_junc_canonical = fusion_stats_map["fusion_canon"];

		fusion_junc_semi_canonical = fusion_stats_map["fusion_semi_canon"];

		fusion_junc_non_canonical = fusion_stats_map["fusion_non_canon"];
	}
	else if (is_paired == 0 && is_fusion == 0)
	{
		//paired = normal_stats_map["paired"] + normal_stats_map["both_unspliced_paired_multiple"] + normal_stats_map["both_unspliced_paired_unique"] + fusion_stats_map["paired"];

		single = normal_stats_map["single"] + normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["both_unspliced_single_unique"];

		//fusion_paired = fusion_stats_map["fusion_pair"];

		//unmapped = unmapped_stats_map["unmapped_read"];

		multiple = normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["single_multiple"];

		unique = normal_stats_map["both_unspliced_single_unique"] + normal_stats_map["single_unique"];

		unspliced = normal_stats_map["both_unspliced_single_multiple"] + normal_stats_map["both_unspliced_single_unique"] + normal_stats_map["spliced"];

		spliced = normal_stats_map["spliced"];

		insertion = normal_stats_map["insertion"];

		deletion = normal_stats_map["deletion"];

		//fusion_reads = fusion_stats_map["fusion_reads"];

		clipped = normal_stats_map["clipped"];

		normal_canonical = normal_stats_map["canoical"];

		normal_semi_canonical = normal_stats_map["semi_canonical"];

		normal_non_canonical = normal_stats_map["non_canonical"];

		//fusion_canonical = fusion_stats_map["fusion_canonical"];

		//fusion_semi_canonical = fusion_stats_map["fusion_semi_canonical"];

		//fusion_non_canonical = fusion_stats_map["fusion_non_canonical"];

		//fusion_multiple = fusion_stats_map["fusion_multiple"];

		//fusion_unique = fusion_stats_map["fusion_unique"];

		normal_junc_canonical = normal_stats_map["canonical_junction"];

		normal_junc_semi_canonical = normal_stats_map["semi_canonical_junction"];

		normal_junc_non_canonical = normal_stats_map["non_canonical_junction"];

		small_del_junction = normal_stats_map["small_deletion_junction"];

		small_ins_junction = normal_stats_map["small_insertion_junction"];

		//fusion_junc_canonical = fusion_stats_map["fusion_canon"];

		//fusion_junc_semi_canonical = fusion_stats_map["fusion_semi_canon"];

		//fusion_junc_non_canonical = fusion_stats_map["fusion_non_canon"];
	}

	candidate_fusion = read_by_line(candidate_fusion_file);

	well_annotated_fusion = read_by_line(well_annotated_fusion_file);

	notwell_annotated_fusion = read_by_line(not_well_annotated_fusion_file);

	circular_RNAs = read_by_line(circular_fusion_file);

	if (is_paired == 0)
	{
		for (int i = 13; i < argc; ++i)
			unmapped += read_by_line(argv[i]);
	}

	ofstream ofs(output_stats);

	ofs << "################Read Alignment Stats####################"<< endl;

	ofs << "By mapping uniqueness"<<endl;

	size_t total = paired + single + fusion_paired + unmapped;
	
	ofs << "total_reads\t"<< total << endl;
	
	ofs << "unique\t" << unique << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total == 0 ? 0 : ((unique*100)/double(total))) << "%" <<endl;
		
	ofs << "multiple\t" << multiple << "\t" << setiosflags(ios::fixed) << setprecision(2) <<  (total == 0 ? 0 : ((multiple*100)/double(total))) << "%" << endl;

	ofs << "unmapped_reads\t"<< unmapped << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total == 0 ? 0 : ((unmapped*100)/double(total))) << "%" << endl;
	
	ofs << endl;
	
	
	ofs << "By pairing types"<<endl;

	ofs << "paired_reads\t"<< paired << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total == 0 ? 0 : ((paired*100)/double(total))) << "%" << endl;
	
	ofs << "fusion_paired_reads\t"<< fusion_paired << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total == 0 ? 0 : ((fusion_paired*100)/double(total))) << "%" << endl;

	ofs << "unpaired_reads\t"<< single << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total == 0 ? 0 : ((single*100)/double(total))) << "%" << endl;

	ofs << endl;
	
	
	ofs << "By alignment types"<<endl;

	ofs << "unspliced\t"<< unspliced << endl;

	ofs << "spliced\t"<< spliced << endl;

	//ofs << "insertion\t"<< insertion << endl;

	//ofs << "deletion\t"<< deletion << endl;

	ofs << "fusion\t"<< fusion_reads << endl;

	ofs << endl;


	/*ofs << "Normal spliced alignment types"<<endl;

	ofs << "normal_canonical\t" << normal_canonical << endl;

	ofs << "normal_semi_canonical\t" << normal_semi_canonical << endl;

	ofs << "normal_non_canonical\t" << normal_non_canonical << endl;

	ofs << endl;

	ofs << "Fusion spliced alignment types"<<endl;

	ofs << "fusion_canonical\t" << fusion_canonical << endl;

	ofs << "fusion_semi_canonical\t" << fusion_semi_canonical << endl;

	ofs << "fusion_non_canonical\t" << fusion_non_canonical << endl;

	ofs << "fusion_multiple\t" << fusion_multiple << endl;

	ofs << "fusion_unique\t" << fusion_unique << endl;

	ofs << endl;*/
	

	ofs << "################Junction Stats####################"<< endl;

	ofs << "Splice junctions"<<endl;
	
	size_t total_junction = normal_junc_canonical + normal_junc_semi_canonical + normal_junc_non_canonical;

	ofs << "total_splice_junction\t" << total_junction << endl;
	
	ofs << "canonical_splice_junction\t"<< normal_junc_canonical << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total_junction == 0 ? 0 : ((normal_junc_canonical*100)/double(total_junction))) << "%" << endl;

	ofs << "semi_canonical_splice_junction\t"<< normal_junc_semi_canonical << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total_junction == 0 ? 0 : ((normal_junc_semi_canonical*100)/double(total_junction))) << "%" << endl;

	ofs << "non_canonical_splice_junction\t"<< normal_junc_non_canonical << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total_junction == 0 ? 0 : ((normal_junc_non_canonical*100)/double(total_junction))) << "%" << endl;

	ofs << endl;
	
	
	ofs << "Indels"<< endl;

	ofs << "small_deletions\t"<< small_del_junction<< endl;

	ofs << "small_insertions\t"<< small_ins_junction<< endl;

	ofs << endl;
	

	ofs << "Fusion junctions"<<endl;
	
	size_t total_fusion = fusion_junc_canonical + fusion_junc_semi_canonical + fusion_junc_non_canonical;

	ofs << "canonical_fusion_junction\t"<< fusion_junc_canonical << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total_fusion == 0 ? 0 : ((fusion_junc_canonical*100)/double(total_fusion))) << "%" << endl;

	ofs << "semi_canonical_fusion_junction\t"<< fusion_junc_semi_canonical << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total_fusion == 0 ? 0 : ((fusion_junc_semi_canonical*100)/double(total_fusion))) << "%" << endl;

	ofs << "non_canonical_fusion_junction\t"<< fusion_junc_non_canonical << "\t" << setiosflags(ios::fixed) << setprecision(2) << (total_fusion == 0 ? 0 : ((fusion_junc_non_canonical*100)/double(total_fusion))) << "%" << endl;

	ofs << endl;
	
	ofs << "Filtered fusions"<<endl;

	ofs << "candidate_fusion\t"<< candidate_fusion<< endl;

	ofs << "well_annotated_fusion\t"<< well_annotated_fusion<< endl;

	ofs << "not_well_annotated_fusion\t"<< notwell_annotated_fusion<< endl;

	ofs << "circular_RNAs\t"<< circular_RNAs<< endl;	
}