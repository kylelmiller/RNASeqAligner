#include "AlignmentHandler.h"
//#include "JunctionAccumulator.h"

int main(int argc, char** argv)
{
	cout << sizeof(unsigned short) << endl;

	if (argc < 5)
		cout << "junction_file is_paired max_mate_dist max_hits filtered_alignment_file max_read_length ";
		cout << "chrom_dir do_filter min_ins max_del min_anchor min_junc_anchor min_mismatch do_clip mate_dist max_anchor_diff chrom_size input_sam1 input_sam2"<<endl;

#ifdef LINUX
	pthread_mutex_init(&inc_num_threads, NULL); //inc_num_threads = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_init(&fusion_lock, NULL); // fusion_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

	char* junction_file = argv[1];

	int is_paired = atoi(argv[2]);

	size_t max_mate_dist = (size_t)atoi(argv[3]);

	size_t max_hits = (size_t)atoi(argv[4]);

	char* filtered_alignment_file = argv[5];

	size_t max_read_length = (size_t) atoi(argv[6]);

	char* chrom_dir = argv[7];

	int do_filter = atoi(argv[8]);

	int min_ins = atoi(argv[9]);

	int max_del = atoi(argv[10]);

	size_t min_anchor = (size_t)atoi(argv[11]);

	size_t min_junc_anchor = (size_t) atoi(argv[12]);

	size_t min_mismatch = (size_t) atoi(argv[13]);

	int add_S = atoi(argv[14]);

	int matedistsd = atoi(argv[15]);

	mate_dist_sd = matedistsd;

	max_anchor_diff = atoi(argv[16]);

	char* chrom_size_file = (argv[17]);

	fusion_region = atoi(argv[18]);

	threads_number = atoi(argv[19]);

	int introndistsd = atoi(argv[20]);

	int min_coverage = atoi(argv[21]);

	fragment_length = strtod(argv[22], NULL);

	fragment_length_sd = strtod(argv[23], NULL);

	avearge_fragment_length = strtod(argv[24], NULL);

	boundary = atoi(argv[25]);

	min_isoform_length = atoi(argv[26]);

	min_encompass_count = atoi(argv[27]);

	char* endptr;

	min_entropy = strtod(argv[28], &endptr);
	
	low_support_threshold = atoi(argv[29]);

	if (is_paired == 0)
	{
		min_encompass_count = 0;

		min_isoform_length = min_isoform_length / 2;
	}

	intron_dist_sd = introndistsd;

	vector<string> sam_files;

	for (int i = 30; i < argc; ++i)
		sam_files.push_back(argv[i]);

	AlignmentHandler test(sam_files, junction_file, is_paired > 0, add_S > 0, max_mate_dist, max_hits, filtered_alignment_file, max_read_length, 
		chrom_dir, do_filter, min_ins, max_del, min_anchor, min_mismatch, min_junc_anchor, chrom_size_file, min_coverage);


	test.FilterAlignment();
}