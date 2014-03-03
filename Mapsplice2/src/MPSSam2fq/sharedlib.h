#ifndef SHAREDLIB_H
#define SHAREDLIB_H

#define _CRT_SECURE_NO_WARNINGS
//#define VS
//#define LINUX

//#define DEBUG

#include <iostream>
#include <vector>

#include <string>

#ifdef VS
#include <hash_map> //vc only
#else
#include <ext/hash_map> //g++ only
#endif

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>
//#include <dirent.h>
#include <iomanip>
#include <map>
#include <set>
#include <queue>
#include <list>

#include <cmath>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <iterator>
#include <time.h>

#include <sstream>

//#include <pthread.h>
//#include <signal.h>
//#include <unistd.h>

//#include "agbl/specialfunctions.h"


using namespace std;

#ifdef VS
using namespace stdext;
#endif

#ifndef VS
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
#endif

#ifndef VS
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;

namespace __gnu_cxx
{
	template<class Traits, class Allocator>
	struct hash<std::basic_string<char, Traits, Allocator> >
	{
		size_t operator()(const std::basic_string<char, Traits, Allocator>& __s) const
		{
			return __stl_hash_string(__s.c_str());
		}
	};
}

#endif

#define IS_PAIRED 0x0001
#define IS_PAIRED_MAPPED 0x0002
#define IS_UNMAPPED 0x0004
#define MATE_UNMAPPED 0x0008
#define IS_REVERSE 0x0010
#define IS_MATE_REVERSE 0x0020
#define IS_FIRST_END 0x040
#define IS_SECOND_END 0x0080
#define IS_PRIMARY 0x0100
#define IS_FAILED_QUAL_CHECK 0x0200
#define IS_PCR_DUP 0x0400

#define THIRTY_TWO 32
const size_t ALL_BITS_ON = static_cast<size_t>(-1);
const size_t LOWER_THIRTY_TWO_MASK = ALL_BITS_ON >> THIRTY_TWO;
const size_t UPPER_THIRTY_TWO_MASK = LOWER_THIRTY_TWO_MASK << THIRTY_TWO;

void readchrom(const char* filename, string& longseq);

char complement(int i);

string revcomp(const string& s);

bool compare_pair_region(const pair<size_t, size_t>& lhs, const pair<size_t, size_t>& rhs);

string	basename2(string filename);

enum FILTERED_TYPE
{
	NOT_FILTERED,
	FILTERED_BY_SMALL_ANCHOR,
	FILTERED_BY_SMALL_DELETION,
	FILTERED_BY_LARGE_MULTIPLE_PAIRED,
	FILTERED_BY_LARGE_MIN_ANCHOR_DIFF,
	FILTERED_BY_UNBALANCED_LEFT_RIGHT_PAIR,
	FILTERED_BY_NOPAIRED,
	FILTERED_BY_INSERTION,
	FILTERED_BY_LARGE_MISMATCH,
	FILTERED_BY_NONCAN_ERROR,
	FILTERED_BY_NONCAN_ENTROPY,
	FILTERED_BY_NONCAN_MULTI,
	FILTERED_BY_NONCAN_LEFT_RIGHT_PAIR,
	FILTERED_BY_CAN_ENTROPY,
	FILTERED_BY_ENTROPY,

	FILTERED_BY_FUSION_SMALL_ANCHOR,
	FILTERED_BY_FUSION_SMALL_DELETION,
	FILTERED_BY_FUSION_LARGE_MULTIPLE_PAIRED,
	FILTERED_BY_FUSION_LARGE_MIN_ANCHOR_DIFF,
	FILTERED_BY_FUSION_UNBALANCED_LEFT_RIGHT_PAIR,
	FILTERED_BY_FUSION_NOPAIRED,
	FILTERED_BY_FUSION_INSERTION,
	FILTERED_BY_FUSION_LARGE_MISMATCH,
	FILTERED_BY_FUSION_NONCAN_ERROR,
	FILTERED_BY_FUSION_NONCAN_ENTROPY,
	FILTERED_BY_FUSION_NONCAN_MULTI,
	FILTERED_BY_FUSION_NONCAN_LEFT_RIGHT_PAIR,
	FILTERED_BY_FUSION_CAN_ENTROPY,
	FILTERED_BY_FUSION_LOW_COVERAGE,
	FILTERED_BY_FUSION_ISOLATED_EXON,
	FILTERED_BY_FUSION_NO_EXON,

	FILTERED_BY_FUSION_LARGE_MIN_MISMATCH,
	FILTERED_BY_FUSION_SMALL_ENTROPY,
};

enum PAIRED_TYPE
{
	NORMAL_PAIRED,
	FUSION_PAIRED,
	SINGLE,
	UNMAPPED,
};

//extern class SamRec;

extern size_t mate_dist_sd;
extern size_t intron_dist_sd;
extern size_t max_anchor_diff;
extern size_t boundary;
extern size_t fusion_region;
extern size_t buf_size;
extern size_t threads_number;
extern double fragment_length;
extern double fragment_length_sd;
extern double avearge_fragment_length;
extern size_t global_do_filter;
extern size_t min_isoform_length;
extern size_t min_encompass_count;

extern double min_entropy;

extern bool disable_unmapped;


extern vector<vector<int> >* graph_ptr;
extern vector<size_t> DFS_stack;
extern vector<int> DFS_in_stack;
extern vector<vector<int> >* stored_path_ptr;

extern void DFS_VISIT(size_t u, int path_type);

extern void DFS(vector<vector<int> >& graph, vector<vector<int> >& stored_path, size_t u);

extern void DFS_VISIT(size_t u, int path_type);

#ifdef LINUX
extern pthread_mutex_t inc_num_threads;
extern pthread_mutex_t fusion_lock;
#endif

extern size_t doner_side_spanning_pairs_count;
extern size_t accetpr_side_spanning_pairs_count;
extern size_t single_spanning_count;
extern size_t spliceway_true_count;
extern size_t m_fusion_encompassing_reads_doner_count;
extern size_t m_fusion_encompassing_reads_acceptor_count;
extern size_t m_min_mate_dist;

extern ofstream* doner_side_spanning_pairs_ofs;
extern ofstream* accetpr_side_spanning_pairs_ofs;
extern ofstream* single_spanning_ofs;
extern ofstream* spliceway_true_ofs;
extern ofstream* m_fusion_encompassing_reads_doner_ofs;
extern ofstream* m_fusion_encompassing_reads_acceptor_ofs;

#endif
//static const string Is(500, 'I');
