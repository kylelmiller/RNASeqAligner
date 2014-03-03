#include "sharedlib.h"

void
readchrom(const char* filename, string& longseq)
{
	cout << " read chrom: " << filename << endl;
	size_t size;  

	ifstream longfile(filename);
	size = longfile.tellg();
	longfile.seekg(0);

	longseq.reserve(size);

	if (longfile.is_open())
	{
		string skipline;
		getline(longfile,skipline);

		while (!longfile.eof() )
		{
			string line;
			getline(longfile,line);

			if (line.empty())
				continue;
			if (line[strlen(line.c_str()) - 1] == '\r')
				line = line.substr(0, line.length() - 1);
			longseq.append(line);
		}
		longfile.close();
	}
	else cout << "Unable to open file";

	cout <<"chrom size:"<< longseq.size() << endl;
}

char
complement(int i) {
	static const int b2c_size = 20;
	static const char b2c[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
	};
	static const char b2cl[] = {
		//A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
		't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
	};
	if (i - 'A' >= 0 && i - 'A' < b2c_size)
		return b2c[i - 'A'];
	else if (i - 'a' >= 0 && i - 'a' < b2c_size)
		return b2cl[i - 'a'];
	else return 'N';
}

string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}

string
basename2(string filename) {

	//cout << "bef: "<<filename<<endl;
	const string s(filename);//filename.substr(0, filename.find_last_of(".")));
	size_t final_slash = s.find_last_of("/");

	if (final_slash == string::npos)
		final_slash = s.find_last_of("\\");
	if (final_slash != string::npos)
	{
		//cout << "aft 1: "<<s.substr(final_slash + 1)<<endl;
		return s.substr(final_slash + 1);
	}
	else
	{
		//cout << "aft 2: "<<s<<endl;
		return s;
	}
}

bool 
compare_pair_region(const pair<size_t, size_t>& lhs, const pair<size_t, size_t>& rhs)
{
	return lhs.first < rhs.first;
}


vector<vector<int> >* graph_ptr;
vector<size_t> DFS_stack;
vector<int> DFS_in_stack;
vector<vector<int> >* stored_path_ptr;

void DFS(vector<vector<int> >& graph, vector<vector<int> >& stored_path, size_t u)
{
	DFS_stack.clear();
	DFS_in_stack.clear();
	graph_ptr = &graph;
	stored_path_ptr = &stored_path;

	DFS_in_stack.resize(graph.size(), 0);

	DFS_VISIT(u, 1);
	//for (size_t i = 0; i < (*graph_ptr)[u].size(); ++i)
	//{
	//	//if ((*graph_ptr)[u][i] != 0)
	//}
}

void DFS_VISIT(size_t u, int path_type)
{
	DFS_in_stack[u] = path_type;

	DFS_stack.push_back(u);

	bool found_new_node = false;

	for (size_t i = 0; i < (*graph_ptr)[u].size(); ++i)
	{
		//adjacent node
		if ((*graph_ptr)[u][i] > 0)
		{
			if (DFS_in_stack[i] == 0)
			{
				found_new_node = true;

				DFS_VISIT(i, (*graph_ptr)[u][i]);
			}
		}
	}

	if (!found_new_node)//node can't find new node, print path
	{
		vector<int> cur_path;

		for(size_t i = 0; i < DFS_stack.size(); ++i)
		{
			if (DFS_in_stack[DFS_stack[i]] == 1)
				cur_path.push_back(static_cast<int> (DFS_stack[i] + 1));
			else if (DFS_in_stack[DFS_stack[i]] == 2)
			{
				if (cur_path.back() > 0)
					cur_path[cur_path.size() - 1] = -cur_path.back();

				cur_path.push_back(- (static_cast<int> (DFS_stack[i] + 1)));
			}
			else
				cout << "graph abnormal"<<endl;
			//cout << DFS_stack[i] <<'\t';
		}

		(*stored_path_ptr).push_back(cur_path);
		//cout << endl;
	}

	//pop stack
	DFS_in_stack[u] = 0;

	DFS_stack.pop_back();
}

size_t mate_dist_sd = 0;

size_t intron_dist_sd = 0;

size_t max_anchor_diff = 40;

size_t boundary = 36;

size_t fusion_region = 5000;

size_t buf_size = 1000000;

size_t threads_number = 8;

double fragment_length = 400;

double fragment_length_sd = 10;

double avearge_fragment_length = 225;

size_t global_do_filter = 0;

size_t m_min_mate_dist = 10000;

bool disable_unmapped = false;

size_t doner_side_spanning_pairs_count = 0;
size_t accetpr_side_spanning_pairs_count = 0;
size_t single_spanning_count = 0;
size_t spliceway_true_count = 0;
size_t m_fusion_encompassing_reads_doner_count = 0;
size_t m_fusion_encompassing_reads_acceptor_count = 0;

size_t min_isoform_length = 100;
size_t min_encompass_count = 1;
double min_entropy = -0.00001;

ofstream* doner_side_spanning_pairs_ofs;
ofstream* accetpr_side_spanning_pairs_ofs;
ofstream* single_spanning_ofs;
ofstream* spliceway_true_ofs;
ofstream* m_fusion_encompassing_reads_doner_ofs;
ofstream* m_fusion_encompassing_reads_acceptor_ofs;

#ifdef LINUX
pthread_mutex_t inc_num_threads; 
pthread_mutex_t fusion_lock;
#endif