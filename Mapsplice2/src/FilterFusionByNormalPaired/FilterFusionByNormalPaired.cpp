#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
//#include <hash_map>
#include <cmath>
#include <iterator>
using namespace std;

int keepcircular = 0;

inline char
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

inline string
revcomp(const string& s) {
	string r;
	transform(s.begin(), s.end(), back_inserter(r), complement);
	reverse(r.begin(), r.end());
	return r;
}


//void seperatemapreads(const char* mapreadfile, map<string, pair<int, string> >& readsmap);

void
countmapreads(const char* filename, set<size_t >& readsmap)
{
	ifstream ifs(filename);
	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "" || line[0] == '@')
				continue;

			char readname[1000];

			size_t junc_id;

			sscanf(line.c_str(), "%llu", &junc_id);

			readsmap.insert(junc_id);

		}

		ifs.close();

	}
	else cout << "Unable to open file: " << filename << endl;
}

void
countrepeatreads(const char* filename, set<size_t >& readsmap)
{
	ifstream ifs(filename);

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "" || line[0] == '@' || line[0] != '>')
				continue;

			char readname[1000];

			size_t junc_id;

			sscanf(line.c_str() + 1, "%llu", &junc_id);

			readsmap.insert(junc_id);

		}

		ifs.close();

	}
	else cout << "Unable to open file: " <<filename<<endl ;
}

void verifyreads(const char* readsfile, set<size_t>& readsmap, set<size_t>& readsmaplong, const char* filtered_junc, const char* filtered_junc_filtered)
{
	ofstream ofs(filtered_junc);

	ofstream ofs_filtered((string(filtered_junc_filtered)).c_str());

	ifstream ifs(readsfile);

	size_t count = 0;

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line;

			getline(ifs, line);

			if (line.empty())
				continue;

			char chromname[1000], strand1, strand2;

			char from[1000], type[1000], genestrand[100];

			char gene1[1000], gene2[1000];

			size_t start, end;

			char skip[1000];

			//chr10~chr10	73579521	73579505	JUNC_5	5	-+	255,0,0	2	
			//22,33,22,154,	0,39,	1.054920	6	GTAG	0	0	0.000000	22	0	6	3	2	2	2	0	0	2	3	30	
			//73579543	73579538	73579505,153M339N21M|	73579521,137M339N96M1540N13M|	0	0	0.465621	0.439167	
			//174	174	246	246	not_matched	not_matched	AGGAGGTCAGCCCTGAGCTGGT	TGAGCTGGTGTGCAGCATGCTGCAC

			sscanf(line.c_str(), "%s\t%llu\t%llu\t%s\t%s\t%c%c", chromname, &start, &end, skip, skip, &strand1, &strand2);

			string chrstr = chromname;

			size_t splitidx = chrstr.find("~");

			string chr1 = chrstr.substr(0, splitidx);

			string chr2 = chrstr.substr(splitidx + 1, chrstr.length() - splitidx - 1);

			bool iscircular = false;

			if (chr1 == chr2 && strand1 == strand2)
			{
				if (strand1 == '+' && start > end)
					iscircular = true;
				else if (strand1 == '-' && end > start)
					iscircular = true;
			}

			++count;

			if (readsmaplong.find(count) != readsmaplong.end())
				ofs_filtered << line << endl;
			else if (readsmap.find(count) == readsmap.end() || (iscircular && (line.find("normal") != string::npos || keepcircular)))
				ofs << line << endl;
			else
				ofs_filtered << line << endl;
		}
		ifs.close();
	}
	else cout << "Unable to open file: "<<readsfile << endl;
}

int
main(int argc, const char** argv)
{
	if (argc < 2)
		return 0;

	set<size_t > readsmap;

	set<size_t > readsmaplong;

	countmapreads(argv[1], readsmap);

	keepcircular = atoi(argv[4]);

	if (argc > 5)
		countrepeatreads(argv[5], readsmap);

	if (argc > 6)
		countrepeatreads(argv[6], readsmaplong);

	const char* fusion_filtered_junc;

	string filtered;

	if (argc > 7)
		fusion_filtered_junc = argv[7];
	else
	{
		filtered = argv[3];

		filtered.append(".filtered");

		fusion_filtered_junc = filtered.c_str();
	}

	verifyreads(argv[2], readsmap, readsmaplong, argv[3], fusion_filtered_junc);

	return 0;
}