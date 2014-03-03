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

size_t mate_dist_sd = 0;

size_t intron_dist_sd = 0;

size_t max_anchor_diff = 40;

size_t boundary = 36;

size_t fusion_region = 5000;

size_t buf_size = 1000000;
