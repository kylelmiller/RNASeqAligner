#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
//#include <hash_map>
#include <cmath>
#include <iterator>
using namespace std;


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


void
countmapreads(const char* filename, const char* normal_alignments, map<string, size_t >& readsmap)
{
	ifstream ifs(filename);

	ofstream ofs((string(normal_alignments) + ".normal").c_str());

	if (ifs.is_open())
	{
		//string skipline;
		//getline(ifs,skipline);

		while (!ifs.eof() )
		{
			string line;
			getline(ifs,line);
			if (line == "" || line[0] == '@' || line.find("ZF:Z:FUS") != string::npos)
				continue;

			char /*chromname[1000],*/ readname[1000]/*, flankseq[10], strand[10],readseq[1000]*/;
			//int prim;

			sscanf(line.c_str(), "%s", /*chromname, */readname/*, &prim, flankseq, strand, readseq*/);

			ofs << line << endl;

			string readnamestr = readname;

			size_t readid;

			if (readnamestr[readnamestr.length() - 1] == '1')
				readid = 1;
			else if (readnamestr[readnamestr.length() - 1] == '2')
				readid = 2;
			else
			{
				cout << line << endl;
				continue;
			}

			readnamestr = readnamestr.substr(0, readnamestr.length() - 1);

			map<string, size_t >::iterator rit = readsmap.find(readnamestr);

			if (rit == readsmap.end())
			{
				readsmap.insert(make_pair(readnamestr, readid));

				//cout << "readid:"<< readid <<endl;
			}
			else
			{
				rit->second = rit->second | readid;

				//cout << "rit->second:"<< rit->second << endl;
			}
		}
		ifs.close();
	}
	else cout << "Unable to open file";
	
	cout<<"total\t"<<readsmap.size()<<endl;

}

void verifyreads(const char* readsfile, map<string, size_t>& readsmap)
{
	ifstream ifs(readsfile);

	//cout << readsmap.begin()->first << endl;

	ofstream ofs(( string(readsfile) + ".unmapped").c_str());

	if (ifs.is_open())
	{
		while (!ifs.eof() )
		{
			string line, line2, line3, line4;

			getline(ifs,line);
			getline(ifs,line2);
			getline(ifs,line3);
			getline(ifs,line4);

			if (line == "")
				continue;

			string readname = line.substr(1, line.length() - 2);

			map<string, size_t>::iterator rit = readsmap.find(readname);

			if (rit == readsmap.end() || rit->second != 3)
			{
				//cout << readname << endl;
				//cout << rit->second << endl;
				ofs << line << endl;
				ofs << line2 << endl;
				ofs << line3 << endl;
				ofs << line4 << endl;
			}
		}
		ifs.close();
	}
	else cout << "Unable to open file";
}

int
main(int argc, const char** argv)
{
	if (argc < 2)
		return 0;

	map<string, size_t > readsmap;

	countmapreads(argv[1], argv[2], readsmap);

	for (size_t i = 3; i < argc; ++i)
	{
		verifyreads(argv[i], readsmap);

	}
	return 0;
}