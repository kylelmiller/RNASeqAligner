#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <errno.h>

using namespace std;

size_t minlen = 10000, maxlen = 0;

size_t minqual = -1;

size_t maxqual = 0;

size_t readlen = 0;

int countall=0;

void convert_ex_paired_sep(char* extendFile, char* extendFile2, int fasta, int min_read_len)
{
	size_t read_id = 0;

	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	ifstream input_fs2(extendFile2);
	if( !input_fs2 ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}   

	string readidline;
	string line, line2, line3, line4;
	string line5, line6, line7, line8;
	
	while(!input_fs.eof() && !input_fs2.eof())
	{
		getline(input_fs, line);
		getline(input_fs, line2);

		if (line.empty() || line2.empty())
			continue;

		++countall;

		//output_fs <<line[0] <<countall<<'~'<<  line.c_str() + 1 << endl << line2 << endl;

		if (!fasta)
		{
			getline(input_fs, line3);	
			getline(input_fs, line4);

			if (line3.empty() || line4.empty())
				continue;

			for (size_t i = 0; i < line4.size(); ++i)
			{
				char qual = line4[i];

				if (minqual > qual)
					minqual = qual;

				if (maxqual < qual)
					maxqual = qual;
			}
			//output_fs << line3 << endl << line4 << endl;
		}

		getline(input_fs2, line5);	
		getline(input_fs2, line6);

		if (line5.empty() || line6.empty())
			continue;

		++countall;

		//output_fs <<line[0]  <<countall<<'~'<< line5.c_str() + 1 << endl << line6 << endl;

		if (!fasta)
		{
			getline(input_fs2, line7);	
			getline(input_fs2, line8);

			if (line7.empty() || line8.empty())
				continue;

			for (size_t i = 0; i < line8.size(); ++i)
			{
				char qual = line8[i];

				if (minqual > qual)
					minqual = qual;

				if (maxqual < qual)
					maxqual = qual;
			}

			//output_fs << line7 << endl << line8 << endl;
		}

		++read_id;

		/*if (line.find(" ") != string::npos || line.find("\t") != string::npos)
		{
			cerr << "read name contain blank space or tab"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		/*if (line.substr(line.length() - 2, 2) != "/1" && line.substr(line.length() - 2, 2) != "/2" && line.substr(line.length() - 1, 1) != "a" && line.substr(line.length() - 1, 1) != "b")
		{
			cerr << "pairend read name not end with /1 or /2"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		/*if (line5.find(" ") != string::npos || line5.find("\t") != string::npos)
		{
			cerr << "read name contain blank space or tab"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		/*if (line5.substr(line5.length() - 2, 2) != "/1" && line5.substr(line5.length() - 2, 2) != "/2" && line5.substr(line5.length() - 1, 1) != "a" && line5.substr(line5.length() - 1, 1) != "b")
		{
			cerr << "pairend read name not end with /1 or /2"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile2 << endl;
			cerr << line5<<endl;

			exit(1);
		}*/		

		/*if (line.substr(0, line.length() - 1) != line5.substr(0, line5.length() - 1))
		{
			cerr << "base name of two ends not consistent"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< " and " << extendFile2 << endl;
			cerr << line << endl << line5<<endl;

			exit(1);
		}*/
		
		string read_name1;
		string read_name2;
		size_t space_index1 = line.find(" ");
		size_t space_index2 = line5.find(" ");
		if(space_index1 == string::npos)
			read_name1 = line;
		else
			read_name1 = line.substr(0, space_index1);
		if(space_index2 == string::npos)
			read_name2 = line5;
		else
			read_name2 = line5.substr(0, space_index2);
		if(read_name1.substr(read_name1.length() - 2, 2) == "/1" && read_name2.substr(read_name2.length() -2 , 2) == "/2")
		{
			if(read_name1.substr(0, read_name1.length() - 1) != read_name2.substr(0, read_name2.length() - 1))
			{
				cerr << "Base name of two ends not consistent"<<endl;
				cerr << "The " << read_id << "th read in "<< extendFile<< " and " << extendFile2 << endl;
				cerr << line << endl << line5<<endl;
				exit(1);	
			}
		}
		else if(read_name1.substr(read_name1.length() - 2, 2) != "/1" && read_name2.substr(read_name2.length() -2 , 2) != "/2")
		{
			if(read_name1 != read_name2)
			{
				cerr << "Base name of two ends not consistent"<<endl;
				cerr << "The " << read_id << "th read in "<< extendFile<< " and " << extendFile2 << endl;
				cerr << line << endl << line5<<endl;
				exit(1);	
			}
		}
		else
		{
			cerr << "Base name of two ends not consistent"<<endl;
			cerr << "The " << read_id << "th read in "<< extendFile<< " and " << extendFile2 << endl;
			cerr << line << endl << line5<<endl;
			exit(1);
		}
		
		if (line2.length() < readlen)
		{
			cerr << "Read length < specified read length"<<endl;
			cerr << "The " << read_id << "th read in " << extendFile << endl;
			cerr << line << endl << line2 << endl << line3 << endl << line4 <<endl;

			exit(2);
		}

		if (line6.length() < readlen)
		{
			cerr << "Read length < specified read length"<<endl;
			cerr << "The " << read_id << "th read in " << extendFile2 << endl;
			cerr << line5 << endl << line6 << endl << line7 << endl << line8<<endl;
			exit(2);
		}

		if (!fasta)
		{
			if (line2.length() != line4.length())
			{
				cerr << "Read length and quality string length not consistent"<<endl;
				cerr << "The " << read_id << "th read in " << extendFile << endl;
				cerr << line << endl << line2 << endl << line3 << endl << line4 <<endl;

				exit(3);
			}

			if (line6.length() != line8.length())
			{
				cerr << "Read length and quality string length not consistent"<<endl;
				cerr << "The " << read_id << "th read in " << extendFile2 << endl;
				cerr << line5 << endl << line6 << endl << line7 << endl << line8<<endl;

				exit(3);
			}
		}

		if (line2.length() > maxlen)
			maxlen = line2.length();

		if (line2.length() < minlen)
			minlen = line2.length();

		if (line6.length() > maxlen)
			maxlen = line6.length();

		if (line6.length() < minlen)
			minlen = line6.length();
	}

	input_fs.close();
	
	//cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

void convert_ex_paired_comb(char* extendFile, int fasta, int min_read_len)
{
	size_t read_id = 0;

	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	string readidline;
	string line, line2, line3, line4;
	string line5, line6, line7, line8;
	
	while(!input_fs.eof())
	{
		
		getline(input_fs, line);
		getline(input_fs, line2);

		if (line.empty() || line2.empty())
			continue;

		++countall;

		//output_fs <<line[0] <<countall<<'~'<<  line.c_str() + 1 << endl << line2 << endl;

		if (!fasta)
		{
			getline(input_fs, line3);	
			getline(input_fs, line4);

			if (line3.empty() || line4.empty())
				continue;

			for (size_t i = 0; i < line4.size(); ++i)
			{
				char qual = line4[i];

				if (minqual > qual)
					minqual = qual;

				if (maxqual < qual)
					maxqual = qual;
			}

			//output_fs << line3 << endl << line4 << endl;
		}

		

		getline(input_fs, line5);	
		getline(input_fs, line6);

		if (line5.empty() || line6.empty())
			continue;

		++countall;

		//output_fs <<line[0]  <<countall<<'~'<< line5.c_str() + 1 << endl << line6 << endl;

		if (!fasta)
		{
			getline(input_fs, line7);	
			getline(input_fs, line8);

			if (line7.empty() || line8.empty())
				continue;

			for (size_t i = 0; i < line8.size(); ++i)
			{
				char qual = line8[i];

				if (minqual > qual)
					minqual = qual;

				if (maxqual < qual)
					maxqual = qual;
			}

			//output_fs << line7 << endl << line8 << endl;
		}

		++read_id;

		/*if (line.find(" ") != string::npos || line.find("\t") != string::npos)
		{
			cerr << "read name contain blank space or tab"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		/*if (line.substr(line.length() - 2, 2) != "/1" && line.substr(line.length() - 2, 2) != "/2" && line.substr(line.length() - 1, 1) != "a" && line.substr(line.length() - 1, 1) != "b")
		{
			cerr << "pairend read name not end with /1 or /2"<<endl;
			cerr << "the " << read_id * 2 - 1 << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		/*if (line5.find(" ") != string::npos || line5.find("\t") != string::npos)
		{
			cerr << "read name contain blank space or tab"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		/*if (line5.substr(line5.length() - 2, 2) != "/1" && line5.substr(line5.length() - 2, 2) != "/2" && line5.substr(line5.length() - 1, 1) != "a" && line5.substr(line5.length() - 1, 1) != "b")
		{
			cerr << "pairend read name not end with /1 or /2"<<endl;
			cerr << "the " <<  read_id * 2 << "th read in "<< extendFile << endl;
			cerr << line5<<endl;

			exit(1);
		}*/

		/*if (line.substr(0, line.length() - 1) != line5.substr(0, line5.length() - 1))
		{
			cerr << "base name of two ends not consistent"<<endl;
			cerr << "the " << read_id * 2 - 1 << "th read in "<< extendFile<< " and " << read_id * 2 << "th read in" << extendFile << endl;
			cerr << line << endl << line5<<endl;

			exit(1);
		}*/

		string read_name1;
		string read_name2;
		size_t space_index1 = line.find(" ");
		size_t space_index2 = line5.find(" ");
		if(space_index1 == string::npos)
			read_name1 = line;
		else
			read_name1 = line.substr(0, space_index1);
		if(space_index2 == string::npos)
			read_name2 = line5;
		else
			read_name2 = line5.substr(0, space_index2);
		if(read_name1.substr(read_name1.length() - 2, 2) == "/1" && read_name2.substr(read_name2.length() -2 , 2) == "/2")
		{
			if(read_name1.substr(0, read_name1.length() - 1) != read_name2.substr(0, read_name2.length() - 1))
			{
				cerr << "Base name of two ends not consistent"<<endl;
				cerr << "The " << read_id * 2 - 1 << "th read in "<< extendFile<< " and " << read_id * 2 << "th read in" << extendFile << endl;
				cerr << line << endl << line5<<endl;
				exit(1);
			}
		}
		else if(read_name1.substr(read_name1.length() - 2, 2) != "/1" && read_name2.substr(read_name2.length() -2 , 2) != "/2")
		{
			if(read_name1 != read_name2)
			{
				cerr << "Base name of two ends not consistent"<<endl;
				cerr << "The " << read_id * 2 - 1 << "th read in "<< extendFile<< " and " << read_id * 2 << "th read in" << extendFile << endl;
				cerr << line << endl << line5<<endl;
				exit(1);
			}
		}
		else
		{
			cerr << "Base name of two ends not consistent"<<endl;
			cerr << "The " << read_id * 2 - 1 << "th read in "<< extendFile<< " and " << read_id * 2 << "th read in" << extendFile << endl;
			cerr << line << endl << line5<<endl;
			exit(1);
		}

		//cout << line2.length() << ':' << readlen << endl;

		if (line2.length() < readlen)
		{
			cerr << "Read length < specified read length"<<endl;
			cerr << "The " << read_id * 2 - 1 << "th read in " << extendFile << endl;
			cerr << line << endl << line2 << endl << line3 << endl << line4 <<endl;

			exit(2);
		}

		if (line6.length() < readlen)
		{
			cerr << "Read length < specified read length"<<endl;
			cerr << "The " << read_id * 2 << "th read in " << extendFile << endl;
			cerr << line5 << endl << line6 << endl << line7 << endl << line8<<endl;

			exit(2);
		}

		if (!fasta)
		{
			if (line2.length() != line4.length())
			{
				cerr << "Read length and quality string length not consistent"<<endl;
				cerr << "The " << read_id * 2 - 1 << "th read in " << extendFile << endl;
				cerr << line << endl << line2 << endl << line3 << endl << line4 <<endl;

				exit(3);
			}

			if (line6.length() != line8.length())
			{
				cerr << "Read length and quality string length not consistent"<<endl;
				cerr << "The " << read_id * 2 << "th read in " << extendFile << endl;
				cerr << line5 << endl << line6 << endl << line7 << endl << line8<<endl;

				exit(3);
			}
		}

		if (line2.length() > maxlen)
			maxlen = line2.length();

		if (line2.length() < minlen)
			minlen = line2.length();

		if (line6.length() > maxlen)
			maxlen = line6.length();

		if (line6.length() < minlen)
			minlen = line6.length();
	}

	input_fs.close();
	
	//cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

void convert_ex_single(char* extendFile, int fasta, int min_read_len)
{
	size_t read_id = 0;

	ifstream input_fs(extendFile);
	if( !input_fs ) 
	{
		fprintf(stderr,"error: open extended file error\n");exit(1);
	}

	string readidline;
	string line, line2, line3, line4;
	
	while(!input_fs.eof())
	{
		

		getline(input_fs, line);
		getline(input_fs, line2);

		if (line.empty() || line2.empty())
			continue;

		++countall;

		//output_fs <<line[0] <<countall<<'~'<<  line.c_str() + 1 << endl << line2 << endl;

		if (!fasta)
		{
			getline(input_fs, line3);	
			getline(input_fs, line4);

			if (line3.empty() || line4.empty())
				continue;

			for (size_t i = 0; i < line4.size(); ++i)
			{
				char qual = line4[i];

				if (minqual > qual)
					minqual = qual;

				if (maxqual < qual)
					maxqual = qual;
			}

			//output_fs << line3 << endl << line4 << endl;
		}

		++read_id;

		/*if (line.find(" ") != string::npos || line.find("\t") != string::npos)
		{
			cerr << "read name contain blank space or tab"<<endl;
			cerr << "the " << read_id << "th read in "<< extendFile<< endl;
			cerr << line<<endl;

			exit(1);
		}*/

		if (line2.length() < readlen)
		{
			cerr << "Read length < specified read length"<<endl;
			cerr << "The " << read_id << "th read in " << extendFile << endl;
			cerr << line << endl << line2 << endl << line3 << endl << line4 <<endl;

			exit(2);
		}

		if (!fasta)
		{
			if (line2.length() != line4.length())
			{
				cerr << "Read length and quality string length not consistent"<<endl;
				cerr << "The " << read_id * 2 - 1 << "th read in " << extendFile << endl;
				cerr << line << endl << line2 << endl << line3 << endl << line4 <<endl;

				exit(3);
			}
		}

		if (line2.length() > maxlen)
			maxlen = line2.length();

		if (line2.length() < minlen)
			minlen = line2.length();
	}

	input_fs.close();
	
	//cout<<countall<<" reads converted for "<<extendFile<<endl; 
}

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr,"error: too few arguments\n");
		fprintf(stderr,"reads_1 reads_2 ... read_len ispaired(1|0) isfasta\n");
		exit(1);
	}

	int min_read_len = atoi(argv[argc - 1]);

	int fasta;

	ifstream ifs(argv[1]);

	string firstline;

	getline(ifs, firstline);

	if (firstline[0] == '>')
	{
		fasta = 1;

		cout << "read_format -f"<<endl;
		
		cerr << "-----[Read Format: FASTA]"<<endl;
	}
	else if (firstline[0] == '@')
	{
		fasta = 0;

		cout << "read_format -q"<<endl;
		
		cerr << "-----[Read Format: FASTQ]" <<endl;
	}
	else
	{
		cerr << "Input query sequence appears to be unknown format"<<endl;

		cerr << firstline<<endl;

		exit(1);
	}

	int pairend = atoi(argv[argc - 2]);

	readlen = atoi(argv[argc - 3]);

	//cout << "readlen:" << readlen << endl;

	if (pairend && argc % 2 == 0)
	{
	  cerr << "-----[Read Type: Pair End]" << endl;  

		for (size_t i = 1; i < argc - 3; i += 2)

		{
			convert_ex_paired_sep(argv[i], argv[i + 1], fasta, min_read_len);
		}
	}
	else if (pairend && argc % 2 != 0)
	{
	  cerr << "-----[Read Type: Pair End]" << endl;  

		for (size_t i = 1; i < argc - 3; i += 1)
		{
			convert_ex_paired_comb(argv[i], fasta, min_read_len);
		}
	}
	else
	{
	  cerr << "-----[Read Type: Single]" << endl;  

		for (size_t i = 1; i < argc - 3; ++i)
			convert_ex_single(argv[i], fasta, min_read_len);
	}

	cerr <<"-----[Total # Reads: "<< countall << "]" <<endl;

	cerr <<"-----[Max Read Length: " << maxlen << "]" << endl;
	
	cerr <<"-----[Min Read Length: " << minlen << "]" << endl;
	
	cout << "num_read " << countall << endl;
	
	cout <<"minlen "<<minlen << endl;
	
	cout <<"maxlen "<<maxlen<< endl;

	if (!fasta)
	{
		cout <<"maxqualscore "<<maxqual<<" minqualscore "<<minqual<<endl;

		//cerr <<"Max quality score: "<<maxqual<<" Min quality score: "<<minqual<<endl;
		
		cerr << "-----[Max Quality Score: " << maxqual << "]"<< endl;
		
		cerr << "-----[Min Quality Score: " << minqual << "]"<< endl;
	}

	if (minqual > 64 && maxqual < 126)
	{
		//cerr << "Quality score scale appears to be Phred+64"<<endl;
		
		cerr << "-----[Quality Score Scale: Phred+64]"<<endl;

		cout << "qual_scale phred64"<< endl;

	}
	else if (minqual > 59 && maxqual < 126)
	{
		//cerr << "Quality score scale appears to be Solexa+64"<<endl;
		
		cerr << "-----[Quality Score Scale: Phred+64]"<<endl;

		cout << "qual_scale solexa64"<< endl;
	}
	else if (minqual > 33 && maxqual < 126)
	{
		//cerr << "Quality score scale appears to be Phred+33"<<endl;
		
		cerr << "-----[Quality Score Scale: Phred+33]"<<endl;

		cout << "qual_scale phred33"<< endl;
	}
	else
	{
		cerr << "-----[Quality Score Scale: Unknown]" <<endl;
		
		cout << "qual_scale unknown" << endl;
	}

	return 0;
}
