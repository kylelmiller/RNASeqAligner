#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wctype.h>
#include <assert.h>

#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <string>
using namespace std;

#define MAX_SIZE	(300 * 1024 * 1024)

int fexist(char* filename)
{
	FILE* fi;
	fi = fopen(filename, "r");
	if (!fi)
		return 0;

	fclose(fi);
	return 1;
}

void bb(char* build, char* target)
{
	char command[1001];
	memset(command, 0, sizeof(command));

	sprintf(command, "%s -q %s %s", build, target, target);

	system(command);
}

int formatwrite(char* outfile, char* data, char* firstline, int ext, int totalsize, int linesize)
{
	int ret = -1;

	FILE *fo;
	char* buffer;
	int copysize;

	buffer = (char*)calloc(linesize + 1, sizeof(char));
	if (!buffer)
	{
		perror("calloc");
		return ret;
	}

	fo = fopen(outfile, "w");
	if (!fo)
	{
		perror("fopen");
		return ret;
	}

	fprintf(fo, "%s_%d\n", firstline, ext);

	while (totalsize > 0)
	{
		if (totalsize > linesize)
		{
			copysize = linesize;
		}
		else
		{
			copysize = totalsize;
			memset(buffer, 0, linesize + 1);
		}
		//printf("%d\n", copysize);

		memcpy(buffer, data, copysize);
		fprintf(fo, "%s\n", buffer);
		data += copysize;
		totalsize -= copysize;
	}
	fclose(fo);

	free(buffer);
	return ret;
}

int readfile(char* infile, char* data, int maxdatasize, char* line, int maxlinesize, int *linedatalen)
{
	int ret = -1;

	char buffer[1001];
	int buflen;
	FILE* fi;
	char* pdata;
	int copylen;

	fi = fopen(infile, "r");
	if (!fi)
	{
		perror("fopen");
		return ret;
	}

	// skip first line
	fgets(line, maxlinesize - 1, fi);

	// clear '\n'
	ret = strlen(line) - 1;
	line[ret] = 0;

	ret = 0;
	pdata = data;

	while (!feof(fi) && (ret < maxdatasize - 1))
	{
		memset(buffer, 0, sizeof(buffer));
		fgets(buffer, sizeof(buffer) - 1, fi);

		buflen = strlen(buffer) - 1;

		if (*linedatalen <= 0)
			*linedatalen = buflen;

		if (buflen <= 0)
			break;

		if (buflen + ret < maxdatasize - 1)
			copylen = buflen;
		else
			copylen = maxdatasize - 1 - ret;

		memcpy(pdata, buffer, copylen);
		pdata += copylen;
		ret += copylen;
	}

	fclose (fi);

	return ret;
}

int getexpectlen(int size, int ext, int loverlap, int roverlap)
{
	assert (ext > 0);

	if (ext == 1)
		return size + roverlap;

	return size + loverlap + roverlap;
}

string basename2(string filename) 
{

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

void bsb(char* buildpath, char* infile, int size, int loverlap, int roverlap, char* outputpath)
{
	char firstline[1001];

	char* wholedata = (char*)calloc(MAX_SIZE, sizeof(char));
	int wholedatalen = 0;

	char* pwb;

	int extension = 1;
	char outfile[256];

	int remaindata = 0;
	int writesize = 0;

	int linedatalen = 0;

	if (!wholedata )
	{
		perror("calloc");
		return;
	}

	pwb = wholedata;	

	memset(firstline, 0, sizeof(firstline));
	memset(wholedata, 0, MAX_SIZE);
	wholedatalen = readfile(infile, wholedata, MAX_SIZE, firstline, sizeof(firstline), &linedatalen);

	remaindata = wholedatalen;
	pwb = wholedata;

	printf("total length: %d\n", wholedatalen);

	string base_file_name = basename2(infile);

	string new_file_name = outputpath;

	new_file_name.append(base_file_name);

	while (remaindata > 0)
	{
		writesize = getexpectlen(size, extension, loverlap, roverlap);

		memset(outfile, 0, sizeof(outfile));
		sprintf(outfile, "%s.%d", new_file_name.c_str(), extension);

		if (remaindata <= writesize)
		{
			writesize = remaindata;
			remaindata = 0;
		}
		else
		{
			remaindata -= writesize - loverlap - roverlap;
		}

		//printf("%d %d %d %d\n", writesize, loverlap, roverlap, remaindata);

		formatwrite(outfile, pwb, firstline, extension, writesize, linedatalen);
		bb(buildpath, outfile);

		pwb += writesize - loverlap - roverlap;
		++extension;

	}

	free(wholedata);

}

void printusage()
{
	printf("thiscommand bowtiepath splitsize leftoverlapsize rightoverlapsize file1 [file2 ...]\n\n");
	return;
}


/*
* thiscommand bowtiepath splitsize leftoverlapsize rightoverlapsize file1 [file2 ...]
*/
int main(int argc, char* argv[])
{
	int i, splitsize, loverlapsize, roverlapsize;

	if (argc <= 6)
	{
		printusage();
		return -1;
	}

	splitsize = atoi(argv[2]);
	loverlapsize = atoi(argv[3]);
	roverlapsize = atoi(argv[4]);

	if (!fexist(argv[1]))
	{
		return -1;
	}

	for (i = 0; i < argc - 6; ++i)
		bsb(argv[1], argv[6 + i], splitsize, loverlapsize, roverlapsize, argv[5]);

	return 0;
}
