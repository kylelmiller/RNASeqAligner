#include <stdlib.h>
#include <string>
#include <string.h>

using namespace std;

class HeaderSection_Info
{
public:
	string headerSectionInfoStr;

	string int_to_str(int numerical)
	{
			char c[100];
			sprintf(c,"%d",numerical);
			string str(c);
			return str;
	}

	HeaderSection_Info(Index_Info* indexInfo)
	{
		headerSectionInfoStr = "@SG\tSN:" + indexInfo->chrNameStr[0] 
			+ "\tLN:" + int_to_str(indexInfo->chromLength[0]) + "\n";
		for(int tmp = 1; tmp < (indexInfo->chrNameStr.size()); tmp++)
		{
			headerSectionInfoStr = headerSectionInfoStr
				+ "@SG\tSN:" + indexInfo->chrNameStr[tmp]
				+ "\tLN:" + int_to_str(indexInfo->chromLength[tmp]) + "\n"; 
		}

		headerSectionInfoStr = headerSectionInfoStr.substr(0, headerSectionInfoStr.length()-1);
	}

};
