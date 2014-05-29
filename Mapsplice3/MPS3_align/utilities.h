
#ifndef __UTILITIES_H_INCLUDED__
#define __UTILITIES_H_INCLUDED__

#include <string>

using namespace std;

class Utilities
{
public:
	Utilities() {}

	static string toUpper(string value)
	{
		for(int i=0; i<value.length();i++)
			value[i] = toupper(value[i]);

		return value;
	}

	static string int_to_str(int numerical)
	{
			char c[100];
			sprintf(c,"%d",numerical);
			string str(c);
			return str;
	}
};

#endif
