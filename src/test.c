#include <ctype.h>


char * struptolow(str) 
char * str;
{
	register char * s = str;
	/* make upper case alphabetics into lower case */
	for(;*s;s++) if(isascii(*s) && isupper(*s) ) *s = tolower(*s);
	return str;
}
