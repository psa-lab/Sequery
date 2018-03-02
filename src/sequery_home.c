/*
** $Log:	sequery_home.c,v $
 * Revision 1.1  92/03/16  19:25:47  mp
 * Initial revision
 * 
*/

static char rcs_id[] = "@(#) $Header: /export/asd/prog/sequery/devel/src/sequery_home.c,v 1.1 92/03/16 19:25:47 mp Exp Locker: mp $\n";
#include <stdio.h>

#ifdef STARTUP_FILE_PROCESSING /* not implemented yet - M Pique */
read_startup()
{
	FILE *startup, *fopen();

	/*
	** First, try to read the configuration file
	** from the current directory.
	*/
	startup = fopen(STARTUPFILE, "r");

	if (!startup)
	{
		char path[80];
		/*
		** If unsuccessful, try to read the configuration file
		** from user's HOME directory.
		*/
#ifdef unix
		struct passwd *getpwuid();
		struct passwd *p_entry;

		if ( (p_entry = getpwuid(getuid())) )
		{
			sprintf(path, "%s/%s", p_entry->pw_dir, STARTUPFILE);
			startup = fopen(path, "r");
		}
#endif unix
	}

	/*
	** get initial commands from the startup file
	*/
	if (startup)
	{
		register int c;

		while ((c = getc(startup)) != EOF)
			lexputchar(c);
		fclose(startup);
	}
}
#endif

/*
** return string naming path of file 'fname' in the
** directory where run-time library files live.
** If the enviroment variable SEQUERY_HOME is set, use that value, otherwise
**	use default standard library 
** Thus to use current directory for testing, setenv SEQUERY_HOME '.'
**
** If fname is empty, only the directory name itself is returned.
**
** Caution: the returned value points to a static area
** which is overwritten with each call.
*/

char *
sequery_home(fname)
char *fname;
{
#ifndef SEQUERY_HOME
#define SEQUERY_HOME "/psa/lprog/XXX"   /* normally defined in makefile using -D... */
#endif
	static char libname[1024];
	static char filename[1024];
	char *getenv();

	if (!*libname)
	{
		if (getenv("SEQUERY_HOME"))
			strcpy(libname, getenv("SEQUERY_HOME"));
		else
			strcpy(libname, SEQUERY_HOME);
	}

	/*
	** Now we know the library name, build the full filename.
	**  In Unix environment, just concatenate them with a slash.
	*/
	if (!*fname)
		return(libname);
	sprintf(filename, "%s/%s", libname, fname);
	return(filename);
}
