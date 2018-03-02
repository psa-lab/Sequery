/* sequery:
 *  search a file of sequences for occurrences of specified patterns.
 *
 * Sequery begins by loading a file of sequences into memory. 
 * It then reads patterns from the standard input, one per line, and
 * reports the sequences (called "matches") that contained the given pattern.
 *
 * Sequery exits when there are no more pattens: at the end-of-file
 *  of the standard input (normally, control-D from a terminal).
 *
 * Sequery expands incoming patterns for user convenience 
 *  before matching against the sequences:
 *   1. Each digit is replaced by a corresponding line from a "wilddef.dat"
 *	file, numbered from line one, surrounded by square brackets [ ... ].
 * 	A digit immediately following a \ is protected from this replacement.
 *      Replacement lines may contain lower-case letters, or other patterns,
 *	but not digits (that is, the line is scanned for digits only once).
 *   2. Each lower-case letter (a-z) is replaced by a corresponding line
 *	from a "sequery.defs" file, the line beginning with that letter as
 *	a key, with the second field in the line replacing the letter.
 *   3. Each "X" character is replaced by a dot (which matches any residue)
 *
 * Sequery expansions are intended to make it easier to prepare input:
 *  for exploring classes of residues similar by hydrophobicity, size, 
 *  charge, turn propensity, etc.
 *   
 *  Expansion example: suppose file "wilddef.dat" contained 
 *	ND      <- line 0
 *	RA      <- line 1
 *   and file "sequery.defs" contained
 *      a [AGS]    <- replacement for a
 *      w [YWF]    <- replacement for w
 *
 *   Then the incoming pattern YW1*wAQ would be expanded
 *     first to YW[RA]*wAQ  by replacing the 1 with [RA]
 *     then to  YW[RA]*[YWF]AQ by replacing the w with [YWF]
 *
 *   As explained below, this pattern would match a Y followed by a W,
 *    followed by zero or more repetitions (the *) of R, or A,
 *    followed by either Y, W, or F, followed by A, followed by Q.
 *   
 *
 *  The resulting pattern is then considered to be a "regular expression"
 *	and is compared against each sequence in turn.  A sequence matches
 *	a regular expression if the sequence contains the regular expression
 *	at least once.
 *  The rules for regular expressions are as documented in the Unix "ed"
 *	manual entry, and are only summarized here. Note "ed" patterns are more
 *	restrictive than full "egrep" and "awk" patterns (no parens or ORs).
 *	"A", "B", and "C" are any single characters (residue types).
 *	Metacharacters are brackets, braces, and commas.
 *
 *	Expression...	Matches...
 *	C		the single nonmetacharacter C
 *	\C		the single character C
 *	^		beginning of sequence
 *	$		end of sequence
 *	.		any character
 *	[ABC]		any character in the class of characters A,B,C
 *	[^ABC]		any character NOT in the class...
 *			(currently, may not include lower-case expansions)
 *	(r)\{m,n\}	m to n repetitions of r. M and N are digits. Note
 *			 that each digit must be preceded by a \ for sequery.
 *			Example: 3 to 7 residues of any kind: .\{\3,\7\}
 *
 *	Note that any database of alphabetic sequences could be searched, 
 *	provided that its file is in the format expected by this program.
 *
 *
 * Arguments (all optional):  
 *
 *  -d DEFINITION_FILE : use custom file for 1-letter amino acid
 *			acceptable substitutions.
 *			Default: $SEQUERY_HOME/lib/sequery.defs
 *
 *  -w WILDDEFS_FILE : use custom file for single-digit 
 *			shorthand substitutions. Default: ./wilddef.dat
 *			Note there is no "global" default file.
 *
 *  -s SEQUENCE_FILE : use custom file of protein/DNA/etc sequences.
 *			Default: $SEQUERY_HOME/lib/pdbseq.asc
 *
 *  -x NUMBER_OF_CONTEXT_RESIDUES : show this many residues on each side
 *			of the match.  Default: 4
 *
 *  -v : (verbose) : give more output, mostly for debugging.
 *  -q : (quiet) :  give no output except error messages.
 * 
 *  -o OUTPUT_FILE : use OUTPUT_FILE for copy of matches.
 *			Default: sequery.match.
 *			- means standard output (normally used in
 *			  combination with -q (quiet) option so that
 *			  nothing but matches goes to standard output).
 *
 *  -? or -h : give version and help info
 *
 * Input is always from stdin.  If stdin is a terminal, rather than a file
 *   or pipe, the user is prompted for lines and is given reports on number
 *   of matches found that are suppressed otherwise unless -v is given.
 *   Also, if stdin is a terminal the wilddef and definition files are
 *   re-read before each query: if stdin is not a terminal, they are
 *   read only once per sequery run.
 *
 * Output is always to stdout, and after each pattern (not after each run
 *  of sequery) a copy of the matches is placed in file "/tmp/sequeryXXXXXX" 
 * (XXX=Process ID,  see mktemp(3)).
 *
 * Example of running sequery: 
 *
 *  Again supposing file "wilddef.dat" contained 
 *	ND
 *	Ra
 *   and file "sequery.defs" contained
 *   a [AGS]
 *
 *   sequery << EOF
 *   aW
 *   Ha01K 
 *   EOF
 *
 * will perform two searches: the first for any sequence that contains
 *   aW -> [AGS]W   i.e., either A,G, or S, followed by W
 *
 * The second for any sequence that contains
 *   Ha01K -> Ha(N|D)(R|G)
 *   
 *          
 * Output is to standard output, which can be redirected into a file. 
 *
 *
 * The -q and -o options work together for performing batch searches.
 *   Suppose file "search.patterns" contains patterns, one per line.
 *   Then the command
 *     sequery -q -o search.matches < search.patterns
 *   will test each pattern and write all matches into search.matches.
 *
 * This program was written in the programming language C under the
 * SunOS Unix operating system.
 *
 * Copyright (C) 1990 Scripps Clinic and Research Foundation 
 *
 * by Michael E. Pique, Michael A. Siani, Leslie A. Kuhn, 
 * Elizabeth D. Getzoff, and John A. Tainer.
 *
 * All Rights Reserved
 * Program may be distributed with permission of the authors.
 *
 * Research Institute of Scripps Clinic
 * Department of Molecular Biology MB-5
 * 10666 N. Torrey Pines Road
 * La Jolla, CA 92037  
 * Telephone: (619) 554-8119
 *
 * $Log:	sequery.c,v $
 * Revision 1.18  92/08/21  21:04:31  mp
 * Write matches into a temporary file instead of directly into a "pipe"
 * so that batch searches that fail to match anything run much faster.
 * Added -q (quiet) option, improved -o (output) file handling.
 * 
 * Revision 1.17  92/08/18  01:47:05  mp
 * Changed remaining references to "Sequery" from old "Searchwild" name.
 * 
 * Revision 1.16  92/08/18  01:38:05  mp
 * Only re-read definition and wilddef file if input is a terminal.
 * This should slightly speed up background searches.
 * 
 * Revision 1.15  92/03/23  18:29:27  mp
 * Adapted to use subroutines in "resnum_subs.c" instead of private copy.
 * Fixed "X" being substitued by ".X" instead of desired ".".
 * 
 * Revision 1.14  92/03/16  19:24:37  mp
 * Incorporated Gary Liao's September 5 1991 changes as well as
 * calls to new "sequery_home" to get directory name for data files.
 * 
 * Revision 1.13  91/06/25  16:24:56  liao
 * If output file specified exists, prompt user for new name, or delete existing file
 * 
 * Revision 1.12  91/05/06  15:28:01  liao
 * Fixed initialization of output pattern in replace_wild()
 * 
 * Revision 1.11  91/05/03  17:07:01  liao
 * Fixed embedded sequence pattern problem. Print names of input files being used at start of program.
 * Added command line argument to change output file from default (sequery.match)
 * Append results to output file rather than overwrite.
 * 
 * Revision 1.10  91/03/27  21:42:39  mp
 * Changed name from old "searchwild" to new "sequery", including changing
 * the name of the file it reads from "searchwild.defs" to "sequery.defs".
 * 
 * Revision 1.9  90/12/03  20:26:26  mp
 * Made default locations for sequence & definition files "/usr/local/lib"
 * instead of personal directories. Improved definition file reading to
 * be more robust & to allow comments (beginning line with #-sign) and
 * blank lines.
 * 
 * Revision 1.8  90/11/27  23:19:03  mp
 * Modified search loop to report all matches of the pattern in each sequence,
 * rather than only the first match.
 * Corrected lower-case conversion of printed sequence (non-upper-case letters
 * such as asterisks) had been mis-printed due to improper use of "tolower".
 * 
 * Revision 1.7  90/10/29  01:57:31  mp
 * Added command-line options to change file names from default values.
 * Improved digit-to-lower-case substitution behavior.
 * 
 * Revision 1.6  90/10/27  00:29:11  mp
 * Replaced calls to "regex(3)" package with "regexp(3)" which
 * allows determining the exact range of residues matched and
 * provides "1-to-n" length restrictions using \{1,n\} patterns
 * (entered as \{\1,\n\} ).
 * 
 * Revision 1.5  90/10/26  20:49:35  mp
 * Divided up recognition & sequence LOCATING portions in preparation
 * for adapting to use regexp(3) instead of regex(3) library routines.
 * 
 * Revision 1.4  90/10/26  12:24:34  mp
 * Corrections to reporting of where the match occurred, changed to accept
 * "origin-count" instead of "count" PDB sequence files.  
 * 
 * Revision 1.3  90/07/01  23:02:15  mp
 * Uses "sort(1)" [inefficiently..] to put matches into usable order.
 * 
 * Revision 1.2  90/06/27  12:19:15  mp
 * Demo version, using C library pattern-matching subroutines.
 * 
 * Revision 1.1  90/06/16  20:56:33  mp
 * Initial revision
 * 
 *
 */
#ifndef lint
static char rcsid [] = 
  "@(#)RELEASE $Header: /export/asd/prog/sequery/devel/src/sequery.c,v 1.18 92/08/21 21:04:31 mp Exp $";
#endif

#include <stdio.h>
#include <strings.h>
#define streq(a,b) (!strcmp((a),(b)))

#include <ctype.h>
 /* define our own "safe" case converters: */
#define lower(c) ( (isascii(c) && isupper(c)) ? tolower(c) : (c) )
#define upper(c) ( (isascii(c) && islower(c)) ? toupper(c) : (c) )

#include "resnum_subs.h" /* defines "seq" structure and access fcns */

char * pgmname;

main(argc, argv)
int argc;
char ** argv;
{
#define MAXSEQLEN 2048
#define PATTERNLEN 1024
#define MAXNSEQ 8000



char pat_in[PATTERNLEN]; /* pattern given as input */
char pat1[PATTERNLEN]; /* pattern after first file expansion (digits) */
char pat2[PATTERNLEN]; /* pattern being searched for */
int pat_len; /* size of pattern being searched for, before second expansion */

/* interface to regexp(3) library routines, modelled after regex(3) */
char * re_compile();
int reg_match();

/* interface to getopt(3) library routines */
extern char *optarg;
extern int optind, opterr;

char * sequery_home();

char seqfilename[1024];
FILE * seqfile; /* file from which sequences are read */
static struct seq seq[MAXNSEQ]; /* in-core array of sequences */
struct seq * seqp; /* pointer to sequence being examined */
int in_core = 1; /* for future "big file" version.... */
int n_seqs; /* number of sequences in core */

int interactive; /* true if input is a terminal, not pipe or file */

char * wilddeffilename = "wilddef.dat";

/* file for building unsorted list of matches: */
static char matchfilename[] = "/tmp/sequerywXXXXXX";
FILE * matchfile;

/* file for building sorted list of matches: */
static char sortfilename[] = "/tmp/sequeryXXXXXX";

char deffilename[1024];

FILE * testfile; /* for testing access to named files */
FILE * outfile;

#define OUTFILE "sequery.match"
char * outfilename = OUTFILE;
char cat_cmd[256];
char shell_cmd[256];
char sort_cmd[256];
char matrix_header[256];

int errflg=0;
int verbose = 0;
int quiet = 0;

int c;

#define CONTEXT 4
static int context_pre = CONTEXT;
static int context_post = CONTEXT;

	pgmname = argv[0];
	interactive = isatty(0);

	/* set defaults or pick up from environment: */
	strcpy(seqfilename, sequery_home("lib/pdbseq.asc"));
	strcpy(deffilename, sequery_home("lib/sequery.defs"));

	/* set from command line options: */
	while (( c = getopt(argc, argv, "s:w:d:x:vqo:h?")) != -1 ) switch(c) {

 case 's':
	strcpy(seqfilename, optarg); break;
 case 'w':
	wilddeffilename = optarg; break;
 case 'd':
	strcpy(deffilename, optarg); break;
 case 'x':
	context_pre = context_post = atoi(optarg); break;
 case 'v':
	verbose = 1; break;
 case 'q':
	quiet = 1; break;
 case 'o':
        outfilename = optarg; break;
 case 'h':
 case '?':
	fprintf(stderr, "%s: version %s of %s\n",
	  pgmname, "$Revision: 1.18 $", "$Date: 92/08/21 21:04:31 $");
	errflg = 1; break;
 default:
	errflg = 1; break;
	}
	
	if(errflg) {
		fprintf(stderr, "%s: usage : \n", pgmname);
		exit(2);
		}

	/* Seems there's no easy way to say you DON'T want to save
	 * the matches into a logfile, which I think OUGHT to be optional,
	 * thus the check against /dev/null. M Pique.
	 */
	if(interactive) while (
	  !streq(outfilename, "-") &&
	 (testfile = fopen(outfilename, "r"))!=NULL) {
		char answer,abuf[80];
		static int count = 0;
		if(count++ > 10) {
			fprintf(stderr, "%s: no legal matchfile\n",
			 pgmname);
			exit(-1);
			}
		fclose(testfile);
		if(streq(outfilename, "/dev/null")) {
			outfilename = NULL;
			break;
			}
		fprintf(stderr,"%s already exists: Overwrite?",outfilename);
		gets(abuf);
		sscanf(abuf,"%c",&answer);
                switch (answer) {
                  case 'y':
                  case 'Y':     
				testfile = fopen(outfilename, "w");
				if(testfile==NULL)  unlink(outfilename);
				else fclose(testfile);
				goto outfile_ok; /* leave 2 levels */
                  case 'n':
                  case 'N':     printf("Output Filename?");
                                gets(abuf);
                                sscanf(abuf,"%s",outfilename);
                                break;
                  default:      break;
                }
	}
outfile_ok:;

	seqfile = fopen(seqfilename, "r");
	if(seqfile==NULL) {
		perror(seqfilename);
		fprintf(stderr,"%s: can't open sequence file %s\n",
		  pgmname, seqfilename);
		exit(-1);
		}
        else {
                if(!quiet)printf("Sequence file: %s\n",seqfilename);
        }

	if(in_core) {
		n_seqs=fget_seq( seq, MAXNSEQ, seqfile);
		if(verbose) printf("read %d sequences from %s\n",
		  n_seqs, seqfilename);
		}

	/* check that the two (optional) shorthand files are present, warn
	 * user if not there or not readable.
	 */
	testfile = fopen(deffilename, "r");
	if(NULL == testfile && !quiet) {
		fprintf(stderr, "%s: no \"%s\" definition file,\n",
		 pgmname, deffilename);
		fprintf(stderr,  " so \"lower case\" abbreviations will not work.\n");
		}
	else {
                if(!quiet) printf("Deffilename: %s\n",deffilename);
		while((NULL!=(fgets(matrix_header,sizeof(matrix_header),testfile)))&&(strncmp(matrix_header,"# Matrix",8)!=0)) ;
		if(strncmp(matrix_header,"# Matrix",8)==0) {
			if(outfilename!=NULL) {
				if(streq(outfilename,"-")) {
					if(!quiet) fprintf(stdout,"%s",matrix_header);
					}
				else {
					outfile=fopen(outfilename,"a");
					if(!quiet) fprintf(outfile,"%s",matrix_header);
					fclose(outfile);
					}
				}
			if(!quiet) printf("%s",matrix_header);
			}
		else {
			if(!quiet) fprintf(stderr,"No Matrix Header in %s\n",deffilename);
		} 
		fclose(testfile);
	}

	testfile = fopen(wilddeffilename, "r");
	if(NULL == testfile && !quiet ) {
		fprintf(stderr, "%s: no \"%s\" definition file,\n",
		 pgmname, wilddeffilename);
		fprintf(stderr,  " so \"digit\" abbreviations will not work.\n");
		}
	else {
                if(!quiet) printf("%s: Wild def filename: %s\n",
		 pgmname, wilddeffilename);
		fclose(testfile);
	}

		
	mktemp(matchfilename); /* modifies matchfilename */
	mktemp(sortfilename); /* modifies sortfilename */

	matchfile = fopen(matchfilename, "w");

	if(outfilename!=NULL && !quiet) printf("Output File: %s\n",outfilename);

	/* build command to strip sort keys and write to stdout */
	if(quiet)
	sprintf(sort_cmd, "sort %s | sed 's/^[^ ]* *[^ ]* *//' > %s", 
	  matchfilename, sortfilename);
	else
	sprintf(sort_cmd, "sort %s | sed 's/^[^ ]* *[^ ]* *//' | tee %s", 
	  matchfilename, sortfilename);

	/* build command to append to sequery.match or other match file */
	if(outfilename!=NULL) { 
		if(streq(outfilename, "-"))
		sprintf(cat_cmd,"cat %s", sortfilename);
		else
		sprintf(cat_cmd,"cat %s >> %s", sortfilename, outfilename);
		}
	/* note from mp: this whole output processing scheme needs to
	 * be re-engineered now that we have some idea what's useful.
	 */


	/* main loop .... */
	while(interactive && !quiet && fprintf(stdout, " > ") , 
	  NULL != gets(pat_in)) {
		int sequences_examined = 0;
		int sequences_matched = 0;
		int matches_found = 0;

		if(strlen(pat_in) == 0) continue;

		pat_len = replace_wild(wilddeffilename, interactive,
		  pat_in, pat1); /* sets pat1 */
		if(replace_defs(deffilename, interactive, 
		  pat1, pat2)==0) continue; /* sets pat2 */
		if(pat_len<=0) continue;
		if(pat_len<=1) {
			if(!quiet)fprintf(stderr, " too short for safety...\n");
			continue;
			}

		if(NULL!= re_compile(pat2)) {
			fprintf(stderr, re_compile(pat2));
			continue;
			}
		else if(!quiet) fprintf(stdout,"%s (length %d) -> %s\n",
		  pat1, pat_len,pat2);
		fflush(stdout);

		seqp = &seq[0];
		if(!in_core) rewind(seqfile);
  
		/* loop over sequences to be examined */
		while ( in_core ?
		  seqp<&seq[n_seqs] : 
		  fget_seq( seqp, 1, seqfile) == 1) {
			int start_index = 0; /* for multiple searches per seq */
			int any_matches_in_this_seq = 0;
			int seq_len;
			int bgn;
			int match_len;
			sequences_examined++;
			seq_len = strlen(seqp->sequence);
			while( start_index<seq_len &&
			 reg_match(seqp->sequence+start_index, 
			  &bgn, &match_len)) {
				int i;
				char bgn_resnum[5],end_resnum[5];
				if(verbose>2)
				  fprintf(stdout, "%s %s %s\n",seqp->name,
				 seqp->chain, seqp->sequence);

				matches_found++;
				any_matches_in_this_seq = 1;

				/* "bgn" is relative to start_index, make
				 * it absolute...
				 */
				bgn += start_index;

				/* ... advance start_index for next search */
				start_index = bgn + 1;

				/* write out match to use as sort key */
				for(i=bgn;i<bgn+match_len;i++)
				  putc(seqp->sequence[i],matchfile);
				/* print protein name (w/ number after)
				 * to use as secondary sort key*/
				 fprintf(matchfile," %s",seqp->name+1);
				 fprintf(matchfile,"%c",seqp->name[0]);

/*
				fprintf(matchfile," %s %s %4s to %4s -> ",
				  seqp->name, seqp->chain, bgn+1, bgn+match_len-1+1);
*/
				(void) get_resnumber(bgn,seqp,bgn_resnum);
				(void) get_resnumber(bgn+match_len-1,seqp,end_resnum);
				fprintf(matchfile," %s %s %4s to %4s -> ",
				  seqp->name, seqp->chain, bgn_resnum, end_resnum);

				/* print part of sequence before match*/
				for(i=bgn-context_pre;i<bgn;i++) 
				  putc(i<0?' ':lower(seqp->sequence[i]),matchfile);
				/* print match */
				for(i=bgn;i<bgn+match_len;i++)
				  putc(seqp->sequence[i],matchfile);

				/* print part of sequence after match*/
				for(i=bgn+match_len;i<bgn+match_len+context_post;i++) 
				  putc(i>=strlen(seqp->sequence)?' ':lower(seqp->sequence[i]),matchfile);
				if(0!=strncmp(pat_in,&seqp->sequence[i],strlen(pat_in)))
				 fprintf(matchfile," matching %s", pat_in);
				fprintf(matchfile,"\n");
				}
			if(any_matches_in_this_seq) sequences_matched++;
			if(in_core) seqp++;
			else free(seqp->sequence); /* MPique: suspect bug here, memory leak */
			}
		if(verbose || (sequences_matched>0 && interactive) && !quiet)
		 fprintf(stdout, "%d match%s in %d out of %d sequences:\n",
		 matches_found, matches_found==1?"":"es",
		 sequences_matched, sequences_examined);

		if(sequences_matched>0) {
			fflush(matchfile);
			/* sort matchfile into sortfile, print to stdout  */
			system(sort_cmd); 


			if(interactive)
			 if(!quiet) fprintf(stdout, 
			 "%d match%s in %d out of %d sequences.\n",
			 matches_found, matches_found==1?"":"es",
			 sequences_matched, sequences_examined);

			/* append to sequery.match log file  */
			if(outfilename!=NULL && (errflg=system(cat_cmd)!=0)) {
				  fprintf(stderr,
				    "Error %d: Can't append %s to %s\n",
				    errflg, matchfilename, outfilename);
				  exit(-1);
				  }

			/* I don't think closing and reopening *ought*
			 * to be necessary, rewinding ought to be sufficient,
			 * but it doesn't work correctly for me that way. mp
			 */
			fclose(matchfile);
			matchfile = fopen(matchfilename, "w");
			}
		} /* end main loop */
	unlink(matchfilename);
	unlink(sortfilename);
	return 0;
	}

 int
replace_wild(filename, interactive, pat_in, pat)
char *filename;
int interactive; /* re-read file each time if true */
char pat_in[], pat[];
{
	/* read wild cards from file "wilddef.dat", 
	 * copies pat_in to pat, replacing "wild card" digits
	 * with correspondingly-numbered lines from file.
	 * Returns pattern length.
	 */

#include <ctype.h>

register char *d, *s;

#define MAXWILDS 12
#define WILDLEN 100
static char wild[MAXWILDS][WILDLEN+1];
FILE * wildfile;
static int nwilds;
char buf[200];

static int firstcall = 1;

if(interactive || firstcall) {
	for(nwilds=0;nwilds<MAXWILDS;nwilds++) strcpy(wild[nwilds],"");

	wildfile = fopen(filename, "r");
	if(wildfile == NULL) nwilds = 0;
	else {
		nwilds=0;
		while(nwilds<MAXWILDS) {
			if(NULL==fgets(buf, WILDLEN, wildfile)) break;
			if(buf[0]=='\n') continue; /* skip empties */
			if(buf[0]=='#') continue; /* skip comments */
			if(1==sscanf(buf, "%s", wild[nwilds])) nwilds++;
			}
		fclose(wildfile);
		}

	firstcall = 0;
	}

	/* expand pattern by wildcards (represented by integers 0-9) */
	s = pat_in; /* src */
	d = pat; /* dest */
	*d = '\0';
	while ( *s ) {
		if(isspace(*s)) {
			s++;
			continue;
			}
		/* replace 'X' (unknown) by '.' (matches any character) */
		if(*s=='X') {
			*d++ = '.';
			s++;
			continue;
			}

		/* peek ahead if we see a backslash - protect the next digit
		 * from our own integer wildcard expansion. 
		 */
		if(*s=='\\' && isdigit(*(s+1))) {
			*d++ = *++s; /* skip backslash, copy digit */
			}
		else if(isdigit(*s)) {
			char * c;
/* GHL 5-30-91
   Start counting wilddef replacement index at 1. i.e. 1 is replaced by
   the first line in wilddef
*/
			if(*s-'1' > nwilds-1) {
				/* error */
				fprintf(stderr," wildcard entry %c not defined\n", *s);
				pat[0] = '\0'; /* error exit */
				return 0;
				}
			/* make a new string, surrounded by [ ] */
			strcat(d,"[");
			c = wild[*s-'1'];
			while(*c) {
				static char a[2];
				a[0] = *c;
				strcat(d,a);
				c++;
				}
			strcat(d,"]");
			d = &pat[strlen(pat)]; /* first unoccupied char */
			}
		else *d++ = *s;
		s++;
		*d='\0'; /* keep dest string terminated, for strcat */
		}
	return strlen(pat);
	}

 int
replace_defs(filename, interactive, pat_in, pat)
char *filename;
int interactive; /* re-read file each time if true */
char pat_in[], pat[];
{
	/* read string substitution definitions from file "sequery.defs", 
	 * copies pat_in to pat, replacing lower-case letters with 
	 * expansions suitable for "ed"-style regular
	 * expressions.
	 *
	 * Returns pattern length.
	 */

register char *d, *s;

#include <ctype.h>

#define MAXDEFNS 40
#define DEFNLEN 100
#define KEYLEN 1
static char key[MAXDEFNS][KEYLEN+1];
static char defn[MAXDEFNS][DEFNLEN+1];
FILE * defnfile;
static int ndefns;
char buf[200];
char output_buf[PATTERNLEN];
int linenumber;
int in_brackets = 0; /* boolean */
static int firstcall = 1;


if(interactive || firstcall ) {
	for(ndefns=0;ndefns<MAXDEFNS;ndefns++) key[ndefns][0] = defn[ndefns][0] = '\0';

	defnfile = fopen(filename, "r");
	if(defnfile == NULL) ndefns = 0;
	else {
		ndefns=0;
		linenumber=0;
		while(ndefns<MAXDEFNS) {
			if(NULL==fgets(buf, DEFNLEN, defnfile)) break;
			linenumber++;
			if(buf[0] == '\n') continue; /* skip empties */
			if(buf[0] == '#') continue; /* skip comments */
			if(2==sscanf(buf, "%1s %s", 
			  key[ndefns], defn[ndefns])) ndefns++;
			else {
				fprintf(stderr,
			"%s: format problem in file %s, line %d :%s",
				pgmname, filename, linenumber, buf);
				fprintf(stderr," line ignored.\n");
				}
			}
		fclose(defnfile);
		}
	firstcall = 0;
	}

	/* count characters in pattern */
	s = pat_in;

	/* expand character-by-character (protected by \ ) */
	s = pat_in; /* src */
	d = output_buf; /* dest */
	while ( *s ) {
		while(isspace(*s)) s++;

		if(*s=='[') in_brackets++;
		if(*s==']') in_brackets--;

		/* peek ahead if we see a backslash - protect the next character
		 * from our own character expansion. 
		 */
		if(*s=='\\' && isdigit(*(s+1))) {
			s++; /* skip backslash */
			*d = *s; /* copy char */
			}
		else if(isascii(*s) && islower(*s)) {
			/* substitute it from file entry */
			int i;
			for(i=0;i<ndefns;i++) if(upper(*s) == upper(key[i][0])) break;
			if(i == ndefns) {
				/* error */
				fprintf(stderr,"%s: %c not defined in %s\n",
				 pgmname, *s, filename);
				pat[0] = '\0'; /* error exit */
				return 0;
				}
			if(in_brackets) {
				/* remove this layer of brackets - messy. mp */
				sprintf(d, "%s", defn[i]+1);
				d = &output_buf[strlen(output_buf)-2];
				}
			else {
				sprintf(d, "%s", defn[i]);
				d = &output_buf[strlen(output_buf)-1];
				}
			}
		else *d = *s;
		d++;
		s++;
		*d = '\0';
		}
	*d = '\0';
	sprintf(pat, "%s", output_buf);
	return strlen(pat);
	}

/* remainder of this source file is interface to REGEXP(3) : M Pique */
#define INIT register char *sp = instring;
#define GETC() *sp++
#define PEEKC() *sp
#define UNGETC(c) --sp
#define RETURN(c) return
#define ERROR(c) re_error(c)

#define ESIZE 2048
static char expbuf[ESIZE];

#include <regexp.h>

 char  *
re_compile(instring)
char * instring;
{
	compile(instring, expbuf, expbuf+ESIZE, '\0');
	return NULL; /* OK */
	}
re_error(code)
 int code;
 {
 char * msg;
 extern char * pgmname;
 switch(code) {
	case 11: msg = "Range endpoint too large."; break;
	case 16: msg = "Bad number."; break;
	case 25: msg = "``\\ digit'' out of range."; break;
	case 36: msg = "Illegal or missing delimiter."; break;
	case 41: msg = "No remembered search string."; break;
	case 42: msg = "\\( \\) imbalance."; break;
	case 43: msg = "Too many \\(."; break;
	case 44: msg = "More than 2 numbers  given  in \\{ \\}."; break;
	case 45: msg = "} expected after \\."; break;
	case 46: msg = "First number exceeds second in \\{ \\}."; break;
	case 49: msg = "[] imbalance."; break;
	case 50: msg = "Regular expression too long."; break;
	}
	fprintf(stderr,"%s: %s\n", pgmname, msg);
	}
 int
 reg_match(string, p_bgn, p_len)
 char * string;
 int * p_bgn, *p_len;
 {
 extern char *loc1, *loc2, *locs;
	if(step(string, expbuf)) {
		/* match */
		*p_bgn = loc1 - string;
		*p_len = loc2 - loc1;
		return loc2-loc1;
		}
	else return 0;
	}
/* End of interface to REGEXP(3). Do not put other source code here. 
 * (as a matter of style...) Mike Pique 
 */
