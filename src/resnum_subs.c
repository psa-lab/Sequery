/*
* $Log:	resnum_subs.c,v $
 * Revision 1.3  92/08/21  21:03:58  mp
 * Moved "include ctype.h" for System V compatibility.
 * 
 * Revision 1.2  92/03/23  18:30:34  mp
 * Added additional "origin" info to allow most PDB sequences (those that
 * do not contain any suffixed residue names) to fit within the numeric
 * origin-count scheme, reducing "sequery" run-time space requirements
 * from 7 meg to 2 meg.
 * 
 * Revision 1.1  92/03/17  19:11:32  mp
 * Initial revision
 * 
*/
#ifndef lint
static char rcs_id[] =
 "@(#) $Header: /export/asd/prog/sequery/devel/src/resnum_subs.c,v 1.3 92/08/21 21:03:58 mp Exp Locker: mp $";
#endif

#include <stdio.h>
#include "resnum_subs.h"

/*Vishal's changes*/
#include <ctype.h>
 
#define MAXSEQLEN 2048
#define MAXNSEQ 8000

  char * 
get_resnumber(num,seqp,c_num)
int num;
struct seq *seqp;
char *c_num; /* modified */
{
	/* writes residue number/name into "c_num" which should be allocated by
	 * caller as a character array big enough to hold largest residue number/name.
	 */
	if (seqp->resnumber==NULL || seqp->resnumber[num]==NULL && seqp->origin_is_numeric ) {
		sprintf(c_num, "%d", num+seqp->origin_n); /* name is merely index in array */
		}
	//these lines modified because doing a strcpy with a NULL source results in a segmentation fault in Linux using gcc
	else if (seqp->resnumber[num] == NULL) {
		c_num = NULL;
	} else {
		strcpy(c_num, seqp->resnumber[num]); /* copy name from char array */
	}
	
	return(c_num);
}



 int
fget_seq( seqp, count, seqfile)
struct seq * seqp;
int count;
FILE * seqfile;
{
	/* read up to "count" sequences from file. Return number read.
	 * Note: this malloc's the sequence strings.
	 */
int k=0; /* count of ones read */
register int /* char */ c;
register int i;
char * malloc();
extern char * pgmname;
int non_standard; /* flag for sequence that contains "weird" residue numbers. */
#include <ctype.h>


while(count--) {
	if(4!=fscanf(seqfile, 
	  "%8s %1s %s %d", 
	  seqp->name, seqp->chain, seqp->origin, &seqp->len) 
	  ) return k;
	(void) struptolow(seqp->name); /* force to lower case */

	/* set "origin_n" if origin is purely numeric */
	seqp->origin_is_numeric = is_numeric(seqp->origin);
	if(seqp->origin_is_numeric) sscanf(seqp->origin, "%d", &seqp->origin_n);
		
	seqp->sequence = (char *) malloc(1+seqp->len); /* needs better checking... */
	seqp->resnumber=(char **)NULL;
	non_standard=0;
	i=0;
	while( i<seqp->len) {
		char *resnum_ptr,*c_ptr;
		int resnum,j;
		int itoa();
		c = getc(seqfile);
		if(c==EOF) {
			free(seqp->sequence); /* incomplete */
			seqp->len = 0;
			return k;
			}
		if( isspace(c) ) continue;
		if ( c == '#' ) {
			/* bypass comments */
			while ( '\n' != (c= getc(seqfile)))  if(c==EOF) return k;
			}
		if(c=='(') {
			if(!non_standard) {
				seqp->resnumber=(char **)calloc(seqp->len+1, sizeof(char *));
				non_standard=1;
				}
			c_ptr=resnum_ptr=(char *)malloc(9);
			for(j=0;(c=getc(seqfile))!=')';j++) *c_ptr++=c;
			*c_ptr='\0';
			seqp->resnumber[i]=resnum_ptr;
			sscanf(resnum_ptr,"%d",&resnum);
			/* check to see if we're really beginning a standard sequence from
			 * the specified origin (MP experiment):
			 */
			if(i==0 && seqp->origin_is_numeric && resnum == seqp->origin_n) {
				/* if so, make us "standard" to save space & speed */
				free(seqp->resnumber[i]);
				free(seqp->resnumber);
				seqp->resnumber = NULL;
				non_standard=0;
				}
			while(isspace(seqp->sequence[i]=getc(seqfile))); /* MP: why this asgnmt? */
			++i;
		  	}
		else if(non_standard) {
			/* ordinary number but not equal to index in array ... */
			seqp->resnumber[i]=(char *)malloc(9);
			sprintf(seqp->resnumber[i],"%d",++resnum);
			seqp->sequence[i++] = c;
			}
		else seqp->sequence[i++] = c;
		}
	seqp->sequence[i] = '\0';
	k++;
	seqp++;
	}
return k;
}


 char *
struptolow(str)
char * str;
{
#include <ctype.h>
register char * s = str;
	/* make upper case alphabetics into lower case */
	for(;*s;s++) if(isascii(*s) && isupper(*s) ) *s = tolower(*s);
	return str;
	}

 int
is_numeric(s)
char * s;
{
#include <ctype.h>

	for(;*s;s++) if(!( isascii(*s) && isdigit(*s) ) ) return 0;
	return 1;
	}
