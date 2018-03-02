/*
* $Log:	matchextractpdb.c,v $
 * Revision 1.2  92/03/23  18:28:57  mp
 * Adapted to use subroutines in "resnum_subs.c" instead of private copy.
 * 
 * Revision 1.1  92/03/17  19:12:42  mp
 * Initial revision
 * 
*/
#ifndef lint
static char rcsid[] =
 "@(#) $Header: /asd/prog/sequery/devel/src/matchextractpdb.c,v 1.2 92/03/23 18:28:57 mp Exp $";
#endif

#include <stdio.h>


#include "resnum_subs.h" /* defines seq structure & access functions. */

#define MAXSEQLEN 2048
#define MAXNSEQ 8000

main(argc, argv)
int argc;
char **argv;
{

/* interface to getopt(3) library routines */
extern char *optarg;
extern int optind, opterr;

struct seq *seqp;
static struct seq seq[MAXNSEQ]; /* in-core array of sequences */
char * sequery_home();
char * get_resnumber();
char * c_ptr;
char buf[128],chain[2],start_res[8],stop_res[8],pdbcode[12];
char pdbfilename[132],baby_pdbfilename[132];
char chain_id,pdb_res_seq[8],c_res[5];
char baby_dir[80];
FILE *seqfile,*matchfile,*pdbfile,*baby_pdbfile;
int start_index,stop_index;
int start_flag,stop_flag,last_res_flag,eof_flag;
int i,j;
int c;
static int context_pre=6;
static int context_post=6;
int errflg=0;

  sprintf(baby_dir,".");
  while (( c = getopt(argc, argv, "f:x:")) != -1 ) switch(c) {

  case 'f':
	 strcpy(baby_dir,optarg); break;
  case 'x':
	 context_pre = context_post = atoi(optarg); break;
  default:
   	 errflg = 1; break;
  }
  
  if(errflg) {
	fprintf(stderr, "usage: matchextractpdb -x extra_residues < matchfile \n");
	exit(-1);
  }

  seqfile = fopen(sequery_home("lib/pdbseq.asc"),"r");
  seqp = &seq[0];
  fget_seq(seqp, MAXNSEQ, seqfile);  
  while(gets(buf)!=NULL) {
    if(buf[0]=='#') {
	printf("%s",buf);
	continue;
    	}
    seqp = &seq[0];
    sscanf(buf,"%s %c %s %*s %s",pdbcode,chain,start_res,stop_res); 
    while((strcmp(pdbcode,seqp->name)!=0)||(strcmp(chain,seqp->chain)!=0)) seqp++;
    for(start_index=0;(strcmp(start_res,c_ptr=get_resnumber(start_index,seqp,c_res))!=0)&&(start_index<seqp->len);start_index++);
    for(stop_index=0;(strcmp(stop_res,c_ptr=get_resnumber(stop_index,seqp,c_res))!=0)&&(stop_index<seqp->len);stop_index++);
    start_index-=context_pre;
    stop_index+=context_post;
    if (start_index<0) start_index=0;
    if (stop_index>seqp->len-1) stop_index=seqp->len-1;
    sprintf(pdbfilename,"/mb/data/pdb/struct/%s.pdb",pdbcode);
    if ((pdbfile=fopen(pdbfilename,"r"))==NULL) {
	fprintf(stderr,"Error opening file: %s\n",pdbfilename);;
	exit(-1);
    	}
    if (chain[0]=='_') chain_id=' '; 
      else chain_id=chain[0];

    sprintf(baby_pdbfilename,"%s/%s.%s.%s.%s.pdb",baby_dir,pdbcode,chain,start_res,stop_res);
    if ((baby_pdbfile=fopen(baby_pdbfilename,"w"))==NULL) {
      	fprintf(stderr,"Error opening file: %s\n",baby_pdbfilename);
	exit(-1);
    	}

    start_flag=stop_flag=last_res_flag=0;
    while((eof_flag=fgets(buf,sizeof(buf),pdbfile)!=NULL)&&(!stop_flag)) {
	if ((strncmp(&buf[0],"ATOM",4)==0)||(strncmp(&buf[0],"HETATM",6)==0)) {
                
/* Remove leading blanks in Residue seq. no. */
		j=0;
                for(i=23;i<28;i++) if (buf[i]!=' ') pdb_res_seq[j++]=buf[i];
		pdb_res_seq[j]='\0';
		if((strcmp(pdb_res_seq,get_resnumber(start_index,seqp,c_res))==0)&&(chain_id==buf[21])) start_flag=1;
		if((start_flag)&&(strcmp(pdb_res_seq,get_resnumber(stop_index,seqp,c_res))==0)) last_res_flag=1;
		if((last_res_flag==1)&&(strcmp(pdb_res_seq,get_resnumber(stop_index,seqp,c_res))!=0)) stop_flag=1;
		if((start_flag)&&(!stop_flag)) fprintf(baby_pdbfile,"%s",buf);
        }
    }
    if (!start_flag) fprintf(stderr,"Error: Can't find start res seq %s in file %s\n",get_resnumber(start_index,seqp,c_res),pdbfilename);
    if ((!stop_flag)&&(last_res_flag==0)) fprintf(stderr,"Error: Can't find last res seq %s in file %s\n",get_resnumber(stop_index,seqp,c_res),pdbfilename);
    fclose(pdbfile);
    fclose(baby_pdbfile);
  }
}
