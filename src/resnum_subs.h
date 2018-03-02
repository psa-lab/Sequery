char * get_resnumber();
int fget_seq();
char * struptolow();

static char resnum_subs_h_rcsid[] = 
 "@(#) $Header";

struct seq {
        char name[12]; /* protein name, as in PDB file name */
        char chain[2]; /* one-character chain ID */
        char origin[6]; /* Residue number + possible insertion suffix. GHL*/
	int  origin_is_numeric, origin_n; /* flag: true if origin is purely numberic,
		which is true in the majority of PDB entries (few begin with letter suffix),
		and if so, what that origin is, normally 1 . */
        int  len; /* number of residues */
        char *sequence; /* array of 1-character residue types */
        char **resnumber; /* Residue number/name - char for non-integer name, e.g. 106B,
			   * or any residue whose number/name is not merely its 1-origin
			   * index in the array, which is indicated by NULL pointer.
			   */
        };

