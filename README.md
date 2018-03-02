#  Sequery

Please note the following:

- Use of this software implies that you agree with the license
  terms in the accompanying file, license.txt

- Installation instructions and usage documentation for Consolv 1.0 
  are available at the following location:
  http://www.bch.msu.edu/~kuhn/software/consolv/doc.html

- A copy of the LICENSE.txt file is available in this directory.

--- 


## Introduction

Sequery is a tool to search the sequences of the protein structures in the [Protein Data Bank](http://www.rcsb.org) (PDB) for a particular pattern of residues, which may include exact matches and acceptable substitutions based on a user-specified amino acid substitution matrix and/or a numerical threshold. Sequery was developed by Michael E. Pique, Michael A. Siani, Leslie A. Kuhn, Elizabeth D. Getzoff, and John A. Tainer, Department of Molecular Biology, The Scripps Research Institute and is distributed here with permission of the authors. Development of Sequery was supported in part by NSF grant BIR9631436.

A useful complement to Sequery is [SSA](https://github.com/psa-lab/SSA) (Superpositional Structure Assignment), which automates the assignment of secondary structure to tetrapeptides identified by a Sequery search. For more information and software availability, see the [SSA repository](https://github.com/psa-lab/SSA).

Publications pertinent to the use of Sequery include the following:

- J. F. Collawn, M. Stangel, L. A. Kuhn, V. Esekogwu, S. Jing, I. S. Trowbridge, and J. A. Tainer (1990) Transferrin Receptor Internalization Sequence YXRF Implicates a Tight Turn as the Structural Recognition Motif for Endocytosis, Cell 63, 1061-1072.
- L. Craig, P. C. Sanschagrin, A. Rozek, S. Lackie, L. A. Kuhn, and J. K. Scott (1998) [The Role of Structure in Antibody Cross-Reactivity Between Peptides and Folded Proteins (pdf)](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/CraigJMB1998.pdf), J. Mol. Biol., 281, 183-201.
- M. A. Kron, S. Kutlesa, A. Hendrick, A. Liu, J. Leykam, S. Cichenowicz, and L. A. Kuhn (2008) [Using Structural Analysis to Generate Parasite-Selective Monoclonal Antibodies (pdf)](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/Kron_et_al.ProtSci08.pdf), Protein Science 17, 983-989.

For other literature references related to Sequery, please see the section on Algorithmic Details at the end of this page.

Usage information for Sequery can be found in the section on Running Sequery.


## Installation

Sequery will work on any UNIX system with a C-compiler and has been tested on Sun hardware running Solaris and SGI hardware running Irix. Additional scripts require awk, C-shell (csh) and Bourne-shells (sh) installed, as well as an accessible copy of the [Protein Data Bank](http://www.rcsb.org). The latest version of Sequery, as well as the most recent version of this documentation, can be found on this informational page for the Protein Structural Analysis and Design Laboratory, Department of Biochemistry, Michigan State University, and is distributed here with permission of the co-authors at The Scripps Research Institute, Michael Pique, Elizabeth Getzoff, and John Tainer.


To install Sequery, perform the following steps:

- Download this GitHub repository and unzip it
- Place the downloaded folder in the directory into which you wish to install the software
- Nagivate into the Sequery directory and execute the folling command:
    + `make install`   
    (Note: This defaults to using the current directory as the root directory for Sequery. If this is not desired, alter the makefile to include the desired location where noted in the makefile. Also, if Sequery is later moved from this location, it will need to be recompiled to work correctly.)
- There are a few additional modifications that must be made in some of the scripts in sequery/share to point to the correct file and directory locations. These locations will be specific to your system. All settings are at the top of the script code. The scripts and variable settings needed are as follows:

| Script Name       | Variable to set at the top of the script |
|-------------------|------------------------------------------|
| `addname`         | `PDBHOME`                                |
| `genpdbselectseq` | `PDBHOME`  and `SEQUERY_HOME`            |
| `minipdbextract`  | `PDBHOME`                                |

The PDBHOME variable should be set to the directory of the Protein Data Bank files (e.g. /usr/pdb)

The SEQUERY_HOME variable should be set to the directory in which Sequery was installed (e.g. /usr/bin/sequery)


#### Troubleshooting Shell Scripts

Located in the `Sequery/share` directory are several shell scripts, all of which assume the C-shell and Bourne-shell are installed in the standard locations, which are /bin/csh and /bin/sh. If you have difficulty running one or more of these scripts, you may edit the first line of each script to point to the correct shell location for your system. Sequery scripts which are dependent on shell locations include:

| Script Name       | Variable to set at the top of the script |
|-------------------|------------------------------------------|
| `addname`         | `/bin/sh`                                |
| `genpdbselectseq` | `/bin/sh`                                |
| `genpdbseq`       | `/bin/sh`                                |
| `minipdbextract`  | `/bin/csh`                               |


Additionally, these scripts assume you have awk, nawk, or gawk (different implementations of a simple column- or field-oriented programming language) running on your system. (See The Free Software Foundation (GNU) web site for more information on obtaining and [using gawk](https://www.gnu.org/software/gawk/).


## Library Files

Several library files are included in the sequery/lib directory for use with Sequery. The use of these files is explained later in Running Sequery.

- pdbselectAug98.list -- The PDB Select list including PDB sequences with <25% pairwise identity from August 1998. Further information on the PDB Select list can be found at PDB Select List Site at EMBL.
- pdbselect98Aug.all.asc -- Sequery sequence file generated from all entries in the Aug 1998 25% threshold PDB Select list. (Includes both crystallographically-solved and NMR-solved protein structures.) In the distribution package, this is the default sequence file.
- pdbselect98Augseq.xtalonly.asc -- Sequery sequence file generated from entries in the Aug 1998 25% threshold PDB Select list whose structures were solved via x-ray crystallography for studies in which one wishes to use restrict the search space to crystallographically solved structures.
- pdbComplete.asc -- Sequery sequence file generated from the complete PDB database.
- sequery.defs -- Default Sequery sequence file read during execution. In the distribution package, this is a copy of pdbselect98Aug.all.asc.

## Running Sequery

#### Overview of Usage

Sequery runs by searching for matching sequences (strings) in a precomputed ascii text file which contains sequences for each of the PDB files (or a specific subset of PDB files). Sequery is run using input supplied on stdin, either by redirection of a file containing the sequence patterns to match, see end of this section for an example, or by prompting the user for each pattern. It is common to run Sequery using tetrapeptide matching patterns since sequences of this length or shorter often have close matches in the PDB, but any length is possible. Execution is stopped via ctrl-D. Command line options are as follows:

- `-s SequenceFile`: The SequenceFile contains the listing of PDB codes and corresponding amino acid sequences. If omitted, Sequery defaults to searching sequery/lib/pdbseq.asc. Sequery uses this list to search for sequence pattern matches. This list can be generated using either genpdbseq or genpdbselectseq. Using genpdbselectseq is the suggested method, including only sequences with less than 25% identity, as this eliminates statistical bias introduced by including structurally related proteins. For more information, see the Additional Scripts section below. Note that any sequence database formatted according to the example below may be used.

Example portion of a sequence file (excerpt from pdbselectAug98.all.asc):

```
1hmc B    4    143 SEYCSHMIGSGHLQSLQRLIDSQMETSCQITFEFVDQEQLKDPVCYLKKA
                   FLLVQDIMEDTMRFRDNTPNAIAIVQLQELSLRLKSCFTKDYEEHDKACV
                   RTFYETPLQLLEKVKNVFNETKNLLDKDWNIFSKNCNNSFAEC
1myp A    1    104 CPEQDKYRTITGMCNNRRSPTLGASNRAFVRWLPAEYEDGFSLPYGWTPG
                   VKRNGFPVALARAVSNEIVRFPTDQLTPDQERSLMFMQWGQLLDHDLDFT
1grx _    1     84 MQTVIFGRSGC(13)YSVRAKDLAEKLSNERDDFQYQYVDIRAEGITKED
                   LQQKAGKPVETVPQIFVDQQHIGGYTDFAAWVKENLDA
                   PEPA
```

The fields are the PDB code, the chain identifier ( `_` indicates that there are no chain ID's in this structure), the first residue number of the chain, the last residue number of the chain, and the sequence of the chain. In some structures, such as 1grx, the sequence field will contain residue number. These indicate an instance of non-sequential number in the sequence, e.g. due to the lack of diffractive density for a mobile loop in the protein, and are used to maintain correct sequence numbering for Sequery output.

- `-d DefinitionFile`: The DefinitionFile is a file containing acceptable amino acid substitutions. If omitted, Sequery defaults to using sequery/lib/sequery.defs. The supplied substitution file with each line corresponding to a set or equivalence class of substitutable amino acids, sequery/lib/sequery.defs, was determined based on the Dayhoff mutation data matrix, although any set of substitutions could be provided in this format. When entering the sequence pattern for a Sequery, an upper-case charater indicates a search for an exact match while a lower-case character indicates that all equivalent residues from this file may be considered as substitutes (e.g. `A` to match alanine only and `a` for all residues equivalent to alanine). Further details can be found below in Sequence Query Patterns.

- `-w WildcardFile`: The WildcardFile contains a listing of user-defined acceptable amino acid substitutions (again, with each line containing a set of amino acids that can substitute for each other), e.g. from acceptable variation observed in a sequence alignment or from mutagenesis studies. When entering sequence patterns during Sequery execution, the user enters the line number within the WildcardFile corresponding to the acceptable amino acids at that position. For example, based on the example WildcardFile (sequery/lib/wilddef.dat), entering a 2AAA would find all patterns starting with tyrosine, phenylalanine, or tryptophan (line 2 in the file), followed by 3 alanines.

- `-o OutputFile`: The OutputFile is the file where the user would like output to be placed. If omitted, output will be written to sequery.match, overwriting any previously existing sequery.match.

- `-x NumberOfContextResidues`: This is the number of residues printed (in lower-case) on either side of the matched sequence pattern (in upper-case). Default is 4.

- `-v verbose` mode: More output (mostly for debugging purposes)'

-`q quiet` mode: Suppresses all output to the screen except for error messages. Pattern matches will still be output to the output file.

- `-? or -h`: Gives version and help information.

Sequery can be run in batch mode via the following (where search.patterns contains one line for each pattern to search):

    sequery -s lib/pdbseq.asc -d lib/sequery.defs \
    -w wilddef.dat -q -o search.matches < search.patterns

#### Sequence Query Patterns

Note: The following assumes the user is running with the default definition and wildcard files.

Sequence query patterns follow the regular expression rules, as documented in the UNIX `ed` man page. Following is a brief overview of pattern rules followed by several examples:

**Rules Summary:**

- Upper-case letter -- match character (single-leter code for amino acid) exactly
- Lower-case letter -- match amino acid according to allowed substitutions in the definition file
- . or X -- match any amino acid
- Digit -- match any of the amino acids in that line number of the wildcard file
- [Upper-case letters] -- match any of these amino acids (inside the brackets) in this position.
- [Lower-case letters] -- match any of these amino acids (inside the brackets), allowing for substitutions according to the definition file
- [^Upper-case letters] -- match any sequence not containing any of these amino acids at this position. (i.e. all other amino acids.)
- [^Lower-case letters] --match any sequence not containing any of these amino acids or their substitutions as allowed for in the definition file at this position. (i.e. match all other amino acids.)
- `r\{m,n\}` -- match m to n repetitions of the character r, where m and n are digits. (Note that each digit must be preceded by a \ for Sequery, since digits without the \ result in wildcard matching). Example: 3 to 7 residues of any kind: `.\{\3,\7\}`; r can be upper-case (exact match) or lower-case (substituted match).

**Query Examples:**

- ACYE -- search for Ala-Cys-Tyr-Glu
- acye -- expands to [A][C][Y][DQEBZ] according to the definitions file -- search for Ala-Cys-Tyr-(Asp or Gln or Glu or Asx or Glx)
- [AC]YEA -- search for (Ala or Cys)-Tyr-Glu-Ala
- [^ACDE]YEA -- search for (not(Ala or Cys or Asp or Glu))-Tyr-Glu-Ala. This is effectively a search for Tyr-Glu-Ala not preceeded by Ala, Cys, Asp, or Glu.
- AC.E or ACXE -- search for Ala-Cys-any single residue-Glu
- AC2E -- expands to AC[YFW]E according to line 2 of the wildcard file.
- A\{\3,\7\} -- expands to search for AAA, AAAA, AAAAA, AAAAAA, or AAAAAAA. (3 to 7 Ala's in a row)

**Example Files**

There are example input, output, and wildcard files included in the sequery/examples directory. The example output was generated with the following command (run in the examples directory):

    
    ../bin/sequery -d ../lib/sequery.defs -s ../lib/pdbseq.asc \
     -w example.wilddef.dat -o example.matches < example.patterns

## Output Files

Sequery will output the results of the pattern queries to the supplied output file (or to sequery.match if no output filename is given). Output will appear as follows:

    1cem _  153 to  156 -> aatdADEDiala matching ADED
    1occ A   93 to   96 -> apdmAFPRmnnm matching 1234

## Diagnostics

Most errors will occur due to improper query patterns. These errors will appear simply as non-run queries. The current version of Sequery shows unpredictable behavior with proteins with residues having negative residue numbers and will occasionally produce segmentation faults if sequence patterns would result in a very large number of matches. (In this case, break the query into two or more subqueries and combine the results.)

## Scripts & Tools

Included in the sequery/share directory are a few additional tools that may be useful and are described below. These scripts must be modified upon installation to point to the correct directories. See the previous section "Installing Sequery" for more information.

- `addname` -- adds protein name from the PDB file to the end of a Sequery output line. Output is to stdout. This is useful for noting functional patterns in the proteins that match.  

Syntax:  `addname SequeryOutputFile [columns]`  

SequeryOutputFile: file generated by sequery
columns: total number of columns of each of line of output (default 80). If not specified, output is truncated to 80 chars per line.

- `genpdbselectseq` -- generates a sequence file which contains the amino acid sequences, along with PDB codes and chain IDs, which is used by Sequery to search for pattern matches, for each of the PDB files in the PDB Select list. For more information about the PDB Select list of representative (<25% sequence identity) protein structures and the latest version, see the PDB Select List Site at EMBL. For more information regarding the sequence file, see above under Running Sequery.

Syntax:  `genpdbselectseq pdb-select-list-file output-sequence-file`  

Explanation of use: The PDB Select list is a list of all proteins/chains in the PDB whose sequence identity is less than a certain percentage. (There are different lists for different identity threshold levels.) Each chain in the list represents a set of related chains. By using the lowest-identity threshold (25%), any structural bias in Sequery analysis is minimized. (This bias arises from the fact that if a sequence query identifies sequences in a series of related proteins whose structures are known, any subsequent structural analysis will contain more bias towards these related structures.)

- `genpdbseq` -- generates a sequence file as a Sequery search database from a PDB file or series of PDB files. Output is directed to stdout. 

Syntax: `genpdbseq PDBFile > SequenceFile`

Examples:

    genpdbseq 2sod.pdb > 2sod.ascseq
    genpdbseq *.pdb > pdb.ascseq


- `minipdbextract` -- generates a PDB formatted file for the residues in each line of Sequery output. It takes the start residue, end residue, and pdbcode from Sequery output, searches the $PDBHOME database for that protein, and extracts coordinate lines from the PDB files. 

Syntax:

    minipdbextract SequeryOutputFile [x]

SequeryOutputFile: file generated by sequery
x (optional): Number of flanking residues to include (default is 0)

## Algorithmic Details

Sequery was developed as a successor to Searchwild, which is described in:

- Collawn JF, Kuhn LA, Liu LF, Tainer JA, Trowbridge IS (1991) Transplanted LDL and mannose-6-phosphate receptor internalization signals promote high-efficiency endocytosis of the transferrin receptor, EMBO J., 10(11) (Nov): 3247-3253.

- Collawn JF, Stangel M, Kuhn LA, Esekogwu V, Jing SQ, Trowbridge IS, Tainer JA (1990), Transferrin receptor internalization sequence YXRF implicates a tight turn as the structural recognition motif for endocytosis Cell, 63(5) (Nov 30): 1061-1072.

Other references related to the use of Sequery include the following:

- Craig L, Sanschagrin PC, Rozek A, Lackie S, Kuhn LA, Scott JK (1998), [The Role of Structure in Antibody Cross-Reactivity Between Peptides and Folded Proteins](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/CraigJMB1998.pdf) (pdf) J. Mol. Biol., 281(1) (Aug 7): 183-201.
- Chang CP, Lazar CS, Walsh BJ, Komuro M, Collawn JF, Kuhn LA, Tainer JA, Trowbridge IS, Farquhar MG, Rosenfeld MG (1993), Ligand-induced internalization of the epidermal growth factor receptor is mediated by multiple endocytic codes analogous to the tyrosine motif found in constitutively internalized receptors J. Biol. Chem., 268(26) (Sep 15): 19312-19320.

## More Information

Scientific inquries concerning Sequery should be directed to Leslie Kuhn at: kuhnlab@msu.edu

## Additional Sequery PDB sequence files

Additional Sequery PDB sequence files are included in the [additional_sequence_files](additional_sequence_files) subdirectory.

### 4,304 non-homologous protein chains

Filename: [additional_sequence_files/cullpdb_pc20_res2.0_R0_25col_4304_chains.txt](additional_sequence_files/cullpdb_pc20_res2.0_R0_25col_4304_chains.txt)

This sequence file contains 4,304 chains from the PDB with pairwise less than 20% sequence identify, resolution of .le. 2.0 Angstroms, and R-factor of .le. 0.25 (starting from a culled PDB list from [Roland Dunbrack's PISCES server](http://dunbrack.fccc.edu/Guoli/PISCES_OptionPage.php)).



















Below is a typical Here is a sample line from my own use of sequery



/psa/share/sequery/bin/sequery-32bit-exe -s /home/kuhn/dehom/cullpdb_pc20_res2.0_R0_25col_4304_chains.txt -w ./wilddef.dat

.  The -s flag denotes the PDB sequence file to be searched for user-specified patterns; it contains many chains from the PDB.  It is in a specific format that can be written for any user-defined set of PDB files to search, too (via genpdbseq*).  The one I have attached, which should be included in the distribution with a little documentation, is a good default set to search.  


