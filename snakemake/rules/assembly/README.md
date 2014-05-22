The assembly package
====================
The assembly package contains all the rules for each implemented assembler and scaffolder.

Prefixes
--------
Each rule is determined whether it has to be executed by prefixes in the filename. 
All prefixes are shown in the following list:  

* wgs.contigs - Create contigs with the WGS assembler
* Allpaths:
	* allpaths.contigs - Create contigs with Allpaths
	* allpaths.scaffolds - Create scaffolds with Allpaths
* pbJelly.scaffolds - Create scaffolds with PBjelly
* scarpa.scaffolds - Create scaffolds with Scarpa

Wgs
------
> Celera Assembler is a de novo whole-genome shotgun (WGS) DNA sequence assembler. 
It reconstructs long sequences of genomic DNA from fragmentary data produced by whole-genome shotgun sequencing. 
Celera Assembler was developed at Celera Genomics starting in 1999. It was released to SourceForge in 2004 as the 
wgs-assembler under the GNU General Public License. The pipeline revised for 454 data was named CABOG.

Allpaths
--------
ALLPATHS-LG is a short read assembler and it works on both small and large (mammalian size) genomes. 
To use it, you should first generate ~100 base Illumina reads from two libraries: 
one with ~180 bp insert size (fragment data), and one with a ~3000 bp insert size (jumping library), both at about 45x coverage. 
Sequence from longer fragments will enable longer-range continuity.

PBjelly
-------
> PBJelly is a highly automated pipeline that aligns long sequencing reads (such as PacBio RS reads or long 454 reads in fasta format) 
to high-confidence draft assembles. PBJelly fills or reduces as many captured gaps as possible to produce upgraded draft genomes. 

Scarpa
------
> Scarpa is a stand-alone scaffolding tool for NGS data with paired end reads. It can be used together with virtually any genome assembler and any NGS read mapper 
that supports SAM format. Other features include support for multiple libraries and an option to estimate insert size distributions from data. 
Scarpa is available free of charge for academic and commercial use under the GNU General Public License (GPL).