# CRISPR_scr
Data analysis scripts related to CRISPR screen. Currently, it has only one functionality for checking the quality of a CRISPR guides library.

**Table of Contents**
- [Usage](#user-content-usage)

##Usage
`usage: guideCount.py [-h] -i INPUT [INPUT ...] (-t TABLE | --bt2 BT2) [-r RANDOMER] [-g GUIDE] [-o OUTPATH] [-p THREAD]`
###Arguments
-i|Illumina fastq file(s)
-t, --table|Table with guide infomation, which will be used for making bowtie2 index in the same directory as this table. The table should have 3 columns: 1) Gene name; 2) Guides name/number; 3) Guide sequence. The first row will be discarded as header which then followed by one guides per row. A fake header is required if there is no header yet. If a Bowtie2 index has already been built before, you can provide it through the "--bt2" option. This option will not be accepted if a Bowtie2 index is provided by the "--bt2" option.
--bt2|Bowtie2 index for expected guides. This option is not accepted if a table with guide designs was provided through the "-t/--table" option.
-r, --randomer|Basepair position for randomer, which will be used to remove PCR duplicates. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".
-g, --guide|Basepair position for 20bp CRISPR guide. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".
-o, --outpath|Path to put outputs. Default is the same directory as the first input fastq.
-p, --thread|Number of parallel search threads for bowtie2 alignment. Default is 8.
-h, --help|show this help message and exit
