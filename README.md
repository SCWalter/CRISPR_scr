# CRISPR_scr
Data analysis scripts related to CRISPR screen. Currently, it has only one functionality for checking the quality of a CRISPR guides library.

**Table of Contents**
- [Usage](#user-content-usage)

##Usage
`guidescount.py [-h] -i INPUT [INPUT ...] (-t TABLE | --bt2 BT2) [-r RANDOMER] [-g GUIDE] [-o OUTPATH] [-p THREAD]`
###Arguments
&nbsp;&nbsp;&nbsp;&nbsp;Argument&nbsp;&nbsp;&nbsp;&nbsp;|Help message
---|---
-i|Illumina fastq file(s). Support multi-file input. An individual result folder and report will be generated for each fastq input.
-t, --table|Table with designed/expected guide infomation, which will be used for making bowtie2 index in the same directory as this table. The table should have 3 columns: 1) Gene name; 2) Guides name/number; 3) Guide sequence. The first row will be discarded as header which then followed by one guides per row. A fake header is required if there is no header yet. If a Bowtie2 index has already been built before, you can provide it through the "--bt2" option. This option will not be accepted if a Bowtie2 index is provided by the "--bt2" option. An error will be raised if both "-t" and "--bt2" options are provided (mutual exclusion).
--bt2|Bowtie2 index for expected guides in CRISPR library design. If provided, this index will be used even "--table" option is also provided. In other words, the "--bt2" option is preferred over "--table" option. In fact, though, an error will be raised if both "-t" and "--bt2" options are provided (mutual exclusion).
-r, --randomer|Basepair position for randomer, which will be used to remove PCR duplicates. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b". For example: "1:9,31:50" stands for a sequence combining nucleotides from position 1 to position 9 with nucleotides from position 31 to position 50.
-g, --guide|Basepair position for 20bp CRISPR guide. As explain above ("--randomer" option), this option also supports both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".
-o, --outpath|Path to put outputs. A single folder named "CRISPR\_guideCounts" will be created in this path. Each fastq input will have a corresponding result folder under the "CRISPR\_guideCounts" folder. One other folder named "View\_Result\_Summary" will also be created under the "CRISPR_guideCounts" folder. It has the HTML report from this pipeline. Check "Outputs" below for details. Default is the same directory as the first input fastq.
-p, --thread|Number of parallel search threads for bowtie2 alignment. Default is 8. This option is the same as the "-p" option for bowtie2 and will be used directly in bowtie2. Check bowtie2 documentation for details.
-h, --help|show this help message and exit
