#!/usr/bin/env python

import argparse, os

parser = argparse.ArgumentParser(description='Check the quality (guides representation) of CRISPR library')
parser.add_argument('-i', '--input', help='illumina fastq file(s)', nargs='+', required=True)
bt2index = parser.add_mutually_exclusive_group(required=True)
bt2index.add_argument('-a', '--fasta', help='sequences of expected guides. Will be used for making bowtie2 index')
bt2index.add_argument('--bt2', help='bowtie2 index for expected guides.')
parser.add_argument("--mm", help="allow 1 mismatch during bowtie2 alignment. Default does not all mismatch.", action="store_true")
parser.add_argument('-o', '--outpath', help='path to put outputs', required=True)

args = parser.parse_args()

if args.bt2:
	bt2index = os.path.basename(args.a).replace(".fasta", "")
	cmd = "bowtie2-build {} {}".format(args.a, bt2index)
	print(cmd)
	os.system(cmd)
else: bt2index = args.bt2
	
for fastq in args.input:
	name = os.path.basename(fastq)
	trimdedup = name.replace(".fastq", "_trimmed_dedup.fasta")
	guideseq = name.replace(".fastq", "_processed_guides.fastq")
	bt2unmap = name.replace(".fastq", "_unmapped.fastq")
	bt2map = name.replace(".fastq", "_mapped.sam")
	bt2stats = name.replace(".fastq", "_stats.txt")
	gcounts = name.replace(".fastq", "_counts.txt")

	### Keep base 1 to 50
	### Keep 100% QS>30
	### Remove duplicates
	cmd_trim = "fastx_trimmer -Q33 -l 50 -i {}".format(fastq)
	cmd_dedup = "fastx_collapser -Q33 -v > {}".format(args.outpath + trimdedup)
	cmd = ' | '.join(cmd_trim, cmd_dedup)
	print(cmd)
	os.system(cmd)
	
	### Make arbitrary fastq since collapser gives only fasta
	### Keep only 20bp guides
	cmd_fastq = "perl fasta_to_fastq.pl {}".format(args.outpath + trimdedup)
	cmd_trim = "fastx_trimmer -Q33 -f 31 > {}".format(args.outpath + guideseq)
	cmd = ' | '.join(cmd_fastq, cmd_trim)
	print(cmd)
	os.system(cmd)
	
	### Align to expected guides
	if args.mm:
		cmd_bt2 = "bowtie2 -p 8 -L 20 -N 0 -x {} {} --un {} -S {} > {} 2>&1".format(bt2index, args.outpath + guideseq, bt2unmap, bt2map, bt2stats)
	else:
		cmd_bt2 = "bowtie2 -p 8 -L 20 -N 0 --no-1mm-upfront -x {} {} --un {} -S {} > {} 2>&1".format(bt2index, args.outpath + guideseq, bt2unmap, bt2map, bt2stats)
	print(cmd)
	os.system(cmd)
	
	### Count mapped reads from a SAM file
	cmd = "samtools view -F 260 {} | cut -f 3 | sort | uniq -c | awk '{printf(\"%s\t%s\n\", $2, $1)}' > {}".format(bt2map, gcounts)
	print(cmd)
	os.system(cmd)
