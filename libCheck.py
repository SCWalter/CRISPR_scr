#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Check the quality (guides representation) of CRISPR library')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('-i', '--input', help='illumina fastq file(s)')
parser.add_argument('-a', '--fasta', help='sequences of expected guides. Will be used for making bowtie2 index')
parser.add_argument('--bt2', help='bowtie2 index for expected guides.')
parser.add_argument("--1mm", help="allow 1 mismatch during bowtie2 alignment. Default does not all mismatch.",
                    action="store_true")

args = parser.parse_args()
print args.accumulate(args.integers)

### Keep base 1 to 50
### Keep 100% QS>30
### Remove duplicates
### Keep only 20bp guides
fastx_trimmer -Q33 \
              -l 50 \
              -i /scratch/users/yunhail/CRISPRscreen07292016/sample1_S1_L001_R1_001.fastq \
    | fastx_collapser -Q33 -v > /scratch/users/yunhail/CRISPRscreen07292016/S1_trimmed_filtered_dedup.fasta
perl /home/yunhail/FAST-iCLIP/bin/fasta_to_fastq.pl /scratch/users/yunhail/CRISPRscreen07292016/S1_trimmed_filtered_dedup.fasta \
    | fastx_trimmer -Q33 \
                    -f 31 > /scratch/users/yunhail/CRISPRscreen07292016/S1_processed_guides.fastq
fastx_trimmer -Q33 \
              -l 50 \
              -i /scratch/users/yunhail/CRISPRscreen07292016/sample2_S2_L001_R1_001.fastq \
    | fastx_collapser -Q33 -v > /scratch/users/yunhail/CRISPRscreen07292016/S2_trimmed_filtered_dedup.fasta
perl /home/yunhail/FAST-iCLIP/bin/fasta_to_fastq.pl /scratch/users/yunhail/CRISPRscreen07292016/S2_trimmed_filtered_dedup.fasta \
    | fastx_trimmer -Q33 \
                    -f 31 > /scratch/users/yunhail/CRISPRscreen07292016/S2_processed_guides.fastq

### Align to expected guides
bowtie2 -p 8 \
        -L 20 \
        -N 0 \
        --no-1mm-upfront \
        -x /scratch/users/yunhail/CRISPRscreen07292016/bowtie2_index/PMJ08012016 \
        S1_processed_guides.fastq \
        --un S1_processed_guides_unmapped.fastq \
        -S S1_processed_guides_mapped.sam > S1_processed_guides_bowtie2_stats.txt 2>&1
bowtie2 -p 8 \
        -L 20 \
        -N 0 \
        --no-1mm-upfront \
        -x /scratch/users/yunhail/CRISPRscreen07292016/bowtie2_index/PMJ08012016 \
        S2_processed_guides.fastq \
        --un S2_processed_guides_unmapped.fastq \
        -S S2_processed_guides_mapped.sam > S2_processed_guides_bowtie2_stats.txt 2>&1

### Count mapped reads from a SAM file
samtools view -F 260 S1_processed_guides_mapped.sam \
    | cut -f 3 \
    | sort \
    | uniq -c \
    | awk '{printf("%s\t%s\n", $2, $1)}' > S1_counts.txt
samtools view -F 260 S2_processed_guides_mapped.sam \
    | cut -f 3 \
    | sort \
    | uniq -c \
    | awk '{printf("%s\t%s\n", $2, $1)}' > S2_counts.txt
