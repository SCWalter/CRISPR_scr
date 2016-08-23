#!/usr/bin/env perl
use strict;
use warnings;

# Usage: perl seq_extract.pl [fastqfile] [start] [length]
# Description: Slice out all or part of the sequence for each read and output only
#	extracted sequences (no name or quality score) for all reads (one per line) to STDOUT.
# [fastqfile]: required
# [start], [length]: use the perl 'substr' grammar to slice sequences. First bp starts at index zero. 
# 	If both [start] and [length] are omitted, returns the whole sequence.
# 	If [start] is negative, starts that far back from the end of the string. 
# 	If [length] is omitted, returns everything through the end of the string. 
# 	If [length] is negative, leaves that many characters off the end of the string. 
# Example: 
# 	De-novo discovery of barcodes (bp5-bp10) and their distribution.
# 	zcat Illumina.fastq.gz | perl seq_extract.pl - 5 6 | sort | uniq -c | sort -nrk1

my $file = $ARGV[0];
open(my $fastq,  "<".$file)  or die "Can't open $file: $!";
my $extract = @ARGV - 1;
my $righttrim = @ARGV - 2;
my ($seq_id, $seq, $subseq, $pqs_id, $pqs);

while (<$fastq>) { # assigns each line in turn to $_
	$seq_id = $_; # skip the first line for name/id
	$seq = <$fastq>;
	chomp $seq;
	if ($extract) {
		if ($righttrim) {
			$subseq = substr $seq, $ARGV[1], $ARGV[2];
		} else {
			$subseq = substr $seq, $ARGV[1];
		}
		print $subseq."\n";
	} else {
		print $seq."\n";
	}
	$pqs_id = <$fastq>; # skip the third line for quality score name/id
	$pqs = <$fastq>; # skip the third line for quality score
}

close $fastq or die "$fastq: $!";