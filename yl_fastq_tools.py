#!/usr/bin/env python
#
# Copyright (c) 2017 Yunhai Luo
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies ofthe Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""The module contains a few tools useful for handling fastq file of NGS data.

The initial aim of this module is to extract randomers (molecular barcodes) from
fastq sequences and use these info after alignment to remove PCR duplicates 
through BT&QT tags in the sam file and picard's MarkDuplicates functionality 
combinatorially.
"""

import argparse
import itertools


def fastq_reader(fastq):
    """Get reads from fastq file.

    A generator function which iterates over a fastq file and returns fastq 
    reads by reading in 4 lines at a time. Every read will be returned  as a 
    list of 4. As for now, the fastq file won't be validated. Thus errors 
    won't be detected or handled properly.

    Arg:
        fastq: Path, including name, of one fastq file.

    Yield:
        A list with 4 elements having the expected order of seqence ID,
        sequence, Phred Quality Scores ID, and Phred Quality Scores. For
        example:
        ['@M02357:256:000000000-ARWLK:1:1101:17612:1950 1:N:0:1',
         'TTAGTTTCGTGTGGAAAGGACGAAACACCGCCTCGAGACTCGTTTACTAGGTTTAAGAGCTATGCT',
         '+',
         '>>11>DDD1C1>11111A11B0000A0FE0AEEA////A11/BEF/1A210DDGB2111BF1F1FG']

    Notes:
        Fastq file won't be validated. Therefore, blank or non-fastq line(s) in
        the file may knock off (shift) the expect frame of fastq reads. Also, if
        there are more than one blank lines at the end of file, all but the last
        blank line will be returned as ''.
    """

    with open(fastq) as fastq_file:
        while True:
            read = [line.rstrip('\n') for line in itertools.islice(fastq_file,
                                                                   4)]
            if not read:
                break
            else:
# TODO Handle fastq errors, including blank lines and unequal quality and read 
#      lengths
                yield read


def _fastq_slicers(idxs_arg):
    """Convert a string input of slicing indexes to a list of slice objects.

    Arg:
        idx_arg: A string defining regions to be sliced out from a sequence. It
            uses the following syntax rules:
                - The first index value is 1;
                - Regions are separated by comma;
                - A single region can be either one single number or a range in 
                    the form of "i:j";
                - One single number stands for a single character of interest;
                - "i:j" stands for a range of interst, staring from position i 
                    and ending at postion j. Both the start and the end are 
                    included. Omitted first index defaults to the start, and 
                    omitted second index defaults to the end;
                - This input will not be validated. Characters other than 
                    number, colon, and comma will be converted to integer and 
                    may lead to unexpected results.
    
    Return:
        slicer_array: A list of slice objects. Each slice object within the 
            array defines a range of interest. The list can be used in for 
            loops. The slice object will be used between square brackets 
            (stringtobesliced[sliceobject]) within each loop.
    """
    
    slicer_array = []
    for idx_slice in idxs_arg.split(','):
        if ':' in idx_slice:
            i, j = idx_slice.split(':')
            if i == '':
                i = 1
            if j == '':
                slicer_array.append(slice(int(i)-1, None))
            else:
                slicer_array.append(slice(int(i)-1, int(j)))
        else:
            slicer_array.append(slice(int(idx_slice)-1, int(idx_slice)))
    return slicer_array


def _barcode_file_map(bcfile_arg):
    """Generate a dictionary to map the barcode to its target output fastq.
    
    Arg:
        bcfile_arg: A string used for demultiplexing which contains barcode 
            sequences and file names of their corresponding fastq output. Use 
            syntax 'barcode1-1,barcode1-2,barcode1-3,...:filename1;barcode2-1,
            barcode2-2,barcode2-3,...:filename2;...' (no white space and no 
            ending semicolon). If the barcode sequence extracted from a read 
            doesn't matches any of the barcodes provided here, this read will be
            written into a file with "_index_undetermined.fastq" postfix.
    Return:
        bc_file_dict: A dictionary mapping one barcode to one target fastq. Keys
            are barcode sequences and values are file objects. Different 
            barcodes can point to one single fastq file. Duplicate barcodes end 
            up with only one file object whose name is the last filename it 
            points to.
    """
    
    bc_file_dict = {}
    for bc_file_tuple in bcfile_arg.split(';'):
        barcodes, fastq_name = bc_file_tuple.split(':')
        fastq_file = open(fastq_name, 'w')
        for barcode_seq in barcodes.split(','):
# TODO Abort and raise exception when finding duplicate barcodes.
            bc_file_dict[barcode_seq] = fastq_file
    return bc_file_dict


def hamming_distance(seq1, seq2):
    """Return the Hamming distance between equal-length sequences
    """
    
    if len(seq1) != len(seq2):
        raise ValueError("Undefined for sequences of unequal length.")
    return sum(bp1 != bp2 for bp1, bp2 in zip(seq1, seq2))


def fastq_writer(fastq, read):
    pass

def process_fastq(fastq, seq_idx = None, bc_idx = None, bc_file = None, 
                  rand_idx = None, default_output = None):
    """Pre-procees a fastq file for downstream analysis
    
    Current functionalities are:
        1. Slice sequencing reads;
        2. Demultiplex libraries according to barcodes;
        3. Extract randomers and save the info for the downstream PCR duplicates
           removal.
    
    Args:
        fastq: Fastq file to be processed; one at a time. Supports for fastqs 
            resulting from other platforms have not been tested. Though pair end
            reads can be processed separately, it is not recommended since it 
            may uncouple the pair.
        seq_idx: A string defining region(s) of target sequence. Sequences 
            outside regions defined by seq_idx will be cropped in output fastq. 
            Use syntax 'start1:stop1,start2:stop2,...' with following rules:
                - The first index value is 1;
                - "start:stop" stands for one region of interst. Both the 
                    "start" and the "stop" are included. Omitted "start" index 
                    defaults to the start of the original sequence, and 
                    omitted "end" index defaults to the end of the original 
                    sequence;
                - This argument will not be validated. Characters other than 
                    number, colon, and comma will be converted to integer and 
                    may lead to unexpected results.
        bc_idx: A string defining region(s) of library barcode sequence, which
            will be used for demultiplexing. Use the same syntax as the 
            "seq_idx" argument.
        bc_file: A string defining a map from barcodes to output file names. Use
            syntax 'barcode1-1,barcode1-2,barcode1-3,...:filename1;barcode2-1,
            barcode2-2,barcode2-3,...:filename2;...' with following rules:
                - Barcodes should not have repeats. File names should not have
                    repeats either. Repeats won't be checked but will cause
                    overwriting and mess up the result.
                - When doing barcode sequence match, degenerated base "N" is 
                    considered as failing the match. Therefore, the barcode 
                    sequence provided here should not have "N".
                - By default, reads failed to match any of the barcodes 
                    provided here will be output to a fastq file named with 
                    "_index_undetermined.fastq" postfix.
        rand_idx: A string defining region(s) of randomer (molecular barcode) 
            sequence, which will be extracted and kept in the corresponding 
            read's name (sequence id). This information will be used later after
            alignment. In the sam file, it should be moved from the read's name 
            to the BC&QT tag, which can be used by the Picard tool during its 
            "MarkDuplicates" process. Use the same syntax as the "seq_idx" 
            argument.
    """
    
    # Prepare indexes and output files
    if seq_idx:
        seq_slicers = _fastq_slicers(seq_idx)
    if bc_idx and bc_file:
        bc_slicers = _fastq_slicers(bc_idx)
        bc_file_map = _barcode_file_map(bc_file)
        if default_output is None:
            default_output = fastq.replace('.fastq', 
                                           '_index_undetermined.fastq')
        undetermined_index = open(default_output, 'w')
    else:
        if default_output is None:
            default_output = fastq.replace('.fastq', '_processed.fastq')
        output = open(default_output, 'w')
    if rand_idx:
        rand_slicers = _fastq_slicers(rand_idx)
        
    # Read and process the fastq file
    for read in fastq_reader(fastq):
        seq_id, seq, _, pqs = read
        if rand_idx:
            rand_seq_list = []
            rand_pqs_list = []
            for slicer in rand_slicers:
                rand_seq_list.append(seq[slicer])
                rand_pqs_list.append(pqs[slicer])
            read[0] = '@{}'.format(':'.join([''.join(rand_seq_list), 
                                             ''.join(rand_pqs_list), 
                                             seq_id[1:]]))
        if seq_idx:
            target_seq_list = []
            target_pqs_list = []
            for slicer in seq_slicers:
                target_seq_list.append(seq[slicer])
                target_pqs_list.append(pqs[slicer])
            read[1] = ''.join(target_seq_list)
            read[3] = ''.join(target_pqs_list)
        if bc_idx:
            bc_seq_list = []
            for slicer in bc_slicers:
                bc_seq_list.append(seq[slicer])
            read_bc = ''.join(bc_seq_list)
            output = undetermined_index
            for barcode in bc_file_map:
                if hamming_distance(barcode, read_bc) <= 0:
                    output = bc_file_map[barcode]
# TODO Handle N in barcodes (both barcode and read_bc) and process reads mapping
#      to more than one barcode
                    break
        read.append('') # So that the last line has EOL after .join
        output.write('\n'.join(read))

    # Close all output files
    if bc_idx:
        for file in bc_file_map.itervalues():
            file.close()
        output = undetermined_index
    output.close()

    return None


def main():
    """Use the main() function to test the major function of this module:

        - Process fastq file.
    """

    parser = argparse.ArgumentParser(description='Pre-procees one fastq file '
                                                 'for downstream analysis.')
    parser.add_argument('-f', '--fastq',
                        help='One fastq file from Illumina sequencing. Fastq '
                             'files coming from other platforms have not been '
                             'tested. Though pair end reads can be processed '
                             'separately, it is not recommended since it may '
                             'uncouple the pair.',
                        required=True)
    parser.add_argument('-s', '--sequence',
                        help='Basepair position for target sequences which '
                             'should be kept for downstream analysis. Support '
                             'both continuous and gapped regions. Use syntax '
                             '\'start1:stop1,start2:stop2,...\'. Both the start'
                             ' and the end are included. If a stop is omitted, '
                             'start-stop pairs after it will be ignored and '
                             'sequencing bases will be kept up to the end. If '
                             'this option is omitted, the whole sequence will '
                             'be kept. However, at least one of the \'-s\', '
                             '\'-i\'+\'-b\' and \'-r\' options must be '
                             'provided.',
                        default=None)
    parser.add_argument('-i', '--bcindex',
                        help='Basepair position for library barcodes which help'
                             ' demultiplex reads. Barcode sequences and '
                             'corresponding output file names must be provided'
                             ' in the \'-b\' option (check it for details). '
                             'Support both continuous and gapped regions. Use '
                             'syntax \'start1:stop1,start2:stop2,...\'. The '
                             'syntax specification is the same as \'-s\' option'
                             ' above (check it for details). At least one of '
                             'the \'-s\', \'-i\'+\'-b\' and \'-r\' options must'
                             ' be provided.',
                        default=None)
    parser.add_argument('-b', '--bcfile',
                        help='Barcode sequences and corresponding output file '
                             'names for demultiplexing. Use syntax '
                             '\'barcode1-1,barcode1-2,barcode1-3,...:filename1;'
                             'barcode2-1,barcode2-2,barcode2-3,...:filename2;'
                             '...\'. No space is allowed. At least one of the '
                             '\'-s\', \'-i\'+\'-b\' and \'-r\' options must be '
                             'provided.',
                        default=None)
    parser.add_argument('-r', '--randomer',
                        help='Basepair position for randomer '
                             '(molecular barcode) which help remove PCR '
                             'duplicates generated during library preparation. '
                             'It will be extracted and saved in corresponding '
                             'reads\' name temporarily. Later in the pipeline '
                             '(in the sam file), it can be moved to the BC&QT '
                             'tag which can be utilized by the Picard tool for '
                             'marking/removing duplicates. Support both '
                             'continuous and gapped regions. Use syntax '
                             '\'start1:stop1,start2:stop2,...\'. The syntax '
                             'specifications are the same as \'-s\' option '
                             'above (check it for details). At least one of the'
                             ' \'-s\', \'-i\'+\'-b\' and \'-r\' options must be'
                             ' provided.',
                        default=None)
    parser.add_argument('-o', '--output',
                        help='Default output file. If \'-i\' and \'-b\' options'
                             ' are provided, the file specified by this option '
                             'will be used for collecting reads failed to match'
                             'any of the barcodes. When this option is omitted,'
                             ' unmatched reads will be outputted to a fastq '
                             'file named with input fastq file name plus a '
                             '"_index_undetermined.fastq" postfix. If \'-i\' '
                             'and \'-b\' options are omitted, all processed '
                             'reads will be outputted to one single fastq file '
                             'specified by this option. When this optioned is '
                             'omitted, all processed reads will be outputted to'
                             ' a fastq file named with input fastq file name '
                             'plus a "_processed.fastq" postfix.',
                        default=None)

    args = parser.parse_args()
    if (args.sequence is None and (args.bcindex is None or args.bcfile is None) 
       and args.randomer is None):
        parser.print_help()
        exit(1)
    else:
        process_fastq(args.fastq, args.sequence, args.bcindex, args.bcfile, 
                      args.randomer, args.output)


if __name__ == '__main__':
    main()

