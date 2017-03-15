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

""" The initial aim of this module is to attach randomers (molecular barcodes)
    extracted during fastq processing (as part of the read's name) to the
    BT&QT tag in the sam file so that the Picard tool can MarkDuplicates using
    the BARCODE_TAG option.
"""
import argparse
import itertools
import pysam


def attach_barcode(sam, output):
    """Attach randomers (molecular barcodes) to the BC&QT tag in the sam file.
    Args:
        sam: A sam file which needs to attach barcodes at BC&QT tags. The
            barcode and its sequencing quality should have been stored in the
            read's name, separated by ":". The "yl_fastq_tools" module can
            extract randomers accordingly. Check it for details.
        output: Name for an unsorted sam file output. Default name is to add a
            "_bcqt" postfix to the name of the input sam file.
    """
    
    if output is None:
        output = sam.replace('.sam', '_bcqt.sam')
    infile = pysam.AlignmentFile(sam, "r")
    outfile = pysam.AlignmentFile(output, "wh", template=infile)
    for read in infile.fetch():
        id_sam = read.query_name
        sep_si = id_sam.index(':')
# TODO Abort and raise exception if randomer info is not kept properly in the 
#      read's name.
        bc_seq = id_sam[0:sep_si]
        sep_qi = sep_si + 1 + len(bc_seq)
        bc_pqs = id_sam[sep_si + 1: sep_qi]
        read.set_tag('BC', bc_seq)
        read.set_tag('QT', bc_pqs)
        read.query_name = id_sam[sep_qi+1:]
        outfile.write(read)
    outfile.close()
    infile.close()


def main():
    """Use the main() function to test this module:
        - Attach randomers (molecular barcodes) to the BC&QT tag of a sam file.
    """

    parser = argparse.ArgumentParser(description='Attach randomers (in read\'s '
                                                 'name) to the BC&QT tag.')
    parser.add_argument('-s', '--sam',
                        help='One aligned sam file which needs to attach '
                             'barcodes at BC&QT tags. The barcode and its '
                             'sequencing quality should have been stored in the'
                             ' read\'s name, separated by ":". Check '
                             '"yl_fastq_tools" module for detailed formats.',
                        required=True)
    parser.add_argument('-o', '--outsam',
                        help='Output unsorted sam file. Default is to add a '
                             '"_bcs" postfix and save in the same directory.',
                        default=None)

    args = parser.parse_args()
    attach_barcode(args.sam, args.outsam)


if __name__ == '__main__':
    main()

