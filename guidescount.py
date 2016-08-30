#!/usr/bin/env python

import argparse
import os
import logging
import datetime

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import jinja2

import yltkv2 as yltk


class CRISPRLibSeqQC:
    """Store info for the result presentation

    Attributes:
        samplename: Name of the fastq file which the following processing
                  info/results belong to.
        outfolder: Directory where you can find processing results for the
                 corresponding (samplename) fastq file.
        bt2stats: A multi-line strings recording the statistic of bowtie2
                alignment.
        picardstats: A multi-line strings recording the statistic of sequencing
                   duplicates.
        filenames: A dict mapping all potentially useful result files.
        figures: A dict mapping 3 summarizing figures.
    """

    def __init__(self, name, out):
        # Create instance variables instead of class variables for class
        # variables are shared across all instances while instance variables are
        # instance specific. Those class variables which are set by immutable
        # value, such as numbers, strings and tuples, behave like instance
        # specific. However, this is not a correct practice since we
        # really/actually need instance variables here. Maybe using objects
        # instead of dictionaries for jinja2 is not a correct practice either.
        self.samplename = name
        self.outfolder = out
        self.bt2stats = ''
        self.picardstats = ''
        self.filenames = {'sam':'', 'guidesCt':'', 'guidespergene':'',
                          'guidemiss':'', 'misspergene':''}
        self.figures = {'numguidespergene':'', 'nummisspergene':'',
                        'coverhist':''}


def main():
    """Pipeline for processing fastq sequencing results of a CRISPR screen
       experiment.

    WARNING: this pipeline was only tested under linux environment. It could
    work on MacOS platform but may not work on Windows (since pysam does not
    work on Windows).
    """

    parser = argparse.ArgumentParser(description='Count number of reads'
                                                 'sequenced for each guide in a'
                                                 'CRISPR screen library.')
    parser.add_argument('-i', '--input', help='Illumina fastq file(s)',
                        nargs='+', required=True)
    bt2index = parser.add_mutually_exclusive_group(required=True)
    bt2index.add_argument('-t', '--table',
                          help='Table with designed/expected guide infomation,'
                          'which will be used for making bowtie2 index in the'
                          'same directory as this table. The table should have'
                          '3 columns: 1) Gene name; 2) Guides name/number; 3)'
                          'Guide sequence. The first row will be discarded as'
                          'header which then followed by one guides per row. A'
                          'fake header is required if there is no header yet.'
                          'If a Bowtie2 index has already been built before,'
                          'you can provide it through the "--bt2" option. This'
                          'option will not be accepted if a Bowtie2 index is'
                          'provided by the "--bt2" option.', default=None)
    bt2index.add_argument('--bt2',
                          help='Bowtie2 index for expected guides. This option'
                               'is not accepted if a table with guide designs'
                               'was provided through the "-t/--table" option.',
                          default=None)
    parser.add_argument('-r', '--randomer',
                        help='Basepair position for randomer, which will be'
                             'used to remove PCR duplicates. Support both'
                             'discrete and continuous regions. Discrete regions'
                             'are separated by ",", and continuous region is'
                             'expressed as "a:b".', default=None)
    parser.add_argument('-g', '--guide',
                        help='Basepair position for 20bp CRISPR guide. Support'
                             'both discrete and continuous regions. Discrete'
                             'regions are separated by ",", and continuous'
                             'region is expressed as "a:b".', default=None)
    parser.add_argument('-o', '--outpath',
                        help='Path to put outputs. Default is the same'
                             'directory as the first input fastq.', default='')
    parser.add_argument('-p', '--thread',
                        help='Number of parallel search threads for bowtie2'
                             'alignment. Default is 8.', default=8)

    args = parser.parse_args()
    if not args.outpath:
        maindir = os.path.join(os.path.dirname(os.path.abspath(args.input[0])),
                               'CRISPR_guideCounts')
    else:
        maindir = os.path.join(args.outpath, 'CRISPR_guideCounts')
    os.system('mkdir -p {}'.format(maindir))
    resultsum = os.path.join(maindir, 'View_Result_Summary')
    os.system('mkdir -p {}'.format(resultsum))
    mainlog = os.path.join(os.path.dirname(maindir), 'CRISPR_guide_counts.log')
    logging.basicConfig(filename=mainlog, level=logging.DEBUG)
    os.system('echo "Please go to {} check for runlog and results."'
              .format(os.path.dirname(maindir)))

    if not args.bt2:
        guidetab_abspath = os.path.abspath(args.table)
        guidetab_dir = os.path.dirname(guidetab_abspath)
        guidetab_base = os.path.basename(guidetab_abspath)
        bt2indexdir = os.path.join(guidetab_dir, 'bowtie2index')
        os.system('mkdir -p {}'.format(bt2indexdir))
        bt2index = os.path.join(bt2indexdir, guidetab_base)
        bt2indexlog = os.path.join(bt2indexdir, '{}_build.log'
                                   .format(guidetab_base))
        fasta = os.path.join(bt2indexdir, guidetab_base + ".fasta")

        tab2fasta = ('tail -n +2 {} | '
                     'awk -F\'\\t\' \'{{print ">"$1":"$2"\\n"$3}}\' '
                     '> {}'.format(guidetab_abspath, fasta))
        logging.info('Generating fasta from the table of guides')
        logging.info(tab2fasta + '\n')
        os.system(tab2fasta)

        bt2build = ('bowtie2-build {} {} 1>{} 2>&1'
                    .format(fasta, bt2index, bt2indexlog))
        logging.info('Making Bowtie2 index from fasta file')
        logging.info(bt2build)
        logging.info('Check %s for progress and stats.%s', bt2indexlog, '\n')
        os.system(bt2build)
    else: bt2index = args.bt2

    lib_results = []
    samplecount = 0

    for fastq in args.input:
        samplecount += 1
        outdir = os.path.join(maindir, 'sample{}'.format(samplecount))
        os.system('mkdir -p {}'.format(outdir))
        outtemp = os.path.join(outdir, 'temp')
        os.system('mkdir -p {}'.format(outtemp))
        name = os.path.basename(fastq)
        logging.info('Start to process %s', name)
        logging.info('--------------------------------------------------------')

        guideseq = os.path.join(outtemp,
                                name.replace('.fastq',
                                             '_processed_guides.fastq'))
        bt2unmap = os.path.join(outtemp,
                                name.replace('.fastq', '_unmapped.fastq'))
        bt2map = os.path.join(outtemp, name.replace('.fastq', '_mapped.sam'))
        bt2stats = os.path.join(outdir, name.replace('.fastq', '_bt2stats.txt'))
        bt2map_bc = os.path.join(outtemp, name.replace('.fastq', '_bc.sam'))
        bt2map_bc_sorted = os.path.join(outtemp,
                                        name.replace('.fastq',
                                                     '_bc_qnsorted.sam'))
        bt2map_samsort_log = os.path.join(outtemp,
                                          name.replace('.fastq',
                                                       '_bc_qnsorted.log'))
        bt2map_dedup = os.path.join(outdir,
                                    name.replace('.fastq',
                                                 '_qnsorted_dedup.sam'))
        dedup_stats = os.path.join(outdir,
                                   name.replace('.fastq', '_picard_dedup.txt'))
        dedup_log = os.path.join(outtemp,
                                 name.replace('.fastq', '_picard_dedup.log'))
        gcounts = os.path.join(outdir, name.replace('.fastq', '_counts.txt'))
        obsguides = os.path.join(outtemp,
                                 name.replace('.fastq', '_observedG.txt'))
        ctbygene = os.path.join(outdir,
                                name.replace('.fastq', '_countsbygene.txt'))
        expguides = os.path.join(outtemp,
                                 name.replace('.fastq', '_expectedG.txt'))
        missguides = os.path.join(outdir, name.replace('.fastq', '_missG.txt'))
        missbygene = os.path.join(outdir,
                                  name.replace('.fastq', '_missbygene.txt'))
        ggfig = os.path.join(resultsum,
                             name.replace('.fastq', '_guidespergene.png'))
        mggfig = os.path.join(resultsum,
                              name.replace('.fastq', '_missguidespergene.png'))
        coverfig = os.path.join(resultsum,
                                name.replace('.fastq', '_coverage.png'))

        # Extract molecular barcode (randomer) to name; Will be used by picard
        # for removing duplicates in sam file
        logging.info('Pre-processing fastq')
        yltk.process_fastq(fastq,
                           bc_bp=args.randomer,
                           slice_bp=args.guide,
                           output=guideseq)

        # Align to expected guides
        bt2align = ('bowtie2 -p {} '
                    '-L 20 -N 0 --no-1mm-upfront -x {} {} '
                    '--un {} '
                    '-S {} '
                    '> {} 2>&1'
                    .format(args.thread,
                            bt2index, guideseq,
                            bt2unmap,
                            bt2map,
                            bt2stats))
        logging.info('Bowtie2 alignment')
        logging.info(bt2align + '\n')
        os.system(bt2align)

        # Attach barcodes to the BT&QT tags of samfile and remove duplicates
        # using picard
        logging.info('Removing duplicates with picard')
        yltk.attach_barcode(bt2map, bt2map_bc)
        picardsort = ('java -Xmx2g -jar $PICARD SortSam '
                      'I={} '
                      'O={} '
                      'SORT_ORDER=queryname '
                      '2>{}'
                      .format(bt2map_bc,
                              bt2map_bc_sorted,
                              bt2map_samsort_log))
        logging.info(picardsort)
        logging.info('Check %s for progress if this process is slow.',
                     bt2map_samsort_log)
        os.system(picardsort)
        dedup = ('java -Xmx2g -jar $PICARD MarkDuplicates '
                 'I={} '
                 'O={} '
                 'M={} '
                 'BARCODE_TAG=BC '
                 'REMOVE_DUPLICATES=True '
                 '2>{}'
                 .format(bt2map_bc_sorted,
                         bt2map_dedup,
                         dedup_stats,
                         dedup_log))
        logging.info(dedup)
        logging.info('Check %s for progress if this process is slow.%s',
                     dedup_log, '\n')
        os.system(dedup)

        # Count mapped reads from a SAM file
        countreads = ('samtools view -F 260 {} | '
                      'cut -f 3 | '
                      'sort | '
                      'uniq -c | '
                      'awk \'{{printf("%s\\t%s\\n", $2, $1)}}\' '
                      '> {}'
                      .format(bt2map_dedup,
                              gcounts))
        logging.info('Counting reads for guides')
        logging.info(countreads)
        os.system(countreads)

        # Find out missing guides and summarize a few stats
        logging.info('Performing statistical summary of results')
        findobs = 'cut -f 1 {} | sort > {}'.format(gcounts, obsguides)
        logging.info(findobs)
        os.system(findobs)
        countsum = ('awk -F\':\' \'{{print $1}}\' {} | '
                    'sort | '
                    'uniq -c | '
                    'awk \'{{print $2"\\t"$1}}\' '
                    '> {}'
                    .format(obsguides,
                            ctbygene))
        logging.info(countsum)
        os.system(countsum)
        findexpect = ('samtools view -H {} | '
                      'grep -P \'@SQ\\tSN:\' | '
                      'awk -F\'\\t\' \'{{print substr($2,4)}}\' | '
                      'sort '
                      '> {}'
                      .format(bt2map_dedup,
                              expguides))
        logging.info(findexpect)
        os.system(findexpect)
        obsvsexp = ('comm -13 --check-order {} {} > {}'
                    .format(obsguides, expguides, missguides))
        logging.info(obsvsexp)
        os.system(obsvsexp)
        missum = ('awk -F\':\' \'{{print $1}}\' {} | '
                  'sort | '
                  'uniq -c | '
                  'awk \'{{print $2"\\t"$1}}\' '
                  '> {}'
                  .format(missguides,
                          missbygene))
        logging.info(missum + '\n')
        os.system(missum)

        logging.info('Generating 3 figures to visualize results')
        # Figures 1 # of guides acquired by each gene
        gghist = (pd.read_table(ctbygene, header=None).iloc[:, 1]
                  .value_counts(sort=False).sort_index())
        barx = np.arange(gghist.size)
        plt.bar(barx, gghist, width=0.8)
        plt.xticks(barx+0.4, gghist.index)
        plt.ylabel('# of genes (with corresponding # of guides)')
        plt.xlabel('# of guides per gene')
        plt.title('Stats on number of guides each gene HAS in the library')
        plt.savefig(ggfig, format='png')

        plt.clf()
        plt.cla()

        # Figures 2 # of guides missed by each gene
        mgghist = (pd.read_table(missbygene, header=None).iloc[:, 1]
                   .value_counts(sort=False).sort_index())
        barx = np.arange(mgghist.size)
        plt.bar(barx, mgghist, width=0.8)
        plt.xticks(barx+0.4, mgghist.index)
        plt.ylabel('# of genes (missing corresponding # of guides)')
        plt.xlabel('# of guides missed')
        plt.title('Stats on number of guides MISSED by each gene in the '
                  'library')
        plt.savefig(mggfig, format='png')

        plt.clf()
        plt.cla()

        # Figures 3 histogram for # of reads per guide
        readhist = pd.read_table(gcounts, header=None).iloc[:, 1].values
        avg_rpg = readhist.sum()/readhist.size
        # If unique values of reads per guide are only a few (<=10), use bar
        # graph to show all real data instead of histogram
        histuniqvalues, histvaluecounts = np.unique(readhist,
                                                    return_counts=True)
        if histuniqvalues.size <= 10:
            barx = np.arange(histuniqvalues.size)
            plt.bar(barx, histvaluecounts, width=1)
            plt.xticks(barx+0.5, histuniqvalues)
            plt.xlabel('reads per guide (real counts)')
        else:
            plt.hist(readhist, bins='auto')
            plt.xlabel('reads per guide (auto bins)')
        plt.axvline(x=avg_rpg, color='r')
        plt.ylabel('# of guides (having certain reads per guide)')
        plt.title('Guide coverage distribution of the library')
        plt.savefig(coverfig, format='png')

        plt.clf()
        plt.cla()

        # Save checking results in the CRISPRLibSeqQC class and report it later
        # using jinja2 templating system
        logging.info('Saving and organizing results for final presentations in '
                     'HTML\n')
        libqc = CRISPRLibSeqQC(name, os.path.relpath(outdir, maindir))
        with open(bt2stats, 'r') as bt2statsfile:
            libqc.bt2stats = (bt2statsfile.read().replace('\n', '<br>')
                              .replace(' ', '&nbsp;'))
        picardstatstemp = []
        with open(dedup_stats, 'r') as picardstatsfile:
            for line in picardstatsfile:
                if line.strip() == ('## METRICS CLASS\t'
                                    'picard.sam.DuplicationMetrics'):
                    break
            for line in picardstatsfile:
                picardstatstemp.append(line.strip().split('\t'))
            for metrics_name, metrics_value in zip(picardstatstemp[0],
                                                   picardstatstemp[1]):
                libqc.picardstats += '{}:&nbsp;{}<br>\n'.format(metrics_name,
                                                                metrics_value)
        libqc.filenames['sam'] = os.path.relpath(bt2map_dedup, outdir)
        libqc.filenames['guidesCt'] = os.path.relpath(gcounts, outdir)
        libqc.filenames['guidespergene'] = os.path.relpath(ctbygene, outdir)
        libqc.filenames['guidemiss'] = os.path.relpath(missguides, outdir)
        libqc.filenames['misspergene'] = os.path.relpath(missbygene, outdir)
        libqc.figures['numguidespergene'] = os.path.relpath(ggfig, resultsum)
        libqc.figures['nummisspergene'] = os.path.relpath(mggfig, resultsum)
        libqc.figures['coverhist'] = os.path.relpath(coverfig, resultsum)
        lib_results.append(libqc)

        logging.info('%s analysis finished', name)
        logging.info('------------------------------------------------------\n')

    # Report all library checking results using jinja2 templating system
    logging.info('Making the final aggregated HTML report')
    shenshe2_searchpath = os.path.dirname(os.path.abspath(__file__))
    shenshe2_loader = jinja2.FileSystemLoader(searchpath=shenshe2_searchpath)
    shenshe2_env = jinja2.Environment(loader=shenshe2_loader)
    reporttemp = shenshe2_env.get_template('report_template.html')
    with open(os.path.join(resultsum, 'Report.html'), 'w') as reporthtml:
        current_date = str(datetime.datetime.today().date())
        reporthtml.write(reporttemp.render(date=current_date,
                                           mainfolder=maindir,
                                           results=lib_results))
    logging.info('Please find your results and the aggregated report at %s',
                 maindir)
    logging.info('CRISPR guide counting finished!')


if __name__ == '__main__':
    main()
