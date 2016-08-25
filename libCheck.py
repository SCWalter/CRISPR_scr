#!/usr/bin/env python

import argparse, os, logging
import yltk

import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import jinja2
import datetime

class CRISPRlibQ:
	### Store info for the result presentation
	def __init__(self, name, out):
	### Create instance variables instead of class variables for class variables are shared across all instances while instance variables are instance specific.
	### Those class variables which are set by immutable value, such as numbers, strings and tuples, behave like instance specific. However,
	### this is not a correct practice since we really/actually need instance variables here.
	### Maybe using objects instead of dictionaries for jinja2 is not a correct practice either...
		self.samplename = name
		self.outfolder = out
		self.bt2stats = ''
		self.picardstats = ''
		self.filenames = {'sam':'', 'guidesCt':'', 'guidespergene':'', 'guidemiss':'', 'misspergene':''}
		self.figures = {'numguidespergene':'', 'nummisspergene':'', 'coverhist':''}
	
parser = argparse.ArgumentParser(description='Check the quality (guides representation) of CRISPR library')
parser.add_argument('-i', '--input', help='Illumina fastq file(s)', nargs='+', required=True)
bt2index = parser.add_mutually_exclusive_group(required=True)
bt2index.add_argument('-t', '--table', help='Table with guide infomation. It will be used for making bowtie2 index in the same directory as this table. The table should have 3 columns: 1) Gene name; 2) Guides name/number; 3) Guide sequence. The first row will be discarded as header which then followed by one guides per row. A fake header is required if there is no header yet.', default=None)
bt2index.add_argument('--bt2', help='Bowtie2 index for expected guides.', default=None)
parser.add_argument('-r', '--randomer', help='Basepair position for randomer, which will be used to remove PCR duplicates. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".', default=None)
parser.add_argument('-g', '--guide', help='Basepair position for 20bp CRISPR guide. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".', default=None)
parser.add_argument('-o', '--outpath', help='Path to put outputs (with slash). Default is the same directory as the first input fastq.', default='')

args = parser.parse_args()
logging.basicConfig(filename='CRISPR_lib_Check.log',level=logging.DEBUG)

if not args.bt2:
	guidetab_abspath = os.path.abspath(args.table)
	guidetab_dir = os.path.dirname(guidetab_abspath)
	guidetab_base = os.path.basename(guidetab_abspath)
	bt2indexdir = os.path.join(guidetab_dir, 'bowtie2index')
	os.system("mkdir -p {}".format(bt2indexdir))
	bt2index = os.path.join(bt2indexdir, guidetab_base)
	fasta = os.path.join(bt2indexdir, guidetab_base + ".fasta")

	tab2fasta = "tail -n +2 {} | awk -F'\\t' '{{print \">\"$1\":\"$2\"\\n\"$3}}' > {}".format(guidetab_abspath, fasta)
	logging.info('Generating fasta from the table of guides')
	logging.info(tab2fasta + '\n')
	os.system(tab2fasta)

	bt2build = "bowtie2-build {} {}".format(fasta, bt2index)
	logging.info('Making Bowtie2 index from fasta file')
	logging.info(bt2build + '\n')
	os.system(bt2build)
else: bt2index = args.bt2

if not args.outpath:
	maindir = os.path.join(os.path.dirname(os.path.abspath(args.input[0])), 'CRISPR_libCheck')
else:
	maindir = os.path.join(args.outpath, 'CRISPR_libCheck')
os.system("mkdir -p {}".format(maindir))
resultsum = os.path.join(maindir, 'View_Result_Summary')
os.system("mkdir -p {}".format(resultsum))
libResults = []
samplecount = 0

for fastq in args.input:
	samplecount += 1
	outdir = os.path.join(maindir, 'sample{}'.format(samplecount))
	os.system("mkdir -p {}".format(outdir))
	outtemp = os.path.join(outdir, 'temp')
	os.system("mkdir -p {}".format(outtemp))
	name = os.path.basename(fastq)
	logging.info('Start to process {}'.format(name))
	logging.info('-------------------------------------------------------------------')

	guideseq = os.path.join(outtemp, name.replace(".fastq", "_processed_guides.fastq"))
	bt2unmap = os.path.join(outtemp, name.replace(".fastq", "_unmapped.fastq"))
	bt2map = os.path.join(outtemp, name.replace(".fastq", "_mapped.sam"))
	bt2stats = os.path.join(outdir, name.replace(".fastq", "_bt2stats.txt"))
	bt2map_bc = os.path.join(outtemp, name.replace(".fastq", "_bc.sam"))
	bt2map_bc_sorted = os.path.join(outtemp, name.replace(".fastq", "_bc_qnsorted.sam"))
	bt2map_dedup = os.path.join(outdir, name.replace(".fastq", "_qnsorted_dedup.sam"))
	dedup_stats = os.path.join(outdir, name.replace(".fastq", "_picard_dedup.txt"))
	gcounts = os.path.join(outdir, name.replace(".fastq", "_counts.txt"))
	obsguides = os.path.join(outtemp, name.replace(".fastq", "_observedG.txt"))
	ctbygene = os.path.join(outdir, name.replace(".fastq", "_countsbygene.txt"))
	expguides = os.path.join(outtemp, name.replace(".fastq", "_expectedG.txt"))
	missguides = os.path.join(outdir, name.replace(".fastq", "_missG.txt"))
	missbygene = os.path.join(outdir, name.replace(".fastq", "_missbygene.txt"))
	ggfig = os.path.join(resultsum, name.replace(".fastq", "_guidespergene.png"))
	mggfig = os.path.join(resultsum, name.replace(".fastq", "_missguidespergene.png"))
	coverfig = os.path.join(resultsum, name.replace(".fastq", "_coverage.png"))

	### Extract molecular barcode (randomer) to name; Will be used by picard for removing duplicates in sam file
	logging.info('Pre-processing fastq')
	yltk.process_fastq(fastq, bc_bp = args.randomer, slice_bp = args.guide, output = guideseq)
	
	### Align to expected guides
	noMMalign = "bowtie2 -p 8 -L 20 -N 0 --no-1mm-upfront -x {} {} --un {} -S {} > {} 2>&1".format(bt2index, guideseq, bt2unmap, bt2map, bt2stats)
	logging.info('Bowtie2 alignment')
	logging.info(noMMalign + '\n')
	os.system(noMMalign)
	
	### Attach barcodes to the BT&QT tags of samfile and remove duplicates using picard
	logging.info('Removing duplicates with picard')
	yltk.attach_barcode(bt2map, bt2map_bc)
	picardsort = "java -Xmx2g -jar $PICARD SortSam I={} O={} SORT_ORDER=queryname".format(bt2map_bc, bt2map_bc_sorted)
	logging.info(picardsort)
	os.system(picardsort)
	dedup = "java -Xmx2g -jar $PICARD MarkDuplicates I={} O={} M={} BARCODE_TAG=BC REMOVE_DUPLICATES=True".format(bt2map_bc_sorted, bt2map_dedup, dedup_stats)
	logging.info(dedup + '\n')
	os.system(dedup)
	
	### Count mapped reads from a SAM file
	countreads = "samtools view -F 260 {} | cut -f 3 | sort | uniq -c | awk '{{printf(\"%s\\t%s\\n\", $2, $1)}}' > {}".format(bt2map_dedup, gcounts)
	logging.info('Counting reads for guides')
	logging.info(countreads)
	os.system(countreads)
	
	### Find out missing guides and summarize a few stats
	logging.info('Performing statistical summary of results')
	findobs = "cut -f 1 {} | sort > {}".format(gcounts, obsguides)
	logging.info(findobs)
	os.system(findobs)
	countsum = "awk -F':' '{{print $1}}' {} | sort | uniq -c | awk '{{print $2\"\\t\"$1}}' > {}".format(obsguides, ctbygene)
	logging.info(countsum)
	os.system(countsum)	
	findexpect = "samtools view -H {} | grep -P '@SQ\\tSN:' | awk -F'\\t' '{{print substr($2,4)}}' | sort > {}".format(bt2map_dedup, expguides)
	logging.info(findexpect)
	os.system(findexpect)
	obsvsexp = "comm -13 --check-order {} {} > {}".format(obsguides, expguides, missguides)
	logging.info(obsvsexp)
	os.system(obsvsexp)
	missum = "awk -F':' '{{print $1}}' {} | sort | uniq -c | awk '{{print $2\"\\t\"$1}}' > {}".format(missguides, missbygene)
	logging.info(missum + '\n')
	os.system(missum)	

	logging.info('Generating 3 figures to visualize results')
	### Figures 1 # of guides acquired by each gene
	gghist = pd.read_table(ctbygene, header=None).iloc[:,1].value_counts(sort=False).sort_index()
	barx = np.arange(gghist.size)
	plt.bar(barx, gghist, width=0.8)
	plt.xticks(barx+0.4, gghist.index)
	plt.ylabel('# of genes (with corresponding # of guides)')
	plt.xlabel('# of guides per gene')
	plt.title('Stats on number of guides each gene HAS in the library')
	plt.savefig(ggfig, format='png')

	plt.clf()
	plt.cla()

	### Figures 2 # of guides missed by each gene
	mgghist = pd.read_table(missbygene, header=None).iloc[:,1].value_counts(sort=False).sort_index()
	barx = np.arange(mgghist.size)
	plt.bar(barx, mgghist, width=0.8)
	plt.xticks(barx+0.4, mgghist.index)
	plt.ylabel('# of genes (missing corresponding # of guides)')
	plt.xlabel('# of guides missed')
	plt.title('Stats on number of guides MISSED by each gene in the library')
	plt.savefig(mggfig, format='png')
	
	plt.clf()
	plt.cla()

	### Figures 3 histogram for # of reads per guide
	readhist = pd.read_table(gcounts, header=None).iloc[:,1].values
	avgRpG = readhist.sum()/readhist.size
	histuniqvalues, histvaluecounts = np.unique(readhist,return_counts=True)
	if histuniqvalues.size <= 10:
		barx = np.arange(histuniqvalues.size)
		plt.bar(barx, histvaluecounts, width=1)
		plt.xticks(barx+0.5, histuniqvalues)
		plt.xlabel('reads per guide (real counts)')
	else:
		plt.hist(readhist, bins='auto')
		plt.xlabel('reads per guide (auto bins)')
	plt.axvline(x=avgRpG, color='r')
	plt.ylabel('# of guides (having certain reads per guide)')
	plt.title('Guide coverage distribution of the library')
	plt.savefig(coverfig, format='png')
	
	plt.clf()
	plt.cla()

	### Save checking results in the CRISPRlibQ class and report it later using jinja2 templating system
	logging.info('Saving and organizing results for final presentations in HTML\n')
	libqc = CRISPRlibQ(name, os.path.relpath(outdir, maindir))
	with open(bt2stats, 'r') as bt2statsfile:
		libqc.bt2stats = bt2statsfile.read().replace('\n', '<br>').replace(' ', '&nbsp;')
	picardstatstemp = []
	with open(dedup_stats, 'r') as picardstatsfile:
		for line in picardstatsfile:
			if line.strip() == '## METRICS CLASS	picard.sam.DuplicationMetrics':
				break
		for line in picardstatsfile:
			picardstatstemp.append(line.strip().split('\t'))
		for n, v in zip(picardstatstemp[0],picardstatstemp[1]):
			libqc.picardstats += '{}:&nbsp;{}<br>\n'.format(n, v)
	libqc.filenames['sam'] = os.path.relpath(bt2map_dedup, outdir)
	libqc.filenames['guidesCt'] = os.path.relpath(gcounts, outdir)
	libqc.filenames['guidespergene'] = os.path.relpath(ctbygene, outdir)
	libqc.filenames['guidemiss'] = os.path.relpath(missguides, outdir)
	libqc.filenames['misspergene'] = os.path.relpath(missbygene, outdir)
	libqc.figures['numguidespergene'] = os.path.relpath(ggfig, resultsum)
	libqc.figures['nummisspergene'] = os.path.relpath(mggfig, resultsum)
	libqc.figures['coverhist'] = os.path.relpath(coverfig, resultsum)
	libResults.append(libqc)

	logging.info('{} analysis finished'.format(name))
	logging.info('-------------------------------------------------------------------\n')

### Report all library checking results using jinja2 templating system
logging.info('Making the final aggregated HTML report')
reporttemp = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=os.path.dirname(os.path.abspath(__file__)))).get_template('report_template.html')
with open(os.path.join(resultsum, 'Report.html'), 'w') as reporthtml:
	reporthtml.write(reporttemp.render(date=str(datetime.datetime.today().date()), mainfolder = maindir, results=libResults))
logging.info('Please find your results and the aggregated report at {}'.format(maindir))
logging.info('CRISPR library check finished!')
