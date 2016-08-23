#!/usr/bin/env python

import argparse
import itertools, pysam

def fastq_reads(fastq):
	f = open(fastq)
	while True:
		read = list(itertools.islice(f, 4))
		if not read:
			f.close()
			return
		else: yield read

def prepslices(idcs_str):
	idcs = []
	for idx_str in idcs_str.split(','):
		if ':' in idx_str:
			a, b = idx_str.split(':')
			idcs.append((int(a)-1,int(b)))
		else: idcs.append((int(idx_str)-1,int(idx_str)))
	return idcs

def multislices(seq, idcs):
	res = ''
	for ia,ib in idcs:
		res += seq[ia:ib]
	return res

def seqmatch(query, ref, queryallowN=False):
	if len(query) != len(ref): return 0
	for q, r in zip(query, ref):
		if q == r: continue
		elif r == 'N': continue
		elif (q == 'N') & queryallowN: continue
		else: return 0
	return 1

def prepdemux(output):
	idxfiledict = []
	for fileidcs in output.split(';'):
		filename, idcs = fileidcs.split(':')
		file = open(filename, 'w')
		for idx in idcs.split(','):
			idxfiledict[idx] = file
	return idxfiledict

def demultiplx(idx_seq, outfiles, fnoidx):
	for fidcs in outfiles:
		for idx in fidcs[1:]:
			if seqmatch(idx_seq, idx): return fidcs[0]
	return fnoidx

def process_fastq(fastq, bc_bp=None, slice_bp=None, idx_bp=None, output=None):
	if bc_bp: bc_idcs = prepslices(bc_bp)
	if slice_bp: slice_idcs = prepslices(slice_bp)
	if idx_bp: 
		idx_idcs = prepslices(idx_bp)
		filedict = prepdemux(output)
		fnoidx = open(fastq.replace('.fastq','_noidx.fastq'), 'w')
	elif not output:
		output = fastq.replace('.fastq','_proced.fastq')
		outfile = open(output,'w')
	else:
		outfile = open(output,'w')

	for read in fastq_reads(fastq):
		seq_id, seq, pqs_id, pqs = read
		seq.rstrip('\n')
		pqs.rstrip('\n')
		if bc_bp:
			bc_seq_id = '@' + ':'.join([multislices(seq, bc_idcs), multislices(pqs, bc_idcs), seq_id[1:]])
			seq_id = bc_seq_id
		if slice_bp:
			seq_proced = multislices(seq, slice_idcs)
			pqs_proced = multislices(pqs, slice_idcs)
		else:
			seq_proced = seq
			pqs_proced = pqs
		if idx_bp:
			idx_seq = multislices(seq, idx_idcs)
			outfile = fnoidx
			for idx in filedict:
				if seqmatch(idx_seq, idx):
					outfile = filedict[idx]
					break
		outfile.write(''.join([seq_id, seq_proced + '\n', pqs_id, pqs_proced + '\n']))

	if idx_bp:
		for file in set(filedict.values()):
			file.close()
	else: outfile.close()

def attach_barcode(sam, output):
	if output is None:
		output = sam.replace('.sam','_bcs.sam')
	infile = pysam.AlignmentFile(sam, "r")
	outfile = pysam.AlignmentFile(output, "wh", template=infile)
	for read in infile.fetch():
		id_sam = read.query_name
		sep_si = id_sam.index(':')
		BC_seq = id_sam[0:sep_si]
		sep_qi = sep_si + 1 + len(BC_seq)
		BC_pqs = id_sam[sep_si + 1: sep_qi]
		read.set_tag('BC', BC_seq)
		read.set_tag('QT', BC_pqs)
		read.query_name = id_sam[sep_qi+1:]
		outfile.write(read)
	outfile.close()
	infile.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Extract barcodes from fastq and attach barcodes back to aligned sam')
	subparsers = parser.add_subparsers(dest = 'subcom')

	parser_extract = subparsers.add_parser('process_fastq', help='Extract barcodes from one fastq file and save them in the read name.')
	parser_extract.add_argument('-f', '--fastq', help='One fastq file to be processed', required=True)
	parser_extract.add_argument('-b', '--barcode', help='Basepair position for molecular barcode (randomer) which help remove duplicates. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".', default=None)
	parser_extract.add_argument('-s', '--slice', help='Basepair position for target sequences to keep. Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".', default=None)
	parser_extract.add_argument('-i', '--index', help='Basepair position for indexes which help demultiplex reads. File names for corresponding reads should be specified in the "outfastq" option (check it for details). Support both discrete and continuous regions. Discrete regions are separated by ",", and continuous region is expressed as "a:b".', default=None)
	parser_extract.add_argument('--outfastq', help='Output fastq file. For single file output, the default is to add a "_proced" postfix and save in the same directory. For demultiplexing, the expected format is "filename1:index1-1,index1-2,index1-3,...;filename2:index2-1,index2-2,index2-3,...;...". "N" is supported in indexes but not sequencing results', default=None)

	parser_attach = subparsers.add_parser('attach_sam', help='Attach barcodes to the BC&QT tag in one aligned sam file.')
	parser_attach.add_argument('-s', '--sam', help="One aligned sam file which needs to attach barcodes at BC&QT tags.", required=True)
	parser_extract.add_argument('--outsam', help='Output unsorted sam file. Default is to add a "_bcs" postfix and save in the same directory.', default=None)

	args = parser.parse_args()
	if args.subcom == 'process_fastq':
		process_fastq(args.fastq, bc_bp=args.b, slice_bp=args.s, idx_bp=args.i, output=args.outfastq)
	elif args.subcom == 'attach_sam':
		attach_barcode(args.sam, args.outsam)
	else:
		parser.print_help()
		exit(1)
