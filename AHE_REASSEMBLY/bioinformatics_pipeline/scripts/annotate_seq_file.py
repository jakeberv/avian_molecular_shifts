import re
import glob
import gzip
import sys
import os
import subprocess
from operator import itemgetter
import argparse

parser = argparse.ArgumentParser(description='annotate')
parser.add_argument('--seq', help="sequence file to run this on")
args = parser.parse_args()

# protein file
ref = '/home/user/Desktop/Prum_data/loci_references/v2/Taeniopygia_guttata.bTaeGut1_v1.p.pep.all.fa'
# nucleotide file
uce = '/home/user/Desktop/Prum_data/loci_references/AHE_renamed.in_orientation.fasta'
seq = args.seq
# one-to-one matching file
mfile = '/home/user/Desktop/Prum_data/loci_references/v2/AHE_renamed.matches.csv'

def get_seq(seqfile):
	if re.search('gz$', seqfile):
		f = gzip.open(seqfile, 'r')
	else:
		f = open(seqfile, 'r')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			id = re.sub('\|.*$', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip().upper()
	f.close()
	return seq


def get_match(mfile):
	m = {}

	mf = open(mfile, 'r')
	for l in mf:
		d = re.split(',', l.rstrip())
		if d[1] != 'NA':
			m[d[0]] = d[1]
	mf.close()

	return m


def do_other(c, s, ref, m, o):
	query = '%s_query.fasta' % run
	subprocess.call("samtools faidx %s %s > %s" % (ref, m[c], query), shell=True)
	target = '%s_target.fasta' % run

        t = open(target, 'w')
        t.write('>%s\n%s\n' % (c, s))
        t.close()
	
	call = subprocess.Popen("exonerate -m protein2genome -q %s -t %s --showtargetgff TRUE -n 1 --showalignment NO --showvulgar 0 | grep 'insert'" % (query, target), shell=True, stdout=subprocess.PIPE)

	lines = [line.rstrip() for line in call.stdout]
	exons = []
	for l in lines:
		d = re.split('\t', l)
		exons.append([int(d[3]), int(d[4])])
	

	if len(exons) > 0:
		# need to sort exons!!!
		if exons[0][1] > exons[0][0]:
			exons = sorted(exons, key=itemgetter(0))
		else:
			print('oh no!!')
			exons = sorted(exons, key=itemgetter(1))

		if exons[0][0] > 1:
			o.write('%s,non-coding,1-%s\n' % (c, exons[0][0] - 1))
		if len(exons) > 1:
			for a, b in zip(exons, exons[1:]):
				if a[1] < b[0]:
					o.write('%s,non-coding,%s-%s\n' % (c, a[1] + 1, b[0] - 1))
		if exons[-1][1] < len(s):
			o.write('%s,non-coding,%s-%s\n' % (c,exons[-1][1] + 1, len(s)))

		for exon in exons:
			o.write('%s,exon,%s-%s\n' % (c, exon[0], exon[1]))
		
	else:
		if re.search('UCE', c.upper()):
			uce_edges(c, s, uce, o)
		else:
			o.write('%s,unknown,1-%s\n' % (c, len(s)))
	
	os.remove(target)
	os.remove(query)

def uce_edges(c, s, uce, o):	
	target = run + '_target.fa'
	t = open(target, 'w')
	t.write('>%s\n%s\n' % (c, s))
	t.close()

	query = run + '_query.fa'
	q = open(query, 'w')
	q.write('>%s_probe\n%s\n' % (c, uce[c]))
	q.close()

	call = subprocess.Popen("exonerate -q %s -t %s --showtargetgff TRUE -n 1 --showalignment NO --showvulgar 0 |grep 'similarity'" % (query, target), shell=True, stdout=subprocess.PIPE)
	lines = [line.rstrip() for line in call.stdout]

        ends = []
        for l in lines:
                d = re.split('\s+', l.rstrip())
                ends += [int(d[3]), int(d[4])]	

        if len(ends) > 0:
                start = min(ends)
                end = max(ends)

                if end < len(s) and start > 1:
                        o.write('%s,uce_flank,1-%s\n' % (c, start - 1))
			o.write('%s,uce_flank,%s-%s\n' % (c, end + 1, len(s)))
                elif end < len(s) and start == 1:
                        o.write('%s,uce_flank,%s-%s\n' % (c,end + 1, len(s)))
                elif end == len(s) and start > 1:
                        o.write('%s,uce_flank,1-%s\n' % (c,start - 1))
                o.write('%s,uce_core,%s-%s\n' % (c, start, end))
        else:
                o.write('%s,unknown,1-%s\n' % (c, len(s)))
        
        os.remove(query)
        os.remove(target)

run = re.sub('\.fa.*$', '', seq) 
m = get_match(mfile)
seq = get_seq(seq)
uce = get_seq(uce)
for c, s in uce.items():
	c = re.sub('_p\d$', '', c)
	c = re.sub('_c\d$', '', c)
	uce[c] = s

o = open(run + '_results', 'w')
for c, s in seq.items():
	if c in m:
		do_other(c, s, ref, m, o)
	else:
		if re.search('UCE', c.upper()):
			uce_edges(c, s, uce, o)
		else:
			pass
