import re
import argparse
import gzip
import os
import copy

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--lineage", required=True, help='lineage to do.')
parser.add_argument("-b", "--basedir", required=True, help='base dir')
parser.add_argument('--qual', type=int, default=20, help='Minimum quality to retain variant for creating final call set.')
parser.add_argument('--dp', type=int, default=10, help='Minimum depth to retain variant for creating final call set.')
args = parser.parse_args()
lineage = args.lineage

prg = os.path.join(args.basedir, 'PRG', '%s.fasta' % lineage)
vcf = os.path.join(args.basedir, 'variants', '%s.qual_filtered%s.cov_filtered%s.vcf.gz' % (lineage, args.qual, args.dp))

def get_seq(f):
	seq = {}
	id = ''
	
	f = open(f, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	mut = {}
	for id, s in seq.items():
		seq[id] = list(s)
		mut[id] = list(s.lower())

	return seq, mut

codes = {'A': {'A': 'A', 'T': 'W', 'C': 'M', 'G': 'R'}, 
	'T': {'A': 'W', 'T': 'T', 'C': 'Y', 'G': 'K'},
	'C': {'A': 'M', 'T': 'Y', 'C': 'C', 'G': 'S'},
	 'G': {'A': 'R', 'T': 'K', 'C': 'S', 'G': 'G'}, 'N': {'N': 'N'}}


seq, mut = get_seq(prg)
inds = []

v = gzip.open(vcf, 'r')
for l in v:
	if re.search('^#CHROM', l):
		d = re.split('\t', l.rstrip())
		inds = d[9:]
		mut2 = {}
		for ind in inds:
			mut2[ind] = mut
	elif not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		pos = int(d[1])
		c = d[0]

		ref = d[3]
		genos = d[9:]
		aa = [d[3]] + re.split(',', d[4])
                alleles = {}
                for ix, a in enumerate(aa):
                        alleles[str(ix)] = a
                alleles['.'] = 'N'

		for ind, geno in zip(inds, genos):
			if re.search('^0/0', geno):
				mut2[ind][c][pos - 1] = ref
			else:
				genox = re.search('^(\S\S\S)', geno).group(1)
				genox = [alleles[x] for x in re.split('/', genox)]
				mut2[ind][c][pos - 1] = codes[genox[0]][genox[1]]
v.close()


outdir = os.path.join(args.basedir, 'variable_PRG')
if not os.path.exists(outdir): 
	os.makedirs(outdir)

for ind in inds:
	out = os.path.join(outdir, '%s.qual_filtered%s.cov_filtered%s.fasta' % (ind, args.qual, args.dp))
	o = open(out, 'w')
	for c in seq:
		mutseq = ''.join(mut2[ind][c])
		# mutseq = re.sub('^N+', '', mutseq)
		# mutseq = re.sub('N+$', '', mutseq)
		o.write('>%s\n%s\n' % (c, mutseq))
	o.close()
