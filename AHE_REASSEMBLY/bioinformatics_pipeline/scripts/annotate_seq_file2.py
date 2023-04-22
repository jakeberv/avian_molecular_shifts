import re
import glob
import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description='annotate')
parser.add_argument('--seq', help="sequence file to run this on")
args = parser.parse_args()

# protein file
ref = '/home/user/Desktop/Prum_data/loci_references/v2/Taeniopygia_guttata.bTaeGut1_v1.exons_pep.fa'
# sequence to annotate
seq = args.seq
# results file
res = re.sub('\.fa.*$', '', seq) + '_results'
res2 = re.sub('\.fa.*$', '', seq) + '_results2'
# one-to-one matching file
mfile = '/home/user/Desktop/Prum_data/loci_references/v2/AHE_renamed.matches.csv'
# sample
sample = re.sub('\.fa.*$', '', seq)
sample = re.sub('^.*\/', '', sample)


def get_match(mfile):
	m = {}

	mf = open(mfile, 'r')
	for l in mf:
		d = re.split(',', l.rstrip())
		if d[1] != 'NA':
			m[d[0]] = re.sub("\.\d+", "", d[1])
	mf.close()
	return m
	
def get_seq(seqfile):
	f = open(seqfile, 'r')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip().upper()
	f.close()
	return seq

def get_results(mfile):
	f = open(mfile, 'r')
	m = {}
	for l in f:
		d = re.split(',', l.rstrip())
		if d[0] not in m:
			m[d[0]] = {}
		ends = [int(x) for x in re.split('-', d[2])]
		m[d[0]][ends[0]] = {'name': None, 'type': d[1], 'end': ends[1]}
	return m
		
ref = get_seq(ref)
exons = {}
for id, s in ref.items():
	ids = re.split('_', id)
	if ids[0] not in exons:
		exons[ids[0]] = {}
	exons[ids[0]][int(ids[1])] = s
m = get_results(res)
seq = get_seq(seq)
refm = get_match(mfile)

for c in m:
        s = seq[c]
	# print out exons
	qfile = '%s_%s_exons' % (c, sample)
	q = open(qfile, 'w')
	for start in m[c]:
		if m[c][start]['type'] == 'exon':
			tmps = s[(start - 1): m[c][start]['end']]
			q.write('>%s\n%s\n' % (start, tmps))
	q.close()
	
	# print out ref
	dbfile = '%s_%s_refexons' % (c, sample)
	r = open(dbfile, 'w')
	refc = refm[c]
	for exon, s in exons[refc].items():
		r.write('>%s\n%s\n' % (exon, s))
	r.close()

	cm = {}
	subprocess.call("makeblastdb -in %s -dbtype prot" % dbfile, shell = True)
	call = subprocess.Popen("blastx -db %s -query %s -outfmt=10 -evalue 1e-10" % (dbfile, qfile), shell=True, stdout=subprocess.PIPE)
        bm = {}
	for line in call.stdout:
		d = re.split(',', line.rstrip())
		d0 = int(d[0])
		if not m[c][d0]['name']:
			m[c][d0]['name'] = int(d[1])
                # if d0 not in bm:
                #        bm[d0] = {'match': int(d[1]), 'val': float(d[10])}
        # for d0 in bm:
        #        match = bm[d0]['match']
        #        eval = bm[d0]['val']

        #       winner = eval
        #        for x in bm:
        #                if bm[x]['match'] == match:
        #                        if bm[x]['val'] < winner:
        #                                winner = bm[x]['val']
        #
        #        if eval == winner:
        #                m[c][d0]['name'] = int(d[1])
                                

	starts = sorted(m[c].keys())
	for ix, start in enumerate(starts):
		if ix == 0 and m[c][start]['type'] == 'non-coding':
			# utr
			if m[c][starts[ix + 1]]['name'] == None:
                                m[c][start]['name'] = None
			elif m[c][starts[ix + 1]]['name'] == 1:
				m[c][start]['name'] = "5utr"
			else:
				m[c][start]['name'] = 'intron%s' % (m[c][starts[ix + 1]]['name'] - 1)
		# last one
		elif ix == (len(starts) - 1) and m[c][start]['type'] == 'non-coding':
			# utr                                                                                  
			if m[c][starts[ix - 1]]['name'] == None:
				m[c][start]['name'] =	None
			elif m[c][starts[ix - 1]]['name'] == max(list(exons[refc].keys())):
				m[c][start]['name'] = "3utr"
			else:
				m[c][start]['name'] =	'intron%s' % m[c][starts[ix - 1]]['name']
		elif m[c][start]['type'] == 'non-coding':
			pre = m[c][starts[ix - 1]]['name']
			post = m[c][starts[ix + 1]]['name']
			if pre == None or post == None:
				m[c][start]['name'] = None
			elif post == pre + 1:
				m[c][start]['name'] = 'intron%s' % pre
			else:
				m[c][start]['name'] = None
	
	# cleanup
        [os.remove(x) for x in glob.glob("%s_%s_*" % (c, sample))]
        
r = open(res2, 'w')
for c in m:
        for start in sorted(m[c].keys()):
                r.write('%s,%s,%s,%s-%s\n' % (c, m[c][start]['type'], m[c][start]['name'], start, m[c][start]['end']))
r.close()
