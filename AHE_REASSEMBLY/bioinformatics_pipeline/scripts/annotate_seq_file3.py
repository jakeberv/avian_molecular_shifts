import os
import glob
import re
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='annotate')
parser.add_argument('--outdir', help="output dir")
parser.add_argument('--file', help="csv file with the info")
parser.add_argument('--best', action = "store_true", default = False, help= "best match only")
parser.add_argument('--haplo', action = "store_true", default = False, help = "use haplotypes instead of diplotypes")
parser.add_argument('--qual', type=int, default=20, help='Minimum quality to retain variant for creating final call set.')
parser.add_argument('--dp', type=int, default=10, help='Minimum depth to retain variant for creating final call set.')
args = parser.parse_args()
best = args.best
haplo = args.haplo

bdir = '/home/user/Desktop/Prum_data/'

def get_seq(seqfile):
    f = open(seqfile, 'r')

    seq = {}
    id = ''
    for l in f:
        if re.search('>', l):
            id = re.search('>(\S+)', l.rstrip()).group(1)
            seq[id] = ''
        else:
            seq[id] += l.rstrip()
    f.close()
    return seq


def get_seq2(seqfile, sample):
    locus = re.sub('.*/', '', seqfile)
    locus = re.sub('.aln.fasta', '', locus)
    
    f = open(seqfile, 'r')

    seq = {}
    id = ''
    for l in f:
        if re.search('>', l):
            id = re.search('>(\S+)', l.rstrip()).group(1)
            seq[id] = ''
        else:
            seq[id] += l.rstrip()
    f.close()

    seq2 = {}
    for id, s in seq.items():
        id = re.split('_', id)
        if id[0] == sample:
            seq2[id[1]] = s
    
    return seq2, locus


def get_results(resfile):
    f = open(resfile, 'r')
    r = {}
    for l in f:
        d = re.split(',', l.rstrip())
        if d[0] not in r:
            r[d[0]] = {}
        lim = [int(x) for x in re.split('-', d[3])]
        # best
        c = d[0]
        name = d[2]
        if name in r[c]:
            if best:
                curlen = r[c][name]['end'][0] - r[c][name]['start'][0]
                if (lim[1] - lim[0]) > curlen:
                    r[d[0]][name] = {'name': d[2], 'type': d[1], 'start': [(lim[0] - 1)], 'end': [lim[1]]}
            else: 
                r[d[0]][name]['start'].append((lim[0] - 1))
                r[d[0]][name]['end'].append(lim[1])
        else:
            r[d[0]][name] = {'name': d[2], 'type': d[1], 'start': [(lim[0] - 1)], 'end': [lim[1]]}
    f.close()
    return r


ff = {}

d = pd.read_csv(args.file)
for ix, row in d.iterrows():
    sample = row['sample']
    lineage = row['lineage']
    subdir = row['run']
    
    resfile = os.path.join(bdir, subdir, 'PRG', '%s_results2' % lineage)
    # this variable can be flagged to only
    # pick best match or all matches
    r = get_results(resfile)

    seq  = {}
    # later add in the flag for diplotype
    if haplo:
        hapdir = os.path.join(bdir, subdir, 'variants', '%s_qual%s_cov%s' % (lineage, args.qual, args.dp))
        hapfiles = glob.glob(hapdir + '/*fasta')
        for hapfile in hapfiles:
            tmpseq, locus = get_seq2(hapfile, sample)
            seq[locus] = tmpseq
    else:
        seqfile = os.path.join(bdir, subdir, 'variable_PRG', '%s.qual_filtered%s.cov_filtered%s.fasta' % (sample, args.qual, args.dp))
        if os.path.isfile(seqfile):
            seq = get_seq(seqfile)
            
    for c in r:
        for start in r[c]:
            if r[c][start]['type'] == 'exon':
                fname = '%s_%s%s.fasta' % (c, r[c][start]['type'], r[c][start]['name'])
            else:
                if r[c][start]['name'] == 'None':
                    fname = '%s_noncodingUnknown.fasta' % c
                else:
                    fname = '%s_%s.fasta' % (c, r[c][start]['name'])
            if fname not in ff:
                ff[fname] = {}

            if haplo:
                tmpseq1 = ''
                tmpseq2 = ''
                for a, b in zip(r[c][start]['start'], r[c][start]['end']):
                    tmpseq1 += seq[c]['1'][a:b]
                    tmpseq2 += seq[c]['2'][a:b]
	        ff[fname][sample + '_1'] = tmpseq1
                ff[fname][sample + '_2'] = tmpseq2
            else:
                tmpseq = ''
                for a, b in zip(r[c][start]['start'], r[c][start]['end']):
                    tmpseq += seq[c][a:b]
                ff[fname][sample] = tmpseq

outdir2 = os.path.join(args.outdir, 'qual%s_cov%s' % (args.qual, args.dp))
if args.haplo:
    outdir2 += '_haplo'
else:
    outdir2 += '_diplo'
if args.best:
    outdir2 += '_bestonly'
if not os.path.isdir(outdir2):
    os.mkdir(outdir2)
    
for fname in ff:
    fname2 = os.path.join(outdir2, fname)
    f = open(fname2, 'w')
    for c, s in ff[fname].items():
        f.write('>%s\n%s\n' % (c, s))
    f.close()
