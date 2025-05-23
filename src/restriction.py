# -*- coding: utf-8 -*-
## yufei.luo@institutimagine.org
## 23/10/2017
## Modified by Antoine (antoine.favier@institutimagine.org) --> convert from python 2 to Python 3
import sys

argvL = sys.argv
apprisF = argvL[1]
inGencode = argvL[2]
outGencode = argvL[3]
keep = argvL[4] # one if keep only one appris principal transcript per gene, all if keep all of them

def getAtt(s):
	return s.split(' ')[1].replace('"','')

def getGene(st):
	for s in st: 
		if s.startswith("gene_id"): 
			gene_id = getAtt(s).split('.')[0]
		elif s.startswith("transcript_id"):
			transcript_id  = getAtt(s).split('.')[0]
	return gene_id, transcript_id

def getTranscript(st): 
	for s in st: 
		if s.startswith("gene_id"): 
			gene_id = getAtt(s).split('.')[0]
		elif s.startswith("gene_name"):
			gene_name = getAtt(s).split('.')[0]
		elif s.startswith("transcript_id"):
			transcript_id  = getAtt(s).split('.')[0]
		elif s.startswith("exon_number"):
			exonNumber = getAtt(s)
	return gene_id+':'+gene_name, transcript_id, exonNumber

A = {} # APPRIS principal isoforms
with open(apprisF, "r") as af:
    for l in af:
        t = l.rstrip().split('\t')
        geneSymbol = t[0]
        gene = t[1]
        transcript = t[2]
        ccds_id = t[3].split('.')[0]
		# if ccds_id == '-': continue # don't keep the transcripts which are not having a ccds_id
        if keep == 'one':
            if gene not in A:
                A[gene] = {transcript:ccds_id}
            elif transcript not in A[gene]:
                A[gene][transcript] = ccds_id
        else:
            if gene not in A: #and gene not in R:
                A[gene] = {transcript:ccds_id}
            elif gene in A and transcript not in A[gene]:
                A[gene][transcript] = ccds_id

# Construct a dictionary of gene - transcript - length
D={} # {gene --> transcript --> transcript length}
with open(inGencode, "r") as f:
	for l in f:
		if l.startswith('#'): continue
		t = l.rstrip().split('\t')
		if "protein_coding" in t[8] and t[2] in ["transcript", "CDS", "start_codon", "stop_codon"]:
			st = t[8].split('; ')
			gene, transcript = getGene(st)
			if gene in A and transcript in A[gene]: # if is principal transcript
				if gene not in D: # normally should be transcript line
					D[gene] = {transcript:{'length':0, 'startCodon':0, 'stopCodon':0}}
				elif transcript not in D[gene]: # normally should be transcript line
					D[gene][transcript] = {'length':0, 'startCodon':0, 'stopCodon':0}
				elif t[2] == "CDS": # no need to do this but to be sure, it's for CCDS line
					length = int(t[4]) - int(t[3]) + 1
					D[gene][transcript]['length'] += length
				elif t[2] == "start_codon":
					D[gene][transcript]['startCodon'] = 1
				elif t[2] == "stop_codon":
					D[gene][transcript]['stopCodon'] = 1

# keep the appris principal transcripts
	for g in D: 
		maxLen = 0
		theTrans = ''
		if keep == "one":
			for t in list(D[g]):
				if D[g][t]['startCodon'] == 0 or D[g][t]['stopCodon'] == 0 or D[g][t]["length"]%3 != 0 or D[g][t]['length'] <= maxLen:
					del D[g][t]
				elif D[g][t]['length'] > maxLen:
					theTrans = t
					maxLen = D[g][t]['length']
		else:
			for t in list(D[g]):
				if D[g][t]['startCodon'] == 0 or D[g][t]['stopCodon'] == 0 or D[g][t]["length"]%3 != 0:
					del D[g][t]

# output the filtered gencode.cds.bed
with open(inGencode, "r") as f, open(outGencode, "w") as og: 
	for l in f:
		if l.startswith('#'): continue
		t = l.rstrip().split('\t')
		if t[0] != "chrM" and t[2] == "CDS": 
			st = t[8].split('; ')
			gene, transcript, exonNumber = getTranscript(st)
			gene_id = gene.split(':')[0]
			if gene_id in D and transcript in D[gene_id] and gene_id in A and transcript in A[gene_id]:
				chrom = t[0].replace('chr', '')
				start = t[3]
				end = t[4]
				strand = t[6]
				og.write("%s\n" % '\t'.join([chrom, start, end, transcript+':'+exonNumber, gene, strand]))
