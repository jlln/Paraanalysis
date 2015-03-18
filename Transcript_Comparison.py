import sys
import pandas
import numpy as np 

#Takes two arguments (two directories) and does set comparisons on the matched reads.


directory1=sys.argv[1]
directory2=sys.argv[2]

ucsc_db=pandas.read_table('UCSC_transcripts.txt')


def get_transcripts(cluster):
	transcripts=eval(cluster['transcripts'])
	list_transcripts.extend(transcripts)


datasets=[]
translation_status=pandas.read_csv(directory1+"_Translation_Status_Totals.csv")
for x in translation_status.columns[1:]:
	datasets.append(x)


transcript_dict={}

list_transcripts=[]
for d in datasets:
	
	data=pandas.read_csv(directory1+d+'matchedreads.csv')
	data.apply(get_transcripts,axis=1)
transcripts=set(list_transcripts)
transcripts_one=transcripts




datasets=[]
translation_status=pandas.read_csv(directory2+"_Translation_Status_Totals.csv")
for x in translation_status.columns[1:]:
	datasets.append(x)





list_transcripts=[]
for d in datasets:
	data=pandas.read_csv(directory2+d+'matchedreads.csv')
	data.apply(get_transcripts,axis=1)
transcripts=set(list_transcripts)
transcripts_two=transcripts




intersection=transcripts_one.intersection(transcripts_two)
union=transcripts_one.union(transcripts_two)
print "Transcripts"
print directory1
print len(transcripts_one)
print 'intersection'
print len(intersection)
print len(transcripts_two)



gene_db=pandas.read_csv('UCSC_genenames.csv')

transcripts_one=list(transcripts_one)
transcripts_two=list(transcripts_two)

genes_one=[]
genes_two=[]

for x in transcripts_one:
	gene=gene_db[gene_db['#kgID']==x]
	genes_one.append(gene['geneSymbol'].values[0])
for x in transcripts_two:
	gene=gene_db[gene_db['#kgID']==x]
	genes_two.append(gene['geneSymbol'].values[0])


genes_one=set(genes_one)
genes_two=set(genes_two)
genes_intersect=genes_one.intersection(genes_two)


print "Genes"
print directory1
print len(genes_one)
print "intersection"
print len(genes_intersect)
print directory2
print len(genes_two)



