import pandas
dataset= "Cuffdiff output in Paralyzer format.csv"
ucsc_db=pandas.read_csv("UCSC_transcripts.txt",delim_whitespace=True,low_memory=False) # From the UCSC Database; a list of transcripts.
ucsc_names=pandas.read_csv('UCSC_genenames.csv')                      #To get meaningful names for the ucsc transcriptIds.



clusters=(pandas.read_csv(dataset))
dataset=dataset[:16]


print ucsc_names.head()
transcripts_list=[] 
        
for cluster in clusters.iterrows():
        c=cluster[1]
        cluster_transcripts=[]
        chromosome=c['Chromosome']
        clusterStart=c['ClusterStart']
        clusterEnd=c['ClusterEnd']
        strand = c['Strand']
        # transcripts_on_read_chromosome=ucsc_db[ucsc_db['chrom']==chromosome]
        # transcripts_on_read_chromosome.apply(match_check,axis=1)
        transcripts = ucsc_db[ucsc_db['chrom']==chromosome]
        transcripts = transcripts[transcripts['strand']==strand]
        transcripts = transcripts[transcripts['txStart']<clusterEnd]
        transcripts=transcripts[transcripts['txEnd']>clusterStart]
        transcripts_list.append(transcripts['#name'].tolist())
        print (transcripts['#name'])

       
clusters['transcripts']=pandas.Series(transcripts_list)
        
print clusters.head()
clusters.to_csv('peis_possible_transcript_matches.csv')
