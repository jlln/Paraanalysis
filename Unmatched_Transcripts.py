#Finds the distance from an unmatched cluster to the next transcript start site on the chromosome.
import pandas
import sys
import operator

ucsc_db=pandas.read_table('UCSC_transcripts.txt')

directory=sys.argv[1]

datasets=[]
translation_status_df=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status_df.columns[1:]:
	datasets.append(x)

list_of_bound_transcripts=[]
def accumulate_transcripts(cluster,list_of_transcripts):
	list_of_transcripts.extend(eval(cluster['transcripts']))

def find_closest_downstream_transcript(cluster,dataset):
	cluster_end=cluster['ClusterEnd']
	Chromosome=cluster['Chromosome']
	transcripts_on_chromosome=ucsc_db[ucsc_db['chrom']==Chromosome]
	downstream_transcripts=transcripts_on_chromosome[transcripts_on_chromosome['txStart']>cluster_end].sort_index(by="txStart")
	
	distances.append(downstream_transcripts.iloc[0]['txStart']-cluster_end)
	cluster_id=cluster['ClusterID']
	cluster_start=cluster['ClusterStart']
	transcript=downstream_transcripts.iloc[0]['#name']
	annotation=pandas.read_csv(directory+transcript+'annotations.csv')
	print annotation

	



for d in datasets:
	clusters=pandas.read_csv(directory+d+'matchedreads.csv')
	distances=[]
	unmatched_clusters=clusters[clusters['transcripts']=="[]"]
	unmatched_clusters.apply(find_closest_downstream_transcript,axis=1,args=([d]))
	df_d={d:distances}
	df=pandas.DataFrame(df_d)
	df.to_csv(directory+d+"unmatched_clusters_distance_to_nearest_downstream_transcript.csv")



# clusters1=pandas.read_csv(directory+datasets[0]+'matchedreads.csv')
# clusters1['Dataset']=datasets[0]


# clusters2=pandas.read_csv(directory+datasets[1]+'matchedreads.csv')
# clusters2['Dataset']=datasets[1]

# clusters=pandas.concat([clusters1,clusters2])
# displacement=[]
# def pair_up(cluster):
# 	chrom=cluster['Chromosome']
# 	start=cluster['ClusterStart']
# 	dataset=cluster['Dataset']
# 	if dataset==datasets[0]:
# 		relevent_clusters=clusters[clusters['Chromosome']==chrom]
# 		relevent_clusters=relevent_clusters[relevent_clusters['Dataset']!=dataset]
# 		start_coords=relevent_clusters['ClusterStart'].values
# 		if len(start_coords)>0:
# 			displacement.append(min(start_coords,key=lambda x:abs(x-start))-start)


# clusters.apply(pair_up,axis=1)

# d={"Displacement":displacement}

# d=pandas.DataFrame(d)

# d.to_csv(directory+'distance '+datasets[0]+' upstream of '+datasets[1])

