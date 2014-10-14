#Finds the distances from each cluster to various features.


import pandas
import sys
import operator
import numpy as np

ucsc_db=pandas.read_table('UCSC_transcripts.txt')

directory=sys.argv[1]

datasets=[]
translation_status_df=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status_df.columns[1:]:
	datasets.append(x)


def get_distances(feature,cluster_start,cluster_end,clusterID,transcript_length):
	
	feature_start=feature['feature_start']
	feature_end=feature['feature_end']
	feature_classification=feature['feature_classification']
	feature_type=feature['feature_type']
	feature_type_index=feature['feature_type_index']
	feature_subtype_index=feature['feature_subtype_index']
	
	abs_distance=min([abs(cluster_start-feature_end),abs(cluster_end-feature_start),abs(cluster_start-feature_start),abs(cluster_end-feature_end)])
	if cluster_start<feature_end:
		if cluster_end>feature_start:
			abs_distance=0
	if feature_start<cluster_end:
		if feature_end>cluster_start:
			abs_distance=0
	if cluster_end>feature_end:
		if cluster_start<feature_start:
			abs_distance=0
	if cluster_end<feature_end:
		if cluster_start>feature_start:
			abs_distance=0
	rel_distance=float(abs_distance)/transcript_length
	# print abs_distance
	df_dataset.append(d)
	df_cluster.append(clusterID)
	df_feature_type.append(feature_type)
	df_feature_subtype.append(feature_classification)
	df_abs_distance.append(abs_distance)
	df_rel_distance.append(rel_distance)
	df_feature_type_index.append(feature_type_index)
	df_feature_subtype_index.append(feature_subtype_index)
	# print feature_classification
	print feature_subtype_index





def get_transcripts_and_distances(cluster):
	clusterID=cluster['ClusterID']
	cluster_start=cluster['ClusterStart']
	cluster_end=cluster['ClusterEnd']
	transcripts=eval(cluster['transcripts'])
	for transcript in transcripts:
		transcript_annotation=pandas.read_csv(directory+transcript+'annotations.csv')
		transcript_length=max(transcript_annotation['feature_end'])-min(transcript_annotation['feature_start'])
		# print transcript
		transcript_annotation=transcript_annotation[transcript_annotation['feature_type']=="Splice Junction "] #Select the type of feature to look at
		if len(transcript_annotation.index)>0:

			transcript_annotation.apply(get_distances,axis=1,args=(cluster_start,cluster_end,clusterID,transcript_length))



df_dataset=[]
df_cluster=[]
df_feature_type=[]
df_feature_subtype=[]
df_feature_type_index=[]
df_feature_subtype_index=[]
df_rel_distance=[]
df_abs_distance=[]

for d in datasets:

	cluster_data=pandas.read_csv(directory+d+'matchedreads.csv')
	cluster_data.apply(get_transcripts_and_distances,axis=1)


dict_df={"Dataset":df_dataset,"Cluster":df_cluster,"FeatureType":df_feature_type,"FeatureSubtype":df_feature_subtype,
"FeatureTypeIndex":df_feature_type_index,"FeatureSubtypeIndex":df_feature_subtype_index,"Absolute Distance":df_abs_distance,"Relative Distance":df_rel_distance}
results=pandas.DataFrame(dict_df)



results=results.sort(['Cluster','Absolute Distance'])


closest_features=[]
clusters=list(set(results['Cluster'].values))
for cluster in clusters:
	closest_features.append((results[results['Cluster']==cluster]).head(1))

closest=pandas.concat(closest_features)

def insert_newlines(string, every=30):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)


def fix_labels(row):
	label=row['feature_classification']
	label=insert_newlines(label,every=30)
	row['label']=label


closest.to_csv(directory+'closest_splicing_junctions.csv')



