# "Where features and clusters come together"


import pandas
import sys
import numpy as np

pandas.set_option('line_width', 160)


directory=sys.argv[1]
datasets=[]

filename_index=pandas.read_csv(directory+"cluster_and_transcript_annotation_file_index.csv")
translation_status=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status.columns[1:]:
	datasets.append(x)



list_of_features=[]
def describe_features_within_cluster(cluster):

	cluster_start=cluster['ClusterStart']
	cluster_end=cluster['ClusterEnd']
	matched_features_list=[]
	for transcript in eval(cluster['transcripts']):
		annotation=pandas.read_csv(directory+transcript+'annotations.csv')
		
		matched_features=annotation[annotation['feature_start']<cluster_end]
		matched_features=matched_features[matched_features['feature_end']>cluster_start]
		matched_features['transcript']=transcript
		matched_features_list.append(matched_features)
	if len(matched_features_list)>0:
		cluster_features=pandas.concat(matched_features_list)
		del cluster_features['Unnamed: 0']
		print "========================================================"
		print cluster['ClusterID']
		cluster_features.rename(columns=lambda x: x.replace('feature_classification', 'feature_subtype'), inplace=True)
		groups=cluster_features.groupby(['feature_start','feature_end','feature_subtype'])
		groups=groups.size().reset_index()
		groups.columns=["Feature Start","Feature End","Feature Subtype", "Feature Occurance Count within Cluster"]
		print groups
		list_of_features.extend(groups['Feature Subtype'])
	else:
		print "****** Fell through no clusters from this dataset in this transcript "
		print cluster['transcripts']



feature_dfs=[]
for d in datasets:
	list_of_features=[]
	matched_reads=pandas.read_csv(directory+d+"matchedreads.csv")
	matched_reads.apply(describe_features_within_cluster,axis=1)
	feature_subtypes=set(list_of_features)
	print feature_subtypes
	df_type_list=[]
	df_type_count=[]
	for f in feature_subtypes:
		df_type_list.append(f)
		df_type_count.append(list_of_features.count(f))
	df_d={"Feature_Type":df_type_list,"%s Count"%d:df_type_count}
	totals=pandas.DataFrame(df_d)
	totals.to_csv(directory+d+"feature_cluster_table.csv")
	feature_dfs.append(totals)
	


t=pandas.merge(left=feature_dfs[0],right=feature_dfs[1],on="Feature_Type",how="outer")
t=t.fillna(0)
tt=t.T
tt.columns=tt.iloc[0,]
df = tt.drop(tt.index[0])
df=df.div(df.sum(axis=1), axis=0)*100
df.to_csv(directory+'feature_binding_proportions.csv')
	
