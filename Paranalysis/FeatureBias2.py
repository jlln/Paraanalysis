import pandas
import sys
import numpy as np


#Come back to this once transcripts are properly identified.


pandas.set_option('line_width', 160)


directory=sys.argv[1]
datasets=[]

filename_index=pandas.read_csv(directory+"cluster_and_transcript_annotation_file_index.csv")
translation_status=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status.columns[1:]:
	datasets.append(x)


def describe_features(transcript,cluster_start,cluster_end,clusterID,matched_features):
	annotation=pandas.read_csv(directory+transcript+'annotations.csv')
	
	up_match=annotation[annotation['feature_start']<cluster_end]
	down_match=up_match[up_match['feature_end']>cluster_start]
	cds_start=annotation[annotation['feature_type']=="5'-UTR"]['feature_end']
	cds_end=annotation[annotation['feature_type']=="3'-UTR"]['feature_start']
	
	if "5'-UTR" in down_match['feature_type'].values:
		down_match=down_match[down_match['feature_type']=="5'-UTR"]
		down_match['translation_status']="mRNA"
	elif "3'-UTR" in down_match['feature_type'].values:
		down_match=down_match[down_match['feature_type']=="3'-UTR"]
		down_match['translation_status']="mRNA"
	#If the cluster didn't match with a UTR alone, then it must be outside the UTRs, even if it matches an exon spanning into the UTR.

	elif down_match['feature_classification'].values[0]=="Intron in a ncRNA":
		down_match['translation_status']="ncRNA"
		if down_match['feature_type_index'].values[0]==1:
			down_match['feature_type']="First ncRNA Intron"
		else:
			down_match['feature_type']="Non-First ncRNA Intron"

	elif down_match['feature_type'].values[0]=="Exon":
		down_match['translation_status']="mRNA"
		down_match['feature_type']="mRNA Exon"
	elif down_match['feature_type'].values[0]=="Intron":
		down_match['translation_status']="mRNA"
		if down_match['feature_type_index'].values[0]==1:
			down_match['feature_type']="First mRNA Intron"
		else:
			down_match['feature_type']="Non-First mRNA Intron"
		
	

	

	matched_features.append(down_match)


overall_counts=[]


def describe_transcripts(cluster):
	clusterID=cluster['ClusterID']
	list_transcripts=[]
	cluster_start=cluster['ClusterStart']
	cluster_end=cluster['ClusterEnd']
	transcripts=eval(cluster['transcripts'])
	matched_features=[]
	for t in transcripts:
		describe_features(t,cluster_start,cluster_end,clusterID,matched_features)
	if len(matched_features)!=0:
		matched_features=pandas.concat(matched_features)
		uniques=matched_features.groupby(['feature_start','feature_end','feature_type','translation_status'])
		collapsed=uniques.size().reset_index()
		overall_counts.append(collapsed)

	# else:
		# print "***** No matched transcripts in this cluster *****"


normal=lambda x: (x/x.sum()) #To give relative counts.

result_list=[]
for d in datasets:
	list_transcripts=[]
	data=pandas.read_csv(directory+d+'matchedreads.csv')
	data.apply(describe_transcripts,axis=1)
	overall_df=pandas.concat(overall_counts)
	print overall_df
	groups=overall_df.groupby(['feature_type','translation_status'])
	result=groups.size().reset_index()
	result['dataset']=d
	result.columns=['feature_type','translation_status','Count','Dataset']
	groups=result.groupby(['feature_type','dataset','translation_status'])


overall_df=pandas.concat(result_list)
overall_df.to_csv(directory+)
	


