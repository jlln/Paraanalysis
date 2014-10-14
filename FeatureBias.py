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


matches_list=[]
def describe_matched_feature(transcripts,cluster_start,cluster_end,cluster_ID,dataset):
	for t in transcripts:
		annotations=pandas.read_csv(directory+t+'annotations.csv')
		if len(annotations.index)>0:
			if "5'-UTR" in annotations['feature_classification'].values:
				translation_status="mRNA"
			else:
				translation_status="ncRNA"
			
			matched_features=annotations[(annotations['feature_start']<cluster_end) & (annotations['feature_end']>cluster_start)]	

			if translation_status=="mRNA":
				if "5'-UTR" in matched_features['feature_classification'].values:							#If a cluster matched with a UTR alone(ie it is outside of the CDS, then exon/intron status of that location is ignored)
					matched_feature=matched_features[matched_features['feature_classification']=="5'-UTR"]
				elif "3'-UTR" in matched_features['feature_classification'].values:
					matched_feature=matched_features[matched_features['feature_classification']=="3'-UTR"]
				else:
					matched_feature=matched_features
				
					if matched_features['feature_type'].values[0]=="Intron":
						if matched_features['feature_type_index'].values[0]==1:
							matched_feature['feature_type']="First Intron"
						elif matched_features['feature_type_index'].values[0]>1:
							matched_feature['feature_type']="Non-First-Intron"





			if translation_status=="ncRNA":
				matched_feature=matched_features
				if matched_features['feature_type'].values[0]=="Intron":
						if matched_features['feature_type_index'].values[0]==1:
							matched_feature['feature_type']="First Intron"
						elif matched_features['feature_type_index'].values[0]>1:
							matched_feature['feature_type']="Non-First-Intron"

		matched_feature['clusterID']=cluster_ID
		matched_feature['Dataset']=dataset
		matched_feature['Translation_Status']=translation_status
		matches_list.append(matched_feature)
					
	



def retrieve_transcripts(cluster,dataset):
	transcripts=eval(cluster['transcripts'])
	cluster_ID=cluster['ClusterID']
	cluster_start=cluster['ClusterStart']
	cluster_end=cluster['ClusterEnd']
	describe_matched_feature(transcripts,cluster_start,cluster_end,cluster_ID,dataset)


normal=lambda x: (x/x.sum())


for d in datasets:
	clusters=pandas.read_csv(directory+d+'matchedreads.csv')
	clusters.apply(retrieve_transcripts,axis=1,args=[d])

result=pandas.concat(matches_list)
grouped=result.groupby(['feature_type','Dataset','Translation_Status'])
result=grouped.size().reset_index()
result.columns=['feature_type','Dataset','Translation_Status','Count']
groups=result.groupby(['Dataset','Translation_Status'])


normalized=groups.transform(normal)




result['Relative Count']=normalized

print result

result.to_csv(directory+'FeatureBindingProportions.csv')




