import pandas
import sys
import operator

ucsc_db=pandas.read_table('UCSC_transcripts.txt')

directory=sys.argv[1]

datasets=[]
translation_status_df=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status_df.columns[1:]:
	datasets.append(x)


def find_closest_junction(junctions_df,cluster_start,cluster_end):
	distance_list=[]
	distance_entries=[]



def get_transcript_data(transcript,cluster_start,cluster_end):
	splice_junctions=transcript[transcript['feature_type']=="Splice Junction "]
	if len(splice_junctions.index)>0:
		find_closest_junction(splice_junctions,cluster_start,cluster_end)
		


def analyse_cluster(cluster):
	transcripts=eval(cluster['transcripts'])
	cluster_start=cluster['ClusterStart']
	cluster_end=cluster['ClusterEnd']
	for transcript in transcripts:
		transcript_data=pandas.read_csv(directory+transcript+'annotations.csv')





for d in datasets:
	cluster_data=pandas.read_csv(directory+d+'matchedreads.csv')
	