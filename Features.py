import pandas
import sys
import time
from operator import itemgetter
from itertools import groupby

directory=sys.argv[1]


filename_index=pandas.read_csv(directory+"cluster_and_transcript_annotation_file_index.csv")

translation_status=pandas.read_csv(directory+"Translation_Status.csv")



	




def translated(annotation_file,transcript_file):
	print annotation_file
	annotation_data=pandas.read_csv(annotation_file).drop('Unnamed: 0',axis=1)
	



	UTR5_end=annotation_data[annotation_data['feature_type']=="5'-UTR"]['feature_end'].values[0]
	UTR5=set(range(0,UTR5_end))
	UTR3_end=annotation_data[annotation_data['feature_type']=="3'-UTR"]['feature_end'].values[0]
	UTR3_start=annotation_data[annotation_data['feature_type']=="3'-UTR"]['feature_start'].values[0]
	
	exon_df=annotation_data[annotation_data['feature_type']=="Exon"]
	exon_nucleotides=[]
	exon_starts=[]
	exon_ends=[]
	
	exon_df['feature_start'].map(lambda x:exon_starts.append(x))
	exon_df['feature_end'].map(lambda x:exon_ends.append(x))
	for n in range(0,len(exon_starts)-1):
		exon_nucleotides.extend(range(exon_starts[n],exon_ends[n]))
	exon_nucleotides=set(exon_nucleotides)

		

	
		
		
	transcript=set(range(0,UTR3_end))
	intron_nucleotides=transcript.difference(exon_nucleotides)

	introns=[]
	for key,group in groupby(enumerate(intron_nucleotides),lambda (index,item):index-item):	
		group=map(itemgetter(1),group)
		introns.append([group[0],group[-1]])
	intron_feature=[]
	intron_start=[]
	intron_end=[]
	for intron in introns:
		intron_start.append(intron[0])
		intron_end.append(intron[1])
		intron_feature.append('Intron')
		
	intron_df_dict={"feature_type":intron_feature,"feature_start":intron_start,"feature_end":intron_end}
	intron_df=pandas.DataFrame(intron_df_dict)
	
	concat=pandas.concat([annotation_data,intron_df]).sort('feature_start',axis=0)
	print concat

	concat.to_csv(annotation_file)
    
	

	



def untranslated(annotation_file,transcript_file):	##Finish Later
	return

def translation_sorter(entry):
	
	transcript=entry['Transcript']
	
	status=entry['Translation_Status']
	# print status
	transcript_file=directory+transcript+".csv"
	# print transcript_file
	annotation_file=directory+transcript+"annotations.csv"
	if status=="Translated":
		
		translated(annotation_file,transcript_file)
	elif status=="Untranslated":
		untranslated(annotation_file,transcript_file)






translation_status.apply(translation_sorter,axis=1)



