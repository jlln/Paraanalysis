import pandas
import numpy as np
import sys


directory=sys.argv[1]


datasets=[]
translation_status_df=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status_df.columns[1:]:
	datasets.append(x)



start_coord=[]
end_coord=[]
translation_status_list=[]
dataset=[]


def intron_list_adder(intron):
	start_coord.append(intron['feature_start'])
	end_coord.append(intron['feature_end'])
	
	dataset.append(d)
	


def get_stats(transcript):

	annotation_file=directory+transcript+'annotations.csv'

	annotation=pandas.read_csv(annotation_file)

	if annotation.iloc[0,3]=="5'-UTR":
		
		translation_status="Protein coding"
		print translation_status
		introns=annotation[annotation['feature_type']=="Intron"]
		if len(introns.index)>0:
			introns.apply(intron_list_adder,axis=1)
			translation_status_list.extend([translation_status]*len(introns.index))



	else:
		
		translation_status="Non coding"
		print translation_status
		introns=annotation[annotation['feature_type']=="Intron"]
		if len(introns.index)>0:
			introns.apply(intron_list_adder,axis=1)
			translation_status_list.extend([translation_status]*len(introns.index))
			


def get_transcripts_in_cluster(cluster):
	transcripts=cluster['transcripts']
	transcripts=eval(transcripts)
	if len(transcripts)!=0:
		for transcript in transcripts:
			get_stats(transcript)


for d in datasets:
	data=pandas.read_csv(directory+d+'matchedreads.csv')
	data.apply(get_transcripts_in_cluster,axis=1)

print len(translation_status_list),len(start_coord),len(end_coord),len(dataset)

result_dict={'feature_start':start_coord,'feature_end':end_coord,'translation_status':translation_status_list,'dataset':dataset}

result=pandas.DataFrame(result_dict)

result_g=result.groupby(['feature_start','feature_end','translation_status','dataset'])
result_g=result_g.size().reset_index()
result_g.columns=['start_coord','end_coord','translation_status','dataset','count']
result_g['length']=result_g['end_coord']-result_g['start_coord']

print result_g.head()
result_g.to_csv(directory+'overall_intron_lengths.csv')



