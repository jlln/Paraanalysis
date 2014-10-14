import sys
import pandas
import numpy as np 
import scipy
import numpy.random as npr

directory=sys.argv[1]

datasets=[]
translation_status=pandas.read_csv(directory+"_Translation_Status_Totals.csv")
for x in translation_status.columns[1:]:
	datasets.append(x)



matched_transcripts=[]


def accumulate_transcripts(cluster,list_of_transcripts):
	list_of_transcripts.extend(eval(cluster['transcripts']))


name=pandas.read_csv('UCSC_genenames.csv')

for d in datasets:
	result_list=[]
	dataset_list=[]
	matched_transcripts=[]
	matched_reads=pandas.read_csv(directory+d+'matchedreads.csv')
	matched_reads.apply(accumulate_transcripts,axis=1,args=([matched_transcripts]))
	matched_transcripts=set(matched_transcripts)
	matched_transcripts=list(matched_transcripts)
	result_list.extend(matched_transcripts)
	dataset_list.extend([d]*len(matched_transcripts))
	print len(result_list),len(dataset_list)


	ucsc_transcripts={"transcript":result_list,"dataset":dataset_list}
	ucsc_transcripts=pandas.DataFrame(ucsc_transcripts)
	print ucsc_transcripts.head()


	print name.head()


	merged=pandas.merge(left=ucsc_transcripts,right=name,left_on="transcript",right_on="#kgID",how="left")


	print merged.head()
	merged.to_csv(d[:-4]+"for_ingenuity.csv")





