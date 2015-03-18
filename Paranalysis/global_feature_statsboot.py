import sys
import pandas
import numpy as np 
import scipy
import numpy.random as npr
import itertools



def bootstrap_median_confidence_95_interval(data):
	medians=[]
	for i in range(1,50000):
		sample=npr.choice(data,size=len(data),replace=True)
		medians.append(int(np.median(sample)))
		
	medians=np.sort(medians)
	upper_ci=medians[46250]
	lower_ci=medians[1250]
	median=np.median(data)
	return upper_ci,lower_ci,median



ucsc_db=pandas.read_table('UCSC_transcripts.txt')


feature_types=[]
feature_starts=[]
feature_ends=[]
chromosomes=[]
translation_statuses=[]




def accumulate_features(transcript):
	print "ParsingTranscript"
	chromosome=transcript['chrom']
	

	transcript_name=transcript['#name']
	
	txStart=transcript['txStart']
	txEnd=transcript['txEnd']
	
	cdsStart=transcript['cdsStart']
	cdsEnd=transcript['cdsEnd']
	Translation_Status=""

	
	if cdsStart!=cdsEnd:
	    Translation_Status="mRNA"
	    translation_statuses.append(Translation_Status)
	    feature_starts.append(txStart)
	    feature_ends.append(cdsStart)
	    feature_types.append("5'-UTR")
	    chromosomes.append(chromosome)
		
	    feature_starts.append(cdsEnd)
	    feature_ends.append(txEnd)
	    feature_types.append("3'-UTR")
	    chromosomes.append(chromosome)
	    translation_statuses.append(Translation_Status)
	else:
	    
	    Translation_Status="ncRNA"
	    

	
	n_exons=transcript['exonCount']

	exon_starts=(transcript['exonStarts']).split(",")[:-1]
	exon_starts[:]=[int(x) for x in exon_starts]
	  

	exon_ends=(transcript['exonEnds']).split(",")[:-1]
	exon_ends[:]=[int(x) for x in exon_ends]
	  

	feature_starts.extend(exon_starts)
	feature_ends.extend(exon_ends)
	feature_types.extend(['Exon']*n_exons)
	translation_statuses.extend([Translation_Status]*n_exons)
	chromosomes.extend([chromosome]*n_exons)

	
	
	'''
	A new, faster way to create explicit introns
	'''
	
	exon_boundaries=zip(exon_starts,exon_ends)
	exon_boundaries=sorted(exon_boundaries)
	intron_count=0
	for x in range(len(exon_boundaries)-1):
		us=exon_boundaries[x]
		ds=exon_boundaries[x+1]
		distance=ds[0]-us[1]
		if distance>0:
			if intron_count==0:
				feature_starts.append(us[1]+1)
				feature_ends.append(ds[0]-1)
				feature_types.append("First Intron")
				chromosomes.append(chromosome)
				translation_statuses.append(Translation_Status)
				intron_count+=1
			else:
				feature_starts.append(us[1]+1)
				feature_ends.append(ds[0]-1)
				feature_types.append("Non-First Intron")
				chromosomes.append(chromosome)
				translation_statuses.append(Translation_Status)
				intron_count+=1
	print "Transcript Done"

		


ucsc_db.apply(accumulate_features,axis=1)


global_feature_dictionary={"feature_type":feature_types,"feature_start":feature_starts,"feature_end":feature_ends,"chromosome":chromosomes,"translation_status":translation_statuses}

global_feature_df=pandas.DataFrame(global_feature_dictionary)

global_feature_grouped=global_feature_df.groupby(['feature_type','feature_start','feature_end','translation_status','chromosome'])
global_feature_df=global_feature_grouped.size().reset_index()
del global_feature_df[global_feature_df.columns[4]]	
global_feature_df['length']=global_feature_df['feature_end']-global_feature_df['feature_start']
global_feature_grouped=global_feature_df.groupby(['feature_type','translation_status'])



feature_list=[]
translation_status_list=[]
lower_ci_list=[]
median_list=[]
upper_ci_list=[]

print "Bootstrapping"
for name, group in global_feature_grouped:
	feature_list.append(name[0])
	translation_status_list.append(name[1])
	lengths=group['length']
	bootstrap_result=bootstrap_median_confidence_95_interval(lengths)
	upper_ci_list.append(bootstrap_result[0])
	lower_ci_list.append(bootstrap_result[1])
	median_list.append(bootstrap_result[2])
result_dict={"Feature_Type":feature_list,"Translation_Status":translation_status_list,"UpperCI":upper_ci_list,"LowerCI":lower_ci_list,"Median":median_list}
result=pandas.DataFrame(result_dict)
result['dataset']="UCSC Transcript Database"
result.to_csv('global_transcriptome_feature_bootstrap.csv')
print result



