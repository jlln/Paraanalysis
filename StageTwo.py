#In which introns become explicit, and spliced regions are known for whence they came and where they went.

import pandas
import sys

import numpy as np



np.set_printoptions(threshold='nan')
directory=sys.argv[1]


filename_index=pandas.read_csv(directory+"cluster_and_transcript_annotation_file_index.csv")

translation_status=pandas.read_csv(directory+"Translation_Status.csv")





def features(transcript):

	cluster_data=pandas.read_csv(transcript['Transcript_Clusters'])
	annotation_file=transcript['Transcript_Annotation']
	annotation_data=pandas.read_csv(annotation_file).drop('Unnamed: 0',axis=1)
	print "===========================original================================="
	print annotation_data
	print "============================================================"
	txEnd=max(annotation_data['feature_end'])
	txStart=min(annotation_data['feature_start'])
	#Make the introns explicit
	exon_df=annotation_data[annotation_data['feature_type']=="Exon"]
	exon_nucleotides=[]
	exon_starts=[]
	exon_ends=[]
	# print exon_df
	exon_df['feature_start'].map(lambda x:exon_starts.append(x))
	exon_df['feature_end'].map(lambda x:exon_ends.append(x))
	
	for e in range(len(exon_starts)):
		
		exon_nucleotides.extend(range(exon_starts[e],exon_ends[e]))
	exon_nucleotides=set(exon_nucleotides)
	transcript_nucleotides=set(range(txStart,txEnd))
	intron_nucleotides=transcript_nucleotides.difference(exon_nucleotides)
	

	
	if len(transcript_nucleotides)!=len(exon_nucleotides)+len(intron_nucleotides):
		print"#########################################"+transcript

	intron_nucleotides=list(intron_nucleotides)
	intron_starts=[]
	intron_ends=[]
	intron_features=[]
	if len(intron_nucleotides)!=0:
		intron_nucleotides=np.sort(np.array(intron_nucleotides))	#have to sort the intron nucleotides otherwise the splitting operation will do subtly weird things.
		introns=np.array_split(intron_nucleotides,np.where(np.diff(intron_nucleotides)!=1)[0]+1)
		
		
		for intron in introns:
			intron_starts.append(intron[0])
			intron_ends.append(intron[-1])
			intron_features.append('Intron')



		intron_dict={"feature_type":intron_features,"feature_start":intron_starts,"feature_end":intron_ends}
		intron_df=pandas.DataFrame(intron_dict)
		
		concat=pandas.concat([annotation_data,intron_df]).sort_index(by='feature_start')
		concat['feature_length']=concat['feature_end']-concat['feature_start']
		annotation_data=concat.copy()
		# print concat
		
	if annotation_data.iloc[0,2]=="5'-UTR":																				
		FiveUTREnd=annotation_data[annotation_data['feature_type']=="5'-UTR"]['feature_end'].values[0]
		ThreeUTR_Start=annotation_data[annotation_data['feature_type']=="3'-UTR"]['feature_start'].values[0]
		ThreeUTR_End=annotation_data[annotation_data['feature_type']=="3'-UTR"]['feature_end'].values[0]

		Five_UTR_Starters=annotation_data[annotation_data['feature_start']<FiveUTREnd]
		Five_UTR_Starters['Start_Region']="5'-UTR"

		Three_UTR_Starters=annotation_data[annotation_data['feature_start']>ThreeUTR_Start]
		Three_UTR_Starters['Start_Region']="3'-UTR"

		Coding_Region_Starters=annotation_data[(annotation_data['feature_start']>FiveUTREnd)]
		Coding_Region_Starters=Coding_Region_Starters[Coding_Region_Starters['feature_start']<ThreeUTR_Start]
		Coding_Region_Starters['Start_Region']="Coding_Region"
		

		fiveutr=annotation_data[annotation_data['feature_type']=="5'-UTR"]
		
		threeutr=annotation_data[annotation_data['feature_type']=="3'-UTR"]



		annotation_data=pandas.concat([Five_UTR_Starters,Three_UTR_Starters,Coding_Region_Starters])
		
		Five_UTR_Enders=annotation_data[annotation_data['feature_end']<FiveUTREnd]
		Five_UTR_Enders['End_Region']="5'-UTR"

		Three_UTR_Enders=annotation_data[annotation_data['feature_end']>ThreeUTR_Start]
		Three_UTR_Enders['End_Region']="3'-UTR"

		Coding_Region_Enders=annotation_data[annotation_data['feature_end']>FiveUTREnd]
		Coding_Region_Enders=Coding_Region_Enders[Coding_Region_Enders['feature_end']<ThreeUTR_Start]
		Coding_Region_Enders['End_Region']="Coding_Region"
		

		
		

		annotation_data=pandas.concat([Five_UTR_Enders,Three_UTR_Enders,Coding_Region_Enders,fiveutr,threeutr])

		annotation_data['feature_classification']=annotation_data['feature_type']+" from "+annotation_data['Start_Region']+" to "+annotation_data['End_Region']
		annotation_data[annotation_data['feature_type']=="3'-UTR"]['feature_classification']=="3'-UTR"
		annotation_data[annotation_data['feature_type']=="5'-UTR"]['feature_classification']=="5'-UTR"
	# if annotation_data.iloc[0,2]=="Untranslated Transcript":
	else:	
		annotation_data['feature_classification']=annotation_data['feature_type']+" in a ncRNA"
		


	intron_data=annotation_data[annotation_data['feature_type']=="Intron"].reset_index()
	if len(intron_data)>0:
		del intron_data['index']
		intron_data['feature_type_index']=intron_data.index+1
		
		
		intron_data_subframes=[]
		intron_subtypes=intron_data['feature_classification'].unique()

		for subtype in intron_subtypes:
			data=intron_data[intron_data['feature_classification']==subtype].reset_index()
			del data['index']
			data['feature_subtype_index']=data.index+1
			intron_data_subframes.append(data)
			
		intron_data=pandas.concat(intron_data_subframes)
			
	

	exon_data=annotation_data[annotation_data['feature_type']=="Exon"].reset_index()
	if len(exon_data)>0:
		del exon_data['index']
		exon_data['feature_type_index']=exon_data.index+1
		
		
		exon_data_subframes=[]
		exon_subtypes=exon_data['feature_classification'].unique()

		for subtype in exon_subtypes:
			data=exon_data[exon_data['feature_classification']==subtype].reset_index()
			del data['index']
			data['feature_subtype_index']=data.index+1
			exon_data_subframes.append(data)
			
		exon_data=pandas.concat(exon_data_subframes)
	
	five_utr_data=annotation_data[annotation_data['feature_type']=="5'-UTR"].reset_index()
	del five_utr_data['index']
	five_utr_data['feature_type_index']=1
	five_utr_data['feature_subtype_index']=1
	five_utr_data['feature_classification']="5'-UTR"
	

	three_utr_data=annotation_data[annotation_data['feature_type']=="3'-UTR"].reset_index()
	del three_utr_data['index']
	three_utr_data['feature_type_index']=1
	three_utr_data['feature_subtype_index']=1
	three_utr_data['feature_classification']="3'-UTR"

	annotation_data=pandas.concat([five_utr_data,three_utr_data,exon_data,intron_data]).sort_index(by='feature_start')
	print "=============================updated==============================="
	print annotation_data

	annotation_data.to_csv(annotation_file)

filename_index.apply(features,axis=1)



