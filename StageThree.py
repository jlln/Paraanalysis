# Inserting boundary features into the annotations.

import pandas
import sys
import numpy as np




directory=sys.argv[1]


filename_index=pandas.read_csv(directory+"cluster_and_transcript_annotation_file_index.csv")

def boundaries(transcript):
	
	annotation_file=transcript['Transcript_Annotation']
	annotation_data=pandas.read_csv(annotation_file).drop('Unnamed: 0',axis=1)
	annotation_data=annotation_data.sort_index(by='feature_start')
	
	if annotation_data.iloc[0,5]=="5'-UTR":
		#mRNAs
		introns_exons=annotation_data[annotation_data['feature_type']!="5'-UTR"]
		introns_exons=introns_exons[introns_exons['feature_type']!="3'-UTR"]
		n_junctions=len(introns_exons.index)-1
		junction_start=[]
		junction_end=[]
		junction_feature_type=[]
		junction_classification=[]
		junction_type_index=[]
		j_type_index=1
		for x in range(n_junctions):
			upside=introns_exons.iloc[x,]
			downside=introns_exons.iloc[x+1,]
			junction_start.append(upside['feature_end'])
			junction_end.append(downside['feature_start'])
			
			junction_feature_type.append("Splice Junction ")
			junction_classification.append("Splice Junction "+upside['feature_classification']+" to "+downside['feature_classification'])
			junction_type_index.append(j_type_index)
			j_type_index+=1
		
		fiveutr=annotation_data[annotation_data['feature_type']=="5'-UTR"]
		threeutr=annotation_data[annotation_data['feature_type']=="3'-UTR"]

		junction_start.append(fiveutr['feature_end'][0])
		junction_end.append(fiveutr['feature_end'][0])
		junction_feature_type.append("Translation_Start")
		junction_classification.append("Translation_Start")
		junction_type_index.append(1)

		junction_start.append(threeutr['feature_start'].values[0])
		junction_end.append(threeutr['feature_start'].values[0])

		junction_feature_type.append("Translation_End")
		junction_classification.append("Translation_End")
		junction_type_index.append(1)
		


		d={"feature_type":junction_feature_type,"feature_start":junction_start,"feature_end":junction_end,"feature_classification":junction_classification,"feature_type_index":junction_type_index}
		junctions_df=pandas.DataFrame(d)



		junction_subtypes=junctions_df['feature_classification'].unique()
		print junction_subtypes
		junctions_df_subframes=[]

		
		for subtype in junction_subtypes:
			data=junctions_df[junctions_df['feature_classification']==subtype].reset_index()
			del data['index']
			data['feature_subtype_index']=data.index+1
			junctions_df_subframes.append(data)
			
		junctions_df=pandas.concat(junctions_df_subframes)
		

		annotation_data=pandas.concat([junctions_df,annotation_data]).sort(['feature_start'])
		

	else:
		print "NoncodingRNA"
		n_junctions=len(annotation_data.index)-1
		junction_start=[]
		junction_end=[]
		junction_feature_type=[]
		junction_classification=[]
		junction_type_index=[]
		j_type_index=1
		for x in range(n_junctions):
			upside=annotation_data.iloc[x,]
			downside=annotation_data.iloc[x+1,]
			junction_start.append(upside['feature_end'])
			junction_end.append(downside['feature_start'])
			
			junction_feature_type.append("Splice Junction ")
			junction_classification.append("Splice Junction "+upside['feature_classification']+" to "+downside['feature_classification'])
			junction_type_index.append(j_type_index)
			j_type_index+=1
		d={"feature_type":junction_feature_type,"feature_start":junction_start,"feature_end":junction_end,"feature_classification":junction_classification,"feature_type_index":junction_type_index}
		junctions_df=pandas.DataFrame(d)




		annotation_data=pandas.concat([junctions_df,annotation_data]).sort(['feature_start','feature_end'])

	annotation_data.to_csv(annotation_file)

filename_index.apply(boundaries,axis=1)


	
	




	


