##This script accepts Paralyzer output files (to be listed as arguments on the command line), and requires a
##copy of UCSC_transcripts.txt to be in the directory.







import pandas
import numpy as np
import os


import sys

ucsc_db=pandas.read_csv("UCSC_transcripts.txt",delim_whitespace=True) # From the UCSC Database; a list of transcripts.
ucsc_names=pandas.read_csv("UCSC_genenames.csv")                      #Extra data
input_filenames=sys.argv[1:]
datasets=[]
dataset_names=[]
for f in input_filenames:
    dataset_names.append(f[:16])


n_transcripts=len(ucsc_db)



def match_check(transcript):

    if transcript['txStart']<clusterEnd:
        if transcript['txEnd']>clusterStart:
            transcripts.append(transcript['#name'])


clusters_dict={}
clusters_by_transcripts_dict={}
matched_transcripts_dict={}

#Establish Output Directory
transcript_output_directory_name='_'.join(dataset_names)+"Transcripts"
transcript_output_directory="./"+transcript_output_directory_name+"/"
if not os.path.exists(transcript_output_directory):
    os.makedirs(transcript_output_directory)

#First identify which transcripts contain each cluster.

for dataset in input_filenames:
    clusters=pandas.read_csv(dataset)
    dataset=dataset[:16]
    datasets.append(dataset)
    clusters_dict[dataset]=clusters
    if os.path.isfile(transcript_output_directory+"%smatchedreads.csv"%(dataset)) is True:               #checking if the dataset has already been processed,because this step is relatively slow
        print "Using previously matched reads for ",dataset
        clusters=pandas.read_csv(transcript_output_directory+"%smatchedreads.csv"%(dataset))
        transcripts_list=clusters['transcripts'].tolist()
        matched_transcripts=[]
        for transcript in transcripts_list:
            for g in transcript[1:-1].split(", "):
                if len(g)>2:
                    matched_transcripts.append(g[1:-1])
        matched_transcripts=frozenset(matched_transcripts)


    else:
        transcripts_list=[]
        print "Matching reads for ",dataset

        for cluster in clusters.iterrows():
            transcripts=[]
            c=cluster[1]
            cluster_transcripts=[]
            chromosome=c['Chromosome']
            clusterStart=c['ClusterStart']
            clusterEnd=c['ClusterEnd']
            transcripts_on_read_chromosome=ucsc_db[ucsc_db['chrom']==chromosome]
            transcripts_on_read_chromosome.apply(match_check,axis=1)
            transcripts_list.append(transcripts)


        clusters['transcripts']=pandas.Series(transcripts_list)

        clusters.to_csv(transcript_output_directory+"%smatchedreads.csv"%(dataset))
        matched_transcripts=[]
        for l in transcripts_list:
            for transcript in l:
                matched_transcripts.append(transcript)
        matched_transcripts=frozenset(matched_transcripts)








#Second, do the reverse, ie identify the clusters that are contained within each transcript, and merge with the ucsc table. Store the resulting dataframe in a dicionary keyed to the dataset.
    clusters_lists=[]

    for transcript in matched_transcripts:
        f=lambda x:transcript in x
        mask=clusters['transcripts'].apply(f)
        clusters_lists.append((clusters[mask]['ClusterID']).values)

    matched_transcripts=list(matched_transcripts)

    index_matched_transcripts=range(len(matched_transcripts))
    read_counts=clusters['ReadCount']



    d={"transcript":matched_transcripts,
"reads":clusters_lists}
    matched_df=pandas.DataFrame(d,index=index_matched_transcripts)




    ucsc_transcript_merge=pandas.merge(matched_df,ucsc_db,left_on='transcript',right_on="#name",how="left")


    label=[dataset]*len(ucsc_transcript_merge)
    ucsc_transcript_merge['dataset']=label


    clusters_by_transcripts_dict[dataset]=ucsc_transcript_merge



#Thirdly, make a table, which, for each transcript, lists the matched clusters from each dataset.




dataset_dfs_to_merge=[]
column_names_in_merged_df=['transcript']
transcripts_with_clusters=clusters_by_transcripts_dict[datasets[0]]
colnames=transcripts_with_clusters.columns
colnames=list(colnames)
colnames[0]=datasets[0]
transcripts_with_clusters.columns=colnames
transcripts_with_clusters=transcripts_with_clusters[['transcript',datasets[0]]]

for dataset in datasets[1:]:
    dataset_clusters=clusters_by_transcripts_dict[dataset][['transcript','reads']]
    dataset_clusters.columns=['transcript',dataset]
    dataset_dfs_to_merge.append(dataset_clusters)
    transcripts_with_clusters=pandas.merge(left=transcripts_with_clusters,right=dataset_clusters,on='transcript',how='outer')

transcripts_with_clusters_df=pandas.merge(left=transcripts_with_clusters,right=ucsc_db,left_on='transcript',right_on='#name',how='left')
named_transcripts_with_clusters_df=pandas.merge(transcripts_with_clusters_df,ucsc_names,left_on="#name",right_on="#kgID",how="left")


print named_transcripts_with_clusters_df.head()
named_transcripts_with_clusters_df.to_csv(transcript_output_directory+'named_transcripts_with_cluster.csv')





#Step Four: Produce a dataframe for each transcript containing clusters. 5 columns:clusterID,conversioneventcount,dataset,cluster_start,cluster_end. Only the nucleotides within clusters are described. Nucleotides are numbered by distance from transcription start.

def retrieve_clusters(transcript):


    transcript_name=transcript['#name']
    print transcript_name
    df_clusterID=[]
    df_dataset=[]
    df_conversioneventcount=[]
    df_cluster_start=[]
    df_cluster_end=[]

    start_coord=transcript['txStart']
    end_coord=transcript['txEnd']
    for dataset in datasets:
        print dataset
        clusters=transcript[dataset]                #retrieve the clusters of that dataset(if any) matched to the transcript.
        clusters_dataframe=clusters_dict[dataset]   #The original clusters dataset, from which the additional cluster data will be retrieved.

        if type(clusters) is not float:             #Pandas considers empty cells to be floats
            for clusterID in clusters:
                cluster=clusters_dataframe[clusters_dataframe['ClusterID']==clusterID]

                conversion_event_count=cluster['ConversionEventCount'].values[0]
                cluster_start=cluster['ClusterStart'].values[0]
                cluster_end=cluster['ClusterEnd'].values[0]
                df_clusterID.append(clusterID)
                df_dataset.append(dataset)
                df_conversioneventcount.append(conversion_event_count)
                df_cluster_start.append(cluster_start-start_coord)
                df_cluster_end.append(cluster_end-start_coord)






    df_dict={'Dataset':df_dataset,
    "ClusterID":df_clusterID,
    "ConversionEventCount":df_conversioneventcount,
    "ClusterStart":df_cluster_start,
    "ClusterEnd":df_cluster_end}

    df=pandas.DataFrame(df_dict)


    output_path=transcript_output_directory+transcript_name+'.csv'

    df.to_csv(output_path)

    transcript_cluster_files.append(transcript_output_directory+transcript_name+'.csv')








transcript_cluster_files=[]
transcripts_with_clusters_df.apply(retrieve_clusters,axis=1)

transcripts_reads_filenames={"File":transcript_cluster_files}
transcripts_reads_filenames=pandas.DataFrame(transcripts_reads_filenames)

transcripts_reads_filenames.to_csv(transcript_output_directory+'transcript_list.csv')



#Step Five: Create initial annotations for each of the transcripts previously identified. Each transcript gets its own annotation csv file, naming is 
                                    #as above, but suffixed with 'annotation'.
def create_annotations(transcript):
    print transcript
    feature_start=[]
    feature_end=[]
    feature_type=[]

    transcript_name=transcript['#name']
    list_of_transcripts.append(transcript_name)
    txStart=transcript['txStart']
    txEnd=transcript['txEnd']
    feature_start.append(txStart)
    cdsStart=transcript['cdsStart']
    cdsEnd=transcript['cdsEnd']

    print transcript_name,txStart,txEnd,cdsStart,cdsEnd
    if cdsStart==cdsEnd:
        feature_end.append(txEnd)
        feature_type.append("Untranslated Transcript")
        translated_untranslated_transcripts.append(transcript_name)
        translated_untranslated_status.append("Untranslated")
        print "Untranslated"
    else:
        feature_end.append(cdsStart)
        feature_type.append("5'-UTR")
        feature_start.append(cdsEnd)
        feature_end.append(txEnd)
        feature_type.append("3'-UTR")
        translated_untranslated_transcripts.append(transcript_name)
        translated_untranslated_status.append("Translated")
        print "Translated"
    n_exons=transcript['exonCount']

    exon_starts=(transcript['exonStarts']).split(",")[:-1]
    exon_starts[:]=[int(x) for x in exon_starts]


    exon_ends=(transcript['exonEnds']).split(",")[:-1]
    exon_ends[:]=[int(x) for x in exon_ends]


    feature_start.extend(exon_starts)
    feature_end.extend(exon_ends)
    feature_type.extend(['Exon']*n_exons)


    d={'feature_start':feature_start,
    'feature_end':feature_end,
    'feature_type':feature_type}
    df=pandas.DataFrame(d)
    print df.head()

    print transcript_name
    outputname=transcript_output_directory+transcript_name+'annotations'+'.csv'
    df.to_csv(outputname)
    annotation_filenames.append(outputname)





annotation_filenames=[]
list_of_transcripts=[]  #for the index in step six.
#Two lists to keep track of transcript translation status: one list of transcript names and one list of translation statuses.
translated_untranslated_transcripts=[]
translated_untranslated_status=[]


transcripts_with_clusters_df.apply(create_annotations,axis=1)



# Step Six: Create csv file indexes to assist downstream scripts.
cluster_and_transcript_annotation_file_index={"Transcript_Clusters":transcript_cluster_files,
                                                "Transcript_Annotation":annotation_filenames}



cluster_and_transcript_annotation_file_index=pandas.DataFrame(cluster_and_transcript_annotation_file_index)
cluster_and_transcript_annotation_file_index.to_csv(transcript_output_directory+"cluster_and_transcript_annotation_file_index.csv")
#Step Seven: ClusterCounts for untranslated and translated transcripts
cluster_counts={}
for d in datasets:
    cluster_counts[d]=[]

for transcript in translated_untranslated_transcripts:
    for d in datasets:
        if type(transcripts_with_clusters_df[transcripts_with_clusters_df["#name"]==transcript][d].values[0]) is float:    #empty cell test; pandas considers empty cells as floats.
            cluster_counts[d].append(0)
        else:

            cluster_counts[d].append(len(transcripts_with_clusters_df[transcripts_with_clusters_df["#name"]==transcript][d].values[0]))


translation={"Transcript":translated_untranslated_transcripts,"Translation_Status":translated_untranslated_status}

translation=pandas.DataFrame(translation)




for d in datasets:
    translation[d]=pandas.Series(cluster_counts[d],index=translation.index)




translation.to_csv(transcript_output_directory+"Translation_Status.csv")

grouped_by_translation_status=translation.groupby('Translation_Status')

translation_status_summary=grouped_by_translation_status.aggregate(np.sum).reindex()
print translation_status_summary

translation_status_summary.to_csv(transcript_output_directory+"_Translation_Status_Totals.csv")



print "Output directory (provide as argument to StageTwo.py):"
print transcript_output_directory
