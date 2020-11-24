#!/usr/bin/python
#parse prokka output to pathway binary matrix and abundance matrix

import os
import pandas as pd 

def create_mapping_dict():
	#generate mapping between prokka directory name and gene id
	mapDict={}
	files=["hmp"]
	files=[x for x in files if 'all_annotations' not in x]
	for fi in files:
		fi2=os.listdir(fi)
		fi2=[x for x in fi2 if 'ffn' in x]
		for fi3 in fi2:
			with open(fi+'/'+fi3) as f:
				for i,line in enumerate(f):
					if i==0:
						if fi!='hmp':
							mapDict[line.rstrip().split('_')[0][1:]]=fi
						if fi=='hmp':
							mapDict[line.rstrip().split('_')[0][1:]]=fi3.split('_PROKKA')[0].split('_prokka')[0].split('.denovo')[0]
					else:
						break
	return mapDict

def get_and_load_annotation_dict(ffn):
	os.system('rm all_annotations 2> ../error.log')
	annotationInfo=[]
	if ffn==False: # if we have a tsv file
		#cat annotation data together
		os.system("cat */*tsv > all_annotations ; grep 'CDS' all_annotations > foo ; mv foo all_annotations ; grep -v 'hypothetical protein' all_annotations > foo ; mv foo all_annotations")
		#open annotation info, get gene ids and annotations
		with open('all_annotations') as f:
			for line in f:
				annotationInfo.append([line.split('\t')[0].split('_')[0],line.rstrip().split('\t')[-1]])
	if ffn==True: # if we only have ffn data, get annotations from ffn file
		files=["hmp"]
		#files=[x for x in files if 'build' not in x]
		for fi in files:
			fi2=os.listdir(fi)
			fi2=[x for x in fi2 if 'ffn' in x]
			if len(fi2)==0:
				continue
			for fi3 in fi2:
				os.system("grep '>' %s/%s >> all_annotations"%(fi,fi3))
		with open('all_annotations') as f:
			for line in f:
				foo=' '.join(line.rstrip().split(' ')[1:])
				if foo!='hypothetical protein':
					annotationInfo.append([line.split(' ')[0].split('_')[0][1:],foo])
	return annotationInfo

mapDict=create_mapping_dict()
annotationInfo=get_and_load_annotation_dict(True)

#get list of all annotations
annotations=list(set([x[1] for x in annotationInfo]))

#convert to dataframe, make annotations index to quick subsetting
annotationInfo=pd.DataFrame(annotationInfo)
annotationInfo.index=annotationInfo.iloc[:,-1]

#get list of all gene ids
ind=list(set(annotationInfo.iloc[:,0]))

#iterate through each annotation, find samples that have it
output=[]
for iii,a in enumerate(annotations):
	print iii
	print a
	#create dataframe of all zeros where index is gene (sample) ids, column name is annotation
	tempDf = pd.DataFrame(0, index=ind, columns=[a])
	#subset big dataframe to find samples with that nnotation
	sub=annotationInfo.loc[a,:]
	try:
		sub.index=sub.iloc[:,0]
	except:
		#this is if the annotation only shows up once, the dataframe becomes a series and has to be dealt with differently
		sub=pd.DataFrame(sub)
		sub=sub[:1]
		sub.index=sub.iloc[0]
		sub.loc[:,a]=[1]*len(sub.index)
		tempDf.update(sub)
		output.append(tempDf)
		continue
	#create a vector of ones with the same length as the number of samples with this annotation
	sub.loc[:,a]=[1]*len(sub.index)
	#get counts of annotations with samples
	sub['temp'] = sub.groupby([sub.index])[a].transform('sum')
	sub.loc[:,a]=sub.loc[:,'temp']
	#remove duplicate rows
	sub=sub.drop_duplicates()
	#update tempDf you made above to include 1's from your vector you just made (everything else will remain zeros)
	tempDf.update(sub)
	#append this tempDf (will become row of binary matrix later on) to list for now
	output.append(tempDf)

#cat everything together
output2=pd.concat(output,axis=1)
#map index to sample names with dictionary
foo=output2.index
foo=[mapDict[x] for x in foo]
output2.index=foo

output2=output2.transpose()
#write to csv, this is now ready to be run through the fisher script.
output2.to_csv('all_ancient_ab_inito_count_data.csv',sep='\t')

output2[output2>1.0]=1.0
output2.to_csv('all_ancient_ab_inito_binary_data.csv',sep='\t')







