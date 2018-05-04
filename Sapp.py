import os
from flask import Flask, render_template, request
import pandas as pd
from numpy import mean, sqrt, square
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, recall_score, precision_score,f1_score, fbeta_score
import csv

__author__ = "SitaSC"

app = Flask(__name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))

@app.route("/")
def index():
	return render_template("uploadButton.html")

@app.route("/upload", methods = ['POST'])
def upload():
	target = os.path.join(APP_ROOT, "files/")
	print(target)


	if not os.path.isdir(target):
		os.mkdir(target)

	for file in request.files.getlist("file"):
		print(file)
		filename = file.filename
		destination = "/".join([target, filename])
		print(destination)
		file.save(destination)
	return render_template("algorithmSelectionPage.html", filename=filename, destination=destination)
			#                   ^^^ should render algorithmSelectionPage.html

#below: onClick with "run" button on algorithmSelectionPage.html
@app.route("/run", methods = ['POST'])

def run():
	target = os.path.join(APP_ROOT,"files")

	
	#print(fasta)
	listDescription = []

	# reads in enzymes from .csv file
	#list_of_enzymes = pd.read_csv("RestrictionEnzymes.csv")
	#enzyme_list = list(list_of_enzymes['Recognition Sequence'])
	enzy_filename = "RestrictionEnzymes.csv"
	destination_enz = "/".join([target,enzy_filename])

	enzymes = pd.read_csv(destination_enz)
	enzymes = enzymes[0:84]
	enzy_pattern = enzymes['Recognition Sequence']
	enzyme_patterns = []

	for i in enzy_pattern:
		if 'S' not in i and 'Y' not in i and 'R' not in i and 'N' not in i and 'V' not in i and 'W' not in i and 'B' not in i and 'D' not in i and 'H' not in i and 'K' not in i and 'M' not in i:
		#if 'B' in i:
			enzyme_patterns.append(i)

	names_of_enzymes = []
	
	df_enzymes = pd.DataFrame()
	for i in enzyme_patterns:
		df_enzymes = df_enzymes.append(enzymes.loc[enzymes['Recognition Sequence'] == i],ignore_index=True)
	names_of_enzymes = list(df_enzymes['Enzymes'])
	print(names_of_enzymes)
	i = 0
	while i < 3:
		k = 0
		if i == 0:
			descriptor = "/numCuts"
		if i == 1:
			descriptor = "/RMS"
		if i == 2:
			descriptor = "/avgMW"
		while k < len(names_of_enzymes):
			enzyme = names_of_enzymes[k]
			#print(enzyme)
			
			stringDescription = enzyme + descriptor
			listDescription.append(stringDescription)
			k += 1
		i += 1
	print(len(listDescription))
	
	j = 0
	#fasta = phylopandas.read_fasta('fasta.fasta')
	seq_filename = "SequenceA-D_DF.csv"
	destination_seq= "/".join([target,seq_filename])
	fasta = pd.read_csv(destination_seq)
	fasta.drop_duplicates('sequence',keep="first")
	print(len(fasta))
	c = 0
	
	#for i in fasta['sequence']:

	for i in fasta['sequence']:
		num_of_cuts = []
		listOfRMS = []
		temp_cutLength_list = []
		listAvgMoMa = []
		tempFragmentMoMa_list = []
		#sequence = fasta.loc[i,'sequence']
		for j in enzyme_patterns:
			num_of_cuts.append(len(i.split(j)))
			split1 = i.split(j)
			for k in range(len(split1)):
				temp_cutLength = len(split1[k])
				temp_cutLength_list.append(temp_cutLength)
				tempFragment = split1[k]
				tempFragmentMoMa = 0
				for u in range(len(tempFragment)):
					if tempFragment[u] == 'A':
						tempFragmentMoMa = tempFragmentMoMa + 507
					if tempFragment[u] == 'C':
						tempFragmentMoMa = tempFragmentMoMa + 483
					if tempFragment[u] == 'G':
						tempFragmentMoMa = tempFragmentMoMa + 523
					if tempFragment[u] == 'T':
						tempFragmentMoMa = tempFragmentMoMa + 483
					tempFragmentMoMa_list.append(tempFragmentMoMa)
			rms = (sqrt(mean(square(temp_cutLength_list))))
			listOfRMS.append(rms)
			mean_moma = mean(tempFragmentMoMa_list)
			listAvgMoMa.append(mean_moma)
			print(listAvgMoMa)
			print (listOfRMS)
			print (temp_cutLength_list)
			#print(listAvgMoMa)


		#with open('/Users/sabihabarlaskar/Desktop/MSCS/PINC/code/enzymeCutting_SC/gui_app/files/number_of_cuts.csv', 'a') as f1:
		#	csvWri = csv.writer(f1, delimiter=',')
		#	csvWri.writerow(num_of_cuts)
			#listData.append(string1)
			#print(listData)
			#print(num_of_cuts)

		#with open('/Users/sabihabarlaskar/Desktop/MSCS/PINC/code/enzymeCutting_SC/gui_app/files/RMS', 'a') as f1:
		#	csvWri = csv.writer(f1, delimiter=',')
			#listData.append(string1)
			#print(listData)
		#	csvWri.writerow(listOfRMS)

		#with open('/Users/sabihabarlaskar/Desktop/MSCS/PINC/code/enzymeCutting_SC/gui_app/files/molarWeight.csv', 'a') as f2:
		#	csvWri2 = csv.writer(f2, delimiter=',')
			#listData.listAvgMoMa(string1)
			#print(listAvgMoMa)
		#	csvWri2.writerow(listAvgMoMa)"""
	
	cuts = "number_of_cuts.csv"
	rms = "RMS.csv"
	moMa = "molarWeight.csv"
	destination_cuts = "/".join([target,cuts])
	#print(destination_cuts)
	destination_rms = "/".join([target,rms])
	destination_moMa = "/".join([target,moMa])
	no_of_cuts = pd.read_csv(destination_cuts,header = None)
	rms_cuts = pd.read_csv(destination_rms,header = None)
	moMa = pd.read_csv(destination_moMa,header = None)
	final_data = pd.concat([no_of_cuts,rms_cuts,moMa],axis = 1 )
	final_file = "final_data.csv"
	destination_data = "/".join([target,final_file])
	final_data.to_csv(destination_data,header=listDescription)

	data = pd.read_csv(destination_data)
	print(data.head())
	data = data.drop('Unnamed: 0',axis=1)
	labels = []
	for i in fasta['description']:
	    labels.append(i.split('.')[0])

	se = pd.Series(labels)
	data['labels'] = se.values
	#data = data.drop('Unnamed: 0',axis=1)
	#datase.values
	data['labels'] = data['labels'].map({'A':1,'A1':1, 'A2':1,'A3':1,'A4':1,'A5':1,'A6':1,'B':2,'C':3,'D':4})
	features = data.drop('labels', axis = 1)
	labels = data['labels']
	print(data.head())
	
	X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.20, random_state=42,stratify=labels)
	clf = DecisionTreeClassifier()
	clf.fit(X_train,y_train)
	pred = clf.predict(X_test)
	print("Accuracy:",accuracy_score(y_test,pred))
	AccuracyDT = accuracy_score(y_test,pred)
	print("Recall:",recall_score(y_test,pred,average=None))
	RecallDT=recall_score(y_test,pred,average=None)
	print("Precision:",precision_score(y_test,pred,average=None))
	#A,B,C,D = fbeta_score(y_test,pred,beta=1,average='weighted')
	f1 = fbeta_score(y_test,pred,beta=1,average=None)
	#print("F1-score of A:",A)
	#print("F1-score of B:",B)
	#print("F1-score of C:",C)
	#print("F1-score of D:",D)
	print(list(pred)[0:10],list(y_test)[0:10])
	return render_template("complete.html",f1=f1)

			
	#print(len(num_of_cuts))
		



if __name__ == "__main__":
	app.run(port=5000)