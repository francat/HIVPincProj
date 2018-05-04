import os
from flask import Flask, render_template, request
import pandas as pd
import phylopandas as pp
from numpy import mean, sqrt, square
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, recall_score, precision_score,f1_score, fbeta_score
import csv

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
		#print("FILE:",file)
		filename = file.filename
		destination = "/".join([target, filename])
		#print(destination)
		file.save(destination)
		#global fasta_path = destination
		global fasta_path
		fasta_path = destination
		#print(fasta_path)
		#print(type(fasta_path))
		global fasta
		fasta = pp.read_fasta(fasta_path)
		#print(type(fasta))
		#print (fasta)

		

	return render_template("algorithmSelectionPage.html")     
			#                   ^^^ should render algorithmSelectionPage.html

#below: onClick with "run" button on algorithmSelectionPage.html
@app.route("/run", methods = ['POST'])

def run():
	#print(fasta_path)
	fasta = pp.read_fasta(fasta_path)
	#print(fasta)
	target = os.path.join(APP_ROOT,"files")
	listDescription = []
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
	#df_enzymes.to_csv("Enzymesused.csv")
	#print(names_of_enzymes)
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
	#print(len(listDescription))
	
	j = 0
	#fasta = phylopandas.read_fasta('fasta.fasta')
	seq_filename = fasta_path
	destination_seq= "/".join([target,seq_filename])
	#fasta = pd.read_csv(destination_seq)
	#fasta = df
	#fasta.drop_duplicates('sequence',keep="first")
	#print(len(fasta))
	c = 0
	
	#for i in fasta['sequence']:
	testdf = pd.read_csv("pink.csv")
	for seq in fasta['sequence']:
		seq = seq.upper()
		num_of_cuts = []
		listOfRMS = []
		temp_cutLength_list = []
		listAvgMoMa = []
		tempFragmentMoMa_list = []
		tempFragmentMoMa = 0
		#inputfeature =[]
		#sequence = fasta.loc[i,'sequence'
		for j in enzyme_patterns:
			num_of_cuts.append(len(seq.split(j)))
			split1 = seq.split(j)
			for k in range(len(split1)):
				temp_cutLength = len(split1[k])
				temp_cutLength_list.append(temp_cutLength)
				tempFragment = split1[k]
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
			listAvgMoMa.append(num_of_cuts)
	print (len(num_of_cuts),len(listOfRMS),len(listAvgMoMa))
	cuts = pd.DataFrame([num_of_cuts])
	RMSdf=pd.DataFrame([listOfRMS])
	MOMadf = pd.DataFrame([listAvgMoMa])
	inputfeature = pd.concat([cuts,RMSdf,MOMadf],axis=1)
	#inputcsv=inputfeature.to_csv(header=stringDescription)
	print(type(inputfeature))
	print(inputfeature.shape)
	print(inputfeature.head)
	#print(inputcsv)
	#testinput= pd.read_csv(inputcsv)
	#featureinput = pd.DataFrame(listofinput)
	#print(len(featureinput))
	#print(testinput)
	#print(pd.testinput.shape)
		

	return render_template("complete.html",f1=listAvgMoMa)

if __name__ == "__main__":
	app.run(port=8000)