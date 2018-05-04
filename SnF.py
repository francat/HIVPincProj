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
from sklearn.model_selection import KFold
from sklearn import metrics
from sklearn.model_selection import cross_val_score
from sklearn import datasets, linear_model, cross_validation, grid_search
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict

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
		filename = file.filename
		destination = "/".join([target, filename])
		file.save(destination)
		global fasta_path
		fasta_path = destination
		global fasta
		fasta = pp.read_fasta(fasta_path)


	return render_template("algorithmSelectionPage.html")       
	#^^^ should render algorithmSelectionPage.html

#below: onClick with "run" button on algorithmSelectionPage.html
@app.route("/run", methods = ['POST'])

def run():
	target = os.path.join(APP_ROOT,"files")
	final_file = "final_data.csv"
	destination_data = "/".join([target,final_file])
	#final_data.to_csv(destination_data,header=listDescription)


## creates feature vectors for input sequenes, input sequences saved as global fasta variable

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
	
	#creates string of descriptors for dataframe
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
			stringDescription = enzyme + descriptor
			listDescription.append(stringDescription)
			k += 1
		i += 1
	
	j = 0
	seq_filename = fasta_path
	destination_seq= "/".join([target,seq_filename])
	c = 0
	
	#for i in fasta['sequence']:
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
				tempFragmentMoMa = 0
				numOfA = tempFragment.count('A')
				numOfC = tempFragment.count('C')
				numOfG = tempFragment.count('G')
				numOfT = tempFragment.count('T')
				#print(numOfG, numOfT, numOfT, numOfA)
				tempFragmentMoMa = numOfA*507 + numOfC*483 + numOfG*523 + numOfT*483
				#print(tempFragmentMoMa)
				tempFragmentMoMa_list.append(tempFragmentMoMa)
			rms = (sqrt(mean(square(temp_cutLength_list))))
			listOfRMS.append(rms)
			mean_moma = mean(tempFragmentMoMa_list)
			listAvgMoMa.append(mean_moma)
	#creates feature vector out of input sequence(s)
		cuts = pd.DataFrame([num_of_cuts])
		#print("CUTS",cuts)
		RMSdf=pd.DataFrame([listOfRMS])
		#print("RMS",RMSdf)
		MOMadf = pd.DataFrame([listAvgMoMa])
	#print("MoMa:",MOMadf)
	test_data = pd.concat([cuts,RMSdf,MOMadf],axis=1,ignore_index=True)
	test_data.columns=listDescription
	#print("Test features:",test_data)
	test_data['labels'] = 2
	#print(test_data)
	#test_data['labels'] = test_data['labels'].map({'A':1,'A1':1, 'A2':1,'A3':1,'A4':1,'A5':1,'A6':1,'B':2,'C':3,'D':4})
	test_features = test_data.drop('labels', axis = 1)
	test_labels = test_data['labels']
	print("Test features",test_features)



	#Reading in the training data
	data = pd.read_csv(destination_data)
	data = data.drop('Unnamed: 0',axis=1)
	seq_filename = "SequenceA-D_DF.csv"
	destination_seq= "/".join([target,seq_filename])
	train_data = pd.read_csv(destination_seq)
	se = create_labels(train_data)
	#print(se.values)
	data['labels'] = se.values
	data['labels'] = data['labels'].map({'A':1,'A1':1, 'A2':1,'A3':1,'A4':1,'A5':1,'A6':1,'B':2,'C':3,'D':4})

	train_features = data.drop('labels', axis = 1)
	train_labels = data['labels']
	#print(data.head())
	X_train, X_test, y_train, y_test = train_test_split(train_features, train_labels, test_size=0.20, random_state=42,stratify=train_labels)
	output=defaultdict(list)

	#CHECK BOX CONDITIONS
	if request.form.get('SVM') and not request.form.get('Naive-Bayes') and not request.form.get('Decision-Tree'):
		clf = SVC(C=100.0,kernel='linear')
		clf_svc = cross_val(clf,X_train,y_train)
		pred = clf_svc.predict(X_test)
		pred_label = predict_labels(pred)
		y_true = predict_labels(y_test)
		acc,f1,prec,recall= predict_scores(pred,y_test)
		output = create_output(acc,f1,prec,recall,output)
		print(output)
		return render_template("complete.html",predictions = pred_label[1:10], actual = y_true[1:10] ,output = output,test_features = X_test,len_train = len(X_train),len_test = len(X_test))#,f1=f1SVM,Recall=RecallSVM, precision_score=precision_scoreSVM)

	if request.form.get('Decision-Tree') and not request.form.get('SVM') and not request.form.get('Naive-Bayes') :
		clf = DecisionTreeClassifier()
		clf_dt = cross_val(clf,train_features,train_labels)
		pred = clf_dt.predict(X_test)
		pred_label = predict_labels(pred)
		y_true = predict_labels(y_test)
		acc,f1,prec,recall= predict_scores(pred,y_test)
		output = create_output(acc,f1,prec,recall,output)
		return render_template("complete.html",predictions = pred_label[1:10], actual = y_true[1:10],output = output,test_features = X_test,len_train = len(X_train),len_test = len(X_test))#,f1=f1SVM,Recall=RecallSVM, precision_score=precision_scoreSVM)


	if request.form.get('Naive-Bayes') and not request.form.get('SVM') and not request.form.get('Decision-Tree'):
		clf = GaussianNB()
		clf_g = cross_val(clf,train_features,train_labels)
		pred = clf_g.predict(X_test)
		pred_label = predict_labels(pred)
		y_true = predict_labels(y_test)
		acc,f1,prec,recall= predict_scores(pred,y_test)
		output = create_output(acc,f1,prec,recall,output)
		print(output)
		return render_template("complete.html",predictions = pred_label[1:10], actual = y_true[1:10],output = output,test_features = X_test,len_train = len(X_train),len_test = len(X_test))#,f1=f1SVM,Recall=RecallSVM, precision_score=precision_scoreSVM)

	

def create_labels(data):
	labels = []
	for i in data['description']:
		labels.append(i.split('.')[0])

	se = pd.Series(labels)
	
	return se

def cross_val(clf,features,labels):
	kf = KFold(n_splits=10)
	kf.get_n_splits(features)
	for train_index, test_index in kf.split(features, labels):
			# specific ".loc" syntax for working with dataframes
			cv_x_train, cv_x_test = features.iloc[train_index], features.iloc[test_index]
			cv_y_train, cv_y_test = labels.iloc[train_index], labels.iloc[test_index]
			clf.fit(cv_x_train, cv_y_train)
	return clf



def predict_labels(pred):
	pred_labels = []
	for i in pred:
			if i == 1:
				pred_labels.append('A')
			elif i == 2:
				pred_labels.append('B')
			elif i == 3:
				pred_labels.append('C')
			else:
				pred_labels.append('D')
	return pred_labels


	
def predict_scores(pred,test_labels):
	acc = accuracy_score(test_labels,pred)
	f1 = f1_score(test_labels,pred,average='weighted')
	prec = precision_score(test_labels,pred,average='weighted')
	rec = recall_score(test_labels,pred,average='weighted')
	return (acc,f1,prec,rec)


def create_output(acc,f1,prec,recall,output):
	output['Accuracy'].append(acc)
	output['f1'].append(f1)
	output['Precision'].append(prec)
	output['Recall'].append(recall)
	return output






	







if __name__ == "__main__":
	app.run()
