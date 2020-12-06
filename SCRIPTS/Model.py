import numpy as np
import pandas as pd
import pickle
from numpy import genfromtxt
from sklearn import linear_model
from sklearn.linear_model import RidgeClassifier

def strip_first_col(fname, delimiter=None):
	with open(fname, 'r') as fin:
		for line in fin:
			try:
				yield line.split(delimiter, 1)[1]
			except IndexError:
				continue

def ModelRidge(args, tsvfile, predfile):
	if args.type == "prot":
		modelfile=args.data_path+"/MODEL/PROTEIN.pkl"
		f = open(modelfile, 'rb')
		clf = pickle.load(f)
		f.close()
	
		np.random.seed(0)
		TSV=open(tsvfile,"r")
		FO=open(predfile,"w")
		
		for line in TSV:
			if not line.startswith("ID"):
				arr=line.strip().split("\t")
				X_test=np.array(arr[1:3117], dtype='float')
				X_test=X_test.reshape(1,-1)
				Y_pred = clf.predict(X_test)
				if hasattr(clf, "predict_proba"):
					prob_pos = clf.predict_proba(X_test)[:, 1]
				else:
					prob_pos = clf.decision_function(X_test)
				FO.write("%s\t%3f\t%1f\n" % (arr[0], float(prob_pos), float(Y_pred)))
				
		TSV.close()
		FO.close()

	elif args.type == "nucl":
		modelfile=args.data_path+"/MODEL/DNA.pkl"
		f = open(modelfile, 'rb')
		clf = pickle.load(f)
		f.close()
		np.random.seed(0)
		TSV=open(tsvfile,"r")
		FO=open(predfile,"w")
		for line in TSV:
			if not line.startswith("ID"):
				arr=line.strip().split("\t")
				X_test=np.array(arr[1:3117], dtype='float')
				X_test=X_test.reshape(1,-1)
				Y_pred = clf.predict(X_test)
				if hasattr(clf, "predict_proba"):
					prob_pos = clf.predict_proba(X_test)[:, 1]
				else:
					prob_pos = clf.decision_function(X_test)

				FO.write("%s\t%3f\t%1f\n" % (arr[0], float(prob_pos), float(Y_pred)))

		TSV.close()
		FO.close()
