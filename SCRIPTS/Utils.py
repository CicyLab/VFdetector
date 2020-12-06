import os
import sys



def TwoDimDict(thedict, key_a, key_b, val):
	if key_a in thedict.keys():
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})
