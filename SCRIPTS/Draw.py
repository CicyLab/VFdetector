import argparse
import os
import sys
import logging
import matplotlib.pyplot as plt
				

def genecov(args):
	GeneLenFile=args.database
	samfile=args.sam
	Gname=args.gid
	outfile=args.outfile

	Dgl={}
	GLF=open(GeneLenFile,"r")
	for line in GLF.readlines():
		array=line.strip().split()
		Dgl.update({array[0]:int(array[1])})

	GLF.close()

	Dgd={}
	SF=open(samfile,"r")
	for line in SF.readlines():
		if not line.startswith("@"):
			array=line.strip().split()
			if array[2]==Gname:
				posstart=array[3]
				matchlen=len(array[9])
				for i in range(matchlen):
					j=int(posstart)+i
					Dgd[j]=Dgd.setdefault(j, 0) + 1

	SF.close()

	x=[]
	y=[]
	total=0
	cov=0
	maxdepth=0
	for i in range(Dgl[Gname]):
		count=Dgd.setdefault(i, 0)
		if count > 0:
			x.append(int(i))
			y.append(int(count))
			total=total+count
			cov = cov + 1
			if count > maxdepth:
				maxdepth=count

	AverageDepth=total/(Dgl[Gname]+1)
	Coverage=cov/(Dgl[Gname]+1)*100
	Ave=[AverageDepth for i in x]

	plt.figure(figsize=(12,8))
	if Dgl[Gname] < 100:
		plt.xticks(range(0, Dgl[Gname], 10))
	elif Dgl[Gname] < 500:
		plt.xticks(range(0, Dgl[Gname], 50))
	elif Dgl[Gname] < 1000:
		plt.xticks(range(0, Dgl[Gname], 100))
	elif Dgl[Gname] < 2000:
		plt.xticks(range(0, Dgl[Gname], 200))
	elif Dgl[Gname] < 5000:
		plt.xticks(range(0, Dgl[Gname], 500))
	else:
		plt.xticks(range(0, Dgl[Gname], 1000))

	if maxdepth < 100:
		plt.yticks(range(0, maxdepth, 10))
	elif maxdepth < 500:
		plt.yticks(range(0, maxdepth, 50))
	elif maxdepth < 1000:
		plt.yticks(range(0, maxdepth, 100))
	elif maxdepth < 2000:
		plt.yticks(range(0, maxdepth, 200))
	elif maxdepth < 5000:
		plt.xticks(range(0, maxdepth, 500))
	else:
		plt.xticks(range(0, maxdepth, 1000))

	plt.xlabel("position (bp)", fontdict={'size': 12})
	plt.ylabel("depth(X)", fontdict={'size': 12})
	plt.scatter(x,y,s=10)
	plt.plot(x, Ave, color='red', label='AverageDepth: ~ '+ str(int(AverageDepth))+" X\n"
	        +"1X Coverage: ~ " + str(int(Coverage)) + "%")
	plt.legend(loc='best', fontsize=12, markerscale=0.5)
	plt.title(Gname,fontsize=16)
	plt.savefig(outfile,dpi=300,bbox_inches='tight')

