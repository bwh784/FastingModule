#!/usr/bin/env python

import getopt,sys, urllib, time
import random as rd
import networkx as nx
from scipy.stats import norm
import numpy as np

def ds(G,node,geneset):
	D=[]
	for i in range(len(geneset)):
		if(nx.has_path(G,node,target=geneset[i])):
			Lt=nx.shortest_path_length(G,node,target=geneset[i])
			D.append(Lt)
	if(len(D)>0):
		Ls=min(D)
		return Ls
	else:
		return -1
	

def distance(G,geneset1,geneset2):
	D=0
	idx=0
	for i in range(len(geneset1)):
		dt=ds(G,geneset1[i],geneset2)
		if(dt>=0):
			idx=idx+1
			D=D+dt
	if(idx>0):
		return float(D)/idx
	else:
		return -1



def main():


	file0=open("/home/rw99/mypython/PHdrugs/FastingFat_genes","r")  # This will be replaced with your module genes/proteins
	file1=open("/home/rw99/mypython/PHdrugs/Human_Interactome_proteins.txt","r")	# This is the protein list in the human interactome
	file2=open("/home/rw99/mypython/PHdrugs/Human_Interactome.txt","r") # This is the human interactome (protein-protein interactions)
	file3=open("/home/rw99/mypython/PHdrugs/Alzheimerdisease_Phenopedia.txt","r") # This is the drug target file
	file4=open("/home/rw99/mypython/PHdrugs/FastingFat_Alzheimerdisease_proximity.txt","w") # This is the output file
	
	##reading the module genes/proteins
	subprot=[]
	for line in file0:
		tmp=line.split()
		subprot.append(tmp[0])
	file0.close()

	
	##reading the proteins in the human interactome
	prot=[]
	for line in file1:
		tmp=line.split()
		prot.append(tmp[0])
	file1.close()	

	##Store the human interactome in a gragh format
 	G=nx.Graph()
	for line in file2:
		tmp=line.split()
		G.add_node(tmp[1])
		G.add_node(tmp[3])
		G.add_edge(tmp[1],tmp[3])	
	file2.close()
	
	#Calculate the network proximity from drug targets to disease modules
	idx=0
	for line in file3:
		tmp=line.split("\t")
		DM=[]
		idx=idx+1
		print idx
		file4.write("%s\t"%tmp[0])
		for i in range(2,len(tmp)):
			if(G.has_node(tmp[i].strip())):
				DM.append(tmp[i].strip())
		#print DG
		d=distance(G,DM,subprot)
		#print d
		file4.write("%s\t"%round(d,2))
		#file4.write("\n")			
 		
		#Evaluate the significance of the proximity by randomly selecting a module of the same size	
		RR=100  # repeat 1000 times
		Rdd=[]
		for k in range(RR):
			i=0
			print k
			randDM=[]
			randprot=[]
			dsize=len(DM)
			ssize=len(subprot)
			while (i<dsize):
				r=rd.randint(0,len(prot)-1)
				randDM.append(prot[r])
				i=i+1
			i=0
			while (i<ssize):
				r=rd.randint(0,len(prot)-1)
				randprot.append(prot[r])
				i=i+1	
			dd=distance(G,randDM,randprot)
			if(dd>0):
				Rdd.append(dd)
				
		u=np.mean(Rdd)
		o=np.std(Rdd)
		p=norm(u,o).cdf(d)
		file4.write("%s\t%s\t%s\t%s\n"%(u,o,p,dsize))		
			

	
if __name__ == "__main__":
    main()