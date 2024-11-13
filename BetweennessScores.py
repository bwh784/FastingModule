# ! //usr/bin/env python
import sys, os, string, math, networkx, numpy
#import pylab, numpy

## examine a network's degree, degree distribution, connectedness

    
def main():
	#if (len(sys.argv)<>3):
		#print 'Please provide two file'
		#sys.exit()
			
	
	file1 = open("FastingFatStrong-Aging_interactions_in_PPI201806_BipartiteNetwork.txt","r")
	file2 = open("FastingFatStrong-Aging_interactions_in_PPI201806_BipartiteNetwork_proteins_venn.txt","r") 
	file3 = open("FastingFatStrong_genes_in_PPI201806_network_proteins_BC_fasting-aging3.txt","w") 

	 
	node=[]
	idx=0
	for line in file1:
		tmp=line.split()
		if(tmp[0].strip() not in node):
			node.append(tmp[0].strip())
		if(tmp[1].strip() not in node):
			node.append(tmp[1].strip())	
	print (len(node))
	file1.close()
	
	G=networkx.Graph()
	file1 = open("FastingFatStrong-Aging_interactions_in_PPI201806_BipartiteNetwork.txt","r")
	
	for i in range(len(node)):
		G.add_node(node[i])
	
	for line in file1:
		tmp=line.split()
		G.add_edge(tmp[0].strip(),tmp[1].strip())		
	file1.close()
	
	aging=[]
	fasting=[]
	
	for line in file2:
		tmp=line.split()
		if(tmp[1]=="aging"):
			aging.append(tmp[0])
		if(tmp[1]=="fasting"):
			fasting.append(tmp[0])	
		if(tmp[1]=="common"):
			#aging.append(tmp[0])
			fasting.append(tmp[0])				
	file2.close()
	
	Node_betweenness={}
	
	idx=0
	for i in range(len(aging)):
		print (i)
		for j in range(len(fasting)):	
			if(networkx.has_path(G,source=aging[i],target=fasting[j]) and not (aging[i]==fasting[j])):
				SP=networkx.all_shortest_paths(G,source=aging[i],target=fasting[j])
				#print (len(list(SP)))
				for p in list(SP):		
					for v in range(1,len(p)-1):
						if(p[v] not in Node_betweenness):
							Node_betweenness[p[v]]=[0 for t in range(1)]
						Node_betweenness[p[v]][0]=Node_betweenness[p[v]][0]+1
					idx=idx+1
	file2 = open("FastingFatStrong-Aging_interactions_in_PPI201806_BipartiteNetwork_proteins_FastingFatStrong.txt","r") 	
	
	for line in file2:
		tmp=line.split()
		if(tmp[0] in Node_betweenness):
			file3.write("%s\t"%line.strip())
			BS=Node_betweenness[tmp[0]][0]/float(idx)
			file3.write("%s\n"%BS)
		else:
			file3.write("%s\t"%line.strip())
			file3.write("%s\n"%'NA')				
	file3.close()
	file2.close()
     
if __name__=='__main__':
	main()
