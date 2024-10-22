#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import statistics 
from collections import OrderedDict

# Store simplified mapping file into a dictionary
def read_parsed_paf(file_in):
	paf={}
	with open(file_in, 'r') as f:   
		for line in f: 
			parts = line.strip().split()

			paf[parts[0]] = {
			"ref_scaffold": str(parts[1]),
			"ref_start":  int(parts[2]),
			"ref_end": int(parts[3]),
			"readname": str(parts[4]),
			"readsize": int(parts[5]),
			"strand":str(parts[6]),
			"number_repet": int(parts[7]),
			"read_start":list([int(i) for i in parts[8].split(',')]),
			"read_end": list([int(i) for i in parts[9].split(',')]),
			"MQ":float(parts[10]),


               	 }
	return paf


def compare_reads(paf1):

	# In the dictionnary
	for k1,v1 in list(paf1.items()):
		
		paf2={k: v for k, v in paf1.items()  if (v["readname"] == v1["readname"] and k != k1) }		

		# If the read name, the scaffold and the is the strand is the same and the multiple mapping position of the reads are within 200 bp, it is considered as a repetition
		for k2,v2 in list(paf2.items()):
			#print (k1,k2,"\n")
			if v1["readname"]==v2["readname"] and v1["ref_scaffold"]==v2["ref_scaffold"] and v1["strand"]==v2["strand"] and ((v2["ref_end"]-200 <= v1["ref_end"] and v1["ref_end"] <= v2["ref_end"]+200) or (v2["ref_start"]-200 <= v1["ref_start"] and v1["ref_start"] <= v2["ref_start"]+200)) :
				read_start=list(OrderedDict.fromkeys(v1["read_start"]+v2["read_start"]))
				read_end=list(OrderedDict.fromkeys(v1["read_end"]+v2["read_end"]))
				v1["ref_start"]= min(v1["ref_start"],v2["ref_start"])
				v1["ref_end"]= max(v1["ref_end"],v2["ref_end"])
				v1["read_start"]=read_start
				v1["read_end"]=read_end
				v1["number_repet"]= len(v1["read_start"])
				
				if v1["number_repet"]>=v2["number_repet"] : 
					v1["MQ"]=v1["MQ"]
				else :
					v1["MQ"]=v2["MQ"]
			# If the read name, scaffold, strand and start position is the same, it is considered as a repetition 	
			elif v1["readname"]==v2["readname"] and v1["ref_scaffold"]==v2["ref_scaffold"] and v1["strand"]==v2["strand"] and v2["ref_end"]==v1["ref_end"] and v1["ref_start"] != v2["ref_start"] :
				read_start=list(OrderedDict.fromkeys(v1["read_start"]+v2["read_start"]))
				read_end=list(OrderedDict.fromkeys(v1["read_end"]+v2["read_end"]))
				v1["ref_start"]=  min(paf[k1]["ref_start"],v2["ref_start"])
				v1["ref_end"]= v2["ref_end"]
				v1["read_start"]= read_start
				v1["read_end"]=read_end
				v1["number_repet"]= len(v1["read_start"])

				if v1["number_repet"]>=v2["number_repet"] : 
					v1["MQ"]=v1["MQ"]
				else :
					v1["MQ"]=v2["MQ"]
			# If the read name, scaffold, strand and end position is the same, it is considered as a repetition 	
			elif v1["readname"]==v2["readname"] and v1["ref_scaffold"]==v2["ref_scaffold"] and v1["strand"]==v2["strand"] and v2["ref_end"] != v1["ref_end"] and v1["ref_start"] == v2["ref_start"] :
				read_start=list(OrderedDict.fromkeys(v1["read_start"]+v2["read_start"]))
				read_end=list(OrderedDict.fromkeys(v1["read_end"]+v2["read_end"]))
				v1["ref_start"]= v1["ref_start"]
				v1["ref_end"]= max(v1["ref_end"],v2["ref_end"])
				v1["read_start"]= read_start
				v1["read_end"]=read_end
				v1["number_repet"]= len(v1["read_start"])	

				if v1["number_repet"]>=v2["number_repet"] : 
					v1["MQ"]=v1["MQ"]
				else :
					v1["MQ"]=v2["MQ"]
				#del paf1[k2]
			#print (paf1[k2])

		#del paf1[k2]
				
		#for key in paf2.keys():
			#if k2 in paf1:		
				#del paf1[key]
				#del paf1[k2]
				
	return paf1


def write_paf1(paf):

	for k1 in sorted(paf.keys()):
		print ("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(paf[k1]["ref_scaffold"],paf[k1]["ref_start"],paf[k1]["ref_end"],paf[k1]["readname"],paf[k1]["readsize"],paf[k1]["number_repet"],np.median([abs(y-x) for x, y in zip(paf[k1]["read_start"], paf[k1]["read_end"])]),paf[k1]["read_start"],paf[k1]["read_end"],paf[k1]["MQ"],paf[k1]["strand"]))


		
def main():

	file_in = sys.argv[1]
	paf=read_parsed_paf(file_in)
	paf2=compare_reads(paf)
	#paf3=remove_dup(paf2)
	write_paf1(paf2)

	
	return 1

if __name__ == '__main__':
    main()
