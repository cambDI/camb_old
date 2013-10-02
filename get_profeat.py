#!/bin/env python
# Get Profeat descriptors for a list of UNIPROT identifiers. 
# Isidro Cortes Sept/2013

import os, sys
import numpy as np
from propy import PyPro
from propy.GetProteinFromUniprot import GetProteinSequence

filename=str(sys.argv[1])
file=open(filename,'r')
lines=[x.strip() for x in file.readlines()]
output=(sys.argv[2])
#for i in lines:
for count,i in enumerate(lines):
	print count
	proteinsequence=GetProteinSequence(i)
	DesObject=PyPro.GetProDes(proteinsequence)

	CTD = DesObject.GetCTD() ##calculate 147 CTD descriptors
	AAcomp = DesObject.GetAAComp() ##calculate 20 amino acid composition descriptors
	paac=DesObject.GetPAAC() ##calculate 30 pseudo amino acid composition descriptors 
	DPC=DesObject.GetDPComp() #dipeptide composition descriptors
	TPC = DesObject.GetTPComp() # tripeptide composition descriptors
	QSO=DesObject.GetQSO()
	allD=dict(CTD.items() + AAcomp.items() + paac.items() + DPC.items()+TPC.items()+QSO.items())
	allD=np.asarray(allD.items())
	allD=np.transpose(allD)
	if count == 0:
		np.savetxt(output,allD,fmt='%s')
	else:
		f_handle=open(output,'a')
		np.savetxt(f_handle,allD[1],fmt='%s',newline=" ")
		f_handle.write("\n")
		f_handle.close()

file.close()
