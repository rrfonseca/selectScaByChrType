#sum number of OK sites per scaffold
#sum number of OK sites in 1000 largest scaffolds

#NOTE: sex sca will show up has having zero sites covered

# Xenodacnis_parina.map1.notRepMask.m10k.notXZ.bed
#s1      0       140
#s1      191     6194
#s1      6237    13912

# Xenodacnis.list.median.m10k.M10_scaffSex.tsv 
#Scaffold        Length  X_Z     Abnormal_sex_linked     pVal    homogametic_median      heterogametic_median
#s1      811559  FALSE   FALSE   0.0537612275072474      0.970514700306293       0.990720483969551
#s2      716767  FALSE   FALSE   0.25314167772848        0.980427716237631       0.989210890801131
#s3      561130  FALSE   FALSE   0.0113849729231821      1.26477134388222        1.05932074807368

import sys

minScaSizeKb=int(sys.argv[1])

inOKSites=sys.argv[2] #"Xenodacnis_parina.map1.notRepMask.m10k.notXZ.bed"
inScaLength=sys.argv[3] #Xenodacnis_parina.fa.oldIds.sizes.chr.bed 

prefix=inOKSites[:-4]
outName="%s_min%sKb.sum.txt"%(prefix,minScaSizeKb) 

out=open(outName,"w")

infos=open(inScaLength)
info=infos.readline()
info=infos.readline()

scaL=[]
scaD={} #[sca]=finalSize

while info:
	x=info.rstrip().split()
	sca=x[0]
	size=int(x[1])
	scaD[sca]=0
	scaL.append((sca,size))

	info=infos.readline()

infos=open(inOKSites)
info=infos.readline()

while info:
	x=info.rstrip().split()
	sca=x[0]
	finalSize=scaD[sca]+(int(x[2])-int(x[1]))

	scaD[sca]=finalSize

	info=infos.readline()

out.write("scaffold\tsizeAssembly\tsizeNotRepeat\tpRepeats\n")
i=0
totalOK=0
totalGenome=0
size1000top=0
nscaffolds=0
totalSizeSubset=0

while i<len(scaL):
	sca=scaL[i][0]
	fullSize=scaL[i][1]
	sizeOKsites=scaD[sca]

	if sizeOKsites>0:
		pRepeats=round((fullSize-sizeOKsites)*100/fullSize*0.01,2)
	
		if sizeOKsites >= minScaSizeKb*1000:
			out.write("%s\t%s\t%s\t%s\n"%(sca,fullSize,sizeOKsites,pRepeats))
			totalOK+=sizeOKsites
			if i<1000:
				size1000top+=sizeOKsites
			nscaffolds+=1
		totalSizeSubset+=fullSize
		
	totalGenome+=fullSize

	i+=1
out.close()

if totalSizeSubset>0:
	pRepeats=(totalSizeSubset-totalOK)*100/totalSizeSubset*0.01 #pRepeats in the subset
	pSubset=(totalSizeSubset)*1000/totalGenome*0.001 #size of subset in the genome
	out=open("%s_min%sKb.summary.log"%(prefix,minScaSizeKb),"w")

	out.write("%s\tsum of all sites in Kb\n"%(int(totalSizeSubset/1000)))
	out.write("%s\tfraction of the genome\n"%(pSubset))
	out.write("%s\ttotal size after filtering in Kb \n"%(int(totalOK/1000)))
	out.write("%s\tfraction removed\n"%(pRepeats))
	out.write("%s\tnumber scaffolds remaining\n"%(nscaffolds))
	out.write("%s\tTotal sites in top 1000 scaffolds in Kb (fraction of the genome: %s)\n"%(int(size1000top/1000),size1000top*1000/totalGenome*0.001))
else:
	print "No information in file %s\n"%(inOKSites)

out.close()
