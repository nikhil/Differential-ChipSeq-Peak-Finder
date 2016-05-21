import codecs
import json
#Open Text Gene Reference File
GeneReferenceFile = codecs.open("GeneReference.csv",'rb')

GeneReferenceFile.next()
GeneReferenceDictionary = {}
PromoterReferenceDictionary = {}
GeneChrList = []
GenePromoterList = []
GeneMap = {}
CurrentChrNum = ''
CurrentGeneNum = 0
GeneList = {}
for Singleline in GeneReferenceFile:
	line = Singleline.split()
	if not CurrentChrNum:
		CurrentChrNum = line[0]
	ChrNum = line[0]
	
	if(CurrentChrNum != ChrNum):
		#Do sorting and add  to List
		SortedGeneChrList = sorted(GeneChrList, key=lambda k: k['Start']) 
		SortedPromoterChrList = sorted(GenePromoterList, key=lambda k: k['Start']) 
		GeneReferenceDictionary[str(CurrentChrNum)] = SortedGeneChrList
		PromoterReferenceDictionary[str(CurrentChrNum)] = SortedPromoterChrList
		CurrentChrNum = ChrNum		
		GeneChrList = []
		GenePromoterList = []
	
	NewGeneElement = {'Start':int(line[2]),'End':int(line[3]),'Type':'Gene','Name':line[4]}
	NewPromoterElement = {}
	if str(line[1]) == '+':
		NewPromoterElement = {'Start':int(line[2]) - 1000, 'End':int(line[2])-1,'Type':'Promoter','Name':line[4]}
	elif str(line[1]) == '-':
		NewPromoterElement = {'Start':int(line[3])+1, 'End':int(line[3])+1000,'Type':'Promoter','Name':line[4]}
	GeneChrList.append(NewGeneElement)
	GenePromoterList.append(NewPromoterElement)

if(len(GeneChrList) > 0):
	SortedGeneChrList = sorted(GeneChrList, key=lambda k: k['Start']) 
	GeneReferenceDictionary[str(CurrentChrNum)] = SortedGeneChrList
if(len(GenePromoterList) > 0):
	SortedPromoterChrList = sorted(GenePromoterList, key=lambda k: k['Start']) 
	PromoterReferenceDictionary[str(CurrentChrNum)] = SortedPromoterChrList


ChipSeqchrOrder = ['chr1', 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX','chrY']


for SingleChr in ChipSeqchrOrder:
	CurrentPosition = 0
	CheckOverlap = 1
	GeneArray = GeneReferenceDictionary[str(SingleChr)]
	while CheckOverlap < len(GeneArray):
		CurrentGene = GeneArray[CurrentPosition]
		if 'OverLap' not in CurrentGene:
			CurrentGene['OverLap'] = ''
		NextGene = GeneArray[CheckOverlap]
		if(NextGene['Start']<CurrentGene['End']):
			CurrentGene['OverLap'] = CurrentGene['OverLap'] + str(NextGene['Start']) + '|' + str(NextGene['End']) + '|'+ NextGene['Name']
			CheckMore = CheckOverlap + 1
			while CheckMore < len(GeneArray):
				AdditionalGene = GeneArray[CheckMore]
				if(AdditionalGene['Start']<CurrentGene['End']):
					CurrentGene['OverLap'] = CurrentGene['OverLap'] + ',' + str(AdditionalGene['Start']) + '|' + str(AdditionalGene['End']) + '|'+ AdditionalGene['Name']
					CheckMore = CheckMore + 1
				else:
					break
		GeneArray[CurrentPosition] = CurrentGene				
		CheckOverlap = CheckOverlap + 1
		CurrentPosition = CurrentPosition + 1
	GeneArray[-1]['OverLap'] = ''
	GeneReferenceDictionary[str(SingleChr)] = GeneArray






GeneDictOutput = open('MouseGene.xls','w')
for chrStr in ChipSeqchrOrder:
	for SingleGene in GeneReferenceDictionary[str(chrStr)]:
		GeneDictOutput.write(chrStr +'\t' + str(SingleGene['Start']) + '\t' + str(SingleGene['End']) + '\t' + str(SingleGene['Name']) + '\t' + str(SingleGene['OverLap']) +'\n')



for SingleChr in ChipSeqchrOrder:
	CurrentPosition = 0
	CheckOverlap = 1
	GeneArray = PromoterReferenceDictionary[str(SingleChr)]
	while CheckOverlap < len(GeneArray):
		CurrentGene = GeneArray[CurrentPosition]
		if 'OverLap' not in CurrentGene:
			CurrentGene['OverLap'] = ''
		NextGene = GeneArray[CheckOverlap]
		if(NextGene['Start']<CurrentGene['End']):
			CurrentGene['OverLap'] = CurrentGene['OverLap'] + str(NextGene['Start']) + '|' + str(NextGene['End']) + '|'+ NextGene['Name']
			CheckMore = CheckOverlap + 1
			while CheckMore < len(GeneArray):
				AdditionalGene = GeneArray[CheckMore]
				if(AdditionalGene['Start']<CurrentGene['End']):
					CurrentGene['OverLap'] = CurrentGene['OverLap'] + ',' + str(AdditionalGene['Start']) + '|' + str(AdditionalGene['End']) + '|'+ AdditionalGene['Name']
					CheckMore = CheckMore + 1
				else:
					break
		GeneArray[CurrentPosition] = CurrentGene				
		CheckOverlap = CheckOverlap + 1
		CurrentPosition = CurrentPosition + 1
	GeneArray[-1]['OverLap'] = ''
	PromoterReferenceDictionary[str(SingleChr)] = GeneArray






PromoterDictOutput = open('MousePromoter.xls','w')
for chrStr in ChipSeqchrOrder:
	for SingleGene in PromoterReferenceDictionary[str(chrStr)]:
		PromoterDictOutput.write(chrStr +'\t' + str(SingleGene['Start']) + '\t' + str(SingleGene['End']) + '\t' + str(SingleGene['Name']) + '\t' + str(SingleGene['OverLap']) +'\n')

GeneDictOutput.close()
PromoterDictOutput.close()