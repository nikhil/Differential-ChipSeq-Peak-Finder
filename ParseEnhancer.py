import requests
import json
from bs4 import BeautifulSoup

response = requests.get('http://enhancer.lbl.gov/cgi-bin/imagedb3.pl?show=1;page=1;search.form=no;form=search;search.org=Mouse;page_size=20000;search.result=yes;action=search')
HtmlBody = response.text
HtmlParser = BeautifulSoup(HtmlBody,'html.parser')
TableBody = HtmlParser.find('tbody')
TableRows = TableBody.find_all('tr')
EnhancerRows = TableRows[2:]
EnhancerData = {}
for SingleEnhancerRow in EnhancerRows:
    SingleEnhancerInfo = SingleEnhancerRow.find_all('td')
    ID = SingleEnhancerInfo[0].get_text().strip(' \t\n\r')
    Position = SingleEnhancerInfo[3].get_text().strip(' \t\n\r')
    BracketingGene = SingleEnhancerInfo[4].get_text().strip(' \t\n\r')
    Expression = SingleEnhancerInfo[5].find('img').get('title').strip(' \t\n\r')
    ChrNameList = Position.split(':')
    ChrName = ChrNameList[0]    
    ChrRange = ChrNameList[1].split('-')
    Start = int(''.join(ChrRange[0].split(',')))
    End = int(''.join(ChrRange[1].split(',')))
    if ChrName not in EnhancerData:
    	EnhancerData[ChrName] = []
    NewEnhancerElem = {'Name':str(ID),'Start':Start, 'End':End, 'BracketingGene' : BracketingGene, 'Expression': Expression}
    EnhancerData[str(ChrName)].append(NewEnhancerElem)

#Properly Order Enhancer Data
for key in EnhancerData:
	UnsortedList = EnhancerData[key]
	SortedEnhancerList = sorted(UnsortedList, key=lambda k: k['Start'])
	EnhancerData[key] = SortedEnhancerList 


ChipSeqchrOrder = ['chr1', 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX','chrY']

for SingleChr in ChipSeqchrOrder:
	CurrentPosition = 0
	CheckOverlap = 1
	GeneArray = []
	if str(SingleChr) in EnhancerData:
		GeneArray = EnhancerData[str(SingleChr)]
	while CheckOverlap < len(GeneArray):
		CurrentGene = GeneArray[CurrentPosition]
		if 'OverLap' not in CurrentGene:
			CurrentGene['OverLap'] = ''
		NextGene = GeneArray[CheckOverlap]
		if(NextGene['Start']<CurrentGene['End']):
			CurrentGene['OverLap'] = CurrentGene['OverLap'] + str(NextGene['Start']) + '|' + str(NextGene['End']) + '|'+ NextGene['Name'] +'|' + NextGene['BracketingGene'] + '|' + NextGene['Expression']
			CheckMore = CheckOverlap + 1
			while CheckMore < len(GeneArray):
				AdditionalGene = GeneArray[CheckMore]
				if(AdditionalGene['Start']<CurrentGene['End']):
					CurrentGene['OverLap'] = CurrentGene['OverLap'] + ',' + str(AdditionalGene['Start']) + '|' + str(AdditionalGene['End']) + '|'+ AdditionalGene['Name'] + '|' + AdditionalGene['BracketingGene'] + '|' + AdditionalGene['Expression']
					CheckMore = CheckMore + 1
				else:
					break
		GeneArray[CurrentPosition] = CurrentGene				
		CheckOverlap = CheckOverlap + 1
		CurrentPosition = CurrentPosition + 1
	if len(GeneArray) > 0:
		GeneArray[-1]['OverLap'] = ''
		EnhancerData[str(SingleChr)] = GeneArray

EnhancerDictOutput = open('MouseEnhancer.xls','w')
for chrStr in ChipSeqchrOrder:
	if str(chrStr) in EnhancerData:
		for SingleGene in EnhancerData[str(chrStr)]:
			EnhancerDictOutput.write(chrStr +'\t' + str(SingleGene['Start']) + '\t' + str(SingleGene['End']) + '\t' + str(SingleGene['Name']) + '\t' + str(SingleGene['BracketingGene']) + '\t' + str(SingleGene['Expression']) + '\t' + str(SingleGene['OverLap']) +'\n')

EnhancerDictOutput.close()




