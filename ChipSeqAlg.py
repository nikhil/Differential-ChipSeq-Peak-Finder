from os.path import isfile, splitext, join
import sys
import csv
import codecs
import json

if(len(sys.argv) == 1 or sys.argv[1] == '-help'):
    print "Finds differential peaks within 2 given Chip-Seq files \nUSAGE: \nStart a webserver for a short period: python ChipSeqAlg.py --website \nStart a webserver as a service: sudo supervisord -c supervisord.conf\nTerminal: python ChipSeqAlg.py [ChipSeqFile1] [ChipSeqFile2]"
    sys.exit()


if( sys.argv[1] == '--website'):
    from flask import Flask, render_template, request
    from werkzeug import secure_filename
            
PeakFile1 = None
PeakFile2 = None
GeneFile = codecs.open('MouseGene.xls','rb')
GeneFile1 = codecs.open('MouseGene.xls','rb')
PromoterFile = codecs.open('MousePromoter.xls','rb')
PromoterFile1 = codecs.open('MousePromoter.xls','rb')
EnhancerFile = codecs.open('MouseEnhancer.xls','rb')
EnhancerFile1 = codecs.open('MouseEnhancer.xls','rb')
#print sys.getsizeof(PeakFile1)
PeakFile1Ready = 0
PeakFile2Ready = 0
DifferentDict = {}
SameDict = {}
SameLines = 1
DifferentLines = 1
SamePartNumber = 1
DifferentPartNumber = 1


def PlotOutput(PlotStart,PlotEnd,PlotChrm):
    SamePeakFile = codecs.open("SamePeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','rb')
    DifferentPeakFile = codecs.open("DifferentPeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','rb')
    
    running = 1 
    FileOccurance = {}
    FileOccurance['1'] = []
    FileOccurance['2'] = []
   

    DifferentPeakFile.next()
    DifferentPeakFile.next()
    DifferentPeakFile.next()

    while running == 1:
        try:
            StringList = DifferentPeakFile.next().split()
            chrName = StringList[0]
            Start = StringList[1]
            End = StringList[2]
            FoundOn = StringList[3]
            
            Middle = (int(Start) + int(End))/2
           


            if(chrName == PlotChrm) and (Middle > PlotStart) and (Middle < PlotEnd):
                NewElem = {'x': Middle, 'y': 1}
                FileOccurance[FoundOn].append(NewElem)
                
           
        except StopIteration:
            DifferentPeakFile.close()
            running = 0
            break
  
    return render_template('Analyze.html', chrNum=PlotChrm, File1=FileOccurance['1'], File2=FileOccurance['2'])   
            
def StartServer():
    
    SamePeakFile = codecs.open("SamePeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','rb')
    DifferentPeakFile = codecs.open("DifferentPeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','rb')
    
    running = 1 
    FileOccurance = {}
    FileOccurance['1'] = {'Genes':{},'Enhancer':{},'Promoter':{}}
    FileOccurance['2'] = {'Genes':{}, 'Enhancer':{}, 'Promoter':{}}

    DifferentPeakFile.next()
    DifferentPeakFile.next()
    DifferentPeakFile.next()

    while running == 1:
        try:
            StringList = DifferentPeakFile.next().split()
            chrName = StringList[0]
            FoundOn = StringList[3]
            PromoterName = StringList[4]
            GeneName = StringList[5]
            EnhancerName = StringList[6]

            if chrName not in FileOccurance[FoundOn]:
                FileOccurance[FoundOn][chrName] = {'Genes':0,'Enhancer':0,'Promoter':0}
            else:
                if PromoterName != 'None':
                    PromoterNum = FileOccurance[FoundOn][chrName]['Promoter']
                    FileOccurance[FoundOn][chrName]['Promoter'] = PromoterNum + 1
                if GeneName != 'None':
                    GeneNum = FileOccurance[FoundOn][chrName]['Genes']
                    FileOccurance[FoundOn][chrName]['Genes'] = GeneNum + 1
                if EnhancerName != 'None':
                    EnhancerNum = FileOccurance[FoundOn][chrName]['Enhancer']
                    FileOccurance[FoundOn][chrName]['Enhancer'] = EnhancerNum + 1

            if GeneName not in FileOccurance[FoundOn]['Genes']:
                FileOccurance[FoundOn]['Genes'][GeneName] = 1
            else:
                if GeneName != 'None':
                    GeneNum = FileOccurance[FoundOn]['Genes'][GeneName]
                    FileOccurance[FoundOn]['Genes'][GeneName] = GeneNum + 1
            if PromoterName not in FileOccurance[FoundOn]['Promoter']:
                FileOccurance[FoundOn]['Promoter'][PromoterName] = 1
            else:
                if PromoterName != 'None':
                    PromoterNum = FileOccurance[FoundOn]['Promoter'][PromoterName]
                    FileOccurance[FoundOn]['Promoter'][PromoterName] = PromoterNum + 1
            if EnhancerName not in FileOccurance[FoundOn]['Enhancer']:
                FileOccurance[FoundOn]['Enhancer'][EnhancerName] = 1
            else:
                if EnhancerName != 'None':
                    EnhancerNum = FileOccurance[FoundOn]['Enhancer'][EnhancerName]
                    FileOccurance[FoundOn]['Enhancer'][EnhancerName] = EnhancerNum + 1

        except StopIteration:
            DifferentPeakFile.close()
            running = 0
            break
    
    FirstFileChrList = []
    FirstFileGeneCount = []
    FirstFilePromoterCount = []
    FirstFileEnhancerCount = []

    SecondFileChrList = []
    SecondFileGeneCount = []
    SecondFilePromoterCount = []
    SecondFileEnhancerCount = []

    for SingleChrName in FileOccurance['1']:
        if SingleChrName != 'Genes' and SingleChrName != 'Promoter' and SingleChrName != 'Enhancer': 
            FirstFileChrList.append(SingleChrName)
            FirstChrInfo = FileOccurance['1'][SingleChrName]
            FirstFileGeneCount.append(FirstChrInfo['Genes'])
            FirstFilePromoterCount.append(FirstChrInfo['Promoter'])
            FirstFileEnhancerCount.append(FirstChrInfo['Enhancer'])

    for SingleChrName in FileOccurance['2']:
        if SingleChrName != 'Genes' and SingleChrName != 'Promoter' and SingleChrName != 'Enhancer':
            SecondFileChrList.append(SingleChrName)
            SecondChrInfo = FileOccurance['2'][SingleChrName]
            SecondFileGeneCount.append(SecondChrInfo['Genes'])
            SecondFilePromoterCount.append(SecondChrInfo['Promoter'])
            SecondFileEnhancerCount.append(SecondChrInfo['Enhancer'])

    FirstFileGeneList = FileOccurance['1']['Genes']
    SecondFileGeneList = FileOccurance['2']['Genes']
    FirstFilePromoterList = FileOccurance['1']['Promoter']
    SecondFilePromoterList = FileOccurance['2']['Promoter']
    FirstFileEnhancerList = FileOccurance['1']['Enhancer']
    SecondFileEnhancerList = FileOccurance['2']['Enhancer']

    #print FirstFileGeneList
    SortedFirstFileGeneList = sorted(FirstFileGeneList.iteritems(), key=lambda (k, v): (-v, k))[:20] 
    SortedSecondFileGeneList = sorted(SecondFileGeneList.iteritems(), key=lambda (k, v): (-v, k))[:20] 
    SortedFirstFilePromoterList = sorted(FirstFilePromoterList.iteritems(), key=lambda (k, v): (-v, k))[:20] 
    SortedSecondFilePromoterList = sorted(SecondFilePromoterList.iteritems(), key=lambda (k, v): (-v, k))[:20] 
    SortedFirstFileEnhancerList = sorted(FirstFileEnhancerList.iteritems(), key=lambda (k, v): (-v, k))[:20] 
    SortedSecondFileEnhancerList = sorted(SecondFileEnhancerList.iteritems(), key=lambda (k, v): (-v, k))[:20] 
    
    SortedFirstFileGeneName =  [tuple[0] for tuple in SortedFirstFileGeneList]
    SortedSecondFileGeneName =  [tuple[0] for tuple in SortedSecondFileGeneList] 
    SortedFirstFilePromoterName = [tuple[0] for tuple in SortedFirstFilePromoterList]
    SortedSecondFilePromoterName = [tuple[0] for tuple in SortedSecondFilePromoterList]
    SortedFirstFileEnhancerName =  [tuple[0] for tuple in SortedFirstFileEnhancerList]
    SortedSecondFileEnhancerName =  [tuple[0] for tuple in SortedSecondFileEnhancerList]

    SortedFirstFileGeneValue =  [tuple[1] for tuple in SortedFirstFileGeneList]
    SortedSecondFileGeneValue =  [tuple[1] for tuple in SortedSecondFileGeneList] 
    SortedFirstFilePromoterValue = [tuple[1] for tuple in SortedFirstFilePromoterList]
    SortedSecondFilePromoterValue = [tuple[1] for tuple in SortedSecondFilePromoterList]
    SortedFirstFileEnhancerValue =  [tuple[1] for tuple in SortedFirstFileEnhancerList]
    SortedSecondFileEnhancerValue =  [tuple[1] for tuple in SortedSecondFileEnhancerList]
    
    return render_template('Results.html', FirstGeneValue=SortedFirstFileGeneValue,FirstGeneKey=SortedFirstFileGeneName,
            SecondGeneValue=SortedSecondFileGeneValue, SecondGeneKey=SortedSecondFileGeneName,FirstPromoterValue=SortedFirstFilePromoterValue,
            FirstPromoterKey=SortedFirstFilePromoterName, SecondPromoterValue=SortedSecondFilePromoterValue, SecondPromoterKey=SortedSecondFilePromoterName,
            FirstEnhancerKey= SortedFirstFileEnhancerName, FirstEnhancerValue=SortedFirstFileEnhancerValue, SecondEnhancerKey=SortedSecondFileEnhancerName,
            SecondEnhancerValue=SortedSecondFileEnhancerValue, FirstChrList=FirstFileChrList, FirstGeneCount= FirstFileGeneCount, FirstPromoterCount = FirstFilePromoterCount,
            FirstEnhancerCount = FirstFileEnhancerCount, SecondChrList=SecondFileChrList, SecondGeneCount = SecondFileGeneCount, SecondPromoterCount = SecondFilePromoterCount,
            SecondEnhancerCount = SecondFileEnhancerCount 
            )

    

def StoreSameFile(InParts):
    global SamePartNumber
    global SameDict
    SameDictOutput = None
    if InParts == 1 or SamePartNumber > 1:
        if SamePartNumber == 1:
            SameDictOutput = open("SamePeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','w')
            SameDictOutput.write('chrName' +'\t' + 'Start1' + '\t' + 'End1' + '\t' + 'Start2' + '\t' + 'End2' + '\t' + 'Promoter1' + '\t' + 'Promoter2' + '\t' + 'Gene1' + '\t' + 'Gene2' + '\t' + 'Enhancer1' + '\t' + 'Enhancer2' + '\n')
        else:
            SameDictOutput = open("SamePeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','a')
        SamePartNumber = SamePartNumber + 1
        print "Saving Same Part " + str(SamePartNumber)
    else:
        SameDictOutput = open("SamePeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','w')
        SameDictOutput.write('chrName' +'\t' + 'Start1' + '\t' + 'End1' + '\t' + 'Start2' + '\t' + 'End2' + '\t' + 'Promoter1' + '\t' + 'Promoter2' + '\t' + 'Gene1' + '\t' + 'Gene2' + '\t' + 'Enhancer1' + '\t' + 'Enhancer2' + '\n')
    #counter = 1;
    for chrNumber in SameDict:
        for ListElem in SameDict[chrNumber]:
            #Put X and Y chromosome names back
            chrName = ''
            if chrNumber == '25':
                chrName = 'X'
            elif chrNumber == '26':
                chrName = 'Y'
            elif chrNumber == '27':
                chrName = 'M'
            else:
                chrName = str(chrNumber)
            #counter = counter +1
            #print counter
            if not ListElem['gene1']:
                Gene1List = 'None'
            else:
                Gene1List = ", ".join(ListElem['gene1'])

            if not ListElem['gene2']:
                Gene2List = 'None'
            else:
                Gene2List = ", ".join(ListElem['gene2'])


            if not ListElem['promoter1']:
                Promoter1List = 'None'
            else:
                Promoter1List = ", ".join(ListElem['promoter1'])

            if not ListElem['promoter2']:
                Promoter2List = 'None'
            else:
                Promoter2List = ", ".join(ListElem['promoter2'])

            if not ListElem['enhancer1']:
                Enhancer1List = 'None'
            else:
                Enhancer1List = ", ".join(ListElem['enhancer1'])

            if not ListElem['enhancer2']:
                Enhancer2List = 'None'
            else:
                Enhancer2List = ", ".join(ListElem['enhancer2'])        
         
            SameDictOutput.write('chr'+chrName +'\t' + str(ListElem['start1']) + '\t' + str(ListElem['end1']) + '\t' + str(ListElem['start2']) + '\t' + str(ListElem['end2']) + '\t' + str(Promoter1List) + '\t' + str(Promoter2List) + '\t' + str(Gene1List) + '\t' + str(Gene2List) + '\t' + str(Enhancer1List) + '\t' + str(Enhancer2List) + '\n')
            #raw_input("enter")
    SameDictOutput.close()
    SameDict = {}

def StoreDifferentFile(InParts):
    global DifferentPartNumber
    global DifferentDict
    DifferentDictOutput = None
    if InParts == 1 or DifferentPartNumber > 1:
        if DifferentPartNumber == 1:
            DifferentDictOutput = open("DifferentPeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','w')
            DifferentDictOutput.write('#1 = ' + InputFile1 + '\n')
            DifferentDictOutput.write('#2 = ' + InputFile2 + '\n')
            DifferentDictOutput.write('chrName' +'\t' + 'Start' + '\t' + 'End' + '\t' + 'FoundOn' + '\t' +'Promoter' + '\t' +'Gene' + '\t' +'Enhancer'+ '\n')            
        else:
            DifferentDictOutput = open("DifferentPeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','a')
        DifferentPartNumber = DifferentPartNumber + 1
        print "Saveing Different Part " + str(DifferentPartNumber)
    else:
        DifferentDictOutput = open("DifferentPeaks "+splitext(InputFile1)[0]+"_And_"+splitext(InputFile2)[0]+'.xls','w')
        DifferentDictOutput.write('#1 = ' + InputFile1 + '\n')
        DifferentDictOutput.write('#2 = ' + InputFile2 + '\n')
        DifferentDictOutput.write('chrName' +'\t' + 'Start' + '\t' + 'End' + '\t'+'FoundOn' + '\t' +'Promoter' + '\t' +'Gene' + '\t' +'Enhancer'+ '\n')
    for chrNumber in DifferentDict:
        for ListElem in DifferentDict[chrNumber]:
            #Put X and Y chromosome names back
            chrName = ''
            if chrNumber == '25':
                chrName = 'X'
            elif chrNumber == '26':
                chrName = 'Y'
            elif chrNumber == '27':
                chrName = 'M'
            else:
                chrName = str(chrNumber)

            if not ListElem['gene']:
                GeneList = 'None'
            else:
                GeneList = ", ".join(ListElem['gene'])

            if not ListElem['promoter']:
                PromoterList = 'None'
            else:
                PromoterList = ", ".join(ListElem['promoter'])

            if not ListElem['enhancer']:
                EnhancerList = 'None'
            else:
                EnhancerList = ", ".join(ListElem['enhancer'])
        
            DifferentDictOutput.write('chr'+chrName +'\t' + str(ListElem['start']) + '\t' + str(ListElem['end'])  + '\t' + str(ListElem['foundOn']) + '\t' + str(PromoterList) + '\t' + str(GeneList) + '\t' + str(EnhancerList) + '\n')

    DifferentDictOutput.close()
    DifferentDict = {}

def StoreDifferent(chrNum, chrStart, chrEnd, foundOn,gene,promoter,enhancer):
    global DifferentLines
    DifferentLines = DifferentLines + 1
    if str(chrNum) not in DifferentDict:
        DifferentDict[str(chrNum)] = []
    DifferentDict[str(chrNum)].append({'start': chrStart, 'end':chrEnd,'foundOn':foundOn, 'gene':gene, 'promoter':promoter, 'enhancer':enhancer})
    if DifferentLines > 2000000:
        DifferentLines = 1
        StoreDifferentFile(1)

def StoreSame(chrNum, chr1Start, chr1End, chr2Start, chr2End, gene1, gene2,promoter1,promoter2,enhancer1, enhancer2):
    global SameLines
    global SameDict
    SameLines = SameLines +1
    #print SameLines
    if str(chrNum) not in SameDict:
        SameDict[str(chrNum)] = []
    SameDict[str(chrNum)].append({'start1': chr1Start, 'end1':chr1End,'start2':chr2Start,'end2':chr2End, 'gene1':gene1,'gene2':gene2, 'promoter1':promoter1,'promoter2':promoter2,'enhancer1':enhancer1,'enhancer2':enhancer2})
    #print len(SameDict)
    if SameLines > 2000000:
        SameLines = 1
        StoreSameFile(1)
    #print sys.getsizeof(SameDict)

def FindOverlap(ChangedPeak, StringList, FilePointer, PeakStart, PeakEnd, PeakChrNum,IsEnhancer):
    PeakDone = 0
    PeakOverlap = []
    while PeakDone == 0 and ChangedPeak == 1:
        GeneStart = int(StringList[1])
        GeneEnd = int(StringList[2])
        GeneName = str(StringList[3])
        if(IsEnhancer == 1):
            GeneName = '(' + str(GeneName) + ' ' + str(StringList[4]) + ' ' + str(StringList[5]) +')' 
       


        if(StringList[0][3:] == 'X'):
            GenechrNum = 25
        elif(StringList[0][3:] == 'Y'):
            GenechrNum = 26
        elif(StringList[0][3:] == 'M'):
            GenechrNum = 27
        else:
            GenechrNum = int(StringList[0][3:])
       
            
        if (GenechrNum == PeakChrNum):
            if (GeneEnd < PeakStart):
                StringList = FilePointer.next().split()
            elif (PeakEnd<GeneStart):
                PeakDone = 1
            else:
            #print GeneStart
            #print PeakEnd
            #print PeakEnd
            #print GeneEnd
                
                if (GeneStart<= PeakStart) and (PeakStart <= GeneEnd):
                #Start of peak in gene
                    PeakOverlap.append(GeneName)
                   
                    if (PeakEnd > GeneEnd):
                        StringList = FilePointer.next().split()
                    else:
                       
                        if len(StringList) == 5:
                            
                            OverLappingGeneStr = StringList[4]
                            OverLappingGeneList = OverLappingGeneStr.split(',')
                            for OverLappingElem in OverLappingGeneList:
                                OverLappingInformation = OverLappingElem.split('|')
                                OverlapStart = int(OverLappingInformation[0])
                                OverlapEnd = int(OverLappingInformation[1])
                                OverlapName = OverLappingInformation[2]
                                if(IsEnhancer == 1):
                                    OverlapName = '(' + str(OverlapName) + ' ' + str(OverLappingInformation[3]) + ' ' + str(OverLappingInformation[4]) +')' 
                               
                                if PeakEnd < OverlapStart:
                                    break
                                elif (OverlapStart<= PeakStart) and (PeakStart <= OverlapEnd):
                                    # Start of peak on overlap
                                    PeakOverlap.append(OverlapName)
                                elif (OverlapStart<= PeakEnd) and (PeakEnd <= OverlapEnd):
                                    # End of peak on overlap
                                    PeakOverlap.append(OverlapName)
                                elif (OverlapStart >= PeakStart) and (PeakEnd >= OverlapEnd):
                                    # Peak Overlaps the Overrlapped gene
                                    PeakOverlap.append(OverlapName)
                        PeakDone = 1                    
              
                elif (GeneStart<= PeakEnd) and (PeakEnd <= GeneEnd):
                    #End of peak in gene
                    PeakOverlap.append(GeneName)
                    if len(StringList) == 5:
                        OverLappingGeneStr = StringList[4]
                        OverLappingGeneList = OverLappingGeneStr.split(',')
                        for OverLappingElem in OverLappingGeneList:
                            OverLappingInformation = OverLappingElem.split('|')
                            OverlapStart = int(OverLappingInformation[0])
                            OverlapEnd = int(OverLappingInformation[1])
                            OverlapName = OverLappingInformation[2]
                            if(IsEnhancer == 1):
                                OverlapName = '(' + str(OverlapName) + ' ' + str(OverLappingInformation[3]) + ' ' + str(OverLappingInformation[4]) +')'
                            if PeakEnd < OverlapStart:
                                break
                            elif (OverlapStart<= PeakStart) and (PeakStart <= OverlapEnd):
                                # Start of peak on overlap
                                PeakOverlap.append(OverlapName)
                            elif (OverlapStart<= PeakEnd) and (PeakEnd <= OverlapEnd):
                                # End of peak on overlap
                                PeakOverlap.append(OverlapName)
                            elif (OverlapStart >= PeakStart) and (PeakEnd >= OverlapEnd):
                                # Peak Overlaps the Overrlapped gene
                                PeakOverlap.append(OverlapName)
                    PeakDone = 1
                    
                elif (GeneStart >= PeakStart) and (PeakEnd >= GeneEnd):
                    #Peak Overlaps gene
                    PeakOverlap.append(GeneName)
                    StringList = FilePointer.next().split()
        else:

            if PeakChrNum == 2 and GenechrNum == 19:
                StringList = FilePointer.next().split()
            elif GenechrNum == 2 and PeakChrNum == 19:
                PeakDone = 1
            elif GenechrNum < PeakChrNum:
                StringList = FilePointer.next().split()
            elif GenechrNum > PeakChrNum:
                PeakDone = 1 



    return {'StringList':StringList, 'FilePointer':FilePointer,'PeakOverlap':PeakOverlap}

def run(ToleranceValue):
    global GeneFile 
    global GeneFile1
    global PeakFile1
    global PeakFile2 
    global PromoterFile 
    global PromoterFile1 
    global EnhancerFile 
    global EnhancerFile1 
    global PeakFile1Ready 
    global PeakFile2Ready 
    global DifferentDict 
    global SameDict 
    global SameLines 
    global DifferentLines 
    global SamePartNumber 
    global DifferentPartNumber 
    GeneFile.seek(0)
    GeneFile1.seek(0)
    PromoterFile.seek(0)
    PromoterFile1.seek(0)
    EnhancerFile.seek(0)
    EnhancerFile1.seek(0)
    PeakFile1Ready = 0
    PeakFile2Ready = 0
    DifferentDict = {}
    SameDict = {}
    SameLines = 1
    DifferentLines = 1
    SamePartNumber = 1
    DifferentPartNumber = 1

    PeakFile1.seek(0)
    PeakFile2.seek(0)
    
    while PeakFile1Ready == 0:
        PeakFile1Line = PeakFile1.next()
        if PeakFile1Line[0] != "#" and PeakFile1Line != "\n":
            PeakFile1Ready = 1
        #print PeakFile1.next()       

    while PeakFile2Ready == 0:
        PeakFile2Line = PeakFile2.next()
        if PeakFile2Line[0] != "#" and PeakFile2Line != "\n":
            PeakFile2Ready = 1
        #print PeakFile2.next()       

    running = 1
    tolerance = int(ToleranceValue)

    PeakFile1StringList = PeakFile1.next().split()
    PeakFile2StringList = PeakFile2.next().split()
    GeneFileStringList = GeneFile.next().split()
    GeneFile1StringList = GeneFile1.next().split()
    PromoterFileStringList = PromoterFile.next().split()
    PromoterFile1StringList = PromoterFile1.next().split()
    EnhancerFileStringList = EnhancerFile.next().split()
    EnhancerFile1StringList = EnhancerFile1.next().split()
    NewPreviousLines1 = []
    Peak1Changed  = 1
    Peak2Changed = 1
    while running == 1:
        try:
            #If chromosome is X or Y
            
            if(PeakFile1StringList[0][3:] == 'X'):
                Peak1chrNum = 25
            elif(PeakFile1StringList[0][3:] == 'Y'):
                Peak1chrNum = 26
            elif(PeakFile1StringList[0][3:] == 'M'):
                Peak1chrNum = 27
            else:
                Peak1chrNum = int(PeakFile1StringList[0][3:])

            if(PeakFile2StringList[0][3:] == 'X'):
                Peak2chrNum = 25
            elif(PeakFile2StringList[0][3:] == 'Y'):
                Peak2chrNum = 26
            elif(PeakFile2StringList[0][3:] == 'M'):
                Peak2chrNum = 27
            else:
                Peak2chrNum = int(PeakFile2StringList[0][3:])
            
            Peak1Start = int(PeakFile1StringList[1])
            Peak1End = int(PeakFile1StringList[2])
            Peak1Center = round((Peak1End + Peak1Start)/2,0)
            Peak2Start = int(PeakFile2StringList[1])
            Peak2End = int(PeakFile2StringList[2])
            Peak2Center = round((Peak2End + Peak2Start)/2,0)

            Peak1Done = 0
            Peak2Done = 0
            Peak1Gene = []
            Peak2Gene = []
            Peak1Promoter = []
            Peak2Promoter = []
            Peak1Enhancer = []
            Peak2Enhancer = []
           
            Peak1OverlapGene = FindOverlap(Peak1Changed, GeneFileStringList, GeneFile, Peak1Start, Peak1End, Peak1chrNum, 0)
            GeneFileStringList = Peak1OverlapGene['StringList']
            GeneFile = Peak1OverlapGene['FilePointer']
            Peak1Gene = Peak1OverlapGene['PeakOverlap']



            Peak1Done = 0
            Peak1OverlapPromoter = FindOverlap(Peak1Changed, PromoterFileStringList, PromoterFile, Peak1Start, Peak1End, Peak1chrNum, 0)
            PromoterFileStringList = Peak1OverlapPromoter['StringList']
            PromoterFile = Peak1OverlapPromoter['FilePointer']
            Peak1Promoter =  Peak1OverlapPromoter['PeakOverlap']
            
            Peak1OverlapEnhancer = FindOverlap(Peak1Changed, EnhancerFileStringList, EnhancerFile, Peak1Start, Peak1End, Peak1chrNum, 1)
            EnhancerFileStringList = Peak1OverlapEnhancer['StringList']
            EnhancerFile = Peak1OverlapEnhancer['FilePointer']
            Peak1Enhancer =  Peak1OverlapEnhancer['PeakOverlap']

            
            Peak2OverlapGene = FindOverlap(Peak2Changed, GeneFile1StringList, GeneFile1, Peak2Start, Peak2End, Peak2chrNum,0)
            GeneFile1StringList = Peak2OverlapGene['StringList']
            GeneFile1 = Peak2OverlapGene['FilePointer']
            Peak2Gene = Peak2OverlapGene['PeakOverlap']


           
            Peak2Done = 0
            Peak2OverlapPromoter = FindOverlap(Peak2Changed, PromoterFile1StringList, PromoterFile1, Peak2Start, Peak2End, Peak2chrNum,0)
            PromoterFile1StringList = Peak2OverlapPromoter['StringList']
            PromoterFile1 = Peak2OverlapPromoter['FilePointer']
            Peak2Promoter = Peak2OverlapPromoter['PeakOverlap']

            Peak2OverlapEnhancer = FindOverlap(Peak2Changed, EnhancerFile1StringList, EnhancerFile1, Peak2Start, Peak2End, Peak2chrNum, 1)
            EnhancerFile1StringList = Peak2OverlapEnhancer['StringList']
            EnhancerFile1 = Peak2OverlapEnhancer['FilePointer']
            Peak2Enhancer =  Peak2OverlapEnhancer['PeakOverlap']


            Peak1Changed = 0
            Peak2Changed = 0     

            if Peak1chrNum == Peak2chrNum:
                PeakRegionDifference = Peak1Center - Peak2Center
                if abs(PeakRegionDifference) > tolerance:
                    if(PeakRegionDifference < 0):
                        StoreDifferent(Peak1chrNum,Peak1Start,Peak1End,1,Peak1Gene,Peak1Promoter,Peak1Enhancer)
                        PeakFile1StringList = PeakFile1.next().split()
                        Peak1Changed = 1
                    elif(PeakRegionDifference > 0):
                        StoreDifferent(Peak2chrNum,Peak2Start,Peak2End,2,Peak2Gene, Peak2Promoter, Peak2Enhancer)
                        PeakFile2StringList = PeakFile2.next().split()
                        Peak2Changed = 1
                elif abs(PeakRegionDifference) <= tolerance:
                    StoreSame(Peak1chrNum,Peak1Start,Peak1End,Peak2Start,Peak2End,Peak1Gene,Peak2Gene, Peak1Promoter, Peak2Promoter, Peak1Enhancer, Peak2Enhancer)
                    PeakFile1StringList = PeakFile1.next().split()
                    PeakFile2StringList = PeakFile2.next().split()
                    Peak1Changed = 1
                    Peak2Changed = 1
            else:
                if Peak1chrNum == 2 and Peak2chrNum == 19:
                    StoreDifferent(Peak2chrNum,Peak2Start,Peak2End,2,Peak2Gene, Peak2Promoter, Peak2Enhancer)
                    PeakFile2StringList = PeakFile2.next().split()
                    Peak2Changed = 1
                elif Peak2chrNum == 2 and Peak1chrNum == 19:
                    StoreDifferent(Peak1chrNum,Peak1Start,Peak1End,1,Peak1Gene, Peak1Promoter, Peak1Enhancer)
                    PeakFile1StringList = PeakFile1.next().split()
                    Peak1Changed = 1
                elif Peak1chrNum < Peak2chrNum:
                    StoreDifferent(Peak1chrNum,Peak1Start,Peak1End,1,Peak1Gene, Peak1Promoter, Peak1Enhancer)
                    PeakFile1StringList = PeakFile1.next().split()
                    Peak1Changed = 1
                elif Peak1chrNum > Peak2chrNum:
                    StoreDifferent(Peak2chrNum,Peak2Start,Peak2End,2,Peak2Gene, Peak2Promoter, Peak2Enhancer)
                    PeakFile2StringList = PeakFile2.next().split()
                    Peak2Changed = 1 
        except StopIteration:
            PeakFile1.close()
            PeakFile2.close()
            running = 0
            break

    StoreSameFile(0)
    StoreDifferentFile(0)
    return StartServer()



def allowed_file(filename):
    ALLOWED_EXTENSIONS = set(['txt', 'bed', 'xls', 'xlsm', 'xlsx', 'csv'])
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

def StartWebsite():
    
    app = Flask(__name__)
    UPLOAD_FOLDER = "upload"
    app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
    @app.route('/FindDiff',methods=['GET'])
    def FindDiff():
        return render_template("DiffIndex.html")
    @app.route('/FindDiff',methods=['POST'])
    def upload():
        global PeakFile1
        global PeakFile2
        global InputFile1
        global InputFile2

        file = request.files['file1']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(join(app.config['UPLOAD_FOLDER'], filename))
            InputFile1 = filename
            PeakFile1 = codecs.open(join(app.config['UPLOAD_FOLDER'], filename),'rb')
            
        file = request.files['file2']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(join(app.config['UPLOAD_FOLDER'], filename))
            InputFile2 = filename
            PeakFile2 = codecs.open(join(app.config['UPLOAD_FOLDER'], filename),'rb')
        ToleranceValue = int(request.form['Tolerance'])
        return run(ToleranceValue)
    @app.route('/',methods=['GET'])
    def Index():
        return render_template('Index.html')
    @app.route('/Analyze',methods=['GET'])
    def AnalyzeIndex():
        return render_template('AnalyzeIndex.html')
    @app.route('/Analyze',methods=['POST'])
    def Analyze():
        global PeakFile1
        global PeakFile2
        global InputFile1
        global InputFile2
        file = request.files['file1']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(join(app.config['UPLOAD_FOLDER'], filename))
            InputFile1 = filename
            PeakFile1 = codecs.open(join(app.config['UPLOAD_FOLDER'], filename),'rb')
            
        file = request.files['file2']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(join(app.config['UPLOAD_FOLDER'], filename))
            InputFile2 = filename
            PeakFile2 = codecs.open(join(app.config['UPLOAD_FOLDER'], filename),'rb')
        ToleranceValue = int(request.form['Tolerance'])
        run(ToleranceValue)
        ChrNum = 'chr'+request.form['ChrNum']
        Start = int(request.form['Start'])
        End = int(request.form['End'])
        return PlotOutput(Start,End,ChrNum)

    #app.debug = True
    app.run()



if( sys.argv[1] == '--website'):
    PeakFile1 = None 
    PeakFile2 = None
    StartWebsite()    
else:
    InputFile1 = sys.argv[1]
    InputFile2 = sys.argv[2]
    PeakFile1 = codecs.open(InputFile1,'rb') 
    PeakFile2 = codecs.open(InputFile2,'rb')
    run(100)



