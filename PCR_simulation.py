# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import random as rd
import math

FilePath = "PROVEDIt_1-5-Person CSVs Filtered-2/PROVEDIt_1-5-Person CSVs Filtered_3500_GF29cycles/PROVEDIt_RD14-0003 GF Known Genotypes.csv"
PROVEDItFilePath = "PROVEDIt_1-5-Person CSVs Filtered-2/PROVEDIt_1-5-Person CSVs Filtered_3500_GF29cycles/1-Person/15 sec/RD14-0003_GF_15sec_GM_SE33F_1P.csv"

rd.seed(1)
KnownGenotypes = []
with open(FilePath, newline = "") as KnownGenotypesFile:
    Reader = csv.reader(KnownGenotypesFile)
    for Row in Reader: 
        KnownGenotypes.append(Row)
Genes = KnownGenotypes[0][2:]
KnownGenotypesData = [x[1:] for x in KnownGenotypes[1:]]
GenesBlackList = ["Yindel","AMEL","DYS391"]
Genes = [x for x in Genes if x not in GenesBlackList]

NoiseProbability = 0.1
GaussianMean = 1
GaussianVariance = 0.5
SuccesfulReplicationProbability = 0.95
StutterProbability = 0.005
CE_Sensitivity = 50000000
Cycles = 29
#StartingCopies = 1
InterestId = 1
InterestGeneId = 0
Step = 0.5
ProveditId = 0
Attempts = 2
Capillarity = None

if len(sys.argv) > 1:
    InterestId = int(sys.argv[1])
    if len(sys.argv) > 2:
        InterestGeneId = int(sys.argv[2])
    else:
        print("Select Gene:")
        for i in range(len(Genes)):
            print(i,Genes[i])
        InterestGeneId = int(input("Gene id:"))
    if len(sys.argv) > 3:
        ProveditId = int(sys.argv[3])
    if len(sys.argv) > 4:
        Capillarity = sys.argv[4]

print("Selected Sample Id",InterestId)        
print("Selected Interest Gene",Genes[InterestGeneId])

class Person():
    def __init__(self,KnownGenotype):
        self.MyGenotype = {k:v.split(",") for (k,v) in KnownGenotype.items()}
        self.SimulatedElectropherogram = {}

def PCR(Cycles, Alleles, StartingCopies):
    DnaSample = {int(Allele):[ StartingCopies # ...===... (full lenght original dna)
                         ,0 # ...===____ and ---===...
                         ,0 # ...---=== and ===___...
                         ,0 # === (the one looking for)
                         ]for Allele in Alleles}
    for Cycle in range(Cycles):
        #print(Cycle)
        AlleleList = list(DnaSample.keys())
        for Allele in AlleleList:
            #Replication
            DnaSample[Allele][3] = ReplicationAttempt(2 * DnaSample[Allele][3] + DnaSample[Allele][2])
            DnaSample[Allele][2] = ReplicationAttempt(DnaSample[Allele][2] + DnaSample[Allele][1])
            DnaSample[Allele][1] = ReplicationAttempt(2 * DnaSample[Allele][0])
            #Stutter test
        DnaSample = Stutter(DnaSample)
        #print(DnaSample)
        #input("input")
    return DnaSample

def ReplicationAttempt(ExpectedOutput):
    SuccessfulReplicates = np.random.binomial(ExpectedOutput,SuccesfulReplicationProbability)
    #print("ExpectedOutput",ExpectedOutput,"SuccessfulReplicates",SuccessfulReplicates)   
    return SuccessfulReplicates

def Stutter(DnaSample):
    AlleleList = list(DnaSample.keys())
    for Allele in AlleleList:
        for i in range(len(DnaSample[Allele])-1):
            DNAtype = i+1
            Stutter = np.random.binomial(DnaSample[Allele][DNAtype],StutterProbability)
            DnaSample[Allele][DNAtype] -= Stutter
            #print("Allele",Allele,"DNAtype",DNAtype,"Stutter",Stutter)
            if Allele - 1 > 0 and Stutter > 0:
                if Allele - 1 not in DnaSample:
                    DnaSample[Allele - 1] = [0,0,0,0]
                DnaSample[Allele - 1][DNAtype] += Stutter
    return DnaSample

def AddNoise(Yvalues):
    Lenght = len(Yvalues)
    Noise = np.random.normal(GaussianMean,math.sqrt(GaussianVariance),Lenght)
    #print(Noise)
    SuccessfulNoise = np.random.binomial(1,NoiseProbability,Lenght)
    #print(SuccessfulNoise)
    OutputValues = [Yvalues[i] + (Noise[i] * SuccessfulNoise[i]) if i%2 == 0 else Yvalues[i] for i in range(Lenght)] 
    return OutputValues

def Plot1(x,Title,yset,Legendset,Xlabel,Ylabel):
    ax = plt.subplot(111)
    #for y in yset:
        #ax.plot(x,y)
    ax.plot(x,yset[1],label=Legendset[1],linewidth=3) 
    ax.plot(x,yset[0],label=Legendset[0],linewidth=1.5)    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax.grid()    
    ax.legend(loc = "upper center", bbox_to_anchor = (0.5, -0.115), fancybox = True)
    ax.set_ylabel(Ylabel)
    ax.set_xlabel(Xlabel)
    plt.title(Title)
    plt.show()

def QualitativeComparison(Values1,Values2):
    Lenght = len(Values1)
    SEvalues = [(Values1[i] - Values2[i]) ** 2 if i%2 == 0 else 0 for i in range(Lenght)]
    MSE = sum(SEvalues)/Lenght
    return MSE
    
ProveditData = []
with open(PROVEDItFilePath, newline = "") as ProveditFile:
    Reader = csv.reader(ProveditFile)
    IgnoreRow = True
    for Row in Reader:
        if IgnoreRow:
            Headers = Row
            Indexes = [i for i,x in enumerate(Headers) if x.split(" ")[0] in ["Allele","Height"]]
            FilteredHeaders = [Row[x] for x in Indexes]
            FilteredHeaders.insert(0,"RowDNACapillarity")
            FilteredHeaders.insert(0,"RowMass")
            FilteredHeaders.insert(0,"RowQuality")
            IgnoreRow = False
        else:
            FirstColumnAnalysis = Row[0].split("-")
            RowSampleId = FirstColumnAnalysis[2].split("d")[0]
            RowGene = Row[1]
            #RowQuality = float(FirstColumnAnalysis[4].split("_")[0].replace("Q",""))
            if int(RowSampleId) == InterestId and RowGene == Genes[InterestGeneId]:
                
                RowQuality = float(FirstColumnAnalysis[4].split("_")[0].replace("Q",""))
                if RowQuality >= 0.1:
                    RowDNACapillarity = FirstColumnAnalysis[2][-1]
                    if Capillarity is None or RowDNACapillarity == Capillarity:
                        #print(Row[0],Row[1])
                        print(len(ProveditData)," ",FirstColumnAnalysis)
                        RowMass = float(FirstColumnAnalysis[3].replace("GF",""))
                        
                        FilteredRow = [Row[x] for x in Indexes]
                        FilteredRow.insert(0,RowDNACapillarity)
                        FilteredRow.insert(0,RowMass)
                        FilteredRow.insert(0,RowQuality)
                        ProveditData.append(FilteredRow)

Samples = {}
for Row in KnownGenotypesData:
    SampleId = Row[0]
    Genotype = {}
    for GeneId in range(len(Genes)):
        Genotype[Genes[GeneId]] = Row[GeneId + 1]   
    Samples[int(SampleId)] = Person(Genotype)
    Genotype.clear()

#print(Samples)
KnownAlleles = Samples[InterestId].MyGenotype[Genes[InterestGeneId]]
print("KnownAlleles",KnownAlleles)

print("ProveditData[ProveditId]",ProveditData[ProveditId])
#PotentialAlleleListoriginal = [x for i,x in enumerate(ProveditData[ProveditId]) if i%2!=0 and i>1 and x not in ["","OL"]]
PotentialAlleleIndex = [i for i,x in enumerate(ProveditData[ProveditId]) if i%2!=0 and i>1 and x not in ["","OL"]]
PotentialAlleleList = [int(float(x)) for i,x in enumerate(ProveditData[ProveditId]) if i%2!=0 and i>1 and x not in ["","OL"]]
#PotentialAlleleValues = [float(x) for i,x in enumerate(ProveditData[ProveditId]) if ProveditData[ProveditId][i-1] in PotentialAlleleListoriginal]
PotentialAlleleValues = [float(ProveditData[ProveditId][i+1]) for i in PotentialAlleleIndex]
ProveditCE = {PotentialAlleleList[i]:PotentialAlleleValues[i] for i in range(len(PotentialAlleleValues))}
print("ProveditCE",ProveditCE)
CurrentBestError = float("inf")
BestStartingCopiues = 0
BestAlleles = []
BestProveditHeight = []
BestSimulatedHeight = []
StartingCopiesHistory = []
MSEhistory = []
for StartingCopies in range(1,40):
    for _ in range(Attempts):
        DnaOutput = PCR(Cycles,KnownAlleles,StartingCopies)
        SimulatedCE = {k:(x/CE_Sensitivity) for (k,[_,_,_,x]) in DnaOutput.items()}
        #print(SimulatedCE)
        MinAllele = min(min(ProveditCE),min(SimulatedCE))
        MaxAllele = max(max(ProveditCE),max(SimulatedCE))
        AlleleList = np.arange(MinAllele,MaxAllele,Step)
        SimulatedHeight = [SimulatedCE[x] if x in SimulatedCE else 0 for x in AlleleList]
        #print("SimulatedHeight",SimulatedHeight)
        SimulatedHeight = AddNoise(SimulatedHeight)
        #print("SimulatedHeight",SimulatedHeight)
        ProveditHeight = [ProveditCE[x] if x in ProveditCE else 0 for x in AlleleList]
        ProveditLegend = "Q" + str(ProveditData[ProveditId][0]) + " Mass" + str(ProveditData[ProveditId][1]) + " Cond " + str(ProveditData[ProveditId][2])
        #x, y = Electropherogram(DnaOutput)
        Error = QualitativeComparison(SimulatedHeight,ProveditHeight)
        if Error < CurrentBestError:
            BestStartingCopies = StartingCopies
            CurrentBestError = Error
            BestAlleles = [x for x in AlleleList]
            BestProveditHeight = [x for x in ProveditHeight]
            BestSimulatedHeight = [x for x in SimulatedHeight]
    StartingCopiesHistory.append(StartingCopies)
    MSEhistory.append(Error)

    
Error = QualitativeComparison(BestSimulatedHeight,BestProveditHeight)
print("Best MSE:",Error)
print("Best copies:",BestStartingCopies)
TitleLegend = "Sample Id " + str(InterestId) + ", Gene " + str(Genes[InterestGeneId]) + ", Known Alleles " + str( KnownAlleles) + "\nMSE: {:.2f}".format(CurrentBestError) 
SimulatedLegend = "Simulated EP. Cycles" + str(Cycles) + " StartingCopies" + str(BestStartingCopies)
ProveditLegend = "PROVEDIt EP. Q" + str(ProveditData[ProveditId][0]) + " Mass" + str(ProveditData[ProveditId][1]) + " Cap " + str(ProveditData[ProveditId][2])
Plot1(BestAlleles,TitleLegend,[BestSimulatedHeight,BestProveditHeight],[SimulatedLegend,ProveditLegend],"Alleles","Fluorescence")
#Plot1(StartingCopiesHistory,"MSE given distinct starting allele copies",[MSEhistory],[""],"Starting copies","MSE")   
#himozy
#mass