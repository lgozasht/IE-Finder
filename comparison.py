#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 13:22:16 2019

@author: lgozasht
"""
import os
from pathlib import Path
from sequenceAnalyzer import FastAreader
from Bio.Seq import Seq
import random
import copy



class RetrieveIntroners():
    """
    retrieve putative introners for each species 
    """
    
    def __init__(self, intronerFile, exintFile):
        
        self.intronerFile = intronerFile
        self.ieList = []
        self.exintFile = exintFile
        self.finalList = []
        
    def retrieveIEs(self):
        myReaderIE = FastAreader(self.intronerFile)
        for header, seq in myReaderIE.readFasta():
            headList = header.split('_')
            #print(headList)
            fam = headList[0]

            #print(fam)
            #print(scaff)
            scaff = headList[1]
            #scaff = head[1:len(head):1]
            #print(scaff)
            coord = headList[2]
            coordList = coord.split("-")
            start = int(coordList[0])
            #print(start)
            stop = int(coordList[1])
            #print(stop)
            ieDic = {'fam':fam, 'scaff':scaff,'start':start,'stop':stop, 'seq': seq}
            self.ieList.append(ieDic)
    
    def retrieveAll(self):
        """
        extract from scott's script output file
        """
        allList = []
        myReaderAll = FastAreader(self.exintFile)
        for header, seq in myReaderAll.readFasta():
            #seq.upper()
            data = header + ';' + seq
            allList.append(data)
            
        for data in allList:
           # print(data)
            intron = data.split(';')[-1]
         #   print(data)
            header = data.split(';')[0]
           # print(header)
            for ieDic in self.ieList:
                ieSeq = ieDic['seq']
                #intronList = intron.split() 
                #intronUpper = intronList.join('').upper()
                intronUpper = copy.deepcopy(intron)
                intronUpper = intronUpper.upper()
                if ieSeq in intronUpper:
                   # print(intron)
                    del ieDic['seq']
                    
                    ieDic['seq'] = intron
                    ieDic['header'] = header
                    self.finalList.append(ieDic)
                    #break
                sequence = Seq(ieDic['seq']) 
                revComp = sequence.reverse_complement()
                strrevComp = str(revComp)
                if strrevComp in intronUpper:
                    
                    #print(intron)
                    del ieDic['seq']
                    ieDic['header'] = header

                    ieDic['seq'] = intron
                    self.finalList.append(ieDic)
                    #break
        #print(self.finalList)

        return self.finalList

               

    
    

        
class csvReader():
    """
    Read NCBI eukaryotes.csv file and extract the assembky reference, fasta link, and annotation link for each species in the eukaryotes file.
            
    Return a dictionary containing the assembky reference, fasta link, and annotation link for each species in the eukaryotes file.
    """


    def __init__(self, dataFile = ''):

        self.dataFile = dataFile

    def csv(self):
            assemblies = []
            with open(self.dataFile,'r') as f:
                for line in f:
                    line = line.rstrip()
                    sp = line.split(',')
                    if '#' in sp[0]:
                        continue
                    
                   


                    elif len(sp) == 16:
                        #print('what')
                        refLink = sp[15]
                        refLink = refLink.strip('\"')
                        refLink = refLink.strip('\"')
                        genLink = sp[14]
                        genLink = genLink.strip('\"')
                        genLink = genLink.strip('\"')

                        try:
                            #print(link)
                            refPreAssembly = refLink.split("/")[-1]
                           # print(preAssembly)
                           # preList = preAssembly.split("_")
                           # print(preList[0])
                           # print(preList[1])

                            refAssembly = refPreAssembly.split("_")[0] + "_" + refPreAssembly.split("_")[1]
                           # print(Assembly)
                            refFasta = refLink + '/*_genomic.fna.gz'
                            refAnnotation = refLink + '/*_genomic.gbff.gz'
                            genPreAssembly = genLink.split("/")[-1]
                           # print(preAssembly)
                           # preList = preAssembly.split("_")
                           # print(preList[0])
                           # print(preList[1])

                            genAssembly = genPreAssembly.split("_")[0] + "_" + genPreAssembly.split("_")[1]
                           # print(Assembly)
                            genFasta = genLink + '/*_genomic.fna.gz'
                            genAnnotation = genLink + '/*_genomic.gbff.gz'

                            assemblyDic = {'Species': sp[0],'RefSeq' : refAssembly, 'refFasta' :refFasta , 'refAnnotation' : refAnnotation, 'GenBank' : genAssembly, 'genFasta' :genFasta , 'genAnnotation' : genAnnotation}
                            assemblies.append(assemblyDic)
                        except IndexError:
                        
                            link = sp[14]
                            link = link.strip('\"')
                            link = link.strip('\"')
                            try:
                                #print(link)
                                preAssembly = link.split("/")[-1]
                               # print(preAssembly)
                               # preList = preAssembly.split("_")
                               # print(preList[0])
                               # print(preList[1])
    
                                Assembly = preAssembly.split("_")[0] + "_" + preAssembly.split("_")[1]
                               # print(Assembly)
                               
                                fasta = link + '/*_genomic.fna.gz'
                                annotation = link + '/*_genomic.gbff.gz'
    
                                assemblyDic = {'Species': sp[0],'GenBank' : Assembly, 'genFasta' : fasta, 'genAnnotation' : annotation}
                                assemblies.append(assemblyDic)
                            except IndexError:
                                pass

            return assemblies
            
class Compare():
    
    """
    Run scott's script to compare intron positions for all putative introners
    """
    
    def __init__(self, eukFile, NAME):
        self.eukFile = eukFile
        self.NAME = NAME
        self.dataDic = {}
        self.relevantClade = ''
        self.relatedSpecies = []

    def dataGen(self):
        with open(self.eukFile,'r') as f:
            for line in f:
                line = line.rstrip()
                sp = line.split(',')
                if '#' in sp[0]:
                    continue
                else:
                    try:
                        clade =  sp[1].split(';')[-1]
                            
                        genLink = sp[14]
                        genLink = genLink.strip('\"')
                        genLink = genLink.strip('\"')
    
                        refPreAssembly = genLink.split("/")[-1]
                        refAssembly = refPreAssembly.split("_")[0] + "_" + refPreAssembly.split("_")[1]
                        if self.NAME == refAssembly:
                            self.relevantClade = clade
                        if clade not in self.dataDic.keys():
                            self.dataDic[clade] = [refAssembly]
                        else:
                            self.dataDic[clade].append(refAssembly)
                    except IndexError:
                        pass
    def randomSelect(self):
        relevantList = self.dataDic[self.relevantClade]
        i = 0
        while i < 5:
            randomSpecies = random.choice(relevantList)
            config = Path('../ieFind*/Data_{0}/my{0}.exons-introns'.format(randomSpecies))
            if config.is_file():
                self.relatedSpecies.append(random.choice(relevantList))
                i += 1
        
    def runScottScript(self):
        for species in self.relatedSpecies:
            os.system('cp ../ieFind*/Data_{0}/my{0}.exons-introns .'.format(species))
        
        
        os.system('perl ./HOMOLOGOUS_INTRON_FILTER.PL -filter Data_{0}/putativeIEs.exons-introns Data_{0}/my{1}.exons-introns Data_{0}/my{2}.exons-introns Data_{0}/my{3}.exons-introns Data_{0}/my{4}.exons-introns Data_{0}/my{5}.exons-introns > Data_{0}/HOMO.out'.format(self.NAME, self.relatedSpecies[0], self.relatedSpecies[1], self.relatedSpecies[2], self.relatedSpecies[3], self.relatedSpecies[4]))
        
        
        
        
        
        

def main():
    myData = csvReader('test.csv')
    genomeData = myData.csv()
    """STEP 1"""
    for assembly in genomeData:
        try:
            NAME = assembly['GenBank']
            print(NAME)
            ieFile = Path("Data_{0}/postProteinAlignmentIEs.fa".format(NAME))
    
    #        print("please")
            if ieFile.is_file():
                
    
                print('Processing {0}'.format(NAME))
    
                os.system('wget {0}'.format(assembly['genAnnotation']))
              #  os.system('wget {0}'.format(assembly['genFasta']))

                os.system('gunzip {0}*'.format(NAME))
                os.system('cp {0}* ./Data_{0}'.format(NAME))
                os.system('gunzip ./Data_{0}/*'.format(NAME))
                os.system('rm -r {0}*'.format(NAME))
#                fastaList = assembly['Fasta'].split("/")
#                fastaGz = fastaList[-2]
#                fasta = fastaGz + '_genomic.fna'
                annotationList = assembly['genAnnotation'].split("/")
                annotationGz = annotationList[-2]
                annotation = annotationGz + '_genomic.gbff'
                print(annotation)
#                print(fasta)
                print(NAME)
                """
                Make exint files for each species
                """
                os.system('perl ./exint-gb.pl Data_{0}/{1} > ./Data_{0}/my{0}.exons-introns'.format(NAME, annotation))
    
    
                myPutativeIntroners = RetrieveIntroners("Data_{0}/postProteinAlignmentIEs.fa".format(NAME),'./Data_{0}/my{0}.exons-introns'.format(NAME))
                myPutativeIntroners.retrieveIEs()
                finalList = myPutativeIntroners.retrieveAll()
                print(finalList)
                with open('Data_{0}/putativeIEs.exons-introns'.format(NAME), 'w') as exintFile:
                    for ieDic in finalList:
                        header = ieDic['header']
                        seq = ieDic['seq']
                        exintFile.write('>{0}\n'.format(header))
                        exintFile.write('{0}\n'.format(seq))

                        
                
                
                """
                compare introns (Scott's script)
                
                What about the fasta. Does it just need to be in the same direction? or what?
                """
                os.system('rm -r Data_{0}/GCA*'.format(NAME))
                os.system('rm -r Data_{0}/GCF*'.format(NAME))
                os.system('cp my* Data_{0}'.format(NAME))
                os.system('cp putative* Data_{0}'.format(NAME))

                #myComparison = Compare(ieList)
                
                
                
    
            else:
                print("no IEs")
                continue
        except KeyError:
            print('No Genbank Assembly')
            
            
            
    """Step 2"""
    
    for assembly in genomeData:
        try:
            NAME = assembly['GenBank']
            print(NAME)
            ieFile = Path("Data_{0}/postProteinAlignmentIEs.fa".format(NAME))
            if ieFile.is_file():

                myComparison = Compare('test.csv', NAME)
                myComparison.dataGen()
                myComparison.randomSelect()
                myComparison.runScottScript()
                os.system('cp my* Data_{0}'.format(NAME))
                os.system('cp putative* Data_{0}'.format(NAME))
                os.system('rm -r my*')
                os.system('rm -r putative*')

            else:
                print("no IEs")
                continue
        except KeyError:
            print('No Genbank Assembly')
         
    
        
    
    
    
    
if __name__ == "__main__":
    main() 