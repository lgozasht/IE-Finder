#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:27:09 2019

@author: lgozasht
"""

import os
from pathlib import Path
from sequenceAnalyzer import FastAreader
import copy

class Compare():
    
    """
    Assess blast results from ProteinBlast2 and remove any found paralogs introner element data set.
    
    Return a dictionary containing all introner elements that pass the paralog filter.
    
    
    """
    
    def __init__(self, ieFile, deDuped, annotation, NAME):
        
        self.ieFile = ieFile
        self.deDuped = deDuped
        self.NAME = NAME
       # self.annotation = 'Data_{0}/{1}'.format(self.NAME, annotation)

        self.geneDupList = []
        self.famDic = {}
        self.finalDic = {}
    def deDupedReader(self):
        print('reading protein file ...')

        
        with open(self.deDuped,'r') as f:
            for line in f:
                line = line.rstrip()
                sp = line.split('\t')
                reference = sp[0]
                if reference not in self.geneDupList:
                    self.geneDupList.append(reference)
    def ieReader(self):
        print('reading IE file ...')

        myReaderIE = FastAreader('Data_{0}/finalIEs.fa'.format(self.NAME))
        for header, sequence in myReaderIE.readFasta():
            fam = header[0:2:1]

            if fam[1].isdigit() == False:
                fam = header[0]
                head = header[1:len(header):1]


            else:
                fam = header[0:2:1]
                head = header[2:len(header):1]
            
                
                
                
            #ieDic = {head : sequence}
            if fam not in self.famDic.keys():
                self.famDic[fam] = {head : sequence}
            else:
                self.famDic[fam][head] = sequence
    """
    check this
    """
            
    
    def compare(self):
        print('Removing putative IEs resulting from gene dups ...')
        finalDic = copy.deepcopy(self.famDic)
        for fam, ieDic in self.famDic.items():
            for head, sequence in ieDic.items():
                for reference in self.geneDupList:
                    if reference.split('_')[0] == fam:
                        if head in reference:
                            del finalDic[fam][head]
                            
        return finalDic
                                    
                
                    


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
                    #print(len(sp))
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
                            refPreAssembly = refLink.split("/")[-1]

                            refAssembly = refPreAssembly.split("_")[0] + "_" + refPreAssembly.split("_")[1]
                            refFasta = refLink + '/*_genomic.fna.gz'
                            refAnnotation = refLink + '/*_genomic.gff.gz'
                            genPreAssembly = genLink.split("/")[-1]

                            genAssembly = genPreAssembly.split("_")[0] + "_" + genPreAssembly.split("_")[1]
                           # print(Assembly)
                            genFasta = genLink + '/*_genomic.fna.gz'
                            genAnnotation = genLink + '/*_genomic.gff.gz'

                            assemblyDic = {'Species': sp[0],'RefSeq' : refAssembly, 'refFasta' :refFasta , 'refAnnotation' : refAnnotation, 'GenBank' : genAssembly, 'genFasta' :genFasta , 'genAnnotation' : genAnnotation}
                            assemblies.append(assemblyDic)
                        except IndexError:
                        
                           # print("Fuck you")
                            link = sp[14]
                            link = link.strip('\"')
                            link = link.strip('\"')
                            try:
                                preAssembly = link.split("/")[-1]
    
                                Assembly = preAssembly.split("_")[0] + "_" + preAssembly.split("_")[1]
                               
                                fasta = link + '/*_genomic.fna.gz'
                                annotation = link + '/*_genomic.gff.gz'
    
                                assemblyDic = {'Species': sp[0],'GenBank' : Assembly, 'genFasta' : fasta, 'genAnnotation' : annotation}
                                assemblies.append(assemblyDic)
                            except IndexError:
                                pass

            return assemblies


def main():
    
    
    myData = csvReader('test.csv')
    genomeData = myData.csv()

    for assembly in genomeData:
        try:
            NAME = assembly['GenBank']
            print(NAME)
            ieFile = Path("Data_{0}/finalIEs.fa".format(NAME))
    
    #        print("please")
            if ieFile.is_file():
                os.system('chmod 775 Data_{0}/total_all-vs-all_deduped.tsv'.format(NAME))
    
                print('Processing {0}'.format(NAME))
             #   os.system('rm -r GCA*')
             #   os.system('rm -r GCF*')
                os.system('wget {0}'.format(assembly['genAnnotation']))
                print(assembly['genAnnotation'])
                os.system('wget {0}'.format(assembly['genFasta']))
                print(assembly['genFasta'])

                os.system('gunzip {0}*'.format(NAME))
                os.system('cp {0}* ./Data_{0}'.format(NAME))
               # os.system('gunzip ./Data_{0}/*'.format(NAME))
               # os.system('rm -r {0}*'.format(NAME))
                fastaList = assembly['genFasta'].split("/")
                fastaGz = fastaList[-2]
                fasta = fastaGz + '_genomic.fna'
                annotationList = assembly['genAnnotation'].split("/")
                annotationGz = annotationList[-2]
                annotation = annotationGz + '_genomic.gff'
                print(annotation)
                print(fasta)
                print(NAME)
                noDupFile = 'Data_{0}/gene_duplications.tsv'.format(NAME)
                DeDuped = Compare(ieFile, noDupFile, annotation, NAME)
                DeDuped.deDupedReader()
                DeDuped.ieReader()
                finalDic = DeDuped.compare()
                
                with open('Data_{0}/postProteinAlignmentIEs.fa'.format(NAME), 'w') as f:
                    for fam, ieDic in finalDic.items():
                        for head, seq in ieDic.items():
                            
                            f.write('>{0}_{1}\n'.format(fam, head))
                            f.write('{0}\n'.format(seq))
                
            else:
                print("no IEs")
                continue
        except KeyError:
            print('no genbank assembly')
            
            
if __name__ == "__main__":
    main() 
