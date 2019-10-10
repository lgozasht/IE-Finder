#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 17:43:55 2019

@author: lgozasht
"""
import os


class csvReader():


    def __init__(self, dataFile = ''):
        """
        tialize blasthit as the blast output file.

        """
        self.dataFile = dataFile

    def csv(self):
            """
            Return results, a dictionary containing data extracted from blasthit for each query.


            Edit this to remove alignments that don't align for 80% of sequence length ie. if queryLength/alignmentLength > .8:......
            """
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
    
                                assemblyDic = {'Species': sp[0],'GenBank' : Assembly, 'refFasta' : fasta, 'refAnnotation' : annotation}
                                assemblies.append(assemblyDic)
                            except IndexError:
                                pass

            return assemblies



def main():
    myData = csvReader('eukaryotes2.csv')
    genomeData = myData.csv()
    """STEP 1"""
    for assembly in genomeData: 
        try:
            NAME = assembly['GenBank']
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
            print(NAME)
            os.system('perl ./exint-gb.pl Data_{0}/{1} > ./Data_{0}/my{0}.exons-introns'.format(NAME, annotation))
            os.system('rm -r Data_{0}/GCA*'.format(NAME))

        except KeyError:
            print('No Genbank Assembly')
            
            
if __name__ == "__main__":
    main() 


