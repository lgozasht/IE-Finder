#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 16:18:19 2019

@author: lgozasht
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 09:40:21 2019

@author: lgozasht
"""

"""
compare proteins of putative introner containing genes


"""
import os
from pathlib import Path
from sequenceAnalyzer import FastAreader
#from Bio.Seq import Seq


class retrieveProteins():
    
    def __init__(self, ieFile, annotation, fasta, NAME):
        self.ieFile = ieFile
        self.annotation = annotation
        self.fasta = fasta
        self.locusList = []
        self.dataList = []
        self.proteinDic = {}
        self.geneDic = {}
        self.NAME = NAME
        self.locusDic = {}
        self.ieHeadDic = {}
        self.proteinList = []
    def readIeFile(self):
        myReaderIE = FastAreader(self.ieFile)
        for header, seq in myReaderIE.readFasta():
            fam = header[0:2:1]

            if fam[1].isdigit() == False:
                head = header[1:len(header):1]
                fam = header[0]

            else:
                fam = header[0:2:1]

                head = header[2:len(header):1]
            if fam not in self.ieHeadDic.keys():
                print("NEW FAM")

                self.ieHeadDic[fam] = [head]
            else:
                print('Append')
                self.ieHeadDic[fam].append(head)
        return self.ieHeadDic
    
    
    def readGFF(self,ieHeadList):

        with open('Data_{0}/{1}'.format(self.NAME,self.annotation),'r') as f:
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)

            for line in f:
                if '#' in line:
                    continue
                else:

                    line = line.rstrip()
                    sp = line.split('\t')
                    #print(sp)

                    if sp[2] == 'CDS':
                        self.dataList.append(sp)
                        for head in ieHeadList:
                            
                            start = head.split('_')[1].split('-')[0]
                            #print(start)
                            stop = str(int(head.split('_')[1].split('-')[1]) + 1)
                            #print(stop)
                            if start in sp and head.split('_')[0] in sp:
                               # print('Line 64')
                                self.locusDic[head] = sp[8].split(';')[1]
                            elif stop in sp and head.split('_')[0] in sp:
                                #print('Line 76')
                                self.locusDic[head] = sp[8].split(';')[1]

                                            

    def findGenes(self):
        print("find genes")
        for head, locus  in self.locusDic.items():
            for sp in self.dataList:
                if locus == sp[8].split(';')[1]:
                    #print('Line 86')
                    
                    self.geneDic[head] = sp[8].split(';')[0].split('-')[-1] #check this again
                    """
                    dictionary in which key points to a list
                    """


    
                
    def getProteins(self, fam):
        print("Get Proteins")

        myReaderIE = FastAreader('Data_{0}/{1}'.format(self.NAME,self.fasta))
        for header, sequence in myReaderIE.readFasta():  
            for head, locus in self.geneDic.items():
#                print('New Protein')
                
                
                
                if locus in header and locus not in self.proteinList:
                    self.proteinList.append(locus)
                    with open('Data_{0}/ieProteins2.fa'.format(self.NAME), 'a') as proteinFile:
                            
                        proteinFile.write('>{1}_{2}_{3}_{0}\n'.format(locus, fam, self.NAME, head))
                        proteinFile.write('{0}\n'.format(sequence))
        return self.proteinDic
        


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
                            #print(link)
                            refPreAssembly = refLink.split("/")[-1]
                           # print(preAssembly)
                           # preList = preAssembly.split("_")
                           # print(preList[0])
                           # print(preList[1])

                            refAssembly = refPreAssembly.split("_")[0] + "_" + refPreAssembly.split("_")[1]
                           # print(Assembly)
                            refFasta = refLink + '/*_genomic.fna.gz'
                            refAnnotation = refLink + '/*_genomic.gff.gz'
                            genPreAssembly = genLink.split("/")[-1]
                           # print(preAssembly)
                           # preList = preAssembly.split("_")
                           # print(preList[0])
                           # print(preList[1])

                            genAssembly = genPreAssembly.split("_")[0] + "_" + genPreAssembly.split("_")[1]
                           # print(Assembly)
                            genFasta = genLink + '/*_genomic.fna.gz'
                            genAnnotation = genLink + '/*_genomic.gff.gz'

                            assemblyDic = {'Species': sp[0],'RefSeq' : refAssembly, 'refFasta' :refFasta , 'refAnnotation' : refAnnotation, 'GenBank' : genAssembly, 'genFasta' :genFasta , 'genAnnotation' : genAnnotation}
                            assemblies.append(assemblyDic)
                        except IndexError:
                        
                            print("Fuck you")
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
                                annotation = link + '/*_genomic.gff.gz'
    
                                assemblyDic = {'Species': sp[0],'GenBank' : Assembly, 'genFasta' : fasta, 'genAnnotation' : annotation}
                                assemblies.append(assemblyDic)
                            except IndexError:
                                pass

            return assemblies
            

class proteinBlast():
    
    def __init__(self, proteinDic, NAME):
        self.proteinDic = proteinDic
        self.NAME = NAME
    def proteinFile(self):
        with open('ieProteins.fa', 'w') as proteinFile:
            for gene, seq in self.proteinDic.items():
                proteinFile.write('>{0}\n'.format(gene))
                proteinFile.write('{0}\n'.format(seq))
                
    def allvallBlast(self):
        os.system("/Users/lgozasht/Desktop/PLEASEDONOTFUCKINGDELETE/Mpusillaproject/ncbi-blast-2.7.1+/bin/makeblastdb -dbtype prot -in Data_{0}/ieProteins2.fa -title proteins -out Data_{0}/proteinsDB".format(self.NAME))
        os.system("/Users/lgozasht/Desktop/PLEASEDONOTFUCKINGDELETE/Mpusillaproject/ncbi-blast-2.7.1+/bin/blastp -db Data_{0}/proteinsDB -query Data_{0}/ieProteins2.fa -outfmt 6 -out Data_{0}/proteins_all-vs-all.tsv".format(self.NAME))
        os.system("awk '$1 != $2' Data_{0}/proteins_all-vs-all.tsv > Data_{0}/proteins_all-vs-all_deduped.tsv".format(self.NAME))
#        
#        os.system("../Tools/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype prot -in Data_{0}/ieProteins2.fa -title proteins -out Data_{0}/proteinsDB".format(self.NAME))
#        os.system("../Tools/ncbi-blast-2.8.1+/bin/blastp -db Data_{0}/proteinsDB -query Data_{0}/ieProteins2.fa -outfmt 6 -out Data_{0}/proteins_all-vs-all.tsv".format(self.NAME))
#        os.system("awk '$1 != $2' Data_{0}/proteins_all-vs-all.tsv > Data_{0}/proteins_all-vs-all_deduped.tsv".format(self.NAME))
        
    def totalFile(self):
        with open('Data_{0}/total_all-vs-all_deduped.tsv'.format(self.NAME), 'a') as totalFile:

            with open('Data_{0}/proteins_all-vs-all_deduped.tsv'.format(self.NAME), 'r') as protFile:

                for line in protFile:
                    line = line.rstrip()
                    totalFile.write('{0}\n'.format(line))
        
    def alignLength(self):
        os.system('rm -r ./Data_{0}/gene_duplications.tsv'.format(self.NAME))
        os.system('rm -r ./Data_{0}/deDuped.tsv'.format(self.NAME))
        with open('Data_{0}/total_all-vs-all_deduped.tsv'.format(self.NAME), 'r') as allvallOut:
            for line in allvallOut:
                line = line.rstrip()
                sp = line.split('\t')
                #query1 = sp[0]

                if float(sp[2]) > 50.0 and int(sp[3]) > 100 and float(sp[10]) < 1*10**(-5):
                    with open('./Data_{0}/gene_duplications.tsv'.format(self.NAME), 'a') as file2:
                        file2.write('{0}\n'.format(line))
                        
                else:
                    with open('./Data_{0}/deDuped.tsv'.format(self.NAME), 'a') as file3:
                        file3.write('{0}\n'.format(line))
                    
                    

        
        

                
        

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
                os.system('rm -r Data_{0}/total_all-vs-all_deduped.tsv'.format(NAME))
                os.system("rm -r Data_{0}/ieProteins2.fa".format(NAME))

    
                print('Processing {0}'.format(NAME))
                os.system('rm -r GCA*')
                os.system('rm -r GCF*')
                os.system('wget {0}'.format(assembly['genAnnotation']))
                print(assembly['genAnnotation'])
                os.system('wget {0}'.format(assembly['genFasta']))
                print(assembly['genFasta'])
                
                proteinFastaList = assembly['genFasta'].split('_')
                del proteinFastaList[-1]
               # del proteinFastaList[-2]
                print(proteinFastaList)
                proteinFasta = '_'.join(proteinFastaList) + 'protein.faa.gz'
                os.system('wget {0}'.format(proteinFasta))
                print(proteinFasta)

                os.system('gunzip {0}*'.format(NAME))
                os.system('cp {0}* ./Data_{0}'.format(NAME))
                os.system('gunzip ./Data_{0}/*'.format(NAME))
                os.system('rm -r {0}*'.format(NAME))
                #fastaList = assembly['genFasta'].split("/")
               # fastaGz = fastaList[-2]
              #  fasta = fastaGz + '_genomic.fna'
                annotationList = assembly['genAnnotation'].split("/")
                annotationGz = annotationList[-2]
                annotation = annotationGz + '_genomic.gff'
                proteinList = proteinFasta.split("/")
                proteinGz = proteinList[-2]
                proteinFile = proteinGz + '_protein.faa'
                print(proteinFile)
                
                
                print(annotation)
            #    print(fasta)
                print(NAME)
                
                myProteinData = retrieveProteins(ieFile, annotation, proteinFile, NAME)
                headDic = myProteinData.readIeFile()
                for fam, headList in headDic.items():
                    print(fam)
                    print(headList)
                    myProteinData.readGFF(headList)
                    myProteinData.findGenes()
                    proteinDic = myProteinData.getProteins(fam)
                
                    myBLAST = proteinBlast(proteinDic, NAME)
                    myBLAST.allvallBlast()
                    myBLAST.totalFile()
                myBLAST.alignLength()
                
                
            else:
                print("no IEs")
                continue
        except KeyError:
            print('no genbank assembly')
            
            
if __name__ == "__main__":
    main() 
