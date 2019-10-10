#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:53:11 2019

@author: lgozasht
"""

"""
Get rid of gene duplicates.

"""

from sequenceAnalyzer import FastAreader
import os


        
class removeDups():
    
    def __init__(self,ieList,NAME,famDic):
        self.ieList = ieList
        self.NAME = NAME
        self.famDic = famDic
    
    def remove(self):
        """
        if ie is not the result of a gene duplicate and ie is not the 
        result of different isoform predictions print it in a final fasta file
        """
        os.system('rm -r ./Data_{0}/final_exon_overlaps.tsv'.format(self.NAME))
        with open('Data_{0}/exons_deduped.tsv'.format(self.NAME),'r') as file1:
            for line in file1:
                line = line.rstrip()
                sp = line.split('\t')
#                query1 = sp[0].split('_')
#                coords1 = query1[1].split('-')
#                starta = coords1[0]
#                stopa = coords1[1]
                print(float(int(sp[3])/100))
                if float(int(sp[3])/(100)) > 0.8:
                    with open('./Data_{0}/final_exon_overlaps.tsv'.format(self.NAME), 'a') as file2:
                        file2.write('{0}\n'.format(line))
        
        """
        remove shits
        """
        removedDic = {}
        os.system('rm -r ./Data_{0}/finalIEs.fa'.format(self.NAME))
        upList = []
        downList = []
        upDownList = []
        downUpList = []
    #    try:
        with open('Data_{0}/final_exon_overlaps.tsv'.format(self.NAME),'r') as file3:
            for line in file3:
                line = line.rstrip()
                sp = line.split('\t')
                query1 = sp[0].split('_')
                query2 = sp[1].split('_')
                node1 = query1[0]+'_'+query1[1]
                node2 = query2[0]+'_'+query2[1]
                nodeList = [node1, node2]
                nodeList.sort()
                if query1[2] == 'up' and query2[2] == 'up':
                    upList.append(nodeList)
                if query1[2] == 'down' and query2[2] == 'down':
                    downList.append(nodeList)
                if query1[2] == 'up' and query2[2] == 'down':
                    upDownList.append(nodeList)
                if query1[2] == 'down' and query2[2] == 'up':
                    downUpList.append(nodeList)
            for up in range(0,len(upList),1):
                upIE = upList[up]
                for down in range(0,len(downList),-1):
                    downIE = downList[down]
                    if upIE[0] == downIE[0] and upIE[1] == downIE[1]:
                        ieList = upList[up]
                        ie = ieList[0]
                        fam = ie[0:2:1]
                        if fam[1].isdigit() == False:
                             fam = str(fam[0])
                        else:
                             fam = str(fam)

                        self.ieList.remove(ie)
                        if str(fam) not in removedDic.keys():
                            removedDic[str(fam)]=1
                        else:
                            removedDic[str(fam)]+=1
            for up in range(0,len(upDownList),1):
                for down in range(0,len(downUpList),-1):
                    if upDownList[up] == downUpList[down]:
                        ieList = upDownList[up]
                        ie = ieList[0]
                        fam = ie[0:2:1]
                        if fam[1].isdigit() == False:
                             fam = str(fam[0])
                        else:
                             fam = str(fam)

                        self.ieList.remove(ie)
                        if str(fam) not in removedDic.keys():
                            removedDic[str(fam)]=1
                        else:
                            removedDic[str(fam)]=+1
                
                
                
#                if upStream == True and downStream == True:
#                            
#                            fam = ieDic['fam']
#                            self.ieList.remove(ieDic)
#                            if str(fam) not in removedDic.keys():
#                                removedDic[str(fam)]=1
#                            else:
#                                removedDic[str(fam)]=+1
#                            break
#
#
#                        
#
#                
#                
#                
#                for ieDic in self.ieList:
#                    upStream = False
#                    downStream = False
#
#                    for line in file3:
#                        line = line.rstrip()
#                        sp = line.split('\t')
#                        query1 = sp[0].split('_')
#                        coords1 = query1[1].split('-')
#                        start = coords1[0]
#                        stop = coords1[1]
#                        print(ieDic['start'])
#                        print(start)
#                        if int(ieDic['start'])==int(start) and int(ieDic['stop'])==int(stop) and query1[2] == 'up':
#                            upStream = True
#                            query2 = sp[1].split('_')
#                            coords2 = query2[1].split('-')
#                            start2 = coords2[0]
#                            stop2= coords2[1]
#
#                        if int(ieDic['start'])==int(start) and int(ieDic['stop'])==int(stop) and query1[2] == 'down':
#                            if upStream == True:
#                                downStream = True
#                            
#                            
#                            
#                        if upStream == True and downStream == True:
#                            
#                            fam = ieDic['fam']
#                            self.ieList.remove(ieDic)
#                            if str(fam) not in removedDic.keys():
#                                removedDic[str(fam)]=1
#                            else:
#                                removedDic[str(fam)]=+1
#                            break
        except FileNotFoundError:
            pass

        
            
        for ieDic in self.ieList:
            fam =  ieDic['fam']
            if self.famDic[str(fam)] - removedDic[str(fam)] > 5:
                print(self.famDic[str(fam)])
                print(removedDic[str(fam)])
                with open('./Data_{0}/finalIEs.fa'.format(self.NAME), 'a') as file4:
                    file4.write(">{0}{1}_{2}-{3}\n".format(ieDic['fam'], ieDic['scaff'], ieDic['start'], ieDic['stop']))
                    file4.write("{0}\n".format(ieDic['seq']))
                            
                

class geneDups():
    
    def __init__(self, ieList, NAME):
        self.ieList = ieList
        self.NAME = NAME
    def makeFile(self):
        with open('Data_{0}/exons.fa'.format(self.NAME), 'a') as file:
            try:
                for ieDic in self.ieList:
                    fam = ieDic['fam']
                    up = ieDic['up']
                    down = ieDic['down']
                    scaff = ieDic['scaff']
                    start = ieDic['start']
                    stop = ieDic['stop']
                    if len(up)> 0 and len(down) > 0:
                        file.write(">{0}{1}_{2}-{3}_up\n".format(fam,scaff,start,stop))
                        file.write("{0}\n".format(up))
                        file.write(">{0}{1}_{2}-{3}_down\n".format(fam,scaff,start,stop))
                        file.write("{0}\n".format(down))
            except KeyError:
                pass
                
    def avaBlast(self):
        os.system("./Tools/ncbi-blast-2.7.1+/bin/makeblastdb -dbtype nucl -in {0}Data_{1}/exons.fa -title introns -out {0}Data_{1}/exonsDB".format(self.NAME))
        os.system("./Tools/ncbi-blast-2.7.1+/bin/blastn -db {0}Data_{1}/exonsDB -query {0}Data_{1}/exons.fa -outfmt 6 -perc_identity 80 -out {0}Data_{1}/exons.tsv".format(self.NAME))
        os.system("awk '$1 != $2' Data_{0}/exons.tsv > Data_{0}/exons_deduped.tsv".format(self.NAME))
        
    
class csvReader():
    
    
    def __init__(self, dataFile = ''):
        """
        Initialize blasthit as the blast output file.
    
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
                        pass
                    else:
                        Assembly = sp[5]
                        link = sp[14]
                        link = link.strip('\"')
                        link = link.strip('\"')
                        Assembly = Assembly.strip('\"')
                        Assembly = Assembly.strip('\"')

                        fasta = link + '/*_genomic.fna.gz'
                        
                        assemblyDic = {'Reference' : Assembly, 'Fasta':fasta}
                        assemblies.append(assemblyDic)            
                
                return assemblies
        
def main():
    myData = csvReader('eukaryotes.csv')
    genomeData = myData.csv()
         
    for assembly in genomeData:
        
        NAME = assembly['Reference']

        if 'GCA' not in assembly:
            pass
        else:
            ieList = []
            print('Removing gene duplicsates from {0}'.format(NAME))            

            os.system('wget {0}'.format(assembly['Fasta']))
            os.system('gunzip {0}*'.format(NAME))
            os.system('cp {0}* ./Data_{0}'.format(NAME))
            os.system('gunzip ./Data_{0}/*'.format(NAME))
            os.system('rm -r {0}*'.format(NAME))
            fastaList = assembly['Fasta'].split("/")
            fastaGz = fastaList[-2]
            fasta = fastaGz + '_genomic.fna'
            print(fasta)
            famDic={}
            
            myReaderIE = FastAreader('Data_{0}/IEfamilies*'.format(NAME))
            for header, seq in myReaderIE.readFasta():
                headList = header.split('_')
                head = headList[0]
                fam = head[0:1:1]
        
                if isinstance(fam,int) == False:
                    fam = int(fam[0])
                else:
                    fam = int(fam)
        
                scaff = headList[2]
                coord = headList[3]
                coordList = coord.split("-")
                start = int(coordList[0])
                stop = int(coordList[1])
                ieDic = {'fam':fam, 'scaff':scaff,'start':start,'stop':stop, 'seq': seq}
                if str(fam) not in famDic.keys():
                    famDic[str(fam)]=1
                else:
                    famDic[str(fam)]+=1
                ieList.append(ieDic)
           # print(ieList)
            
            for i in range(0,len(ieList)-1,1):
                ieDic1 = ieList[i]
                start1 = ieDic1['start']
                stop1 = ieDic1['stop']
                scaff1 = ieDic1['scaff']
                for k in range(0,len(ieList),-1):
                    ieDic2 = ieList[k]
                    start2 = ieDic2['start']
                    stop2 = ieDic2['stop']
                    scaff2 = ieDic2['scaff']
                    if scaff1==scaff2:
                        if start2 < start1 and stop1 < stop2:
                            ieList.remove(ieDic1)
                            print(ieDic)
            

            
            ieDicSeqs = []
            myReaderGenome = FastAreader('Data_{0}/{1}'.format(NAME,fasta))
            for header, sequence in myReaderGenome.readFasta():
                portion =  header.split(" ")
                head = portion[0]
                for ieDic in ieList:
                    if ieDic['scaff']==head:
                        start = ieDic['start']
                        stop = ieDic['stop']
                        up = sequence[int(start)-100:int(start):1]
                        down =sequence[int(stop):int(stop)+100:1]
                        ieDic['up'] = up
                        ieDic['down'] = down

                        ieDicSeqs.append(ieDic)
                     
            myFolder = geneDups(ieList,NAME)
            myFolder.makeFile()
            myFolder.avaBlast
            
            myRemoveDups = removeDups(ieList, NAME, famDic)
            myRemoveDups.remove()
            
            
            
            
            
            
            
            
            
            
            
    
if __name__ == "__main__":
    main() 
    