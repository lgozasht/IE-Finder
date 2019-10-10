#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 11:52:20 2018

@author: Landen Gozashti under the guidance of Russel Corbett-Detig
"""


"""

This script identifies introner elements in genomes and only requires an annotation file and assembly fasta file.

This program is written in an object oriented fashion.

"""

import matplotlib
matplotlib.use('Agg')

from sequenceAnalyzer import FastAreader
import os
from Bio.Seq import Seq

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from networkx.algorithms import community
#plt.use(‘agg’)


class IntronRecognition():
    """
    
    Extract all introns less than 300bp using the exon coordinates passed by 
    
    """
    
    def __init__(self, cdsDic):
        
        self.cdsDic = cdsDic

    def introns(self):
        intronList = []


        for gene in self.cdsDic:
            cdsList = self.cdsDic[gene]
            if len(cdsList) > 1:
                starts = [i[2] for i in cdsList]
                stops = [i[3] for i in cdsList]
                scaffoldList = [i[0] for i in cdsList]
                strandList = [i[1] for i in cdsList]
                if '-' in strandList:
                    
                    stops.reverse()
                    starts.reverse()
                scaffold = scaffoldList[0]
                starts.pop(0)
                i = 0
                while i < len(starts):
                    intronStart = stops[i]
                    intronStop = starts[i]
                    if (int(intronStop) - int(intronStart)) < 2000: #change back to 300
                        intron = [scaffold, intronStart, intronStop]
                        intronList.append(intron)
                    i += 1
        return intronList
        
class RetrieveIntrons():
    
    """
    
    Retrieve introns from the genome assembly using the intron coordinates parsed by IntronRecognition.
    Return a list of introns. Each "intron" in this list appears as a string structured:
        "scaffold_intronstart_intronstop,intronsequence"
    
    """
    
    def __init__(self, header, sequence, introns):
        self.sequence = sequence
        self.introns = introns
        self.header = header
        
    def retrieve(self):
        intronSeqs = []
        for intron in self.introns:
            if intron[0] == self.header:
                intronStart = intron[1]
                intronStop = intron[2]
                scaffold = intron[0]
                if int(intronStop) - int(intronStart) < 100:
                    seq = self.sequence[(int(intronStart) - 0): (int(intronStop)-1 + 0): 1]
                    realstop = int(intronStop) - 1
                    intronSeqs.append('{0}_{1}-{2},{3}'.format(scaffold, intronStart, realstop, seq))
                else:
                    seq = self.sequence[(int(intronStart)): (int(intronStop)-1): 1]
                    realstop = int(intronStop) - 1
                    intronSeqs.append('{0}_{1}-{2},{3}'.format(scaffold, intronStart, realstop, seq))
        return intronSeqs
                 
class MakeFasta():
    
    """
    Write a fasta file containing all extracted intron sequences. Provide a header for each intron with the following format:
        >"First letter of genus followed by species_scaffold_start-stop"
        
    
    """
    
    def __init__(self, intronList, PATH, NAME):
        self.intronList = intronList
        self.PATH = PATH
        self.NAME = NAME
    
    def fasta(self):
        #os.system('mkdir {0}'.format(self.NAME))
        #os.system('mkdir ./Shit_I_actually_care_about/related_species')
        os.system('mkdir ./Data_{1}'.format(1,self.NAME))

        with open ('./Data_{1}/Reads.fa'.format( self.PATH, self.NAME), 'w') as file:
            for item in self.intronList:
                for intron in item:
                    header, sequence = intron.split(',')
                    if len(sequence) > 0: #Throw out short ass introns which are probably not IEs
                        if sequence[0] == 'C':
                            seq = Seq(sequence)
                            revcomp = seq.reverse_complement()
                            file.write('>{0}_{1}\n'.format(self.NAME, header))
                            file.write('{0}\n'.format(revcomp))
                        else:                                
                           file.write('>{0}_{1}\n'.format(self.NAME, header))
                           file.write('{0}\n'.format(sequence))
                    
class FindIEs():
    """
    
    Count the number of putative insertions
    
    """
    
    
    def __init__(self, alignFile):
        self.alignFile = alignFile
        
    def find(self):        
        count = 0
        #ieList = []
        with open(self.alignFile) as f:
            for line in f:
#                line = line.rstrip()
#                sp = line.split('\t')
                count += 1
#                intron = sp[0]
#                if sp[0] == intron:
#                    count += 1
#                else:
#                    count = 0
#                if count == 10:
#                    ieList.append(intron)
#                    count = 0
        return count
                
        
                               
                
class GeneDataDic():
    """
    Parse the assembly annotation file.
    Retrieve exon coordinates for each gene: the scaffold, strand, start, and stop for each exon. 
    Return a dictionary with genes as keys and lists of exon coordinates within each gene as values.
    """
    
    def __init__(self, annotation = ''):
        self.annotation = annotation
        
    
    def genedatadic(self):
      
        cdsCoordinates = {}
        
        with open(self.annotation,'r') as f:
            next(f) 
            next(f)     
            next(f)          
            next(f)
            next(f)
            next(f)
            next(f)
            next(f)

            for line in f:
                line = line.rstrip()
                sp = line.split('\t')
                if '#' in line:
                    pass
                else:
                    try:
                        codon = sp[2]
                    except IndexError:
                        pass
                                        
                    if codon == 'exon':
                        try:
                            portion =  sp[8].split(";")
                            gene = portion[1]
                            #print(gene)
                            cdsStart = sp[3]
                            cdsStop = sp[4]
                            scaffold = sp[0]
                            strand = sp[6]
                            cdsList = [scaffold, strand, cdsStart, cdsStop]
                        except IndexError:
                            pass
                        if gene not in cdsCoordinates:
                            cdsCoordinates[gene]=[cdsList]
                                       
                        else:
                            cdsCoordinates[gene].append(cdsList)
                
        return cdsCoordinates
       
class Extraction():
    """
    
    Extract the BLAST results for each intron.  Return a list of sequence strings exrracteed
    
    """
    def __init__(self, Data, header, sequence):
        self.Data = Data
        self.header = header
        self.sequence = sequence
   
    def extract(self):
        resultsList = []
        
        for query in self.Data:
           
            if self.header == query['scaffold']:
                scaffold = query['scaffold']
                if query['stop'] > query['start']: #if the intron is on the coding strand
                    stop = query['stop']
                    start =  query['start']
                    extraction = self.sequence[(query['start']):(query['stop']): 1]
#                    print(start)
#                    print(stop)
#                    print(extraction)
                    resultsList.append('{0}_{1}-{2},{3}'.format(scaffold, start, stop, extraction))

                else: #if the ibtntron is on the reverse complenebt
                    extraction = self.sequence[(query['stop']):(query['start']): 1]
                    start = query['stop']
                    stop =  query['start']
                    seq = Seq(extraction)
                    revcomp = seq.reverse_complement()
                    resultsList.append('{0}_{1}-{2},{3}'.format(scaffold, start, stop, revcomp))
                    
                
        return resultsList

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
                        #print('what')
                        Assembly = sp[5]
                        link = sp[14]
                        link = link.strip('\"')
                        link = link.strip('\"')
                        Assembly = Assembly.strip('\"')
                        Assembly = Assembly.strip('\"')

                        print(link)
                        fasta = link + '/*_genomic.fna.gz'
                        annotation = link + '/*_genomic.gff.gz'
                        
                        assemblyDic = {'Reference' : Assembly, 'Fasta' : fasta, 'Annotation' : annotation}
                        assemblies.append(assemblyDic)
            
                
            return assemblies

class Graph():
    """
    Employs a graph theory algorithm to cluster introns. 
    Plots a scatter plot for cluster visualization.
    Rerurns communities formed from intron clusters.
    
    """
    def __init__(self, minimapOut, NAME):
        self.minimapOut = minimapOut
        self.NAME = NAME
    def graph(self):
       # from random import random
        IEfamilies = []
        IEs = []
        try:
            try:
                
                with open (self.minimapOut,'r') as file1:

                    for line in file1:

                        line = line.rstrip()
                        sp = line.split('\t')
                        query1 = sp[0].split('_')
                        coords1 = query1[3].split('-')
                        starta = coords1[0]
                        stopa = coords1[1]
                        length1 = int(stopa) - int(starta)
#                        query2 = sp[1].split('_')
#                        coords2 = query2[3].split('-')
#                        startb = coords2[0]
#                        stopb = coords2[1]
#                        length2 = int(stopb) - int(startb)
#                            print(float(int(sp[3])/length1))
                        if length1 > 0:
                            if float(int(sp[3])/length1) > 0.8:
                                with open('./Data_{0}/final_introns.tsv'.format(self.NAME), 'a') as file2:
                                    file2.write('{0}\n'.format(line))
#                            else:
#                                
#                                if float(int(sp[3])/length1) > 0.7:
#                                    with open('./Data_{0}/final_introns.tsv'.format(self.NAME), 'a') as file2:
#                                        file2.write('{0}\n'.format(line))
#                        else:
#                            print(float(int(sp[3])/length2))
#                            if length2 > 180:
#                                if float(int(sp[3])/length2) > 0.6:
#                                    with open('./Data_{0}/final_introns.tsv'.format(self.NAME), 'a') as file2:
#                                        file2.write('{0}\n'.format(line))       
#                            else:
#                                if float(int(sp[3])/length2) > 0.7:
#                                    with open('./Data_{0}/final_introns.tsv'.format(self.NAME), 'a') as file2:
#                                        file2.write('{0}\n'.format(line))
                
                data = pd.read_csv('./Data_{0}/final_introns.tsv'.format(self.NAME), sep='\t', header=None)
                df = nx.from_pandas_edgelist(data, source=0, target=1, edge_attr=True)
                
                
                communities_generator = community.girvan_newman(df)
                top_level_communities = next(communities_generator)
                IEfam = sorted(map(sorted, top_level_communities))
                for fam in IEfam:
                    if len(fam) > 5:
                        IEfamilies.append(fam)
                    else:
                        IEs = IEs + fam
                for fam in IEs:
                    df.remove_node(fam)
                print(len(IEfamilies))
                pos = nx.spring_layout(df)
                plt.figure(figsize=(20,20))
               # colors = [(random(), random(), random()) for _i in range(10)]
                nx.draw_networkx(df, pos, with_labels=False) 
                #plt.title('{0}'.format(self.NAME))
                
                #print(IEfamilies)
                colorListTwo = ['Red',
                                'Lime',
                                'Blue',
                                'Yellow',
                                'Cyan',
                                'Magenta',
                                'Maroon',
                                'Olive',
                                'Green',
                                'Purple',
                                'Teal',
                                'Navy', 
                                'darkgoldenrod', 
                                'cadetblue', 
                                'aqua', 
                                'darkseagreen',
                                'lightblue', 
                                'deeppink',
                                'mediumslateblue',
                                'gold',
                                'coral',
                                'slategrey',
                                'indigo',
                                'lawngreen',
                                'beige',
                                'rosybrown',
                                'mediumspringgreen',
                                'mediumblue',
                                'orchid',
                                'aliceblue',
                                'darkmagenta',
                                'darkseagreen',
                                'limegreen',
                                'powderblue',
                                'whitesmoke',
                                'navajowhite',
                                'oldlace',
                                'mediumpurple',
                                'mediumturquoise',
                                'peru',
                                'plum',
                                'wheat',
                                'thistle',
                                'sienna',
                                'linen',
                                'lemonchiffon',
                                'khaki',
                                'darkorchid',
                                'cornsilk',
                                'cadetblue'
                                ]
                count = 0
                for fam in IEfamilies:
                   # print(fam)
                    nx.draw_networkx_nodes(df, pos, with_labels=False, nodelist = fam, node_color = colorListTwo[count])
                    count += 1
                    
                    
                os.system("rm -r ./Data_{0}/final_introns.tsv".format(self.NAME))
                plt.savefig('./Data_{0}/Blast_fig.png'.format(self.NAME))
            except EmptyDataError:
                pass
        except NameError:
            IEfamilies = 'None'

        return IEfamilies


class GeneDups():
    """
    Remove introns from duplicate genes from putative IE file
    """
    def __init__(self, intronList, header, sequence):
        
        self.intronList = intronList
        self.sequence = sequence
        self.header = header
        self.flankList = []
    def flanks(self):
        
        for intron in self.intronList:
            if intron[0] == self.header:
                if int(intron[1]) < int(intron[2]): 
                    rightFlank = self.sequence[int(intron[2]) : int(intron[2]) + 20: 1]
                    leftFlank = self.sequence[int(intron[1]) - 20 : int(intron[1]) : 1]
                    self.flankList.append([intron[0], intron[1], intron[2], leftFlank, rightFlank])
                    #print(leftFlank)
                else:
                    rightFlank = self.sequence[int(intron[1]) : int(intron[1]) + 20: 1]
                    leftFlank = self.sequence[int(intron[2]) - 20 : int(intron[2]) : 1]
                  
                    self.flankList.append([intron[0], intron[1], intron[2], leftFlank, rightFlank])
        
    def prune(self):
        
        """
        Remove gene duplicates based on flanking regions of introns
        """
        noDups = []
        for myIntron in self.flankList:
            leftFlank = myIntron[3]
            rightFlank = myIntron[4]
            dupCount = 0
            for intron in self.flankList:
                if intron[3] == leftFlank and intron[4] == rightFlank:
                    dupCount += 1
#                    print('Removed {0}'.format(myIntron))
#                    print(leftFlank)
#                    print(intron[3])
            if dupCount == 1:
                newIntron = [myIntron[0], myIntron[1], myIntron[2]]
                noDups.append(newIntron)
               
                
            
                    
        return noDups
                    
                
            
                

def main():
    """
    Introners are for the boys.
    """

    myData = csvReader('algae.csv')
    genomeData = myData.csv()
         
    for assembly in genomeData:
        
        PATH = './'
        
        NAME = assembly['Reference']
        if 'GCA' not in NAME:
            pass
        else:
            print('Downloading files for {0} assembly'.format(NAME))
            #print(assembly['Reference'])
            
            os.system('mkdir Data_{0}'.format(NAME))
            os.system('rm -r ./Data_{0}/blastOutIntrons.fa'.format(NAME))
            print(assembly['Fasta'])
            os.system('wget {0}'.format(assembly['Fasta']))
            print(assembly['Annotation'])
            os.system('wget {0}'.format(assembly['Annotation']))
            os.system('gunzip {0}*'.format(NAME))
            os.system('cp {0}* ./Data_{0}'.format(NAME))
            os.system('gunzip ./Data_{0}/*'.format(NAME))
            os.system('rm -r {0}*'.format(NAME))
            
            annotationList = assembly['Annotation'].split("/")
            annotationGz = annotationList[-2]
            annotation = annotationGz + '_genomic.gff'
            print(annotation)
            
            fastaList = assembly['Fasta'].split("/")
            fastaGz = fastaList[-2]
            fasta = fastaGz + '_genomic.fna'
            print(fasta)
            
            print('Finding introner elements in {0}'.format(NAME))
            
            mygeneData = GeneDataDic('{0}Data_{1}/{2}'.format(PATH, NAME, annotation)) 
            cdsData = mygeneData.genedatadic()
            
            comparison = IntronRecognition(cdsData)
            intronList = comparison.introns()
            #Get rid of gene duplicates
            ###########################
            intronSeqs = []
            noDupList = []
                portion =  header.split(" ")
                head = portion[0]
                myDups = GeneDups(intronList, head, sequence)
                myDups.flanks()
                newList = myDups.prune() 
                noDupList = noDupList + newList
            #print(noDupList)
            
            ###########################
            
            
            
            print('Extracting Introns')
                        myReaderGenome = FastAreader('{0}Data_{1}/{2}'.format(PATH, NAME, fasta))
            for header, sequence in myReaderGenome.readFasta():

            
            for header, sequence in myReaderGenome.readFasta():
                portion =  header.split(" ")
                head = portion[0]
                MyIntrons = RetrieveIntrons(head, sequence, noDupList) #changed this from intronList
                intronSeqs.append(MyIntrons.retrieve())    
            finalIntronList = list(filter(None, intronSeqs))
            MyReads = MakeFasta(finalIntronList, PATH, NAME)    
            MyReads.fasta()
            ################################################################
            
            #print('Performing all vs all alignment with minimap2')
           # os.system("./Tools/minimap2/minimap2 -X -N 1000 {0}Data_{1}/Reads.fa {0}Data_{1}/Reads.fa | awk '$10>50' > {0}Data_{1}/overlaps.paf".format(PATH, NAME))
           # #os.system("./Tools/minimap2/minimap2 -X -N 1000 {0}Data_{1}/Reads.fa {0}Data_{1}/Reads.fa > {0}Data_{1}/overlaps.paf".format(PATH, NAME))
            ###############################################################
            print("Performing all-v-all BLAST")
            
            os.system("./Tools/ncbi-blast-2.7.1+/bin/makeblastdb -dbtype nucl -in {0}Data_{1}/Reads.fa -title introns -out {0}Data_{1}/intronsDB".format(PATH, NAME))
            os.system("./Tools/ncbi-blast-2.7.1+/bin/blastn -db {0}Data_{1}/intronsDB -query {0}Data_{1}/Reads.fa -outfmt 6 -perc_identity 80 -out {0}Data_{1}/all-vs-all.tsv".format(PATH,NAME))
            os.system("awk '$1 != $2 && awk $4 > 30' {0}Data_{1}/all-vs-all.tsv > {0}Data_{1}/all-vs-all_deduped.tsv".format(PATH,NAME))
            

            print('Clustering introns from minimap output')
            #Data = Graph('./Data_{0}/overlaps.paf'.format(NAME), NAME)
            Data = Graph('./Data_{0}/all-vs-all_deduped.tsv'.format(NAME), NAME)
            IEfamilies = Data.graph()
           # myReaderReads = FastAreader('./Data_{0}/Reads.fa'.format(NAME))
            count = 1
            with open('./Data_{0}/IEfamilies.fa'.format(NAME), 'w') as file:
            
                for family in IEfamilies:
                    if len(family) > 5:
                        #print(family)
                        #print(len(family))
                        for header, genomeSeq in myReaderGenome.readFasta():

                            for ie in family:
                                portion =  header.split(" ")
                                head = portion[0]
                                ieLabelList = ie.split('_')
                                scaff = ieLabelList[2]
                                coords = ieLabelList[3].split('-')
                                start = coords[0]
                                stop = coords[1]
                                if head == scaff:
                                    sequence = genomeSeq[int(start):int(stop):1]
                                    if sequence[0] == 'C': #If intron was found on the noncoding strand
                                        seq = Seq(sequence)
                                        revcomp = seq.reverse_complement() #Return reverse complement so that all introns are in the same orientation
        
                                        file.write('>{1}{0}\n'.format(ie, count))
                                        file.write('{0}\n'.format(revcomp))
                                    else:                                
                                        file.write('>{1}{0}\n'.format(ie, count))
                                        file.write('{0}\n'.format(sequence))
                        count += 1
    
            #Running minimap2 on Blastn results 
#            print('Running BLAST on putative introners')
#            os.system('./Tools/ncbi-blast-2.7.1+/bin/blastn -query {0}Data_{1}/IEfamilies.fa -subject {0}Data_{1}/{2} -perc_identity 85 -outfmt 6 >{0}Data_{1}/blasthit.txt'.format(PATH, NAME, fasta)) 
#        
#            
#            data = DataDic('{0}Data_{1}/blasthit.txt'.format(PATH, NAME))
#            
#            blastOut = data.datadic()
#            blastOutIntrons = []
#            blastOutDups = []
#           # print(blastOut)
#            for header, sequence in myReaderGenome.readFasta():
#                portion =  header.split(" ")
#                head = portion[0]
#                extractions = Extraction(blastOut, head, sequence)
#                blastOutDups.append(extractions.extract())
#            
#             #Check with Russ, we could accidently remove insertions here
#           # print(blastOutDups)
#            for result in blastOutDups: #Remove duplicates
#                if result is not '':
#                    for elem in result:
#                        if elem not in blastOutIntrons:
#                            blastOutIntrons.append(elem)
#                        else:
#                            print('Removed {0}'.format(elem))
#        
            print('Writing final IE fasta file')
#           
           # os.system('./bin/fastx_collapser < ./Data_{0}/blastOutIntrons.fa > ./Data_{0}/uniqueIEs.fa'.format(NAME))
           
            os.system('mv blastOutIntrons.fa . ./Data_{0}'.format(NAME))
            os.system("rm -r {0}Data_{1}/all-vs-all.tsv".format(PATH,NAME))
           # os.system("rm -r {0}Data_{1}/all-vs-all_deduped.tsv".format(PATH,NAME))
            os.system("gzip {0}Data_{1}/all-vs-all_deduped.tsv".format(PATH,NAME))
            os.system("rm -r {0}Data_{1}/intron*".format(PATH,NAME))
            os.system('rm -r ./Data_{0}/{0}*'.format(NAME))
            os.system('rm -r ./Data_{0}/o*'.format(NAME))

            print('-------------------------------wow----wow-----wee----wow-----')
            print('Just took a fat dub')
                
            
            
           # myIEs = FindIEs('{0}Data_{1}/blasthit.txt'.format(PATH, NAME))
          #  ieList = myIEs.find()
          
            #print('Identified {0} putative insertions in {1}'.format(ieList, NAME))
    
    
if __name__ == "__main__":
    main() 
    
    
    
    
