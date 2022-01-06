import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from os import chdir, chdir, mkdir
from os.path import exists
import multiprocessing as mp
import statistics
import numpy
import scipy.stats
import sys
import shutil



RunName=sys.argv[1]

samplesFile="allPortugalSamples.txt"
sampleExtension=".fna"
samplesDir="/mnt/storage5/anton/Portugal_fastas/"
gffDir="/mnt/storage5/anton/prokka_genus_wo_hypotheticaln/"
minIdentity=95
maxLenDifference=0.05

wd="/mnt/storage5/anton/entropy/Portugal_Kp_all/"
tempData=wd+"/"+RunName+"Temp/"
outDir=wd#+"/"+RunName+"/"



#clear all temp folders
if exists(tempData):
    shutil.rmtree(tempData)
mkdir(tempData)
for directory in ["allignedFasta/","blastResults/","unallignedFasta/","genes/"]:
    mkdir(tempData+directory)

def runBlast(sample):
    global samplesDir
    global sampleExtension
    global tempData
    sample=sample.replace(sampleExtension,"")
    chdir(tempData+"/genes/")
    subprocess.call("blastn -query " +samplesDir+sample+sampleExtension+" -max_target_seqs 1000000000 -task 'megablast' -db "+tempData+"/genes/"+refSample+"genes -num_threads 4 -evalue 1.0E-4 -word_size 11 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident qseq sstrand \" > "+tempData+"/blastResults/"+sample+".csv", shell=True)
    #subprocess.call("tblastx -query " +samplesDir+sample+sampleExtension+" -max_target_seqs 1000000000 -db "+tempData+"/genes/"+refSample+"genes -num_threads 4 -evalue 1.0E-4 -word_size 3 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident qseq sstrand \" > "+tempData+"/blastResults/"+sample+".csv", shell=True)

def runMafft(gene):
    global tempData
    gene=gene+".fasta"
    subprocess.call("mafft --quiet --auto "+tempData+"/unallignedFasta/"+gene+" > "+tempData+"/allignedFasta/"+gene, shell=True)


samples=[]
with open(outDir+samplesFile) as file:
    for line in file:
        samples.append(line.strip())
refSample=samples[0] #first sample in list is the reference sample. Shouldn't make much difference which to pick.
refSample="1408962.4"

refFasta=SeqIO.to_dict(SeqIO.parse(samplesDir+refSample+sampleExtension, "fasta"))
genes={}
sequencesToWrite={}
for line in open(gffDir+refSample+".gff"):
    if line[0]!="#":
        values=line.split("\t")
        if len(values)>2 and values[2]=="CDS":
            id=values[8].split(";")[0].replace("ID=","")
            genes[id]=[]
            if values[6]=="+":
                sequence=refFasta[values[0]][int(values[3]):int(values[4])]
            else:
                sequence=refFasta[values[0]][int(values[3]):int(values[4])].reverse_complement()
            if len(sequence.seq)>100:
                if not id in sequencesToWrite:
                    sequencesToWrite[id]=""
                if len(sequence.seq)>len(str(sequencesToWrite[id])):
                    sequencesToWrite[id]=str(sequence.seq)

with open(tempData+"/genes/"+refSample+"genes.fasta", "w") as geneFile:
    for id in sequencesToWrite:
        geneFile.write(">"+id+"\n")
        geneFile.write(sequencesToWrite[id]+"\n")

chdir(tempData+"/genes/")
subprocess.call("makeblastdb -in "+ tempData+"/genes/"+refSample+"genes.fasta -blastdb_version 4 -title "+refSample+"genes -out "+refSample+"genes -dbtype nucl", shell=True)



#run blast of individual samples vs all genes of reference
if __name__ == '__main__':

    pool=mp.Pool(10)
    results=pool.map(runBlast, samples)
    #result=runBlast(samples[0])
    pool.close()
    pool.join()

    #collect gene lengths
    print("colleting genes from blast results")
    geneLengths=dict.fromkeys(genes)
    geneSequences={} #key=sample name, value={refGene, sequence} 
    for sample in samples:
        sampleGenes={}
        sample=sample.replace(sampleExtension,"")
        geneSequences[sample]={}
        for line in open(tempData+"/blastResults/"+sample+".csv"):
            values=line.strip().split("\t")
            allignmentLength=int(values[2])-int(values[1])
            if (values[3] not in sampleGenes or sampleGenes[values[3]]<allignmentLength) and float(values[5])>minIdentity: #value[5] is identity %
                sampleGenes[values[3]]=allignmentLength
                if values[7]=="minus":
                    values[6]=str(Seq(values[6]).reverse_complement())
                geneSequences[sample][values[3]]=values[6]
        for sampleGene in sampleGenes.keys():
            genes[sampleGene].append(sampleGenes[sampleGene])

    print("writing genes that passed quality checks")

    genesToAllign=set()
    with open(outDir+"/samplesPerGene.tsv","w") as tempOutput:
        for gene in genes.keys():
            complete=0
            if len(genes[gene])>0:
                geneMedianLB=statistics.median(genes[gene])*(1-maxLenDifference)
                geneMedianUB=statistics.median(genes[gene])*(1+maxLenDifference)
                for sample in genes[gene]:
                    complete+=1 if (sample>geneMedianLB and sample<geneMedianUB) else 0
            if complete==len(samples):
                genesToAllign.add(gene)
                with open(tempData+"/unallignedFasta/"+gene+".fasta", "w") as allignedFile:
                    for sample in geneSequences.keys():
                        allignedFile.write(">"+sample+"\n")
                        allignedFile.write(geneSequences[sample][gene]+"\n")
                #extract alligned sequences
            
            tempOutput.write(gene+"\t"+str(len(genes[gene]))+"\t"+str(complete)+"\n")

    #sys.exit()
    pool=mp.Pool(30)
    results=pool.map(runMafft, list(genesToAllign))
    pool.close()
    pool.join()

    print("Calculating Shannon's entropy")
    tempOutput=open(outDir+"/entropy.tsv","w")

    
    genesToAllign=list(genesToAllign)
    geneEntropy=[]
    for gene in genesToAllign:
        print(genesToAllign.index(gene)/len(genesToAllign))
        data=[]
        geneLength=0
        totalSequences=0
        fastaFile = open(tempData+"/allignedFasta/"+gene+".fasta",'r')
        for record in SeqIO.parse(fastaFile,'fasta'):
            totalSequences+=1
            if len(data)==0:
                data=numpy.array([list(str(record.seq))],dtype=numpy.unicode) 
                geneLength=len(str(record.seq))
            else:
                data=numpy.concatenate((data,numpy.array([list(str(record.seq))],dtype=numpy.unicode)), axis=0)

        totalEntropy=0
        for i in range(0,geneLength):
            values, counts = numpy.unique(data[:,i], return_counts=True)
            #print(values)
            totalEntropy+=scipy.stats.entropy(counts/totalSequences)
        print(gene)
        geneEntropy.append([gene,totalEntropy/geneLength])
        tempOutput.write(gene+"\t"+str(totalEntropy/geneLength)+"\t"+str(geneLength)+"\n")
    tempOutput.close()
    #skip first five genes to avoid forcing bias on results
    geneEntropy=sorted(geneEntropy, reverse=True, key=lambda value: value[1])

    if len(geneEntropy)>5:
        with open(outDir+"/genesForIqtree.txt", "w") as output:
            i=0
            while i<100 and i<len(geneEntropy) and geneEntropy[i][1]>0:
                if i>10:
                    output.write(geneEntropy[i][0]+"\n")
                i+=1
    
    import merge_MSA
    merge_MSA.run(outDir, tempData+"/allignedFasta/")
    chdir(outDir)

    subprocess.call("iqtree -redo -nt AUTO -pre "+outDir+"/"+RunName+" -s "+outDir+"/alligned.fasta -spp "+outDir+"/nexus.txt -bb 1000", shell=True, executable="/bin/bash")
