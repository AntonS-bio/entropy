#!/bin/python
from Bio import SeqIO, AlignIO


def run(outdir, samplesDir):
    #sampleSet=sys.argv[1]

    #wd="/mnt/storage5/anton/entropy/"


    #if sampleSet=="Portugal":
    #    outputDir=wd+"/IQTree/P1/"
    #else:
    #    outputDir=wd+"/IQTree/ST258/"


    #outputDir=wd+"/IQTree/ST258pathogen/"
    #fastaDir=wd+"/allignedFasta/"

    #top 49 from busco pangenome
    genes=open(outdir+"genesForIqtree.txt").read().splitlines()

    sequences={}
    nexFile=open(outdir+"nexus.txt", "w") #write nex partition file along side collecting sequences
    nexFile.write("#nexus\n")
    nexFile.write("begin sets;\n") 
    totalLengths=0
    partitionCounter=1
    for gene in genes:
        msa_fasta=open(samplesDir+gene+".fasta")
        for record in SeqIO.parse(msa_fasta,'fasta'):
            if record.id in sequences:
                sequences[record.id]=sequences.get(record.id)+str(record.seq)[1:len(str(record.seq))-1]
            else:
                sequences[record.id]=str(record.seq)[1:len(str(record.seq))-1]
        msa_fasta.close()
        nexFile.write("    charset part"+str(partitionCounter)+" = "+outdir+"alligned.fasta: "+str(totalLengths+1)+"-"+str(len(sequences.get(record.id)))+";\n")
        partitionCounter+=1
        totalLengths=len(sequences.get(record.id))
    nexFile.write("end;\n")
    nexFile.close()


    mergedFasta=open(outdir+"alligned.fasta", "w")
    for key in sequences.keys():
        mergedFasta.write(">"+key+"\n")
        mergedFasta.write(sequences.get(key)+"\n")
    mergedFasta.close()


    mergedFasta=open(outdir+"alligned.fasta", "r")
    records = AlignIO.parse(mergedFasta, "fasta")

    mergedFasta.close()

# wd="/mnt/storage5/anton/PATRIC/Tree/"
# tempData="/mnt/storage5/anton/PATRIC/Tree/PatricIsoniazidTemp/"
# outDir=wd#+"/"+RunName+"/"
# run(outDir, tempData+"/allignedFasta/")
