# entropy
Compare genes across samples and generate phylogenetic tree based on most diverse genes.

the inputs are listed at the top:
1. file (txt or similar) containing list of fasta files (assemblies) to be used
2. extension of fasta files (fna or fasta)
3. directory with the assemblies
4. directory with corresponding gff files
5. minimum identity (%) between two genes to consider them homologues
6. minimum lenght difference (% of total length) to consider two genes as homologues

Both 5 and 6 must be met.

The most diverse (based on Shannon entropy) genes are taken for construction of phylogenetic tree.

The required programs are mafft, IQTREE and BLAST
