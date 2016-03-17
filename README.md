# ChimSim
**ChimSim: Chimeric transcripts simulator**

From an annotation, a genome index and numbers of chimeric transcripts wanted from 5 different classes 
(read-through, intra-chromosomal, inverted, interstrand, inter-chromosomal), makes a fasta file
with the wanted numbers of chimeric transcripts from the 5 classes

To clone the repository do:
git clone https://github.com/Chimera-tools/ChimSim.git

Running Bash/make_chimtr_from_annot_better.sh will provide usage
Usage: make_chimtr_from_annot_better.sh annot.gtf genome.gem n1 n2 n3 n4 n5 [list_gn_bt.txt]
where
- annot.gtf is a gene annotation (with at least exon and gene rows and where gene id and transcript id are in column 10 and 12) (mandatory)
- genome.gem is a genome index produced by gem (mandatory)
- n1 is the number of chimeric transcripts to make between non overlapping genes that are on the same chr, same str, ok gx order, with no intervening gene in the middle, and less distant than 100kb = read through (mandatory)
- n2 is the number of chimeric transcripts to make from non overlapping genes on same chr, same str, ok gx order, with no intervening gene in the middle and more distant than 100kb = intrachromosomal (mandatory)
- n3 is the number of chimeric transcripts to make from non overlapping genes on same chr, same str, ko gx order = inverted (mandatory)
- n4 is the number of chimeric transcripts to make from non overlapping genes on same chr, diff str = interstrand (mandatory)
- n5 is the number of chimeric transcripts to make from non overlapping genes on diff chr = interchromosomal (mandatory)
- list_gn_bt.txt is a list of gene biotypes to consider from the annotation (optional)

Will produce in the working directory:
- a single fasta file with n1+n2+n3+n4+n5 chimeric transcripts randomly selected from the gene pairs of the 5 classes
- an Aux directory with some intermediate files

NOTE1: cannot be run twice in the same directory without loosing previous outputs since uses fixed names for outputs
NOTE2: needs gem-retriever binay of the same gemtools version as the one used to generate gem genome index used as input
