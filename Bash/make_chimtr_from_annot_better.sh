#!/bin/bash

# make_chimtr_from_annot_better.sh
##################################
# Given as input:
#################
# - a gene annotation (1st arg, gtf file with at least exon and gene rows and where the gene id is in column 10 and transcript id is in column 12), 
# - a genome index produced by gem (2nd argument)
# - 5 numbers (3th to 7th arguments) for the number of chimeric transcripts to make from spliced transcripts of gene pairs that are
#   - not overlapping, on same chr, same str, ok gx order, with no gene in the middle and less distant than 100kb = readthrough (n1)
#   - not overlapping, on same chr, same str, ok gx order, with no gene in the middle and more distant than 100kb = intrachromosomal (n2)
#   - not overlapping, on same chr, same str, ko gx order = inverted (n3)
#   - not overlapping, on same chr, diff str = interstrand (n4)
#   - diff chr = interchromosomal (n5)
# - an optional 1 colum file with the gene biotypes to consider for generating the chimeric transcripts (but not for defining read through) (8th argument) (if not provided then no filter is done)
# Provides as output in the working directory:
##############################################
# - a single fasta file with n1+n2+n3+n4+n5 chimeric transcripts randomly selected from the gene pairs of the 5 classes
# - an Aux directory with some intermediate files
# !!! be careful: the names of the 5 classes is hard-coded here and has to match what is produced by the scripts !!!
# !!! make_gene_pairs_from_annot_better.awk and make_samechrstrok_gnpairs_from_annot.awk !!!
# !!! be careful: cannot be run twice in the same directory without loosing previous outputs since uses fixed names for outputs !!!

# Improvements:
###############
# - change the class names to match what is given in chimpipe's output = readthrough, intrachromosomal, inverted, interstrand, interchromosomal
#   and not the very long names we have here, but note that this depends on the nomenclature of some files produced by 2 other awk scripts 
#   make_gene_pairs_from_annot_better.awk and make_samechrstrok_gnpairs_from_annot.awk, so need to be careful
# - have a max threshold of distance between gene middles (like 100kb), that could be provided as input
#   to consider a gene pair as read through. For the moment it is always 100kb. This is now possible since the script make_samechrstrok_gnpairs_from_annot.awk
#   have this parameter


# example
#########
# cd /no_backup/rg/sdjebali/Chimeras/Benchmark/Data/Make_Chimeras_from_Annot/Test
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
# genome=/users/rg/projects/references/Genome/H.sapiens/hg19/gemtools1.7.1-i3/Homo_sapiens.GRCh37.chromosomes.chr.M.gem
# time make_chimtr_from_annot.sh $annot $genome 50 50 50 50 50 ../wanted_gnbt.txt 2> make_chimtr_from_annot.err
# about 30 minutes

# inputs
# protein_coding
# 1 (1 fields)

# Check if the arguments are correct
####################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ] || [ ! -n "$7" ]
then
    echo "" >&2
    echo "Usage: make_chimtr_from_annot_better.sh annot.gtf genome.gem n1 n2 n3 n4 n5 [list_gn_bt.txt]" >&2
    echo "where" >&2
    echo "- annot.gtf is a gene annotation (with at least exon and gene rows and where gene id and transcript id are in column 10 and 12) (mandatory)" >&2
    echo "- genome.gem is a genome index produced by gem (mandatory)" >&2
    echo "- n1 is the number of chimeric transcripts to make between non overlapping genes that are on the same chr, same str, ok gx order, with no intervening gene in the middle, and less distant than 100kb = read through (mandatory)" >&2
    echo "- n2 is the number of chimeric transcripts to make from non overlapping genes on same chr, same str, ok gx order, with no intervening gene in the middle and more distant than 100kb = intrachromosomal (mandatory)" >&2
    echo "- n3 is the number of chimeric transcripts to make from non overlapping genes on same chr, same str, ko gx order = inverted (mandatory)" >&2
    echo "- n4 is the number of chimeric transcripts to make from non overlapping genes on same chr, diff str = interstrand (mandatory)" >&2
    echo "- n5 is the number of chimeric transcripts to make from non overlapping genes on diff chr = interchromosomal (mandatory)" >&2
    echo "- list_gn_bt.txt is a list of gene biotypes to consider from the annotation (optional)" >&2
    echo "Will produce in the working directory:" >&2
    echo "- a single fasta file with n1+n2+n3+n4+n5 chimeric transcripts randomly selected from the gene pairs of the 5 classes" >&2 
    echo "- an Aux directory with some intermediate files" >&2
    echo "NOTE1: cannot be run twice in the same directory without loosing previous outputs since uses fixed names for outputs" >&2 
    echo "NOTE2: needs gem-retriever binay of the same gemtools version as the one used to generate gem genome index used as input" >&2 
    echo "" >&2
    exit 1
fi


# Assign variables
##################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
anntmp=`basename $annot`
anntmp2=${anntmp%.gtf}
annbase=${anntmp2%.gff}
genome=$2
n1=$3
n2=$4
n3=$5
n4=$6
n5=$7

# Programs
##########
GFF2GFF=$rootDir/../Awk/gff2gff.awk
MAKEGNPAIR1=$rootDir/../Awk/make_samechrstrok_gnpairs_from_annot.awk
MAKEGNPAIR2=$rootDir/../Awk/make_gene_pairs_from_annot_better.awk
MAKEINTRONS=$rootDir/../Awk/make_introns.awk
RECONSTRUCT=$rootDir/../Awk/reconstruct_chimtr.awk
RETRIEVER=$rootDir/../bin/gem-retriever 


# Make a file with the 5 last categories and the numbers of chimeric transcripts wanted for them (useful for the loops)
#######################################################################################################################
# !!! be careful: the names of the 3 last classes have to match what is produced by the script make_gene_pairs_from_annot_better.awk !!!
printf nonoverlap_samechr_samestr_okgxorder_read-through"\t"$n1"\n" > gnpaircategory_nbwanted.tsv
printf nonoverlap_samechr_samestr_okgxorder_not_read-through"\t"$n2"\n" >> gnpaircategory_nbwanted.tsv
printf nonoverlap_samechr_samestr_kogxorder"\t"$n3"\n" >> gnpaircategory_nbwanted.tsv
printf nonoverlap_samechr_diffstr"\t"$n4"\n" >> gnpaircategory_nbwanted.tsv
printf diffchr"\t"$n5"\n" >> gnpaircategory_nbwanted.tsv
# nonoverlap_samechr_samestr_okgxorder_read-through 50
# 5 (1 fields)


# Create an Aux directory for the aux files that we want to keep
################################################################
mkdir -p Aux


# Start with the 7 actions
##########################
# 1. filter the annotation to only get the rows corresponding to a list of biotypes provided by the user
########################################################################################################
# and corresponding to genes with at least one spliced transcript. Sort the file by transcript and then beg and end
###################################################################################################################
# for many downstream steps
###########################
# !!! here the transcript type is now considered and should be in the list of biotypes the user wants !!!
# !!! be careful some genes could be retained but without any transcript in the file so need to adjust for this !!!
if [ -n "$8" ] 
then
awk -v fileRef1=$8 -v fileRef2=$annot 'BEGIN{OFS="\t"; while (getline < fileRef1 >0){ok["\""$1"\""]=1} while (getline < fileRef2 >0){if($3=="exon"){split($12,a,"\""); nbex[a[2]]++}}} $3=="transcript"{n=0; split($10,a,"\""); split($12,b,"\""); split($0,c,"\t"); split(c[9],d,"; "); k=1; while(d[k]!=""){split(d[k],e," "); if(((e[1]=="gene_type")&&(ok[e[2]]==1))||((e[1]=="transcript_type")&&(ok[e[2]]==1))){n++} k++} if(n==2){print a[2], b[2], nbex[b[2]]}}' $annot > $annbase.filt.gnid_trid_nbex.tsv
awk -v fileRef=$annbase.filt.gnid_trid_nbex.tsv 'BEGIN{while (getline < fileRef >0){if($3>=2){ok[$1]=1; ok[$2]=1}}} {split($10,a,"\""); split($12,b,"\""); if((ok[a[2]]==1)&&(ok[b[2]]==1)){print}}' $annot | sort -k12,12 -k4,4n -k5,5n | awk -f $GFF2GFF > $annbase.filt.splicedgn.gtf
# awk -v fileRef1=$8 -v fileRef2=$annot 'BEGIN{while (getline < fileRef1 >0){ok["\""$1"\"\;"]=1;} while (getline < fileRef2 >0){if((ok2[$10]!=1)&&($3=="exon")){nbex[$10,$12]++; if(nbex[$10,$12]>=2){ok2[$10]=1;}}}} (ok2[$10]==1){ok3=0; ok4=0; k=9; while($k!=""){if($k=="gene_type"){if(ok[$(k+1)]==1){print}} k+=2}}' $annot | sort -k12,12 -k4,4n -k5,5n | awk -f $GFF2GFF > $annbase.filt.splicedgn.gtf
else
awk -v fileRef2=$annot 'BEGIN{while (getline < fileRef2 >0){if((ok2[$10]!=1)&&($3=="exon")){nbex[$10,$12]++; if(nbex[$10,$12]>=2){ok2[$10]=1;}}}} (ok2[$10]==1)' $annot | sort -k12,12 -k4,4n -k5,5n | awk -f $GFF2GFF > $annbase.filt.splicedgn.gtf
fi
# chrX	HAVANA	gene	99883667	99894988	.	-	.	gene_id "ENSG00000000003.10"; transcript_id "ENSG00000000003.10"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "TSPAN6"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "TSPAN6"; level 2; havana_gene "OTTHUMG00000022002.1";
# 242 (26 fields)
# 16504 (28 fields)
# 130007 (30 fields)
# 53809 (32 fields)
# 495999 (34 fields)
# 336699 (36 fields)
# 384954 (38 fields)
# 360265 (40 fields)
# 508411 (42 fields)
# 94980 (44 fields)
# 9233 (46 fields)
# 848 (48 fields)
# 116 (50 fields) *** real	3m0.973s

# 2. Generate all the possible non overlapping gene pairs from the first 2 categories from the whole annotation and then split into two groups
##############################################################################################################################################
# with only the spliced gene pairs belonging to the biotypes wanted by the user
###############################################################################
# !!! read through or not is made from the complete annotation and not the biotypes chosen by the user !!!
awk '$3=="gene"' $annot | sort -k1,1 -k4,4n -k5,5n > $annbase.gn.sorted.gff
awk -v dist=100000 -v fileRef=$annbase.gn.sorted.gff -f  $MAKEGNPAIR1 $annbase.gn.sorted.gff 
# chr1  11869   14412   chr1    29554   31109   ENSG00000223972.4:ENSG00000243485.2     0       +       +       read-through
# 29440465 (11 fields)  *** 3 minutes and file is 3G so need to gzip it
gzip -f genepairs_nonoverlap_samechr_samestr_okgxorder.bedpe
zcat genepairs_nonoverlap_samechr_samestr_okgxorder.bedpe.gz | awk -v fileRef=$annbase.filt.splicedgn.gtf 'BEGIN{OFS="\t"; while (getline < fileRef >0){split($10,a,"\""); ok[a[2]]=1}} {split($7,a,":"); if((ok[a[1]]==1)&&(ok[a[2]]==1)){print a[1]"\t"a[2] > "genepairs_nonoverlap_samechr_samestr_okgxorder_"$NF".tsv"}}' 
# ENSG00000187634.6     ENSG00000187961.9
# 8297 (2 fields)
# ENSG00000269308.1     ENSG00000187634.6
# 4702399 (2 fields) *** 2 minutes or so

# 3. generate all the possible non overlapping spliced gene pairs from the last 3 categories and from the biotypes wanted by the user
######################################################################################################################################
#    1) genepairs_nonoverlap_samechr_samestr_kogxorder.tsv 
#    2) genepairs_nonoverlap_samechr_diffstr.tsv 
#    3) genepairs_diffchr.tsv
awk -v fileRef=$annbase.filt.splicedgn.gtf -f $MAKEGNPAIR2 $annbase.filt.splicedgn.gtf
rm genepairs_nonoverlap_samechr_samestr_okgxorder.tsv
# real	9m46.792s  *** 2 same chr files of same size (eg v19 pcg 165M) and the diff chr one much bigger (eg v19 pcg 5.9G)

# 4. gzip all the files produced until now to save space
########################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
gzip -f genepairs_$cat.tsv
done
# ENSG00000187634.6       ENSG00000269308.1
# 730M in total when zipped, instead of 6.4G

# 5. randomly sample n1 elements from 1, n2 elements from 2, n3 elements from 3, n4 elements from 4 and n5 elements from 5, as the user wants
#############################################################################################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
zcat genepairs_$cat.tsv.gz | shuf | head -n $nb > genepairs_$cat\_$nb\_sampled.tsv
done
# ENSG00000130032.11    ENSG00000147378.7
# 50 (2 fields)
# ENSG00000188368.5	ENSG00000183207.8
# 150 (2 fields)
# ENSG00000123836.10	ENSG00000203857.5
# 50 (2 fields)
# ENSG00000115946.3	ENSG00000144161.8
# 50 (2 fields)
# ENSG00000116871.11	ENSG00000165507.8
# 50 (2 fields) *** real	1m32.214s

# 6. for each randomly selected gene pair, randomly pick 1 spliced transcript for each gene and then 1 donor from 1st tr and
############################################################################################################################
#    one acceptor from 2nd transcript, also in a random way
###########################################################
# Make gene pair file but with spliced transcript list information
###################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=$annbase.filt.splicedgn.gtf 'BEGIN{OFS="\t"; while (getline < fileRef >0){if($3=="exon"){split($10,a,"\""); split($12,b,"\""); nbex[a[2],b[2]]++; if(nbex[a[2],b[2]]==2){trlist[a[2]]=(trlist[a[2]])(b[2])(",")}}}} {print $1, $2, trlist[$1], trlist[$2]}' genepairs_$cat\_$nb\_sampled.tsv > Aux/genepairs_$cat\_$nb\_sampled_splicedtrlist_eachgn.tsv
done
# ENSG00000130032.11    ENSG00000147378.7       ENST00000370353.3,ENST00000370354.1,ENST00000448324.1,ENST00000448726.1,ENST00000538575.1,      ENST00000370350.3,ENST00000417321.1,
# 50 (4 fields)
# ENSG00000169313.9     ENSG00000162244.6       ENST00000302632.3,ENST00000468596.1,    ENST00000294189.6,ENST00000466397.1,ENST00000475248.1,ENST00000479017.1,ENST00000480306.1,ENST00000481629.1,ENST00000486565.1,ENST00000492277.1,ENST00000495383.1,
# 100 (4 fields)
# ENSG00000144476.5     ENSG00000198522.9       ENST00000272928.3,ENST00000447924.1,    ENST00000264718.3,ENST00000407583.3,ENST00000424214.1,ENST00000436280.1,ENST00000458167.2,ENST00000461249.1,ENST00000478484.1,ENST00000481754.1,ENST00000503738.1,ENST00000515877.1,ENST00000610189.1,
# 50 (4 fields)
# ENSG00000216937.7     ENSG00000204172.7       ENST00000277657.6,ENST00000362006.5,ENST00000435402.1,ENST00000471905.1,ENST00000476558.1,ENST00000489718.1,ENST00000493685.1,ENST00000535327.1,ENST00000537047.1,ENST00000539197.1,ENST00000545067.1,     ENST00000355232.3,ENST00000413193.2,ENST00000452145.2,
# 50 (4 fields)
# ENSG00000124374.8     ENSG00000126803.8       ENST00000244221.8,      ENST00000394709.1,ENST00000554883.1,
# 50 (4 fields)  *** real	0m19.871s  
# Then the gn pair files but with 1 transcript randomly selected from each list
###############################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
cat Aux/genepairs_$cat\_$nb\_sampled_splicedtrlist_eachgn.tsv | while read g1 g2 trl1 trl2
do
tr1=`echo $trl1 | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
tr2=`echo $trl2 | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
echo $g1 $g2 $tr1 $tr2
done > genepairs_$cat\_$nb\_sampled_rndtr_eachgn.tsv
done
# ENSG00000130032.11 ENSG00000147378.7 ENST00000370354.1 ENST00000370350.3
# 50 (4 fields)
# ENSG00000169313.9 ENSG00000162244.6 ENST00000302632.3 ENST00000479017.1
# 100 (4 fields)
# ENSG00000144476.5 ENSG00000198522.9 ENST00000447924.1 ENST00000461249.1
# 50 (4 fields)
# ENSG00000216937.7 ENSG00000204172.7 ENST00000476558.1 ENST00000452145.2
# 50 (4 fields)
# ENSG00000124374.8 ENSG00000126803.8 ENST00000244221.8 ENST00000554883.1
# 50 (4 fields) *** real	0m3.427s
# Then the donor list of the first tr and the acceptor list of the second tr in addition to these 4 columns
###########################################################################################################
# !!! note: they have to exist since here we only considered spliced tr !!!
# first make all introns from the genes which biotypes are provided by the user
###############################################################################
# !!! note: the annot file is already sorted by transcript id and then beg and end so no need to redo here !!!
awk -v fldgn=10 -v fldtr=12 -f $MAKEINTRONS $annbase.filt.splicedgn.gtf > $annbase.filt.splicedgn.introns.gff
# chr7	HAVANA	intron	127228620	127229136	.	+	.	gene_id "ENSG00000004059.6"; transcript_id "ENST00000000233.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "ARF5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "ARF5-001"; level 2; tag "basic"; tag "appris_principal"; tag "CCDS"; ccdsid "CCDS34745.1"; havana_gene "OTTHUMG00000023246.5"; havana_transcript "OTTHUMT00000059567.2";
# 155 (28 fields)
# 250303 (30 fields)
# 132014 (32 fields)
# 134495 (34 fields)
# 146352 (36 fields)
# 216592 (38 fields)
# 40781 (40 fields)
# 4018 (42 fields)
# 373 (44 fields)
# 53 (46 fields)  *** real	0m9.160s
# then make the list of donors and acceptors of the 1st and 2nd tr respectively
###############################################################################
# !!! be careful introns boundaries are not exonic while here we want donor and acceptor coord to be exonic !!!
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=$annbase.filt.splicedgn.introns.gff 'BEGIN{OFS="\t"; while (getline < fileRef >0){split($12,a,"\""); chr[a[2]]=$1; str[a[2]]=$7; donlist[a[2]]=(donlist[a[2]])($7=="+" ? ($4-1) : ($5+1))(","); acclist[a[2]]=(acclist[a[2]])($7=="+" ? ($5+1) : ($4-1))(",");}} {print $1, $2, $3, $4, chr[$3], str[$3], chr[$4], str[$4], (donlist[$3]!="" ? donlist[$3] : "."), (acclist[$4]!="" ? acclist[$4] : ".")}' genepairs_$cat\_$nb\_sampled_rndtr_eachgn.tsv > Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_donlist_acclist.tsv
done
# ENSG00000130032.11    ENSG00000147378.7       ENST00000370354.1       ENST00000370350.3       chrX    +       chrX    +       150864016,150867293,150868773,  150885744,150889867,150890375,150891100,
# 50 (10 fields)
# ENSG00000169313.9     ENSG00000162244.6       ENST00000302632.3       ENST00000479017.1       chr3    -       chr3    -       151058384,151102480,    52028142,52029122,52029485,
# 100 (10 fields)
# ENSG00000144476.5     ENSG00000198522.9       ENST00000447924.1       ENST00000461249.1       chr2    +       chr2    +       237476484,      27852726,27855500,27857713,27858007,27861077,27861752,27862910,27864107,27865296,27870703,
# 50 (10 fields)
# ENSG00000216937.7     ENSG00000204172.7       ENST00000476558.1       ENST00000452145.2       chr10   +       chr10   -       32735325,32740849,32742364,32745262,32751631,32751977,32760158,32761470,32762951,32806903,32807433,32832313,32833229,32854548,32856819,32860821,   47193533,47194230,47197573,47200261,47207846,47210555,47212922,
# 50 (10 fields)
# ENSG00000124374.8     ENSG00000126803.8       ENST00000244221.8       ENST00000554883.1       chr2    -       chr14   +       71416975,71429582,71454058,     65008055,
# 50 (10 fields)  *** real	0m15.431s  
# finally randomly select one donor and one acceptor from those lists
#####################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
cat Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_donlist_acclist.tsv | while read a b c d e f g h i j
do
don=`echo $i | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
acc=`echo $j | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
printf $a"\t"$b"\t"$c"\t"$d"\t"$e"\t"$f"\t"$g"\t"$h"\t"$don"\t"$acc"\n"
done > genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv
done
# ENSG00000130032.11    ENSG00000147378.7       ENST00000370354.1       ENST00000370350.3       chrX    +       chrX    +       150868773       150890375
# 50 (10 fields)
# ENSG00000169313.9     ENSG00000162244.6       ENST00000302632.3       ENST00000479017.1       chr3    -       chr3    -       151058384       52029122
# 100 (10 fields)
# ENSG00000144476.5     ENSG00000198522.9       ENST00000447924.1       ENST00000461249.1       chr2    +       chr2    +       237476484       27855500
# 50 (10 fields)
# ENSG00000216937.7     ENSG00000204172.7       ENST00000476558.1       ENST00000452145.2       chr10   +       chr10   -       32742364        47197573
# 50 (10 fields)
# ENSG00000124374.8     ENSG00000126803.8       ENST00000244221.8       ENST00000554883.1       chr2    -       chr14   +       71454058        65008055
# 50 (10 fields)   *** real	0m4.671s


# 7. make the complete fasta sequence of all these transcripts while recording in each header all info about genes, 
###################################################################################################################
#    transcripts and donor and acceptor coordinates
###################################################
# !!! note: the annot file is already sorted by transcript id and then beg and end so no need to redo here !!!
# - For the donor side
#   * if tr is on + take exons from 1st to the one with gend=don coord (in THIS order)
#     otherwise take exons from last to the one with gbeg=don coord (in THIS order)
# - For the acc side
#   * if tr is on + take exons from the one with gbeg=acc coord to the last (in THIS order)
#     otherwise take exons from the one with gend=acc coord to the 1st (in THIS order)
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=$annbase.filt.splicedgn.gtf -f $RECONSTRUCT genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv > genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv
done
# ENSG00000130032.11    ENSG00000147378.7       ENST00000370354.1       ENST00000370350.3       chrX    +       chrX    +       150868773       150890375       chrX_150863951_150864016_+,chrX_150867256_150867293_+,chrX_150868444_150868773_+,  chrX_150890375_150890453_+,chrX_150891100_150891666_+,
# 50 (12 fields)
# ENSG00000169313.9     ENSG00000162244.6       ENST00000302632.3       ENST00000479017.1       chr3    -       chr3    -       151058384       52029122        chr3_151102480_151102600_-,chr3_151058384_151058548_-,     chr3_52029058_52029122_-,chr3_52027645_52028142_-,
# 100 (12 fields)
# ENSG00000144476.5     ENSG00000198522.9       ENST00000447924.1       ENST00000461249.1       chr2    +       chr2    +       237476484       27855500        chr2_237476430_237476484_+,        chr2_27855500_27855537_+,chr2_27857713_27857791_+,chr2_27858007_27858101_+,chr2_27861077_27861122_+,chr2_27861752_27861898_+,chr2_27862910_27862992_+,chr2_27864107_27864146_+,chr2_27865296_27865386_+,chr2_27870703_27870761_+,
# 50 (12 fields)
# ENSG00000216937.7     ENSG00000204172.7       ENST00000476558.1       ENST00000452145.2       chr10   +       chr10   -       32742364        47197573        chr10_32735068_32735325_+,chr10_32740520_32740849_+,chr10_32742272_32742364_+,     chr10_47197538_47197573_-,chr10_47194179_47194230_-,chr10_47191844_47193533_-,
# 50 (12 fields)
# ENSG00000124374.8     ENSG00000126803.8       ENST00000244221.8       ENST00000554883.1       chr2    -       chr14   +       71454058        65008055        chr2_71454058_71454213_-,  chr14_65008055_65008293_+,
# 50 (12 fields)  *** real	0m18.986s

# Extract the exons using gem retriever for all categories
##########################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk 'BEGIN{OFS="\t"} {split($11,a,","); k=1; while(a[k]!=""){split(a[k],b,"_"); print b[1], b[4], b[2], b[3]; k++} split($12,a,","); k=1; while(a[k]!=""){split(a[k],b,"_"); print b[1], b[4], b[2], b[3]; k++}}' genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv 
done | sort | uniq > genepairs_distinct_exon_coord.tsv 
# chr10 +       100003848       100004654
# 1940 (4 fields)  *** real	0m0.056s 
cat genepairs_distinct_exon_coord.tsv | $RETRIEVER $genome > genepairs_distinct_exon_coord.seq
# AGAGAAAGCGGTTGGAAGCCAAGCAACGGGAAGACATCTGGGAAGGCAGAGACCAGTCTACAGTTTGAACATCACTCAATGAAAGGGATAATTCCATGAATCAGAAAATGTTTCCATAGCCTTCAGATAAGATGATCCTTCCAGAGCTCTATGTACATGCAGATGTGCATGTTAAAGAGATAAAGTGATCGAGACAAGGACTGACTGGGTATAGAAGGAAGACAGACTCCTGTCTTCACTCCTAAATGCAGTTCTTTGGAATCACCCTACTGTGGTGGGCGTAGTAGGGAGCCATCAGCTAGGAAGAAACGTGGGAGATGTGAATTCCAAGAGTTGCCTGGACAGGGCAAGTCATGTTAGCGTGGGTCACACTTCCAAGATATTTAAAGCAAATACAAAACAGAACAGAGGATTCAAACCGCAAGTATGGGAGATTTAGGCCCTGCAGAGGCAGACCATTCCTTAGTATCTCACAAAGCAGAGTAATACTGGAGGCAGAGTAGGGGGTGGTTGGAGAGCAGTTAGTACAAAGAGGCAGAACAGTGTCTGGTTTACTTGGCATACACAGAATCTGCACTGCCGGTTCCAGAACTGCAAAGTTGGTGAACTACAGGAGATGTGGGTATTTAGACTCCAAAGTTTATACTGAGCTCAGTGCCTGGGACCGCTCCAGCTGCACTGCAGCCAGGGGGGATACCAGCTTCATTCCTCAGAGTAAACCAACCTTGGGAGCTGTGTGCCAGTTCCAGTGTAAGCAGTAATACCATGTTGTAAGAGGGAACATTAAAGCCCATTTTATGACATATA
# 1940 (1 fields)  *** real	0m19.801s
paste genepairs_distinct_exon_coord.tsv genepairs_distinct_exon_coord.seq | awk '{print $1"_"$3"_"$4"_"$2, $5}' > genepairs_distinct_exon_coord_seq.txt
# chr10_100003848_100004654_+ AGAGAAAGCGGTTGGAAGCCAAGCAACGGGAAGACATCTGGGAAGGCAGAGACCAGTCTACAGTTTGAACATCACTCAATGAAAGGGATAATTCCATGAATCAGAAAATGTTTCCATAGCCTTCAGATAAGATGATCCTTCCAGAGCTCTATGTACATGCAGATGTGCATGTTAAAGAGATAAAGTGATCGAGACAAGGACTGACTGGGTATAGAAGGAAGACAGACTCCTGTCTTCACTCCTAAATGCAGTTCTTTGGAATCACCCTACTGTGGTGGGCGTAGTAGGGAGCCATCAGCTAGGAAGAAACGTGGGAGATGTGAATTCCAAGAGTTGCCTGGACAGGGCAAGTCATGTTAGCGTGGGTCACACTTCCAAGATATTTAAAGCAAATACAAAACAGAACAGAGGATTCAAACCGCAAGTATGGGAGATTTAGGCCCTGCAGAGGCAGACCATTCCTTAGTATCTCACAAAGCAGAGTAATACTGGAGGCAGAGTAGGGGGTGGTTGGAGAGCAGTTAGTACAAAGAGGCAGAACAGTGTCTGGTTTACTTGGCATACACAGAATCTGCACTGCCGGTTCCAGAACTGCAAAGTTGGTGAACTACAGGAGATGTGGGTATTTAGACTCCAAAGTTTATACTGAGCTCAGTGCCTGGGACCGCTCCAGCTGCACTGCAGCCAGGGGGGATACCAGCTTCATTCCTCAGAGTAAACCAACCTTGGGAGCTGTGTGCCAGTTCCAGTGTAAGCAGTAATACCATGTTGTAAGAGGGAACATTAAAGCCCATTTTATGACATATA
# 1940 (2 fields)  *** real	0m0.063s
# And then concatenate the exons in the order they appear in the donor and acceptor list
########################################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=genepairs_distinct_exon_coord_seq.txt 'BEGIN{OFS="\t"; while (getline < fileRef >0){seqex[$1]=$2}} {donseq=""; split($11,a,","); k=1; while(a[k]!=""){donseq=(donseq)(seqex[a[k]]); k++} accseq=""; split($12,a,","); k=1; while(a[k]!=""){accseq=(accseq)(seqex[a[k]]); k++} print $0, donseq, accseq}' genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv > Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist_donseq_accseq.tsv
done
# ENSG00000130032.11    ENSG00000147378.7       ENST00000370354.1       ENST00000370350.3       chrX    +       chrX    +       150868773       150890375       chrX_150863951_150864016_+,chrX_150867256_150867293_+,chrX_150868444_150868773_+,  chrX_150890375_150890453_+,chrX_150891100_150891666_+,  GGC...AAA       CAAC...AAA
# 50 (14 fields)
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337	49502580	chr19_42812926_42814337_+,	chr19_49502580_49502621_+,	CCTAG...CCAG	CCACAACCAAAGTCCCGGAGATCCGTGATGTAACAAGGATTG
# 150 (14 fields)
# ENSG00000123836.10	ENSG00000203857.5	ENST00000411990.2	ENST00000235547.6	chr1	+	chr1	+	207235423	120056457	chr1_207222801_207222962_+,chr1_207228046_207228147_+,chr1_207235298_207235423_+,	chr1_120056457_120057681_+,	TTT...AAG	GTAC...TGTG
# 50 (14 fields)
# ENSG00000115946.3	ENSG00000144161.8	ENST00000488728.1	ENST00000409573.2	chr2	+	chr2	-	68400549	112991813	chr2_68400239_68400549_+,	chr2_112991697_112991813_-,chr2_112990837_112990948_-,chr2_112989415_112989524_-,chr2_112988480_112988527_-,chr2_112969102_112974045_-,	TAAG...CTTGG	GATG...TAAA
# 50 (14 fields)
# ENSG00000116871.11	ENSG00000165507.8	ENST00000373151.2	ENST00000496638.1	chr1	+	chr10	-	36644197	45467123	chr1_36621803_36622064_+,chr1_36636572_36636916_+,chr1_36637114_36637182_+,chr1_36638065_36638228_+,chr1_36638965_36639079_+,chr1_36640499_36640609_+,chr1_36641800_36642182_+,chr1_36642298_36642443_+,chr1_36643474_36643802_+,chr1_36644020_36644197_+,	chr10_45466429_45467123_-,	GGG...ACAA	GAGT...GAGG
# 50 (14 fields) *** real	0m0.057s

# Make a proper fasta file with all info in header and then a nt sequence correctly formatted (60nt per row)
############################################################################################################
# !!! include the coord of the junction in transcript space !!!
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v cat=$cat '{print ">"$3"-"$4, cat, $1, $2, $5, $6, $7, $8, $9, $10, length($13); s=($13)($14); n=length(s); n2=int(n/60); for(i=0; i<=(n2-1); i++){print substr(s,i*60+1,60)} if(n>n2*60){print substr(s,n2*60+1,n-n2*60)}}' Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist_donseq_accseq.tsv > genepairs_$cat\_$nb.fasta
done
# >ENST00000370354.1-ENST00000370350.3 nonoverlap_samechr_samestr_okgxorder_read-through ENSG00000130032.11 ENSG00000147378.7 chrX + chrX + 150868773 150890375 434
# 1350 (1 fields)
# 50 (11 fields)
# >ENST00000354289.4-ENST00000476862.1 nonoverlap_samechr_samestr_okgxorder ENSG00000105835.7 ENSG00000085563.10 chr7 - chr7 - 105917469 87230394 251
# 4553 (1 fields)
# 150 (11 fields)
# >ENST00000514078.1-ENST00000504807.1 nonoverlap_samechr_samestr_kogxorder ENSG00000169777.5 ENSG00000183258.7 chr5 - chr5 - 9712284 176939877 207
# 1843 (1 fields)
# 50 (11 fields)
# >ENST00000422000.1-ENST00000555554.1 nonoverlap_samechr_diffstr ENSG00000120647.5 ENSG00000118307.14 chr12 + chr12 - 514730 25267804 191
# 2581 (1 fields)
# 50 (11 fields)
# >ENST00000407793.2-ENST00000546806.1 diffchr ENSG00000013573.12 ENSG00000196511.9 chr12 + chr7 - 31247568 144380071 1665
# 1689 (1 fields)
# 50 (11 fields) *** real	0m0.065s  *** seems to be fine  (here output from another test)
# Makes a single fasta file with all these transcripts (the chimeric transcript class is in the header)
#######################################################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
cat genepairs_$cat\_$nb.fasta
done > genepairs_allcat.fasta
# >ENST00000370354.1-ENST00000370350.3 nonoverlap_samechr_samestr_okgxorder_read-through ENSG00000130032.11 ENSG00000147378.7 chrX + chrX + 150868773 150890375 434
# 8199 (1 fields)
# 300 (11 fields) *** final output with all needed information (also from other test than the previous steps)


# Clean
#######
rm $annbase.filt.splicedgn.gtf
rm $annbase.filt.gnid_trid_nbex.tsv
rm $annbase.filt.splicedgn.introns.gff
rm $annbase.gn.sorted.gff
rm genepairs_nonoverlap_samechr_samestr_okgxorder.bedpe.gz
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
rm genepairs_$cat\_$nb\_sampled.tsv
rm genepairs_$cat.tsv.gz
rm genepairs_$cat\_$nb\_sampled_rndtr_eachgn.tsv
rm genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv
rm genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv
rm genepairs_$cat\_$nb.fasta
done
rm genepairs_distinct_exon_coord.tsv genepairs_distinct_exon_coord.seq
