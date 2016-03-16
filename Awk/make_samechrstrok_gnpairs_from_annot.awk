
# ~/Awk/make_samechrstrok_gnpairs_from_annot.awk

# This script takes as input a gene annotation file with !!!at least gene!!! rows, !!!sorted!!! by chr, start, end coord
# and produces a bedpe file of non overlapping gene pairs that are on the same chr, same strand, and expected gx order 
# with an additional flag saying whether they are read through or not. Here read-through means that they is no intervening gene
# between the two genes AND the two genes are less distant than 100kb (middle to middle) (can be specified with -v dist=x). 
# In general for this it is better to use all long genes, not only pcg genes

# example1
# cd /no_backup/rg/sdjebali/Chimeras/Benchmark/Data/Make_Chimeras_from_Annot/Test
# awk -v dist=150000 -v fileRef=gencv19_long_chr22_sorted_100firstrows.gff -f ~sdjebali/Awk/make_samechrstrok_gnpairs_from_annot.awk gencv19_long_chr22_sorted_100firstrows.gff 
# instantly, see input and output below
# example2
# cd /no_backup/rg/sdjebali/Chimeras/Benchmark/Data/Make_Chimeras_from_Annot/Test
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
# time awk -v fileRef=$annot -f ~sdjebali/Awk/make_samechrstrok_gnpairs_from_annot.awk $annot
# chr1  11869   14412   chr1    29554   31109   ENSG00000223972.4:ENSG00000243485.2     0       +       +       read-through
# 29440465 (11 fields) *** real    3m22.192s *** quite long

# input (which must contain gene rows)
# chr22 HAVANA  exon    16062157        16062316        .       +       .       gene_id "ENSG00000233866.1"; transcript_id "ENST00000424770.1"; gene_type "lincRNA"; gene_status "KNOWN"; gene_name "LA16c-4G1.3"; transcript_type "lincRNA"; transcript_status "KNOWN"; transcript_name "LA16c-4G1.3-001"; exon_number 1; exon_id "ENSE00001660730.1"; level 2; tag "basic"; havana_gene "OTTHUMG00000140195.1"; havana_transcript "OTTHUMT00000276574.1";
# 10 (28 fields)
# 10 (30 fields)
# 10 (32 fields)
# 31 (34 fields)
# 26 (36 fields)
# 10 (38 fields)
# 2 (40 fields)
# 1 (44 fields)

# output
# chr22 16062157        16063236        chr22   16122720        16123768        ENSG00000233866.1:ENSG00000215270.3     0       +       +       read-through
# 39 (11 fields)  


BEGIN{
    if(dist=="")
    {
	dist=100000;
    }
    OFS="\t";
    while (getline < fileRef >0) 
    {
	nbstr=2;
	str[1]="+";
	str[2]="-";
	if($3=="gene")
	{
	    seenchr[$1]++;
	    if(seenchr[$1]==1)
	    {
		nbchr++;
		chr[nbchr]=$1;
	    }
	}
    }
}

$3=="gene"{
    split($10,a,"\"");
    nb[$1,$7]++; 
    gene[$1,$7,nb[$1,$7]]=a[2]; 
    gbeg[$1,$7,nb[$1,$7]]=$4; 
    gend[$1,$7,nb[$1,$7]]=$5;
}


END{
    # for each chromosome k
    for(k=1; k<=nbchr; k++)
    {
	# for each strand i
	for(i=1; i<=nbstr; i++)
	{
	    # for each gene on chr k and strand i in our previous file
	    for(l=1; l<=(nb[chr[k],str[i]]-1); l++)
	    {
		# we look for another gene that is not equal to it nor for which the pair was already reported
		for(m=(l+1); m<=nb[chr[k],str[i]]; m++)
		{
		    # we also require that the two genes are not overlapping
		    if(foverlap(gbeg[chr[k],str[i],l],gend[chr[k],str[i],l], gbeg[chr[k],str[i],m], gend[chr[k],str[i],m])==0)
		    {
			mid1=(gbeg[chr[k],str[i],l]+gend[chr[k],str[i],l])/2; 
			mid2=(gbeg[chr[k],str[i],m]+gend[chr[k],str[i],m])/2; 
			# and further require that the two genes are in the expected genomic order, which means that the 1st one is the most 5' (careful of the strand)
			# here the 1st gene is on 5' so write 1st and then 2nd
			if(((i==1)&&(mid1<mid2))||((i==2)&&(mid1>mid2)))
			{
			    print chr[k], gbeg[chr[k],str[i],l], gend[chr[k],str[i],l], chr[k], gbeg[chr[k],str[i],m], gend[chr[k],str[i],m], gene[chr[k],str[i],l]":"gene[chr[k],str[i],m], 0, str[i], str[i], (((m==(l+1))&&((abs(mid2-mid1))<=dist)) ? "read-through" : "not_read-through") > "genepairs_nonoverlap_samechr_samestr_okgxorder.bedpe";
			}
			# here the 1st gene is on 3' so write 2nd and then 1st
			else
			{
			    if(((i==1)&&(mid1>mid2))||((i==2)&&(mid1<mid2)))
			    {
				print chr[k], gbeg[chr[k],str[i],m], gend[chr[k],str[i],m], chr[k], gbeg[chr[k],str[i],l], gend[chr[k],str[i],l], gene[chr[k],str[i],m]":"gene[chr[k],str[i],l], 0, str[i], str[i], (((m==(l+1))&&((abs(mid2-mid1))<=dist)) ? "read-through" : "not_read-through") > "genepairs_nonoverlap_samechr_samestr_okgxorder.bedpe";
			    }
			}
		    }
		}
	    }
	}
    }
}

function foverlap(beg1,end1,beg2,end2)
{
  return ((end1>=beg2)&&(beg1<=end2)) ? 1 : 0;
}

function abs(x)
{
    return (x>=0 ? x : -1*x);
}
