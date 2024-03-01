#try lamarkian again


#get reads
rsync -avtr /work/07707/dgallery/dom_shared/Lamarck/*.gz $SCRATCH/lamk/lamk_new

>gzips
for file in *.gz
do echo "gunzip ${file}" >>gzips
done
ls6_launcher_creator.py -n gzips -j gzips -q normal -N 4 -w 8 -a $allo -e $email -t 00:10:00
sbatch gzips.slurm

>concatLanes
for file in *L001_R1_001.fastq
 do LANE4_FQ=$file
 LANE5_FQ=${file/_L001_/_L002_}
 CAT_FQ=${file/_L001_R1_001.fastq/.fq}
 echo "cat $LANE4_FQ $LANE5_FQ > $CAT_FQ" >> concatLanes
done

#check the job file
#head concatLanes
#wc -l concatLanes

#launch and submit
ls6_launcher_creator.py -n concatLanes -j concatLanes -q normal -N 5 -w 10 -a $allo -e $email -t 01:00:00
sbatch concatLanes.slurm


######################
### TRIM ADAPTERS ###
######################

>trimse
for file in *.fq
do echo "cutadapt \
-a GATCGGAAGAGCA \
-a AGATCGGAAGAGC \
--minimum-length 20 \
-q 20 \
-o ${file/.fq/}.trim \
$file > ${file/.fq}_trimlog.txt" >> trimse
done

ls6_launcher_creator.py -n trimse -j trimse -q normal -N 5 -w 5 -a $allo -e $email -t 00:50:00

export REF=/work/05410/eabbott/ls6/genomes/Amil_v2.01/Amil.v2.01.chrs.fasta

>mapse
for file in *.trim
do echo "bowtie2 -x $REF -U $file --local -p 4 -S ${file/.trim/}.sam">> mapse
done

ls6_launcher_creator.py -n mapse -j mapse -q normal -N 5 -w 6 -a $allo -e $email -t 01:00:00 

>s2b
for file in *.sam
do echo "samtools sort -O bam -o ${file/.sam}.bam ${file} && samtools index ${file/.sam}.bam" >>s2b
done
ls6_launcher_creator.py -n s2b -j s2b -q normal -N 2 -w 6 -a $allo -e $email -t 00:45:00

conda activate wgs

########################
###quality assessment###
########################
#removing bams with log(coverage)<3SD
# also calculating minimum number of individuals(MI) a locus must be seen in (genotyping rate cutoff)
# if you are mapping to a real genome, replace chr1 on line 146 by a name of a nice long contig (a few megabases). Look this up in a header of any of the *.sam files.
#find the biggest contig
cat Amil.v2.01.chrs.fasta | bioawk -c fastx '{ print length($seq), $name }' | sort -k1,1rn | head -1

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 1000" 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo "ls *.bam > bams && angsd -b bams -r chr1 -GL 1 $FILTERS $TODO -P 1 -out dd && Rscript ~/bin/plotQC.R prefix=dd >qualranks">a0
ls6_launcher_creator.py -j a0 -n a0 -a $allo -e $email -t 00:30:00 -N 1 -w 6 -q normal 

#scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/lamk/lamkn/testang/dd\*.gz .

scp eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/lamk/lamkn/dd.pdf .


######################################
############# GET COUNTS #############
######################################

MY_GFF="/work/05410/eabbott/ls6/genomes/Amil_v2.01/Amil.coding.gff3"; GENE_ID="ID"

#choose the genome you're working with

#run featurecounts
echo "featureCounts -a $MY_GFF -t gene -g $GENE_ID -o feature_counts_out.tsv -T 48 --primary *.bam" > runFeatureCounts
ls6_launcher_creator.py -n runFeatureCounts -j runFeatureCounts -t 00:15:00 -q development -a $allo -e $email -N 2 -w 2

#scp to local
/Users/evelynabbott/Dropbox/project/Projects/lamk/lamk_new
scp -r eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/lamk/lamkn/feature_counts_out.tsv .


#######################################
######### GET WINDOW COUNTS ###########
#######################################

#window generation is shown in picoMethyl_data_processing_pipeline.txt
#window files to use here are:
geneWindowFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/geneBoundaries.bed
exonWindowFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/cdsBoundaries.bed
promoterWindowFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/promoterBoundaries.bed
tssWindowFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/tssBoundaries.bed
window1KbFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/windowBoundaries_1kb.bed


#run bedtools
echo "bedtools multicov -bams *.bam -bed $geneWindowFile > mr_gene_counts.tsv0
bedtools multicov -bams *.bam -bed $exonWindowFile > mr_exon_counts.tsv0
bedtools multicov -bams *.bam -bed $promoterWindowFile > mr_promoter_counts.tsv0
bedtools multicov -bams *.bam -bed $tssWindowFile > mr_tss_counts.tsv0
bedtools multicov -bams *.bam -bed $window1KbFile > mr_1kbWindow_counts.tsv0" > runBedtools


module load bedtools
ls6_launcher_creator.py -n runBedtools -j runBedtools -q normal -N 1 -w 5 -a $allo -e $email -t 08:00:00


#format output
samples=`ls *.bam | sed 's/.bam//' | tr "\n" "\t"`
echo -e "chr\tstart\tend\tname\t$samples" > header.txt

for file in *_counts.tsv0
do echo "cat header.txt $file > ${file/_counts.tsv0/}_multicov.tsv"
done


#----- split by chromosome for parallelizing small windows -----#
window1KbFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/windowBoundaries_1kb.bed
window500bpFile=/work/05410/eabbott/ls6/genomes/Amil_v2.01/windowBoundaries_500bp.bed


WINDOW_FILE=/work/05410/eabbott/ls6/genomes/Amil_v2.01/geneBoundaries.bed
NAME=geneBoundaries.bed

#get chrs
cut -f 1 $WINDOW_FILE | sort | uniq | grep "^chr" > chrs.txt


#build chr beds
while read chr
do echo "${chr}..."
grep -w "^${chr}" $WINDOW_FILE > ${chr}_sub_${NAME}
done < chrs.txt
grep -v "^chr" $WINDOW_FILE > chrUn_sub_${NAME}

#run multicov for each
>paraMulticov
for file in *_sub_${NAME}
do echo "bedtools multicov -bams *.bam -bed $file > mr_${file}_counts.tsv0" >>paraMulticov
done

ls6_launcher_creator.py -n paraMulticov -j paraMulticov -q normal -N 1 -w 15 -a $allo -e $email -t 10:00:00
sbatch paraMulticov.slurm


#returns 15 *counts.tsv0 files. One for each chromosome and 1 for scaffolds
#add headers to each
samples=`ls *.bam | sed 's/.bam//' | tr "\n" "\t"`
echo -e "chr\tstart\tend\tname\t$samples" > ${NAME}.tsv
cat *_sub_${NAME}*tsv0  > ${NAME}.tsv










