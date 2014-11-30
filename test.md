http://creskolab.uoregon.edu/stacks/comp/process_radtags.php
https://genomicislands.wordpress.com/2013/11/20/how-to-perform-a-rad-seq-de-novo-assembly/
http://creskolab.uoregon.edu/stacks/pe_tut.php
http://www.simison.com/brian/Illumina_Rad_process_notes.html

————————————————————————————————————————
Task 1: We need to rename all the reads to the sample that the barcode/lane matches to. 
Then we can put all the lanes together and know each sample name. The sequences are 
cleaned using the Illumina quality tags.

A) gzipped files (.gz) need to be unzipped. Do this by typing in the command line while in
directory of files that needed to be unzipped:

for f in *.gz
do
gzip -d $f
done

which should be the same as:

for f in *.gz; do gzip -d $f; done

B) All files are separated into number of the wells. These need to be concatenated into
1 single file for R1 (read 1) and R2 (backwards read without barcode, but is paired to R1
by number). To concatenate everything in the right order, make sure that the numbers in
the directory are correctly sorted. Example to type in comment line while in directory of 
files needed to be concatenated:
cat Ragweed6_NoIndex_L004_R1_*.fastq > Ragweed6_R1.fastq

C) process_radtags.sh (http://creskolab.uoregon.edu/stacks/comp/process_radtags.php)
Two ways to do this (I used option2). Either run every single file (001,002….) through 
process.radtags, or concatenate all the files as shown in B). 
Option1:

```
******script process.radtags.plate$.sh********
#/bin/bash

# note by simon. if you have trouble decoding a bash script, you can run it in debug mode
# by boing /bin/bash -x
. /etc/profile
module load stacks

#In this example, we are dealing with paired-end reads. Normally Stacks should be able to
#deal with this using the -P flag, but I can’t get it working. Instead, we define the 
#pairs with $1, which requires an argument input in the command line, e.g. 001, for well 
#1. This mean that the script has to be run for every single well, and it will create a 
#subfolder for every single well (see below)

pair_1=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_NoIndex_L004_R1_$1.fastq.gz
pair_2=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_NoIndex_L004_R2_$1.fastq.gz

#process_radtags is only able to process barcodes of the same length. We have 6, 7, and 8 bp length barcodes in the same lane and need to feed this into the process using:
barcode_file=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_6.txt #path to barcode file
barcode_file2=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_7.txt
barcode_file3=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_8.txt
out_dir=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6_output/$1 #path to output

#make output directory, this is a subfolder for every single well
mkdir $out_dir

#run process_radtags, in this example we have paired-end, Illumina HiSeq data and use the 
#-P flag to show this to the program. Also, data are gzipped, and the -i flag is specified:
#-i gzfastq

process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file -e pstI -c -q -r -i gzfastq
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file2 -e pstI -c -q -r -i gzfastq
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file3 -e pstI -c -q -r -i gzfastq
*****end script*****
``

Option 2:

```
******script process.radtags.plate*all.sh********
#/bin/bash
#see option 1 for better descriptions
. /etc/profile
module load stacks
pair_1=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_R1.fastq
pair_2=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6/Project_Kristin_23_04_14/Sample_Ragweed6/Ragweed6_R2.fastq
barcode_file=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_6.txt #path to barcode file
barcode_file2=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_7.txt
barcode_file3=/nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes97to192_8.txt
out_dir=/nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/plate6_output/all #path to output
mkdir $out_dir
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file -e pstI -c -q -r
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file2 -e pstI -c -q -r
process_radtags -1 $pair_1 -2 $pair_2 -o $out_dir -b $barcode_file3 -e pstI -c -q -r
*****end script*****

``

to run the script: type sh process.radtags.plate[whichever you wanna run][all].sh 
[well# whichever you wanna run, in case of option1, e.g. 004]

These process_radtag output files are in ~/ragweed/GBS/raw_common/output_all/plate*_pr 
(pr for process radtags

D) Use rename.sh to rename the barcodes. Input the file that contains both the barcode 
and the sample name and this shell script will change the names.

******script rename.sh********
#!bin/bash

#awk -F specifies the field separator, in this case .
barcode=$(echo $1 | awk -F"." '{print $2}')
sample=$(echo $1 | awk -F"." '{print $1} ')

echo $barcode
echo $sample

mv sample_$barcode.1.fq $sample.R1.fq
mv sample_$barcode.2.fq $sample.R2.fq
********end script***********

Then, type in command line (in directory where files have to be renamed):

for f in `cat /nfs/home/hpcsci/lotteanv/ragweed/barcodes/barcodes_plate1.txt`; 
do sh ~/scripts/rename.sh $f; done

NB: this is to rename file names of plate1

Good practise to copy files into another folder so original process.radtags output files 
are saved. Renamed files (including .rem.fq files) are in 
~/ragweed/GBS/raw_common/output_all/plate*_rn

E) Next, we move all fq files into one folder. Create a new folder and move the 
renamed files into it. E.g. cp -r *R1.fq ../allind

F) We need to remove reads that contain the adapter sequence or part of the adapter sequence.
Stacks is able to do this, but does a shit job. Kay made a perl script which is better, 
adapter_removal.pl

**********start script**************
#!bin/perl
use warnings;
use strict;
#get file
my $fq_1=$ARGV[0];

open FQ1, $fq_1;

open FQ1_clean, ">$fq_1\_clean";
my $i=0;
my $header1 = ();
my $read1 = ();
my $qual1 = ();
my $badreadno1 =0;
while (<FQ1>) {
    my $line1 = $_;
    chomp $line1;
        if($i==0){
                if($line1=~m/^\@/){
                        $header1=$line1;
                        $i=1;
                }elsif($line1=~m/^\+$/){
                        $i=2;
                }
        }elsif($i==1){
                $read1=$line1;
                $i=0;
        }elsif($i==2){
                $qual1=$line1;
                $i=0;
                #remove forward adapter contamination whole or parts
                if($read1=~m/CTGCAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG/ || $read1=~m/CTGCAAGATCGGAAGAGCG/ || $read1=~m/CGGTTCAGCAGGAATGCCGAG/ || $read1=~m/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT/ || $read1=~m/AGATCGGAAGAGCGTCGTGT/ || $read1=~m/AGGGAAAGAGTGT/){
                        $badreadno1 = $badreadno1 +1;
                }else{
                        print FQ1_clean "$header1\n$read1\n+\n$qual1\n";
                }
        }

}
close FQ1_clean;
close FQ1;

print "The number of reads with adapter contamination is $badreadno1\n";

********end script***********
                                                                                
In allind type in command line:

for i in *.fq; do perl ~/scripts/adapter_removal.pl $i ; done

When working, this will give output lines on the screen:

>The number of reads with adapter contamination is 14329
>The number of reads with adapter contamination is 10637
>The number of reads with adapter contamination is 163986
>The number of reads with adapter contamination is 129353
>The number of reads with adapter contamination is 35568
>The number of reads with adapter contamination is 28743

#This will additionally change the *.fq to *.fq_clean, so you know that the file has been 
#looked at. Check and compare file size to make sure files were reduced in size.

#Copy the cleaned, renamed files into another new folder (clean):
cp -R *fq_clean ../clean
#The output files have the ending '.fq_clean' Other programs do not like this. Change the 
#endings by using the command:

rename .fq_clean .fq *fq_clean

#E) Trimming reads.  NOTE: This step has to be skipped when using STACKS, as it is not
#able to handle reads of different lengths. bwa is though.
#Sickle uses sliding windows along with quality and length thresholds 
#to determine when quality is sufficiently low to trim both ends of reads.
#the following option in Sickle takes 2 input files, produces 2 output files and a
#orphan file:
#sickle pe -f input_file1.fastq -r input_file2.fastq -t sanger \
#-o trimmed_output_file1.fastq -p trimmed_output_file2.fastq \
#-s trimmed_singles_file.fastq
#pe stands for paired-end reads, -t is quality encoding (phred+33 is sanger)

for i in *R1.fq; do sickle pe -t sanger -f $i -r ${i/R1.fq/}\R2.fq \
-o ../trimmer/${i/R1.fq/}\R1_tr.fq -p ../trimmer/${i/R1.fq/}\R2_tr.fq \
-s ../trimmer/${i/R1.fq/}\S_tr.fq -q 20 -l 50 \
>> ../trimmer/${i/.R1.fq/}\_tr.log 2>&1 $i; done


#F) It is really important that all the reads are of good quality because STACKS aligns 
#them without caring about the quality score. This step might be skipped when trimmer is
#already used, but I prefer to use both, as the trimmer only trims off from the ends
# until the quality is good. But maybe the average read qual is still crap
#load fastx
. /etc/profile
module load fastax
#run quality filter on .fq.clean and creates output files _q10.fq and .fq.log (log file 
#for filter). Or _tfq20.fq for both trimmed and filtered reads (in dir trimfilt)

for i in *.fq; do  fastq_quality_filter -q 10 -p 90 -Q 33 -v -i $i \
-o ../filter/q10/${i/.fq/}\_q10.fq >> $i_q10.log 2>&1 $i; done

for i in *.fq; do  fastq_quality_filter -q 20 -p 90 -Q 33 -v -i $i \
-o ../filter/q20/${i/.fq/}\_q20.fq >> $i_q20.log 2>&1 $i; done

#to see the amount of reads that are filtered out, I can used the produced log files
#first make a list with the names
ls -laht | awk '{print $NF}' | grep ‘_q10.fq$’ | cut -d. -f1  | sort | uniq > list_sample
awk '{print $0".R1"}' list_sample > list_sampleR1R2
awk '{print $0".R2"}' list_sample >> list_sampleR1R2

#then make a new text file containing this information
echo -e "Sample\nInput\nOutput\nDiscarded\n" > qualfiterlog.txt; 
while read i;  
do echo "$i" >> qualfiterlog.txt; 
cat ${i}.fq.clean.log | awk '$1 ~ /^Input:/ { print $2}' >> qualfiterlog.txt;  
cat ${i}.fq.clean.log | awk '$1 ~ /^Output:/ { print $2}' >> qualfiterlog.txt; 
cat ${i}.fq.clean.log | awk '$1 ~ /^discarded/ { print $2}' >> qualfiterlog.txt; 
done < list_sampleR1R2
#remove whitespaces and empty lines respectively:
sed -i -e 's/^[ \t]*//' -e 's/[ \t]*$//' qualfiterlog.txt
sed -i '/^$/d' qualfiterlog.txt
#put all this in 4 columns
awk 'ORS=NR%4?" ":"\n"' qualfiterlog.txt >> qualfiterlog2.txt

#copy to local computer (not ssh into mcc)

scp -pr lotteanv@msgln4.its.monash.edu:~/ragweed/GBS/raw_common/output_all/clean/qualfiterlog2.txt Documents

#G) Re-pair reads. During filtering some reads have been filtered out, which means that 
#gaps/differences exist between R1 and R2. This can be fixed by re-pairing the R1 and R2
#files. NOTE however that below script does not do anything with "orphans", non-paired reads

Copy all _f.fq files to pairing:
cp -R *_f.fq ../pairing

#saved in ~scripts/fix_fqpair_all.sh
#let's go to your directory where the files live
cd ~/ragweed/GBS/raw_common/output_all/pairing

#throw all basenames of files to be processed into list. here we look at all files ending 
#with .fq - change your grep if you want to be more specific. make sure you dont capture 
#out output files into list_sample. also this assumes the base file name is before the 
#first dot.
ls -laht | awk '{print $NF}' | grep '_q10.fq$' | cut -d. -f1 | sort | uniq > list_sample

#do shit to the list. 
while read i
do
perl ~/scripts/FQ_pair_no.pl ${i}.R1_q10.fq ${i}.R2_q10.fq >> $i_q10.log 2>&1
done < list_sample

mv paired* ../paired

---
Task 2: Align paired files to crappy reference (WGS) using BWA (NB: I chopped everything up in
one-liners but these can be put together in a script)
. /etc/profile
module load bwa

A) Rename sequence identifier of fq files. 

During process.radtags, stacks renamed R1 and R2 sequence identifiers to _1 and _2. When running bwa mem it expects same file names or an error will appear, e.g. [mem_sam_pe] paired reads have different names: 
"5_1202_2403_86748_1", "5_1202_2403_86748_2" 
```sed 's/_1$//g' filenameR1 > newfilenameR1```
```sed 's/_2$//g' filenameR2 > newfilenameR2```

```for i in *.R1_f.fq; do basename=${i/.R1_f.fq}; sed 's/_1$//g' $i > ${basename}.R1.fq; done```
```for i in *.R2_f.fq; do basename=${i/.R2_f.fq}; sed 's/_2$//g' $i > ${basename}.R2.fq; done```

OR (as I have encountered problems with above commandline:

```
ls -laht | awk '{print $NF}' | grep 'R1’ | sort | uniq > list_allR1
ls -laht | awk '{print $NF}' | grep 'R2' | sort | uniq > list_allR2
while read i; do sed  ’s/_1$//g' $i > ${i}_s; done < list_allR1
while read i; do sed  's/_2$//g' $i > ${i}_s; done < list_allR2
```

Copy *fq_s to align folder

rename .fq_s .fq *fq_s

B) (optional) In case of fragmented reference, it is better to make a pseudo-scaffold. This is 
# a concatenated file with all the contigs togethers, separated with 30 A's. 
perl ~/scripts/scaffolds.pl soaprunk61.contig
#output files are soaprunk61.contig.pseudo  & scaffold_order.soaprunk61.contig
#REMEMBER to translate the SNP table back later to the original genome locations
# with scaffold2contig.v2.pl

C) GATK is not able to handle N's, count and replace these with A's
Count number of N’s in the file as GATK can’t handle N’s
#note: whenever counting lines make sure to account for header lines. 
#Use grep -v '>'  first bascially because you don’t want to count characters in the 
#header and you should never ever assume the header won't contain some chars. 
#grep -v is "inverse grep": that is "look for lines without this"
#But without accounting for this:

###Option1 (—> gives string of A’s and counts the characters)
tr -cd N < soaprunk61.contig.pseudo | wc -c 
 
#OR 
###Option2(—> slower option, gives char in lines and then counts the number of lines. 
# wc -c won’t work here as there will be twice as many characters, A and newline)
fgrep -o N soaprunk61.contig.pseudo | wc -l 

# In case of N's replace with A's with sed, not done here as wasn't necessary

# can also be used to see how many A's were added in the scaffolds.pl script
# by counting the difference between soaprunk61.contig and soaprunk61.contig.pseudo
echo '712102436 - 421159706' | bc

D) Index .contig (or in case of B, contig.pseudo) file for samtools, bwa and gatk

#Check maximum line lengths in file with:
awk '{ if (x < length()) x = length() } END { print x }' <referencegenome>
#Check for anything shorter that this maximum length (here 70)
awk 'length($0) < 70' <referencegenome>
#In case not all lines are same length:

#delete white lines
sed '/^$/d' <filename> > <newfilename>

# -a specifies the indexing algorithm bwa uses. There is another option (IS), for smaller
# genomes (<2GB), check bwa man page
bwa index -a bwtsw <referencegenome>

samtools faidx <referencegenome>

#NB: Note by Kay (SNP calling blog) --> if you want to use Picard downstream use the -M option in bwa. This is NOT 
#implemented in below script! From bwa man: Mark shorter split hits as secondary 


E) Run bwa (in align folder) 
# this includes conversion from .sam to .bam and sorting
#The listreadgroups list used later is without the "paired_" prefix, so this has to be 
#removed from the files. Also, this list has 1.* instead of pl1-* etc
for f in paired_pl*; do mv "$f" "${f#paired_pl}"; done
#rename pl1- 1- pl1-*
#etc

#Make a list with the just the basenames (cut off last 6 characters)
ls -laht | awk '{print $NF}' | grep '.fq$' | sed 's/\(.*\)....../\1/' | sort | uniq > list_sample

#while read i; do sh ~/scripts/SNP_calling/bwa_cra_lotte.sh $i; done < list_sample
#The following is implemented in bwapseudo.job in align directory:
#while read i; do sh ~/scripts/SNP_calling/bwa_cra_pseudo.sh $i; done < list_sample

Example output:
[M::main_mem] read 206186 sequences (20000042 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 30466, 47, 0)
[M::mem_pestat] skip orientation FF as there are not enough pairs
[M::mem_pestat] analyzing insert size distribution for orientation FR...
[M::mem_pestat] (25, 50, 75) percentile: (114, 166, 211)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 405)
[M::mem_pestat] mean and std.dev: (169.49, 71.38)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 502)
[M::mem_pestat] analyzing insert size distribution for orientation RF...
[M::mem_pestat] (25, 50, 75) percentile: (39, 39, 39)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (39, 39)
[M::mem_pestat] mean and std.dev: (39.00, 0.00)
[M::mem_pestat] low and high boundaries for proper pairs: (39, 39)
[M::mem_pestat] skip orientation RR as there are not enough pairs
[M::mem_pestat] skip orientation RF

#this should be run in parallel and an example job is found in bwajob.job (align dir)

# to go through conversion of sam to bam as one-liners (this is however implemented in 
# bwa_cra_pseudo.sh script, not in bwa_cra_lotte.sh script):
#for i in *.sam; do basename=${i/.sam};  samtools view -u -t  $SAM_INDEX  -S ${basename}.sam -o ${basename}.bam; done
#for i in *.bam; do basename=${i/.bam}; samtools sort ${basename}.bam ${basename}.sort; done
#for i in *sort.bam; do basename=${i/.sort.bam}; samtools index ${basename}.sort.bam; done

#check number of mapped and unmapped reads
#the file produced will show the number of mapped reads, unmapped reads reads where
#both read pairs are mapped and the bam file size (useful info on 
#http://left.subtree.org/2012/04/13/counting-the-number-of-reads-in-a-bam-file/)
#the du command gives size of file.
echo -e 'sample\nmapped_reads\nunmapped_reads\nboth_reads_mapped\nbam_size' > mappedreads.txt; 
while read i;  do echo "$i" >> mappedreads.txt; samtools view -c -F4 ${i}.sort.bam >> mappedreads.txt;  
samtools view -c -f 4 ${i}.sort.bam >>mappedreads.txt; samtools view -c -f1 -F12 ${i}.sort.bam >> mappedreads.txt; 
du ${i}.sort.bam | awk '{print  $1}' >> mappedreads.txt; done < list_sample

#and transform the rows to columns:
awk 'ORS=NR%5?" ":"\n"' mappedreads.txt > mappedreads.txt

--------------------------------------------------------------
Task 3: Use GATK to call variants
#This task is divided into Task 3.1 using UnifiedGenotyper; and 3.2 using Haplotypecaller

Make list to add Read Groups (which will run in  gatk_cra_pseudo.sh, in step 3.1A) 
#gatk is delicate and can't work without specified read groups. I made a list with
#what I think are appropriate read-groups, with library being the illumina plate, platform
#illumina, platform unit the barcodes and sample name the name of each individual.

****************************************
Task 3.1 Calling SNPs and genotypes with UnifiedGenotyper
A) Re-align around the indels. 
#Indels can cause alignment mismatches, especially near the end of an alignment and for 
#larger indels. This happens because in the alignment algorithm introducing an indel is 
#more 'expensive' than introducing a SNP. Therefore, the global alignment of all samples
#needs to be considered. Reasons to re-align: a) known sites of indels; b) indels in the 
#original alignment (as can be seen in the CIGAR string); c) sites showing evidence of 
#indel (clustering of SNPs)
#It is good to include sites of known indels & sites of known variation in the realignment. 
#RealignTargetCreator creates targets of interest for the realignment and outputs an
#.intervals file. This tool can use a vcf file of known indels. Not necessary if not
#available, but will speed up process.
#The IndelRealigner is a heavier computational process which is doing the actual realignment
#It changes the CIGAR string of re-aligned reads and maintains the original CIGAR string
#with an OC tag. It is thus easy to assess afterwards how many reads have been re-aligned
#by using grep etc.
#A full Smith-Waterman realignment is recommended when: a) older data; b) shorter reads
#(~36bp); c) no known indels; d) if you want to reduce false-positives
#Indel realignment is very important before proceeding with the Base Quality Score
#Recalibration

#########SCRIPT SNP_calling/gatk_cra_pseudo.sh################
#/bin/bash

REF=/nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa

RUN=$(echo $1 | awk -F"." '{print $1}')
LIB=$(echo $1 | awk -F"." '{print $3}')
#SAMPLE=$(echo $1 | awk -F"." '{print $2}')
ID=$(echo $1 | awk -F"." '{print $1"-"$2}')

echo $RUN $LIB $SAMPLE $ID

java -jar $PICARD/AddOrReplaceReadGroups.jar I=$ID.sort.bam O=$ID.rg.bam LB=$LIB PL=illumina PU=$RUN SM=$ID VALIDATION_STRINGENCY=SILENT  >> $ID.log 2>&1

#index
samtools sort $ID.rg.bam $ID.rg.sort
samtools index $ID.rg.sort.bam

#ID indels
java -Xmx15g -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar -T RealignerTargetCreator -S LENIENT -R $REF -nt 10 -o $ID.bam.list -I $ID.rg.sort.bam > $ID.realign1.log

# realign around indels
java -Xmx10g -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar -T IndelRealigner -R  $REF -targetIntervals $ID.bam.list -I $ID.rg.sort.bam -o $ID.realign.bam > $ID.realign2.log
#############END SCRIPT########

#implemented in gatkpseudo.job in bwa_genome_pseudo dir
#while read f; do bash ~/scripts/SNP_calling/gatk_cra_pseudo.sh $f; done < listreadgroups.txt

#Check if all files are present. If not (which happened to me), they can be checked by 
#using the following command (where list_sample2 contains ALL sample names:
ls -laht | awk '{print $9}' | grep 'rg.bam' | cut -f1 -d. > list_got
#this will print list_got with all *rg.bam files in the working directory
diff <(sort list_got) <(sort list_sample2)

#ND --> USE GREP TO CHECK HOW MANY READS HAVE BEEN RE-ALIGNED

C) Base recallibration (look at bwa_red_gatk2.sh, started writing in gatk_cra2.sh)
# this step is only possible after a preliminary run that identifies SNPs and creates a vcf file

D) SNP and genotype calling with UnifiedGenotyper
#Make a list of all the realigned bam files to be the input file for the unified genotyper
ls -laht | awk '{print $9}' | grep 'realign.bam' > bam.list

############START SCRIPT################
#!/bin/bash

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-R /nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa \
-T UnifiedGenotyper \
-I bam.list \
-o snps.raw.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 20.0 \
############END SCRIPT###############

sh gatk_cra3.sh
#when logging a job for 288 samples, this took 2.5 days

E) filter vcf file
#this script will locate called SNPs and genotypes and filters on set quality and depth
#it will produce an 'N' for every position (SNP/individual) for which no genotype was 
#called by UG

#call vcf file and get all the SNPs (grep -v will get everything BUT the command)
cat snps.raw.vcf | grep -v 'INDEL' | /nfs/home/hpcsci/lotteanv/scripts/vcf2vertical_dep_GATK-UG.pl > snptableUG.tab

#get the number of filtered SNPs:
cat snptableUG.tab | grep Scaf* | wc -l

#vcfdepth_lotte.pl is the exact same as vcf2vertical_dep_GATK-UG.pl but will produce
#'depth.txt', which gives the depth of the # of reads for each called genotype

#Explore the depth distribution with R
R
>d<-read.table('depth.txt')
>head(d) #R called my column V1
>hist(d$V1)

F) get summary statistics and filter SNP table based on minor allele frequency,
heterozygosity and missing data using snp_coverage.pl
perl /nfs/home/hpcsci/lotteanv/scripts/snp_coverage.pl snptableUG.tab 
#Minor allele frequency: Low frequency SNPs could be due to errors and are not useful for
#outlier tests and several other tests of selection (although they are for site frequency
#spectrum tests)
#Heterozygosity: High or fixed heterozygosity could indicate parology.
#Missing data: SNPs that have low coverage are removed
#explore data in R using the summary output:
d<-read.table('snptableUG.tab.summary',col.name=c("scaf_no","position","allele_tot","genotype_tot","top_o","sorted","hash_0","mj","sort_1","hash_1","mn","het"))
hist(d$genotype_tot,xlim=c(0,50),breaks=1000)
hist(d$het)
hist(d$mn) 

#once decided on cut-off values (I used mn_cut = 0.05, het < 0.7 and geno_cut = 
#<half of total number of individuals>)
#the snptableUG.tab.table is the new filtered output file

#check quality per individual
#add header names to snptableUG.r.tab.table by deleting # from header in snptableUG.tab to
#new file	
R
df<-read.table("snptable.tab.r192.table", header=T)
head(df)
lst<-colnames(df[,-(1:2)]) # discards the first 2 columns in the lst headers
countN<-sapply(lst,FUN=function(x,df){sum(df[,x]=="N",na.rm=T)},df) 
#counts the total number of N’s for each individual
percN<-sapply(lst,FUN=function(x,df){sum(df[,x]=="N",na.rm=T)}/nrow(df),df) 
#calculates the percentage of N’s per ind divided by total number of SNPs
hist(percN)
write.table(percN,"percN_filt192i.txt",quote=F,sep="\t")
write.table(countN,"countN_filt192i.txt",quote=F,sep="\t")

#copy these files to local computer --> open new terminal but DO NOT ssh into MCC


G) Put scaffolds back into original contig positions
perl /nfs/home/hpcsci/lotteanv/scripts/scaffold2contig.v3.pl snptableUG.tab.table \
	~/ragweed/WGS/soap_assembly/scaffold_order.soaprunk61.contig

#check how this all went
#this will get every SNP for each unique Scaffold
-cat snptableUG.tab.table | grep Scaffold* | cut -c 1-13 | sort | uniq -c | wc -l
353 # this means that 35% of the arbitrary scaffolds made contain SNPs
#After putting back to original contigs
awk 'NR!=1{print $1}' snptableUG.tab.table.contig | sort | uniq -c | wc -l
5244
#that is possible as the Scaffolds were arbitrary, so it makes sense that the SNPs are on 
#more contigs than scaffolds.


********************************************
Task 3.2. Calling haplotypes with HaplotypeCaller
# This is a fairly new tool and works at the moment only with diploid organisms
# Approximately 5 times slower than UnifiedGenotyper
A) Haplotype calling
HaplotypeCaller is approximately 5 times slower than UnifiedGenotyper, depending on
parameters set by user. 
This tool is able to detect large variance and looks at regions of variation instead of 
independent loci. It does this through local denovo assembly. It uses all the reads
localised to a region, so it does not assume the exact alignment made by the alignment
software. For this reason, indel re-alignment is not necessary (but it should not matter
if this step is preformed before using this tool). HC does not use prior known data,
like known indels etc.
This tool works better than UG because it takes multiple loci into account, e.g UG
calls several loci where it's only 1 deletion (same reason for ID realignment). HC also
calls variants on regions that are too big for ID realigner.
HC brings in ALL mate pairs, including un-mapped mates.

#make list of input bam files (I WANT TO CHECK IF RESULTS ARE DIFFERENT WITH AND WITHOUT 
#LOCAL INDEL REALIGNMENT

#with local realignment (implemented in hapl.job):

#############START SCRIPT hapl_cra.sh############
#!/bin/bash

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-R /nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa \
-T HaplotypeCaller \
-I bam.list \
-o hapl.raw.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 20.0 \
-minPruning 3
#############END SCRIPT#################

sed -e 's/$/.rg.sort.bam/' list_sample2 > rgsortbam.list

#############START SCRIPT############
#!/bin/bash

java -jar /opt/sw/gatk-3.1.1/GenomeAnalysisTK.jar \
-R /nfs/home/hpcsci/lotteanv/ragweed/WGS/soap_assembly/soaprunk61.contig.pseudo.fa \
-T HaplotypeCaller \
-I rgsortbam.list \
-o hapl.raw.vcf \
-stand_call_conf 50.0 \
-stand_emit_conf 20.0 \
-minPruning 3
#############END SCRIPT#################
Pruning: the amount of minimum reads that have to show the variation to be included.

**************************************
Task 3.3. Calling haplotypes with Beagle
#needs a vcf as input

————————————————————————————————————————
Task 4: Prepare the SNP-tables for the different programs

STRUCTURE
#Can't deal with linked loci, so need to select at random one SNP per contig
perl /nfs/home/hpcsci/lotteanv/scripts/random_one_per_locus_combined.pl snptableUG.tab.table.contig

#Check if number of selected contigs is the same as number of unique contigs
awk 'NR!=1{print $1}' snptableUG.tab.table.contig | sort | uniq -c | wc -l
5244
awk 'NR!=1{print $1}' snptableUG.tab.table.contig.random | wc -l
10

BAYENV

LFMM



--------------------------------------------

Task X: Run "denovo.pl" This program creates stacks within each individual, then puts these stacks together to form catalog stacks (loci).

A) Create an environment to work in
in output_all
$ mkdir stacks
$ cd stacks
$ mkdir raw samples stacks paired assembled
$ ls
assembled  paired  raw	samples  stacks
Copy all the cleaned files into samples

B) Create a mysql database to import the results into. On the MCC, this requires using credentials to get access to mysql:
$ mysql --user=stacksuser --password=st4cks852 [command]

To create the database:

$ mysql --user=stacksuser --password=st4cks852 -e "CREATE DATABASE 125cra_radtags"

The database table definitions need to be send to the server to create all necessary components of the database, telling the database how to link tables.

$ mysql --user=stacksuser --password=st4cks852 temp_radtags < /opt/sw/stacks/1.13/share/stacks/sql/stacks.sql

To check that this worked: 
$ mysql --user=stacksuser --password=st4cks852 
$ use temp_radtags;
$ show tables;

this should show a list of tables:
+------------------------+
| alleles                | 
| batches                | 
| catalog_alleles        | 
| catalog_annotations    | 
| catalog_genotypes      | 
| catalog_snps           | 
| catalog_tags           | 
| chr_index              | 
| fst                    | 
| genotype_corrections   | 
| markers                | 
| matches                | 
| pileup                 | 
| populations            | 
| ref_radome             | 
| samples                | 
| sequence               | 
| sequence_blast         | 
| snps                   | 
| sumstats               | 
| unique_tags            | 
+------------------------+
21 rows in set (0.00 sec)

To exit mysql:
$exit

C) Next, we can run denovo_map.pl

To print all the files in the required format for the denovo_map.pl type:
$ ls -l *.fq | awk '{print "-s", $9, "\\"}'
copy paste into the command and delete the last \

Remember that if you somehow cause this script to mess up, you will have to go back and re-create your mysql database and re-import the definitions again. This is because the database won't write over itself. Specify an output folder and create it before attempting to run denovo_map.pl

Parameters used:
-m:	Minimum number of raw reads required to make a stack, so the minimum depth
-M:	The maximum number of nucleotide mismatches expected between haplotypes at a locus within an individual. With high expected overall heterozygozity, this value can be set higher than default 3
-n:	The maximum number of nucleotide mismatches expected between any two haplotypes in the population. Follows same reasoning as -M
-T: 	The number of thread cores to run Stacks on. A higher number of cores used, means a lower time to complete analysis. Find out how many cores can be used

Multiple batches can be run in the one created database (as mentioned above). This can be useful when running the assembly with different parameters. Design the batches (make a list of all the different parameters to run) and for every batch only change the output file -o batch name -D and batch number -b.

denovo_map.pl -m 3 -M 3 -n 2 -T 2 -o /nfs/home/hpcsci/lotteanv/ragweed/GBS/raw_common/output_all/stacks/stacks/batch1 -b 1 -B 125cra_radtags -D “125cra_batch1” \
-s AA-1-10.R1_f.fq \
-s AA-1-10.R2_f.fq \
-s AA-1-19.R1_f.fq \
-s AA-1-19.R2_f.fq \
-s AA-1-20.R1_f.fq \
-s AA-1-20.R2_f.fq \
-s AA-1-21.R1_f.fq \
-s AA-1-21.R2_f.fq \
-s AA-1-24.R1_f.fq \
-s AA-1-24.R2_f.fq \
-s AA-10-25.R1_f.fq \
-s AA-10-25.R2_f.fq \
-s AA-10-26.R1_f.fq \
-s AA-10-26.R2_f.fq \
-s AA-10-30.R1_f.fq \
-s AA-10-30.R2_f.fq \
-s AA-11-1.R1_f.fq \
-s AA-11-1.R2_f.fq \
-s AA-11-13.R1_f.fq \
-s AA-11-13.R2_f.fq \
-s AA-11-22.R1_f.fq \
-s AA-11-22.R2_f.fq \
-s AA-11-25.R1_f.fq \
-s AA-11-25.R2_f.fq \
-s AA-11-4.R1_f.fq \
-s AA-11-4.R2_f.fq \
-s AA-12-10.R1_f.fq \
-s AA-12-10.R2_f.fq \
-s AA-12-11.R1_f.fq \
-s AA-12-11.R2_f.fq \
-s AA-12-13.R1_f.fq \
-s AA-12-13.R2_f.fq \
-s AA-12-8.R1_f.fq \
-s AA-12-8.R2_f.fq \
-s AA-13-12.R1_f.fq \
-s AA-13-12.R2_f.fq \
-s AA-13-13.R1_f.fq \
-s AA-13-13.R2_f.fq \
-s AA-13-16.R1_f.fq \
-s AA-13-16.R2_f.fq \
-s AA-13-17.R1_f.fq \
-s AA-13-17.R2_f.fq \
-s AA-14-10.R1_f.fq \
-s AA-14-10.R2_f.fq \
-s AA-14-11.R1_f.fq \
-s AA-14-11.R2_f.fq \
etc

To browse the database in a firefox window, make sure you are logged in through the username@msgln6.its.monash.edu.au -Y and vlm001 -Y (-Y will give access to windows)
To open firefox:
$ firefox &
In the firefox window, go to: https://localhost/stacks/catalog.php?db=temp_radtags&id=1

--------------------------------
Task 3: Match the output from the denovo output to the Loblolly genome using the program BWA.

A) Convert stacks output to fasta format. First line starts with a ">" and has information about the run. Then there is a hard return to make a second line. The second line has the sequence. There are no quality scores associated with FASTA.

batch1: perl /data/GBS/catalogue2consensusfasta.pl clean_pineGBS/pineGBSstacks/batch_1.catalog.tags.tsv > batch1_pineGBSfull_consensus.fasta

output: batch1_pineGBSfull_consensus.fasta

B) Input the consensus file that you created into BWA to align it to the loblolly genome.

# open byobu and change the working directory to clean_pineGBS/pineGBSstacks/

batch1: bwa mem -t 2 /data/Pine_genome/ptaeda.v1.01.fa.masked_cutN.noblank batch1_pineGBSfull_consensus.fasta > batch1_pineGBSfull_consensus.sam
