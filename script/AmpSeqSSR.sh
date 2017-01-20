# AmpSeq-SSR 
# email: weixiong.zhang@wustl.edu


# Program requirement
# 1. reference sequences of each amplicon
# 2. primer sequences of amplicons
# 3. sequenceing data
# 4. bowtie2
# 5. perl 
# 6. blast tool suite



if [ $# != 5 ] ; then 
echo "USAGE: sh AmpSeqSSR.sh <sequencing reads fasta> <reference fasta> <primer fasta> <directory of scripts> <sample name>\n" 
exit 1; 
fi 




#!/usr/bin/bash

DAT=$1
REF=$2
PRIMER=$3
SCRIPT_DIR=$4
sam=$5




# Step 1: processing amplicon reference
echo "Step 1: processing amplicon reference\n"
# 1.1 extract SSRs within reference
perl $SCRIPT_DIR/ssr_f.pl $REF ampliconRefenece 5 2:4,3:4,4:4,5:4

# 1.2 mask repeat sequences with N in reference.
perl $SCRIPT_DIR/maskReference.pl ampliconRefenece.ssr $REF > ampliconRefenece_ssr2N.fa

# 1.3 format reference
formatdb -i ampliconRefenece_ssr2N.fa -p F 
bowtie2-build $REF ampliconRefenece


# Step 2: assign the reads to amplicons
echo "Step 2: Assigning reads to amplicons\n"
# 2.1 map reads to reference by bowtie2
bowtie2 -f -x ampliconRefenece -U $DAT -S $sam.sam --un $sam.unmapped.fa 1>$sam.bwt.log 2>$sam.bwt.err

# 2.2 map unmapped reads to masked reference
perl  $SCRIPT_DIR/bath.pl $sam.unmapped.fa ampliconRefenece_ssr2N.fa $sam 

# 2.3 combined with mapped reads
ls $sam.sam|perl -e '$line=<>;$sam = (split /\./,$line)[0];open IN, "$sam.sam" or die "$sam.sam not exists!\n"; open OT, ">>  $sam.Run100.combin.fa"; while(<IN>){chomp; next if /^\@/;($id,$flag,$aid,$seq)=(split /\s+/,$_)[0,1,2,9]; print OT ">$id\_$aid\n$seq\n" if $flag !=4;} close IN; close OT;'

# Step 3: call genotype of each read
echo "Step 3: gentoype calling\n"
perl $SCRIPT_DIR/getSubseqBetweenSSR.pl ampliconRefenece.ssr $REF $sam.vs.singleFlank.bst $sam.Run100.combin.fa 0 $PRIMER 10

perl $SCRIPT_DIR/extract4SSR.pl ampliconRefenece.ssr $sam.vs.singleFlank.bst $sam.Run100.combin.fa $sam $REF $PRIMER 10
#ATTENTION!!!

# Step 4: tally reads for each alleles
echo "reads tallying\n"
perl $SCRIPT_DIR/process.pl $sam.stat $sam.stat


# Step 5: summarize results
echo "Step 5: results summarizing\n"
perl $SCRIPT_DIR/alleleProfile.pl $sam.stat > $sam.allele_profile

cat $sam.stat | perl -e 'print "ampliconID\t#totalreads\t#naReads\tmajorSSR\tmajorAbd\tsecondSSR\tsecondAbd\tstutter\t#readsOfOtherSSR\n"; while(<>){chomp;($aid,$gt,$pos,$sid,$gtr)=(split /\s+/,$_); $sid=~/\_x(\d+)/;$h{$aid}{$gtr}+=$1; }  foreach my$k(keys %h){ @tmp=(sort {$h{$k}{$b}<=>$h{$k}{$a}} keys %{$h{$k}}); die ("wrong for $k\n")  if @tmp ==0; my$na=0;   my($most,$mostabd)=($tmp[0],$h{$k}{$tmp[0]}); $na =$mostabd if $most=~/NA/;  my($sed,$sedabd)=("-",0); ($sed,$sedabd)=($tmp[1],$h{$k}{$tmp[1]}) if @tmp>=2; $na +=$sedabd if $sed =~/NA/; my$other=0; for(my$i=2;$i<=$#tmp;$i++){$other +=$h{$k}{$tmp[$i]}; $na +=$h{$k}{$tmp[$i]} if $tmp[$i]=~/NA/; $ot .="$tmp[$i],$h{$k}{$tmp[$i]};"; } $total =$mostabd+$sedabd+$other;  my$ratio=($sedabd/$mostabd); print "$k\t$total\t$na\t$most\t$mostabd\t$sed\t$sedabd\t$ratio\t$other\t$ot\n";$ot="";} ' >$sam.stat.xls
perl $SCRIPT_DIR/addAmpliconLen.pl -f $sam.stat -c 1 -m 1 -r $REF -p $PRIMER -o $sam -s $sam.Run100.combin.fa


echo "Finished!!!"

