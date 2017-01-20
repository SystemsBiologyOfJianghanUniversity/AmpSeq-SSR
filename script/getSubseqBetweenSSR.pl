#!/usr/bin/perl -w
use strict;

die("Usage: perl $0 <single.ssr> <single.fa> <outF> <reads> <1:kept flank seq> <primer> <flankLen>
       Used to mapped reads on flank seq of each single amplicon region(ssr).
       single.srr: ssr list for each amplicon region, generate by ssr.pl
       single.fa : seq in fasta formate for each amplicon (ssr).
       reads     : query seq
       kept flank: if kept the temp file for flank seq of each amplicon.
       primer    : primer file used for amplify the target region;
       flankLen  : the minium length of flank sequence for each SSR; this should be KEPT in line with extract4SSR.pl;

       Revised:
            we ignored all SSR located on the primer regions by add line 44, considering that such SSR will be digested in library construction. @ 20160720.
            we still cannot process SSR that locate at both end of reads.

       outF could be furthur processed with extract4SSR.pl.^_ ^\n") unless @ARGV ==7;

my($ssrF,$singleF,$outF,$reads,$kept,$primer,$flankLen)=@ARGV;
my$pre=(split /\./, (split /\//,$reads)[-1])[0];

print `mkdir flank ` if $kept==1;

my%primerLen=();
open IN, " $primer " or die  "cannot open $primer $!";
$/="\>";<IN>;
while(<IN>){
    chomp;
    my($pid,$pseq)=(split /\n+/,$_);
#$primerLen{$pid}=length$pseq;
    $primerLen{$pid}=rindex($pseq,"T");#### all primer will be digested at T base, so the last T will be used as end of flank seq;
    $primerLen{$pid}+=1;
}
close IN;
$/="\n";

my%h=();
my(%pos,%RepeatUnit)=();
#open IN, "single.ssr" or die  "cannot open single.ssr $!";
open IN, " $ssrF " or die  "cannot open $ssrF $!";

while(<IN>){
    chomp;
    my($aid,$repeat,$num,$s,$e,$len)=(split /\s+/,$_)[0,3,4,5,6,7];
#    next if ($s<=$primerLen{"$aid\_F"} or $e>=($len-$primerLen{"$aid\_R"}));
    $RepeatUnit{$aid}{$s}=$repeat;
    $repeat=$repeat x $num; 
    $h{$aid} .= $repeat; 
    $h{$aid} .="_"; 
    $pos{$aid}{$s}=$e;
}

#foreach my$kk(keys %{$pos{"AMPL1119031"}}){print "$kk\t".$pos{"AMPL1119031"}{$kk}."\n";} exit;
#print `rm -rf shann.vs.singleFlank.bst.W9`;
print ` rm -rf $outF `;
#open IN, " single.fa"; 
open IN, " $singleF "  or die  "cannot open $singleF $!";
$/=">"; <IN>;

while(<IN>){
    chomp;
    my($id,$seq)=(split /\n+/,$_); 
    $seq=lc$seq;
    die ("No repeat for $id\n") unless exists $h{$id};
    foreach my$key(keys %{$pos{$id}}){
        my$sub=substr($seq,$key-1,($pos{$id}{$key}-$key+1));
        my$slen=length$sub;
        my$N="N" x $slen;
        my$left=substr($seq,0,($key-1));
        my$right=substr($seq,$pos{$id}{$key});

        my$repeat_tmp=$RepeatUnit{$id}{$key};
        if($left=~/(($repeat_tmp){1,}\w{0,5})$/){
            my$tmpLen=length$1;
            $left=~s/\w{$tmpLen}$//;
            $left.= "N" x $tmpLen;
        }

        if($right=~/^(\w{0,5}($repeat_tmp){1,})/){
            my$tmpLen=length$1;
            $right=~s/^\w{$tmpLen}//;
            my$tmpN ="N" x $tmpLen;
            $right="$tmpN$right";
        }
        $seq="$left$N$right";
#        $seq=~s/$sub/$N/;

    }
    if($seq=~/^N+/ or $seq=~/N+$/){
        print STDERR "No flank seq found for $id \n";
        next;
    }


=cut
    my@t=(split /\_/,$h{$id});

    my$flag=0;
    foreach my$k(@t){
        if(($seq=~/^$k/) or ($seq=~/$k$/)){
            $flag=1;
            last;
        }
        $seq=~s/$k/N/;
    }

    if($flag==1){
        print STDERR "$id need to be check\n";
        next;
    }
print "$seq\n" if $id eq "AMPL1075249";
=cut
#    print ">$id\n$seq\n" if $id eq "AMPL1119031"; 
    my@t=(split /N+/,$seq);
    my$pFlen = $primerLen{"$id\_F"};
    my$pRlen = $primerLen{"$id\_R"};

    my$flankFlag=0;

#    my($tF,$tR,$pFlen,$pRlen)=("","",0,0);

    $flankFlag=1 if ((length$t[0])-$pFlen)<$flankLen;
    $flankFlag=1 if ((length$t[-1])-$pRlen)<$flankLen;
    
    if($flankFlag==1){
        print STDERR "$id need to be check\n";
print "(length$t[0])-$pFlen)<$flankLen\n(length$t[-1])-$pRlen)<$flankLen\n";
    }  
   
    open OT, ">  tmp.$pre.fa " or die  "cannot write tmp.$pre.fa $!";
    my $i=0;

############## if the flank of most left is shorter than 20bp, it will be ignored;
print "$t[0])-$pFlen)>=$flankLen\n" if $id=~/AMPL1118676/;
    if(((length$t[0])-$pFlen)>=$flankLen){
        $i++;
        print OT ">$id\_$i\n$t[0]\n";
    }
    shift @t;

    foreach my$seq(@t){ 
#        $i++;
        next if ((length$seq)<$flankLen);  
        $i++;
        print OT ">$id\_$i\n$seq\n";
    }
    close OT;
    print `cp tmp.$pre.fa ./flank/$id.fa` if $kept ==1;; 
#exit if $id eq "AMPL1075249";
#next;
#    print `grep -A 1  $id shann > tmp_reads.fa `;
    print `grep -A 1  $id $reads > tmp_reads.$pre.fa `;

    print `formatdb -i tmp.$pre.fa -p F`;
    print `blastall  -i  tmp_reads.$pre.fa -d  tmp.$pre.fa -o tmp.vs.tmp_reads.$pre.bst -p blastn -W 9 -m 8 -F F -a 2 `;
#    print `cat tmp.vs.tmp_reads.bst >> shann.vs.singleFlank.bst.W9`; 
    print `cat tmp.vs.tmp_reads.$pre.bst >> $outF`;

}
