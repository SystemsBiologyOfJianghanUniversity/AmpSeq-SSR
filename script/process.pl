#!/usr/bin/perl -w 
use strict;

die("Usage: perl $0 <stat> <otPre> \n
        Used to stat the ssr type and abd(reads num).
        stat: ssr type for each reads. it is generate by       extract4SSR.pl.
        otPre: prefix for output.\n
        
        Revised: to adapte new output format by extractSSR.pl, lines 29  has been changed at 20160701.") unless @ARGV ==2;

my($inputF,$otP)=@ARGV;

#open IN, " cc.shann " or die  "cannot open cc.shann $!";
open IN, " $inputF " or die  "cannot open$inputF $!";
open OT,"> $otP.stat " or die "cannot write $otP.stat $!";

my%h=();
while(<IN>){
    chomp;
#seq21_x1_AMPL1122747 assigned to amplicon AMPL1122747, which inlude two type ssr,GC and ACG;
#two type ssr, GC4 and ACG6, detected in seq21_x1_AMPL1122747, and no complex SSR. if two ssr
# is away 20bp from each other, then they were seem as complex ssr. complex ssr were seperated 
#by comma in the last column.
#AMPL1122747     GC,ACG  seq21_x1_AMPL1122747            GC4,|   ACG6,|
#AMPL1122071     TATC    seq23_x1_AMPL1122071            TATC4,|
#AMPL1121343     TC      seq24_x1_AMPL1121343            TC5,|

#my($aid,$ssr,$seqID,$type)=(split /\s+/,$_,4)[0,1,2,3];
    my($aid,$ssr,$seqID,$type)=(split /\s+/,$_,5)[0,1,3,4];
#print STDERR "$type |$_\n";
    next unless $type=~/\w+\d+/; 
    $type=~s/\,//g;
    $type=~s/\t+//g;
    $type=~s/\|//g;
    $seqID=~/\_x(\d+)/;
    $h{$aid}{$type}+=$1;
}  
foreach my$k(keys %h){
    foreach (keys %{$h{$k}}){ 
        print OT "$k\t$_\t$h{$k}{$_}\n";
    } 
}
close IN;
close OT;

%h=();
open IN , " $otP.stat "  or die  "cannot open $otP.stat $!";
open OT, "> $otP.stat.stat " or die "cannot write $otP.stat.stat $!";
while(<IN>){
    chomp;
    my($id,$ssr,$abd)=(split /\s+/,$_)[0,1,2];
    $h{$id}{$abd}{$ssr} ="";

}

foreach my$k(keys %h){ 
    my$n=0;  
    my$ot="$k\t";
    foreach my$k2(sort {$b<=>$a} keys %{$h{$k}}){
        foreach (sort keys %{$h{$k}{$k2}}){
            $n++; 
            if($n<=2){
                $ot .="$_\t$k2\t";
                next;
            } 
            $ot .="$_,$k2;"; 
        }
    } 
    $ot=~s/\,$//;
    print OT "$ot\n";
}
close IN;
close OT;
