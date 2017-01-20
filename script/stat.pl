#!/usr/bin/perl -w
use strict;

my@sam=();
my(%ssr,%ra)=();
foreach my$f(@ARGV){
    my$sam=(split /\./,$f)[0];  
    push @sam,$sam;
    open IN, "$f" or die  "cannot open $f $!";; 
    while(<IN>){
        chomp;
        my($id,$ssr,$abd1,$abd2)=(split /\s+/,$_)[0,1,2,4]; 
        $ssr{$id}{$sam}=$ssr;
        my$total=$abd1;
        if($abd2){
            $ra{$id}{$sam}=($abd2/$abd1); 
            $total+=$abd2;
        }else{
            $ra{$id}{$sam}=1;
        }
        
        while($_=~s/\,(\d+)\;//){
            $total+=$1;
        }
        $ra{$id}{$sam}=$total;
    }
    close IN; 
} 
my$title=join"\t",@sam; 
print "\t$title\n"; 
print STDERR "\t$title\n";
foreach my$k(sort keys %ssr){ 
    print "$k"; 
    foreach (@sam){
        print "\t$ssr{$k}{$_}" if exists $ssr{$k}{$_}; 
        print "\tNA" unless exists $ssr{$k}{$_};
    }
    print "\n"; 
}
foreach my$k(sort keys %ra){
    print STDERR "$k"; 
    foreach (@sam){
        print STDERR "\t$ra{$k}{$_}" if exists $ra{$k}{$_};
        print STDERR"\tNA" unless exists $ra{$k}{$_};
    } 
    print STDERR "\n"; 
}
