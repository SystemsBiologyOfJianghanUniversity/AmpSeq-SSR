#!/usr/bin/perl -w
use strict;

=cut
seq7_x1 chr4    97.50   40      1       0       1       40      40691   40730   5e-15   71.9
seq7_x1 chr4    100.00  21      0       0       1       21      48959   48979   4e-06   42.1
seq7_x1 chr4    100.00  21      0       0       1       21      51325   51305   4e-06   42.1
seq7_x1 chr4    100.00  20      0       0       4       23      859     840     2e-05   40.1
seq7_x1 chr4    100.00  20      0       0       2       21      10446   10465   2e-05   40.1
seq7_x1 chr4    100.00  20      0       0       2       21      31224   31243   2e-05   40.1
seq7_x1 chr4    100.00  20      0       0       2       21      49706   49725   2e-05   40.1
seq7_x1 chr4    100.00  19      0       0       4       22      3786    3804    7e-05   38.2
seq7_x1 chr4    100.00  19      0       0       4       22      11505   11487   7e-05   38.2
seq7_x1 chr4    100.00  19      0       0       4       22      32650   32632   7e-05   38.2
seq7_x1 chr4    100.00  19      0       0       3       21      40061   40079   7e-05   38.2
seq7_x1 chr4    100.00  19      0       0       3       21      58302   58284   7e-05   38.2
seq7_x1 chr4    100.00  18      0       0       4       21      616     599     3e-04   36.2
seq7_x1 chr4    100.00  18      0       0       4       21      618     601     3e-04   36.2
seq7_x1 chr4    100.00  18      0       0       4       21      620     603     3e-04   36.2
seq7_x1 chr4    100.00  18      0       0       4       21      622     605     3e-04   36.2
seq7_x1 chr4    100.00  18      0       0       4       21      861     844     3e-04   36.2
seq7_x1 chr4    100.00  18      0       0       4       21      863     846     3e-04   36.2

=cut 

#open IN, " 9031.vs.single.bst ";
#open IN, " 9031_prefix2.vs.single.bst ";
#open IN, " 9031_prefix3.vs.single.bst ";
my($bst)=@ARGV;

open IN, " $bst ";
my$l="";
while(<IN>){
	chomp;
	my($tmps,$tmpe)=(split /\s+/,$_)[6,7];
	next if abs($tmpe-$tmps)<30;
	$l=$_;
	last;
}

if($l eq ""){
    print STDERR "No hits longer than 30bp were found, please check!\n";
    exit;
}

my($id,$chr,$qs,$qe,$ss)=(split /\s+/,$l)[0,1,6,7,8];
my$best= $qe-$qs+1;
my$flag=0;
my$len=0;
my%h=();

my($tid,$tchr,$tqs,$tqe,$tss)=();
while(<IN>){
	chomp;
	($tid,$tchr,$tqs,$tqe,$tss)=(split /\s+/,$_)[0,1,6,7,8];
	$len=$tqe-$tqs+1;
    next if $len<30;

	if($tid eq $id){
#		next if $flag==1;
		next if $tchr eq $chr;
        next if ($len-$best)<10;
        $chr=$tchr;
        $best= $len;

#		next if abs($best-$len)>=10;
#		$flag=1;
#		next;
	}

	$h{"$id\_$chr"}="";

	$id=$tid;
#	$flag=0;
	$chr=$tchr;
	$best= $len;
}
$h{"$id\_$chr"}="";
close IN;


#open IN, " 9031.vs.single.bst ";
#open IN, " 9031_prefix3.vs.single.bst ";
open IN, " $bst ";
while(<IN>){
	chomp;
	($id,$chr)=(split /\s+/,$_)[0,1];
	next unless exists $h{"$id\_$chr"};
	print "$_\n";
	delete $h{"$id\_$chr"};
}	
