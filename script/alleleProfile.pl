#!/usr/bin/perl -w
use strict;

my $f = shift @ARGV; #input file
my %h=();
my %h_ref=();

my $sam =(split /\./,$f)[0]; 
open IN, "$f" or die "cannot open $f\n";
while(<IN>)
{
	chomp;
	my($aid,$refGt,$pos,$sid,$gt) =(split /\s+/,$_);
	my@pos =(split /\,/,$pos);
	my@gt =(split /\,/,$gt); 
	my@REFGT = (split /\,/,$refGt);
	$sid =~/\_x(\d+)/; 
	my$abd =$1;
	for(my$i=0;$i<=$#gt;$i++)
	{
		
		unless($gt[$i] eq "NA0")
		{
			$h{$aid}{$pos[$i]}{$gt[$i]} +=$abd;
			$h_ref{$aid}{$pos[$i]}=$REFGT[$i];
		}
	}
} 

close(IN);

foreach my $aid(keys %h)
{
	foreach my $pos(keys %{$h{$aid}})
	{
		my @tmp = (sort {$h{$aid}{$pos}{$b}<=>$h{$aid}{$pos}{$a}} keys %{$h{$aid}{$pos}});
		my $stutter = 0;
		if(@tmp>1)
		{
			$stutter = $h{$aid}{$pos}{$tmp[1]}/$h{$aid}{$pos}{$tmp[0]};
		}
		for(my $i =0;$i<scalar(@tmp);$i++)
		{
			if($i==0)
			{
				print "$aid\t$pos\t$h_ref{$aid}{$pos}\t$stutter\t$tmp[$i]:$h{$aid}{$pos}{$tmp[$i]}";
			}else
			{
				print "\t$tmp[$i]:$h{$aid}{$pos}{$tmp[$i]}";
			}
		}
		print "\n";

	}

}


