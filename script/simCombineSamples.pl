#!/usr/bin/perl -w
# created by Lun @210.42.74.54 on 2016/12/26
# take two samples as parents 
# only consider ssrs with more than THRESHOLD_MAIN_COVERAGE reads supporting the main genotype and stutter ratio < THRESHOLD_STUTTER
# exclude impure ssrs


die "USAGE: perl $0 <coverage_threshold> <stutter threshold> <genotype_profile1> <genotype_profile2> " unless @ARGV ==4;

my ($THRESHOLD_MAIN_COVERAGE,$THRESHOLD_STUTTER,@files) = @ARGV;

use strict;
my %h = ();
my %h_motif=();
my %h_sample = ();

foreach my $f(@files)
{
	open IN, "$f" or die "cannot open $f\n";
	while(<IN>)
	{
		chomp;
		my ($aid,$pos,$motif,@geno) = (split/\s+/,$_);
		$motif = uc($motif);
		if(length($motif)>=5)
		{
			next;
		}
		my $rc = 0;
		my ($mainNum,$secondNum)=(0,0);
		my $ssr = "$aid\_$pos";
		
		$h_motif{$ssr}=$motif;
		
		my $res = "";
		for(my $i = 0;$i < @geno;$i++)
		{
			my ($g,$num) = (split /:/,$geno[$i]);
			if($i==0)
			{
				last if($g=~m/\_/); #if main genotype is impure, then discard this SSR
				$mainNum = $num;
			}
			if($i==1)
			{
				$secondNum = $num;
				last if($g=~m/\_/); #if second genotype is impure, then discard this SSR
			}
			if($i > 1)
			{
				next if($g=~m/\_/); #if genotype is impure, then remove those reads
			}
			
			$rc += $num;
			
			$g =~s/[ATCG]{1,}//ig;
			my $tem = "$g," x $num;
			if(exists $h{$ssr})
			{
				
				$h{$ssr}="$h{$ssr}$tem";
			}else
			{
				$h{$ssr}=$tem;
			}
		}

		if($mainNum >= $THRESHOLD_MAIN_COVERAGE and $secondNum/$mainNum < $THRESHOLD_STUTTER)
		{
			if(exists $h_sample{$ssr})
			{
				$h_sample{$ssr}++;
			}else
			{
				$h_sample{$ssr}=1;
			}
		}
	}
	
	close(IN);
}


foreach my $ssr (keys %h_sample)
{
	next unless $h_sample{$ssr}==2;
	$h{$ssr}=~s/\,$//;
	$h_motif{$ssr} = uc($h_motif{$ssr});
	print "$ssr\t$h{$ssr}\t$h_motif{$ssr}\n";
}
