#!/usr/bin/perl -w
use strict;


die "USAGE perl $0 <reference ssr file> <reference fasta>\n" unless @ARGV >= 2;


my ($ssrFile,$refFile) = @ARGV; 


open SSR_IN, "$ssrFile" or die "cannot open SSR file!\n";
my %h=();
while(<SSR_IN>)
{
	chomp; 
	my ($id,$s,$e)=(split /\s+/,$_)[0,5,6]; 
	
	$s -=1 ;
	$h{$id}{$s}=($e-$s) ;
}  

open IN, "$refFile" or die "cannot open $refFile!\n";
$/=">"; 
<IN>; 
while(<IN>)
{
	chomp; 
	my ($id,$seq)=(split /\n+/,$_)[0,1]; 
	unless (exists $h{$id}){ print ">$_"; next ; } 
	foreach my$s(keys %{$h{$id}})
	{ 
		my $len =$h{$id}{$s}; 
		$seq=~/(\w{$s})(\w{$len})(\w*)/; 
		$len=length$2; 
		my $n="N" x $len; 
		$seq="$1$n$3";
	} 
	print ">$id\n$seq\n";
}

close(SSR_IN);
