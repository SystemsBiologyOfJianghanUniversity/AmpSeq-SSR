#!usr/bin/perl -w
use strict;
use File::Basename;

die("Usage: perl $0 <fa> <ref> <pre> \n
        Used to locate reads to different ref and compare with plugin\n")  unless @ARGV ==3;


my$dir=dirname($0);
my($fa,$ref,$pre)=@ARGV;

my$fList="";
for(my$i=1;$i<=100;$i++){
	
	my$subS=40*($i-1);
	my$short=40*$i;
	$short=80 if $short<80;
	
	my%h=();
	for(my$j=1;$j<$i;$j++){
		last if $i==1;
#		open IN," 9031_prefix$j.vs.single.bst.f";
        open IN," $pre\_prefix$j.vs.single.bst.f";
		while(<IN>){
			chomp;
			$_=(split /\s+/,$_)[0];
			$h{$_} =""; 
		}
		close IN;
	}

#	open IN, "9031.fa " ;
#    open IN," IonXpress_068_rawlib.fa " or die  "cannot open IonXpress_068_rawlib.fa $!";
     open IN," $fa  " or die  "cannot open $fa $!";
	open OT,">$pre\_prefix$i.fa";
	$/=">"; <IN>; 
	while(<IN>){
		chomp;
		my($id,$seq)=(split /\n+/,$_);
		next if exists $h{$id};
		next if (length$seq)<$short;
		my$sub=substr($seq, $subS,40);
		print OT ">$id\n$sub\n";
	}
	close IN;
	close OT;
	$/="\n";
	%h=();

	print `blastall -i $pre\_prefix$i.fa -d $ref -o $pre\_prefix$i.vs.single.bst -p blastn -m 8 -F F -a 6 -e 1e-2 -b 10 `;
	print `perl $dir/cc.pl $pre\_prefix$i.vs.single.bst >$pre\_prefix$i.vs.single.bst.f`;
 
	open IN, " $pre\_prefix$i.vs.single.bst.f";
	while(<IN>){
		chomp;
		my($id,$aid)=(split /\s+/,$_)[0,1];
		$h{$id}=$aid;
	 }
	close IN;
	
	$fList.="$pre\_aid$i.fa ";
	open IN , " $fa"; 
	open OT, "> $pre\_aid$i.fa";
	$/=">"; <IN>;
	while(<IN>){
		chomp;
		my($id,$seq)=(split /\n+/,$_); 
		next unless exists $h{$id}; 
		print OT ">$id\_$h{$id}\n$seq\n";
	}
	close IN;close OT;

	my(%m,%c,%co)=();
	%h=();
=cut
	open IN, "$plugin " or die  "cannot open $plugin $!";
	$/=">"; <IN>; 
	while(<IN>){
		chomp;
		my($id, $seq)=(split /\n+/,$_)[0,1]; 
		$id=(split /\_/,$id)[-1];
		$m{$id}="";
		$c{$id} +=1;
		$h{$id}{$seq}= "";
	}
=cut	
	print `cat $fList >$pre.Run$i.combin.fa`;
	my$f="$pre.method2Unique$i"."run.fa";
	my$st="$pre.abd.stat4method$i"."run";
	open IN, " $pre.Run$i.combin.fa "; 
	open OT, "> $f";
	
	$/=">"; <IN>;
	while(<IN>){
		chomp;
		my($id, $seq)=(split /\n+/,$_)[0,1];
		$_=~/\_x(\d+)/; 
		my$n=$1; 
		$id=(split /\_/,$id)[-1]; 
		$co{$id} +=$n; 
		$m{$id}="";
		print OT ">$_" unless exists $h{$id}{$seq}; 
	}
	
	open STAT, "> $st";
	foreach (keys %m){ 
		$c{$_} =0 unless exists $c{$_}; 
		$co{$_}=0 unless exists $co{$_}; 
		print STAT  "$_\t$c{$_}\t$co{$_}\n";
	}
    $/="\n";
}

`rm -rf   $pre\_prefix* $pre\_aid* $pre.method2Unique* $pre.abd.stat4method* $pre.Run[0-9].combin.fa $pre.Run[0-9][0-9].combin.fa`;
