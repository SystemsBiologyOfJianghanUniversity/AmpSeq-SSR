#!usr/bin/perl -w
use strict;
use Getopt::Long;

my$usage="Usage: perl $0 -f < stat > -c <method> -m <method2> -r <ref> -p <primer> -o <prefix> -s <reads> -v \n
        Used to add the amplicon length.
        stat: 
        ref:  reference of each amplicon, [single.fa].
        primer:  the primer sequence used for amplicom sequecing.

        method2: if the flank sequece of  both 5' and 3' end shoulde be cheked. If it is required, 
                we check the che uniqueness of 30bp flank sequence nearby each primer by blast against
                amplicon reference, and if both end are unique, it will be used as query seq to dertermine
                if each reads have this end sequences, i.e, sequeced fully.
                are unique
        method: how to calculated the amplicon length, two methods are available:
            1. we dicard all reads that not same as the major SSSR of this amplicon, and then 
               discard  20% longest and shortest reads, the average length of remain reads were
                amplicon length.
            2. We identified the most abundant length and set that length to be the amplicon's
               length. If there two or more equal abundant length, then their average length
               was set to be the amplicon's length.

            Attention: both methods were adatped to fully sequenced reads only. If the SSR 
                genotype included NA, then it deem as partially sequenced reads and be excluded
                from calculattion of amplicon length.
        \n";

my($stat,$method,$primer,$ref,$method2,$prefix,$reads,$verbose)=();
GetOptions(
        'v+' => \$verbose,
        'p=s' => \$primer,
        'r=s' => \$ref,
        'c=i' => \$method,
        'f=s' => \$stat,
        'm=s' => \$method2,
        'o=s' => \$prefix,
        's=s' => \$reads,
);

die("$usage") if $verbose;
die("$usage") unless defined$stat ;

$method =1 unless (defined $method);
$method2=1 unless (defined $method);


my$pre=(split /\./,$stat)[0];

my(%h,%l)=();
open IN , "$stat.xls" or die  "cannot open $stat.xls $!"; 
while(<IN>){
    chomp;
    my($id,$gt)=(split /\s+/,$_)[0,3]; 
    next if $gt =~/NA/; 
    $h{$id}=$gt; 
    $l{$id}=$_;
}
close IN;

my%len=();
open IN, " $pre.Run100.combin.fa ";
$/=">"; <IN>; 
while(<IN>){
    chomp; 
    my($id,$seq)=(split /\n+/,$_); 
    $len{$id}=length$seq;
}
close IN; 

if($method2 ==1){
    %len=();
    my(%primerLen,%flankLen,%uniqueFlankID)=();
    
    open IN , " $primer " or die  "cannot open $primer $!";
    $/=">"; <IN>;
    while(<IN>){
        chomp;
        my($pid,$pseq)=(split /\n+/,$_)[0,1];
        $primerLen{$pid}=length$pseq;
    }
    close IN;
    $/="\n";

    open IN, " $ref " or die  "cannot open $ref $!";
    open OT, "> $prefix\_bothEnd.fa " or die  "cannot write  $prefix\_bothEnd.fa $!";
    $/=">"; <IN>;
    while(<IN>){
        chomp;
        my($rid,$rseq)=(split /\n+/,$_)[0,1];
        my$sub1=substr($rseq, $primerLen{"$rid\_F"},30);
        my$sub2=substr($rseq, ((length$rseq)-$primerLen{"$rid\_R"}-30),30);
        print OT ">$rid\_5end\n$sub1\n>$rid\_3end\n$sub2\n";
    }
    $/="\n";
    close IN;
    close OT;

    my$refName=(split /\//,$ref)[-1];
    print `formatdb -i $ref -p F `  unless -f "$ref.nsq";
    print `blastall -i $prefix\_bothEnd.fa -d $ref -o $prefix\_bothEnd.vs.$refName.bst -p blastn  -m 8 -F F -b 10`;

    open IN , " $prefix\_bothEnd.vs.$refName.bst " or die  "cannot open $prefix\_bothEnd.vs.$refName.bst $!";
    while(<IN>){
        chomp;
        my($qid,$sid,$identity,$mlen)=(split /\s+/,$_)[0,1,2,3];
        next if $mlen<=25;
        next if $identity <= 0.9;
        next unless $qid=~/$sid/;

        $uniqueFlankID{$sid} +=1;
    }
    close IN;

    foreach (keys %uniqueFlankID){
        delete  $uniqueFlankID{$_} if $uniqueFlankID{$_} !=2;
    }


    print `formatdb  -i $prefix\_bothEnd.fa -p F` unless -f "$prefix\_bothEnd.fa.nsq";
    print `blastall  -i $reads -d $prefix\_bothEnd.fa -o $reads.vs.$prefix\_bothEnd.bst -p blastn -m 8 -F F -b 20 -a 3`;
    
    my%hash=();
    my$lastqid="";
    
    open IN  , " $reads.vs.$prefix\_bothEnd.bst " or die  "cannot write $reads.vs.$prefix\_bothEnd.bst $!" ;
    while(<IN>){
        chomp; 
        my($qid,$sid,$identity,$mlen,$qs,$qe)=(split /\s/,$_)[0,1,2,3,6,7];
        my$AmpliconID =(split /\_/,$qid)[-1];
        
        next unless exists $uniqueFlankID{$AmpliconID}; ## check if flank seqs at both end are unique;
        next if $mlen<=25;
        next if $identity <= 0.9;
        next unless $qid=~/$AmpliconID/;
        next unless exists $uniqueFlankID{$AmpliconID};
        
        if($lastqid eq ""){
            $lastqid =$qid;
            $hash{$sid}{$mlen}{$identity}{$qs}=$qe;
            next;
        }

        if($lastqid eq $qid){
            $hash{$sid}{$mlen}{$identity}{$qs}=$qe;
            next;
        }

        $lastqid=~/\_(\D+\d+)$/;
        my$aID=$1;

        my$flag=0;
        $flag++ unless exists $hash{"$aID\_3end"};
        $flag++ unless exists $hash{"$aID\_5end"};

        if($flag !=0){
            %hash=();
            $hash{$sid}{$mlen}{$identity}{$qs}=$qe;
            $lastqid=$qid;
            next;
        }

        #two suitation should be checked, reads is same or not same as referece;
        my($minS1,$MaxS1,$minE1,$MaxE1)=(-1,-1,-1,-1);
        my$Mlen =(sort {$b<=>$a} keys %{$hash{"$aID\_5end"}})[0];
        my$Iden =(sort {$b<=>$a} keys %{$hash{"$aID\_5end"}{$Mlen}})[0];
        ($minS1,$MaxS1)=(sort  {$a<=$b}  keys %{$hash{"$aID\_5end"}{$Mlen}{$Iden}})[0,-1];
        $minE1=$hash{"$aID\_5end"}{$Mlen}{$Iden}{$minS1} if $minS1 !=-1;
        $MaxE1=$hash{"$aID\_5end"}{$Mlen}{$Iden}{$MaxS1} if $MaxS1 !=-1;

         my($minS2,$MaxS2,$minE2,$MaxE2)=(-1,-1,-1,-1);
         $Mlen= (sort {$b<=>$a} keys %{$hash{"$aID\_3end"}})[0];
         $Iden= (sort {$b<=>$a} keys %{$hash{"$aID\_3end"}{$Mlen}})[0];
         ($minS2,$MaxS2)=(sort  {$a<=$b}  keys %{$hash{"$aID\_3end"}{$Mlen}{$Iden}})[0,-1];
         $minE2=$hash{"$aID\_3end"}{$Mlen}{$Iden}{$minS2} if $minS2 !=-1;
         $MaxE2=$hash{"$aID\_3end"}{$Mlen}{$Iden}{$MaxS2} if $MaxS2 !=-1;

         if($minS1 <$minS2){### indicated that reads is the same direction with reference;
            if( ($minS1 !=-1) and  ($MaxS2 !=-1)){
                 $len{$lastqid}=$MaxE2-$minS1+1;
             }elsif(($minS1 !=-1) and  ($MaxS1 !=-1)){
                 $len{$lastqid}=$MaxE1-$minS1+1;
             }else{die"No 3' end location found.\n";}
         }else{
             if(($minS2 !=-1) and ($MaxS1 !=-1)){
                 $len{$lastqid}=$MaxE1-$minS2+1;
             }elsif(($minS2 !=-1) and ($minS1 != -1)){
                 $len{$lastqid}=$minE1-$minS2+1;
             }else{die"No 5' end location found.\n$minS1,$MaxS1,$minE1,$MaxE1\n$minS2,$MaxS2,$minE2,$MaxE2\n";}
         }

         $len{$lastqid} += $primerLen{"$aID\_F"}+$primerLen{"$aID\_R"};
         %hash=();
         $hash{$sid}{$mlen}{$identity}{$qs}=$qe;
         $lastqid=$qid;
    }
   close IN;
   $/="\n"; 

#   print `rm -rf $prefix\_bothEnd.fa $prefix\_bothEnd.fa.nsq $prefix\_bothEnd.fa.nin $prefix\_bothEnd.fa.nhr $prefix\_bothEnd.vs.$refName.bst $reads.vs.$prefix\_bothEnd.bst`;
}


$/="\n";  
my%k=();
open IN , " $stat " or die  "cannot open $stat $!";; 
while(<IN>){
    chomp;  
    my($aid,$sid,$gt)=(split /\s+/,$_)[0,3,4];
    next unless exists $h{$aid};
    if($method !=1){
        next if $gt ne $h{$aid};
    } 
    $k{$aid}{$sid}="";
}


open OT, "> $pre.stat.addLen.xls " or die  "cannot write $pre.addLen.xls $!";
print OT "ampliconID\t#totalreads\taverLen\t#naReads\tmajorSSR\tmajorAbd\tsecondSSR\tsecondAbd\tstutter\t#readsOfOtherSSR\n";
foreach my$aid( keys %k){

    my@len=();
    my@ids=keys %{$k{$aid}}; 
    my%count=();

    foreach(@ids){
        next unless exists $len{$_};
        my$len=$len{$_}; 
        $_=~/\_x(\d+)/; 
        my$abd=$1;
        $count{$len} +=$abd;

        for(my$i=1;$i<=$abd;$i++){
            push @len,$len;
        }
    }
    
    next if @len ==0; ### no reads have length infor for this amplicon;


    my$total=@len; 
    my$stat=0.2*$total ; 
    $stat=int$stat; 
    my$end=0.8*$total;  
    $end=int$end; 
    @len=(sort {$a<=>$b} @len);
    my$aver=0;
    for(my$j=$stat;$j<=$end;$j++){ 
        $aver+=$len[$j];
    }
   
    $aver=$aver/($end-$stat+1);
    $aver=int$aver; 

    if($method==1){
        my@lenths=(sort {$count{$b}<=>$count{$a}} keys %count);
        my$nn=0;
        $aver=0;
        for(my$ii=0;$ii<=$#lenths;$ii++){
            last if $count{$lenths[$ii]}<$count{$lenths[0]};
            $aver +=$lenths[$ii] if ($count{$lenths[0]}==$count{$lenths[$ii]});
            $nn++;
        }
        $aver=int($aver/$nn) if $nn !=0;
    }

    (my$nul1,$total,my@tmp)=(split /\s+/,$l{$aid}); 
    pop@tmp; my$line=join "\t",@tmp;
    print OT "$aid\t$total\t$aver\t$line\n"; 
}  

