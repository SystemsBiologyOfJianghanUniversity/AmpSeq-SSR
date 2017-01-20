#!/usr/bin/perl -w
use strict;

die ("Usage: perl $0 <ref.ssr> <query2ref.bst> <query>  <otPre> <refF> <Primer> <flankLen>\n
        Used to extract ssr based on ssr repeat
        ref.srr: ssr list generate by ssr.pl for reference. One line for one ssr amplicon region.
        query2ref.bst: result that query seq mapped on flank seq of ssr amplicon, this is in m8 formate.
        query: this is the reads seq that used for alignment against refsrr flank seq.
        refF : ref file in fasta file, which generate ref.ssr.
        Primer: the primer seq designed for amplicon by ampliSeq.
        flankLen: the minium length of flank sequence for each SSR; this should be KEPT in line with getSubseqBetweenSSR.pl; 

        Attention:
            1.the best hits determine if the reads is in same direction with reference ssr flank Seq.
              all other local hit of reads should be in the same direction ( by \$direction).\n
            2. if two repeat of a complex SSR are same, a WRONG SSR will be called, such as:
                AG[4]aaagaAG[17] will called as AG[17]aaagaAG[17];
               another bad case is: for ref like this, tatc[3]TAtatc[5], a read with 5 repeat tatc will
               be call to tatc[2];

        Revised:
            we ignored all SSR located on the primer regions by add line 90, considering that such SSR will be digested in library construction. @ 20160720.
            we still cannot process SSR that locate at both end of reads.

        Revised:
            we reformated the otput to present each SSR include reads, if NO ssr found it marked as 0 (such as AT0), if not found due to only part of amplicon were sequence, marker as NA0; at 20160704;

        Revised:
            removed the \",\\t\" in the output. 20160702 


        Revised:
            1. the position on ref for each SSR were added, now the formate of output is as follows:
                #amplicon ID    SSRtypeOnRef            ssrPosOnRef             readsID                         SSR1                    ssr2    SSR3
                AMPL1122869     ATCT,TTCT,ATCT,ATCT     1_1,1_2,2_1,3_1,        seq555_x303_AMPL1122869         ATCT6,  TTCT9,  |       ATCT5_| ATCT4_|
                ssrPosOnRef: 1_1 means this is the fisrt ssr (ATCT6) in a compond SSR(SSR1) on reads started from left (5' end) of reference; 1_2 means the second one (TTCT9) in compond SSR;

            2. to be sure all ssr were called started from the most left, the keys(pos of ssr) should be sorted at lines 302; at 20160630;


        Revised: the direction have take account in when check if two hits are overlap for query, and also for subject seq;

                for non-complex SSR, we concentate all ssr with >=4 ssr repeat units, lines 385-392; at 20160428;
        
        The \$otPre.stat could be process by process.pl in forthur.") unless @ARGV ==7;

my($refSsr,$bst,$reads,$otPre,$refF,$primerF,$flankLen)=@ARGV;

########## recorde length of primer
my%pLen=();
open IN, " $primerF " or die  "cannot open $primerF $!";
$/=">";<IN>;
while(<IN>){
    chomp;
    my($id,$seq)=(split /\n+/,$_)[0,-1];
#    $pLen{$id}=length $seq;
    $pLen{$id}=rindex($seq,"T");#### all primer will be digested at T base, so the last T will be used as end of flank seq;
    $pLen{$id}+=1;
}
close IN;
$/="\n";

##########  record length of reads;
my%readLen=();
open IN, " $reads "  or die "cannot open $reads $!";
$/=">"; <IN>;
while(<IN>){
    chomp;
    my($id,$seq)=(split /\n+/,$_)[0,1];
    my$reLength=length$seq;
    $readLen{$id}=$reLength;
}
close IN;
$/="\n";

############  get ssr repeat unit along the 2000 loci;
###  record pos for each ssr repeat;
my%repeat=();
#open IN , " single.ssr " or die "cannot open IN sing.ssr $!";
open IN, " $refSsr " or die  "cannot open $refSsr $!";
while(<IN>){
    chomp;
#AMPL1122307     1       3       gga     5       142     156     273
#AMPL1122307     2       3       gag     6       167     184     273
#AMPL1122307     3       3       gga     7       190     210     273
#AMPL1122307     4       3       gac     6       209     226     273


    my($id,$repeat,$s,$e,$rlen)=(split /\s+/,$_)[0,3,5,6,7];
#    next if ($s<=$pLen{"$id\_F"} or $e>=($rlen-$pLen{"$id\_R"}));# filter all ssr located in primer region;
    $repeat{$id}{$s}="$e,$repeat";
}
close IN;

print STDERR "Finish reading ssr pos ....at".`date`;

open IN, " $refF " or die  "cannot open $refF $!";
$/=">";<IN>;
while(<IN>){
    chomp;
    my($id,$seq)=(split /\n+/,$_)[0,1];
    foreach my$pos(keys %{$repeat{$id}}){
        my($posEnd,$unit)=(split /\,/,$repeat{$id}{$pos})[0,1];

        my($tmpS,$tmpE)=(0,0);
        my$left=substr($seq,0,($pos-1));
        if($left=~/(($unit){1,}\w{0,5})$/){
#print "OK:$1 $2\n$unit \n$left\n" if $id=~/AMPL1118676/;
            $left=~/(($unit){1,}\w{0,5})$/;
#print "ok $1\n|$2\n$left\n$unit\n" if $id=~/AMPL1118676/;
            $tmpS=$pos-(length$1);
        }

        my$right=substr($seq,$posEnd);
        if($right=~/^(\w{0,5}($unit){1,})/){
            $tmpE=$posEnd+(length$1);
#print "OK\n" if $id=~/AMPL1118531/;
        }
        delete$repeat{$id}{$pos} if (($tmpS !=0) or ($tmpE !=0));
#print "$pos $posEnd $tmpS $tmpE $unit\n$left\n$right\n" if $id=~/AMPL1118531/;
        my$pos1=$pos;
        $pos1=$tmpS if $tmpS !=0;
        my$posEnd1=$posEnd;
        $posEnd1=$tmpE if $tmpE !=0;
        $repeat{$id}{$pos}="$pos1,$posEnd1,$unit";
#        $pos=$tmpS if $tmpS !=0;
#        $posEnd=$tmpE if $tmpE !=0;
#        $repeat{$id}{$pos}="$posEnd,$unit";
    }
}
close IN;
$/="\n";


### record ssr repeat unit along 5'->3' of a loci;
my%pos2ssr=();
foreach my$aID(keys %repeat){
    my$i=0;
    my$j=0;
    my$s=0;

    my($firstS_old,@allS)=(sort {$a<=>$b} keys %{$repeat{$aID}});
    (my$firstS,$s)=(split /\,/,$repeat{$aID}{$firstS_old})[0,1];

    if(($firstS-$pLen{"$aID\_F"})>=$flankLen){
        $i++;
        $j=1;
        $pos2ssr{"$aID\_$i"}{$j}=(split /\,/,$repeat{$aID}{$firstS_old})[-1];
    }
#print "start: $firstS $s $i:$j\n";

    foreach my$start_old(@allS){
        my($start,$ssrUnit_t)=(split /\,/,$repeat{$aID}{$start_old})[0,-1];
        if(($start-$s)>=$flankLen){## if two ssr is less than 20 bp away from each other, they will deem as complex SSR;
            $i++;
            $j=1;
            $pos2ssr{"$aID\_$i"}{$j}=$ssrUnit_t;
            $s=(split /\,/,$repeat{$aID}{$start_old})[1];
#print "complex: $start $s $i:$j\n";
            next;
        }
        $j++;
#        $i=1 if $i==0; #### prevent that one ssr start at primer region, the next one start from the first one,  so the $i will be 0; this lead to disagreement between the sequenal No of flakseq and SSR, so this should be comment;
        $pos2ssr{"$aID\_$i"}{$j}=$ssrUnit_t;
        $s=(split /\,/,$repeat{$aID}{$start_old})[1];
#print"NotComp: $start $s $i:$j\n";
    }
}

### used for get all ssr infor for all SSR in each reads;
my%pos2SSR=();
foreach my$key(keys %pos2ssr){
    my($ampliconID_tmp,$No)=(split /\_/,$key)[0,1];
    foreach (keys %{$pos2ssr{$key}}){
        $pos2SSR{$ampliconID_tmp}{$No}{$_}=$pos2ssr{$key}{$_};
    }
}

print STDERR "Finish record ssr unit ....at".`date`;

#foreach my$k(keys %pos2ssr){ foreach (keys %{$pos2ssr{$k}}){  print "$k\t$_\t$pos2ssr{$k}{$_}\n" if $k=~/AMPL1118676/;;}}#exit;
#####################  record the pos a flanks seq hits on  reads;
my%extract=();
my%mLen=();
my%qid2aid=();
my%direction=();
my%checkOverlap4Sub=();
#open IN, " shann.vs.singleFlank.bst.W9 " or die  "cannot open shann.vs.singleFlank.bst.W9 $!";
#open IN, " shann.vs.singleFlank.bst.W11 " or die  "cannot open shann.vs.singleFlank.bst.W11 $!";
#open IN, "test.bst " or die "cannto open test.bst $!";
open IN, " $bst " or die  "cannot open $bst $!";
while(<IN>){
    chomp;
    my($qid,$sid,$mlen,$qs,$qe,$ss,$se)=(split /\s+/,$_)[0,1,3,6,7,8,9];
    unless( exists $direction{$qid}){
        $direction{$qid}=0 if $ss<$se;#forward alignment;
        $direction{$qid}=1 if $ss>$se; # reverse alignment;
    }else{
        next if (($direction{$qid}==0) and  ($ss>$se));
        next if (($direction{$qid}==1) and  ($ss<$se));
    }

    if($ss>$se){
        die("No length found for $qid\n") unless exists $readLen{$qid};
        my$change_qs=$qs;
        $qs=$readLen{$qid}-$qe+1;
        $qe=$readLen{$qid}-$change_qs+1;
    }
#print "$qs $qe\n" if $qid eq "seq602_x288_AMPL1122869";
#seq47_x1_AMPL1142081    AMPL1142081_2   96.79   156     3       2       31      185     1       155     1e-72    254

    my$flag=0;
    foreach my$key(keys %{$extract{$qid}}){
        ###preven this case: one region of reads mapped on flank seq 1, but its subregion also could map on another flank.;
        ###seq565_x1_AMPL1141927      AMPL1141927_3   100.00  79      0       0       47      125     1       79      1e-43    157
        ###seq565_x1_AMPL1141927      AMPL1141927_2   100.00  20      0       0       20      39      1       20      2e-08   40.1
        ###seq565_x1_AMPL1141927      AMPL1141927_1   100.00  9       0       0       32      40      114     106     0.081   18.3
        my($tmp_S,$tmp_E)=(sort {$a<=>$b} ($key,$extract{$qid}{$key}))[0,-1];
        $flag=1 if ($qs>=$tmp_S and $qs<=$tmp_E);
        $flag=1 if ($qe>=$tmp_S and $qe<=$tmp_E);
        if(exists $checkOverlap4Sub{"$qid\_$key"}{$sid}){
            ($tmp_S,$tmp_E)=(split /\_/,$checkOverlap4Sub{"$qid\_$key"}{$sid})[0,1];
            $flag=1 if ($ss>=$tmp_S and $ss<=$tmp_E);
            $flag=1 if ($se>=$tmp_S and $se<=$tmp_E);
            if($ss>$se){
                $flag=1 if ($ss>=$tmp_E and $ss<=$tmp_S);
                $flag=1 if ($se>=$tmp_E and $se<=$tmp_S);
            }
        }

#        if($direction{$qid}==0){
#            $flag=1 if ($qs>=$key and $qs<=$extract{$qid}{$key});
#            $flag=1 if ($qe>=$key and $qe<=$extract{$qid}{$key});
#        }elsif($direction{$qid}==1){ # if it hits on the reverse complemtary strand, then qs if greater than qe;
#            $flag=1 if ($qs<=$key and $qs>=$extract{$qid}{$key});
#            $flag=1 if ($qe<=$key and $qe>=$extract{$qid}{$key});
#        }else{
#            die "No direction infor for $qid\n";
#        }
#print "$key\t$extract{$qid}{$key}\t$flag\t$sid\n";
    }
#print "##############\n";

    next if $flag==1;
    
    ($qs,$qe)=(sort {$a<=>$b} ($qs,$qe))[0,-1];
    unless(exists $mLen{$qid}{$sid}){

        $mLen{$qid}{$sid}="$mlen\_$qs";
        $extract{$qid}{$qs}=$qe;
        $qid2aid{$qid}{$qs}=$sid;
        $checkOverlap4Sub{"$qid\_$qs"}{$sid}="$ss\_$se";
        
        next;
    }

    my($tmp_mlen,$tmp_qs)=(split /\_/,$mLen{$qid}{$sid})[0,1];
#    next if $mLen{$qid}{$sid}>=$mlen;  
    next if $tmp_mlen>=$mlen;

    $mLen{$qid}{$sid}="$mlen\_$qs";
    $extract{$qid}{$qs}=$qe;
    $qid2aid{$qid}{$qs}=$sid;
    $checkOverlap4Sub{"$qid\_$qs"}{$sid}="$ss\_$se";

    delete$extract{$qid}{$tmp_qs};
    delete$qid2aid{$qid}{$tmp_qs};
    delete$checkOverlap4Sub{"$qid\_$tmp_qs"}{$sid};
}

#foreach my$k( keys %extract){my$flag=0; $flag++ if (($k eq "seq555_x303_AMPL1122869") or ($k eq "seq602_x288_AMPL1122869")); next if $flag==0;foreach (keys %{$extract{$k}}){ print "$k\t$_\t$extract{$k}{$_}\t$qid2aid{$k}{$_}\n";}} exit;
print STDERR "Finish analysis bst result ....at ".`date`;

############  reads reads  seq ##

#open IN, " shann " or die  "cannot open shann $!";
my%ssr_No=();###record sequence No of ssr on ref;
my$ssrPos_list="";
open IN," $reads"  or die "cannot open $reads $!";
open OT, "> $otPre.stat " or die "cannot write $otPre.stat $!";
$/=">";<IN>;

while(<IN>){
    my($id,$seq)=(split/\n+/,$_)[0,1];

    $ssrPos_list="";

    next unless exists $extract{$id};

    unless (exists $direction{$id}){
        print STDERR "No direction infor found for $id\n";
        next;
    }
    if($direction{$id}==1){
        $seq=~tr/ATCG/TAGC/;
        $seq=reverse $seq;
    }

    next unless exists $extract{$id};
    my@tmp=(keys %{$extract{$id}});
    next if @tmp==1;#only one flank seq hit;
#next unless $id eq "seq409_x4_AMPL1075307";
#print "$id\t@tmp\n"; exit;

    my($lastS,$lastE,$nextS,$nextE)=(0,0,0,0);

    my$ot="$id\t";
    my%aIDlist=();
    my$ssrUnitList="";

    foreach my$k(sort {$a<=>$b} keys %{$extract{$id}}){
#print STDERR  "extract: $k\n";
        if($lastE==0){
            ($lastS,$lastE)=($k,$extract{$id}{$k});
            next;
        }
        ($nextS,$nextE)=($k,$extract{$id}{$k});

        my$sub=substr($seq,$lastE-1,($nextS-$lastE+1));

        die("Wrong no amplicon flank seq for $id") unless exists $qid2aid{$id}{$lastS};
        my$amplID=$qid2aid{$id}{$lastS};
        
        $amplID=~/(\w+)\_(\d+)/;
        $aIDlist{$1}="";
        my$sequen_No=$2;

        my@ssrUnit=();
        %ssr_No=();
        my$number=-1;
        foreach my$n(sort {$a<=>$b} keys %{$pos2ssr{$amplID}}){
            $number++;
            my$tmp_rpu=$pos2ssr{$amplID}{$n};
            $tmp_rpu=uc$tmp_rpu;
            push@ssrUnit,$tmp_rpu;
            $ssrUnitList.="$tmp_rpu,";
            $ssr_No{$number}="$sequen_No\_$n";  
        }

        $ot.="\t";
        my$return=&process($sub,@ssrUnit);
        my($return0,$return1)=(split /\t+/,$return,2)[0,1];
        $ot.=$return1;
        $ssrPos_list.="$return0";
#        print "$ot\n$sub\n@ssrUnit\n";
=cut
        foreach my$n(keys %{$pos2ssr{$amplID}}){
            my$ssrUnit=$pos2ssr{$amplID}{$n};
            $ssrUnit=uc$ssrUnit;
            $ssrUnitList.="$ssrUnit,";
#print STDERR "$ssrUnit\t$n\t$amplID\n";
            unless ($sub=~/$ssrUnit/){
                print STDERR "$amplID|\t$n|\t$ssrUnit|\t$id\t|$sub|\n";
                $ot.="$ssrUnit,NA";
                next;
            }

            my$max=0;
            my$sub2="";
            my$index=index($sub,$ssrUnit);
#print STDERR "$ssrUnit\t$sub\n";
            while(1){### this is to prevent this case: AAG-AAGAAAGAAAGAAAGAAAGA, aaga was repeat unit;
                $sub2=substr($sub,$index);
                $sub2=~/(($ssrUnit)+)/ig;
                my$repeatN=(length$1)/length($ssrUnit);
                $max=$repeatN if $max<$repeatN;
                $index +=((length$1)-(length$ssrUnit)+1);
                $index=index($sub,$ssrUnit,$index);
                last if $index ==-1;
            }

#            while($sub=~/(($ssrUnit)+)/ig){
#                my$repeatN=(length$1)/length($ssrUnit);
#                $max=$repeatN if $max<$repeatN;
#print STDERR "$ssrUnit\t$sub\t$1\n";
#            }
            $ot.= "\t$ssrUnit$max,";
        }
=cut
        $ot=~s/\|$//;$ot.= "|";

        ($lastS,$lastE)=($nextS,$nextE);
    }

    my@shan=keys %aIDlist; print STDERR "$id\t@shan\n" if @shan>1;
    $ssrUnitList=~s/\,$//;
    print OT "$shan[0]\t$ssrUnitList\t$ssrPos_list\t$ot\n";
} 
close OT; 
$/="\n";


my%tmp=();
print `mv $otPre.stat $otPre.stat_shan`;
open IN,  " $otPre.stat_shan " or die  "cannot open $otPre.stat_shan $!";
open OT,"> $otPre.stat "  or die "cannt not write  $otPre.stat $!";

while(<IN>){
    chomp;
    my$line=$_;
#AMPL1118531     CT,TC,GC        1_1,1_2,2_1,    seq7_x4175_AMPL1118531          CT8,    TC4|    GC6|
    (my$t_aid,my$gt,my$Poslist,my$sid,$_)=(split /\s+/,$_,5);
    next if $_ eq ""; 
    
    $_=~s/,*\|/,/g;
    $_=~s/\s+//g;

    my@gtList=(split /\,/,$_);
    my@POSlist=(split /\,/,$Poslist);

    unless( @gtList == @POSlist){
        print  STDERR ("#No enough SSR and pos found for $_|\t$line\n");
        next;
    }


    %tmp=();
    for(my$ii=0;$ii<=$#POSlist;$ii++){
        my($NO1,$NO2)=(split /\_/,$POSlist[$ii])[0,1]; 
        $tmp{$t_aid}{$NO1}{$NO2}=$gtList[$ii];
    }

    my($ot_gt_ref,$ot_gt)=("","");
    $Poslist="";
    foreach my$keyy1( sort {$a<=>$b} keys %{$pos2SSR{$t_aid}}){
        foreach my$keyy2(sort {$a<=>$b} keys %{$pos2SSR{$t_aid}{$keyy1}}){
            $ot_gt_ref.="$pos2SSR{$t_aid}{$keyy1}{$keyy2},";
            $ot_gt    .="$tmp{$t_aid}{$keyy1}{$keyy2}," if exists $tmp{$t_aid}{$keyy1}{$keyy2}; 
            $ot_gt    .="NA0," unless exists $tmp{$t_aid}{$keyy1}{$keyy2};
            $Poslist  .="$keyy1\_$keyy2,";
        }
    }
    print OT "$t_aid\t$ot_gt_ref\t$Poslist\t$sid\t$ot_gt\n";
}
close IN;
close OT;

print STDERR "All finished!!!! at ".`date`;




sub process{
    my($Sub,@units)=@_;
    $Sub=uc$Sub;
    my%ssrInfor=();

########################### record all pos and lenth for each ssr unit it occured in subseq;
    my$MaxRPN=0;
    foreach my$SsrUnit(@units){
        $SsrUnit=uc$SsrUnit;
        my$sub2="";
        my$index=index($Sub,$SsrUnit);
        next if $index ==-1;
        next if exists $ssrInfor{$SsrUnit};
    
        while(1){### this is to prevent this case: AAG-AAGAAAGAAAGAAAGAAAGA, aaga was repeat unit;
            $sub2=substr($Sub,$index);
            $sub2=~/(($SsrUnit)+)/ig;
            my$rn=((length$1)/(length$SsrUnit));
            my$end=$index+(length$1)-1;
            $ssrInfor{$SsrUnit}{$index}="$rn,$end";
            $MaxRPN=$rn if $MaxRPN<$rn;

            $index +=((length$1)-(length$SsrUnit)+1);
            $index=index($Sub,$SsrUnit,$index);
            last if $index ==-1;
        }
    }

    my$to="";
    my$PosList="";
    if($MaxRPN==0){
        foreach my$keyPos(sort {$a<=>$b} keys %ssr_No){
            $PosList .="$ssr_No{$keyPos},";
        }
        $to=join "0,",@units;
        return "$PosList\t$to"."0";
    }



#my$nn=0;foreach my$k(@units){$nn++; foreach (sort {$a<=>$b} keys %{$ssrInfor{$k}}){print "$nn:$k\t$_\t$ssrInfor{$k}{$_}\n";}}exit;
    
#    return "$units[0]$MaxRPN," if(@units==1);  #no complex SSR;
    if(@units==1){
        foreach my$INDEX(sort {$a<=>$b} keys %{$ssrInfor{$units[0]}}){
            my$RPN=(split /\,/,$ssrInfor{$units[0]}{$INDEX})[0];
            $to.="$units[0]$RPN\_" if $RPN>=4;
        }
        $PosList=$ssr_No{"0"};
        $to="$units[0]$MaxRPN," if $to eq "";
        $to=~s/\_+$//;
        return "$PosList,\t$to"; 
    }

##################################################################
########        get all combination of possible pos;        ######
##################################################################
###  get the first unit of complex SSR;
    my%keptPos=();
    unless (exists $ssrInfor{$units[0]}){
        $keptPos{0}{"0_0"}="";
    }else{
        foreach my$P (keys %{$ssrInfor{$units[0]}}){
            my($rn,$end)=(split /\,/,$ssrInfor{$units[0]}{$P})[0,1];
            $keptPos{0}{"$rn\_$end"}="";;
        }
    }
    $PosList=$ssr_No{"0"};

### get all combination of possible pos;
    for (my$i=1;$i<=$#units;$i++){
        my$j=$i-1;
        my$rpu=$units[$i];
         $PosList.=",$ssr_No{$i}";

        unless (exists $ssrInfor{$rpu}){ #this units does NOT include in this seq;
            foreach (keys %{$keptPos{$j}}){
                my($rpnList,$rpEnd)=(split /\_/,$_)[0,1];
                $keptPos{$i}{"$rpnList,0\_$rpEnd"}="";
            }
            next;
        }

        foreach my$P (keys  %{$ssrInfor{$rpu}}){
            my($rpn,$end)=(split /\,/,$ssrInfor{$rpu}{$P})[0,1];
            foreach my$old( keys %{$keptPos{$j}}){
                my($rpnList,$rpEnd)=(split /\_/,$old)[0,1];
                $keptPos{$i}{"$rpnList,$rpn\_$end"}="" if $P>=$rpEnd; ## one bp was allowed for two units to be overlaped;
            }
        }
    
        unless (exists $keptPos{$i}){#No possible pos for this unit;
            foreach ( keys %{$keptPos{$j}}){
                my($rpnList,$rpEnd)=(split /\_/,$_)[0,1];
                $keptPos{$i}{"$rpnList,0\_$rpEnd"}="";
            }
        }
    }
    
#foreach (keys %{$keptPos{5}}){print "$_\n";} exit;
###stat  all the combination;
    my@tmp=();
    my%comb=();
    my$last=$#units;
    my($total,$posNum,$ssrType)=();
 
    foreach (keys %{$keptPos{$last}}){
        @tmp=(split  /\,/,(split /\_/,$_)[0]);
        ($total,$posNum ,$ssrType)=(0,0,"");
        for(my$i=0;$i<=$#tmp;$i++){
            $total +=$tmp[$i];
            $posNum+=1 if $tmp[$i] ==0;
            $ssrType .= "$units[$i]$tmp[$i],";
        }
        push @{$comb{$posNum}{$total}},$ssrType;
    }
    
#########
#foreach my$k(keys %comb){print "#####$k\n";foreach my$k2(keys %{$comb{$k}}){ my@l=@{$comb{$k}{$k2}}; foreach (@l){print "$k,$k2,$_\n";}}};#exit;

    my$best1= (sort {$a<=>$b} keys %comb)[0];
    my$best2= (sort {$b<=>$a} keys %{$comb{$best1}})[0];
  
#    @{$comb{$best1}{$best2}}[0]=~s/\,$//g;
    @{$comb{$best1}{$best2}}[0]=~s/\,+$//;
    @{$comb{$best1}{$best2}}[0]=~s/\,/,\t/g; 
    return "$PosList,\t".@{$comb{$best1}{$best2}}[0];
}
