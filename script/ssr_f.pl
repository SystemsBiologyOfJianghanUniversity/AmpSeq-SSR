#!/usr/bin/perl -w
use strict;
    die("Usage: perl $0 <fa> <otPre> <cutoff4rpn> <rpnList>\n
            Used to call SSR based on regular expression.
            fa        : sequece file in fasta formate.
            otPre     :
            cutoff4rpn: the max length length of repeat unit.
            rpnList   : in this formate of Len1:num1,Len2:num2, Len is length of repeat Unit, and number is the least count of repeat in sequence.  
            conf      : configure file about least repeat number for different repeat.\n
            
            Attention: it is argue for this case: CTCCTCCTCCTCCTCCTCTCTCTCTCTCTCTCTCTCTC
            ") unless @ARGV ==4;
            
my($fa,$otPre,$cutoff4rpn,$rpnList)=@ARGV;
open IN, " $fa " or die  "cannot open $fa ";
open OT, "> $otPre.ssr" or die  "cannot write $otPre.ssr $!\n";
my%rpu=();
my(%rptRegion,%ssr,%unitList)=();
&getRepeatUnit();

print STDERR"Get all repeat unit\n";
$/=">"; <IN>;
while(<IN>){
    chomp;
    my($id,@seq)=(split /\n+/,$_);
    my$seq=join "",@seq;
    (%rptRegion,%ssr,%unitList)=();

print "match: ".`date`;
    foreach my$unit(keys %rpu){
        my$Num=$rpu{$unit};
        next unless $seq=~/(($unit){$Num,})/i;
        while($seq=~/(($unit){$Num,})/ig){
            my($left,$match,$right)=($`,$1,$&);
            
            my$start=(length$left)+1;
            my$rpLen=length$match;
            my$end=$start+$rpLen-1;
            my$rpuL=length$unit;
            my$rpuN=$rpLen/$rpuL;

            $ssr{$rpuN}{$start}="$end\t$rpuL\t$unit";
        }
    }
print "select".`date`;
    my%keptPos=();
    my$j=0;
    my$seqLen=length$seq;
    foreach my$rpuNum(sort {$b<=>$a} keys %ssr){
        foreach my$Start(sort {$b<=>$a} keys %{$ssr{$rpuNum}}){ #two ssr should be as far as possible from each other;
            my($End,$rpuLen,$Unit)=(split /\t+/,$ssr{$rpuNum}{$Start})[0,1,2];
            $Unit=lc$Unit;

            my$flag=0;
            foreach my$ps(keys %keptPos){
##### TO get all ssr unit, the overlap length for two nearby ssr was no more than the longer repeat: gcgcgcgcgctctctctctctctcacacacatatatatatatatat;
                my($pe,$p_rpuLen)=(split /\_/, $keptPos{$ps})[0,1];
                $p_rpuLen=$rpuLen if $p_rpuLen<$rpuLen;

                $flag=1 if ($Start>=($ps+$p_rpuLen) and $Start<=($pe-$p_rpuLen));
                $flag=1 if ($End  >=($ps+$p_rpuLen) and $End  <=($pe-$p_rpuLen));
                $flag=1 if ($Start>=($ps-$p_rpuLen) and $End <=($pe+$p_rpuLen));
            }
            next if $flag ==1;
            
            $j++;
            print OT "$id\t$j\t$rpuLen\t$Unit\t$rpuNum\t$Start\t$End\t$seqLen\n";
            $keptPos{$Start}="$End\_$rpuLen";
        }
    }
}


sub getRepeatUnit{
    my%RPN=();
    $rpnList=~s/[,:]/ /g;
    my@t=(split /\s+/,$rpnList);

    for(my$i=0;$i<=$#t;$i+=2){
        $RPN{$t[$i]}=$t[$i+1];
    }

    my$rpu="";
    my$count=0;
    my$rpuLen=0;
    foreach my$f(qw/A T C G/){
        foreach my$s(qw/A T C G/){
            $rpu="$f$s";
            $rpuLen=length$rpu;
            $rpu{$rpu}= $RPN{$rpuLen} if $rpuLen<=$cutoff4rpn;
            foreach my$t(qw/A T C G/){
                $rpu="$f$s$t";
                $rpuLen=length$rpu;
                $rpu{$rpu}= $RPN{$rpuLen} if $rpuLen<=$cutoff4rpn;
                foreach my$F(qw/A T C G/){
                    $rpu="$f$s$t$F";
                    $rpuLen=length$rpu;
                    $rpu{$rpu}= $RPN{$rpuLen} if $rpuLen<=$cutoff4rpn;
                    foreach my$v(qw/A T C G/){
                        $rpu="$f$s$t$F$v";
                        $rpuLen=length$rpu;
                        $rpu{$rpu}= $RPN{$rpuLen} if $rpuLen<=$cutoff4rpn;
                    }
                }
            }
        }
    }

    foreach my$key (keys %rpu){
        delete  $rpu{$key} if ($key=~/^A+$/ or $key=~/^T+$/ or $key=~/^C+$/ or $key=~/^G+$/);
    }
}

