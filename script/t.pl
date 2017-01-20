#!/usr/bin/perl -w
use strict;

my%com=();
my%tcom=();
my@sam=();
foreach my$f(@ARGV){
    my$sam=(split /\./,$f)[0]; 
    push @sam,$sam;
        
    my%tmp=();
    open IN, "$sam.stat " or die "cannot open $sam.stat $!";
    while(<IN>){
        chomp;
        my$seqid=(split /\s+/,$_)[2];
        $tmp{$seqid}=""; 
    }
    close IN;

    open IN ," $f " or die  "cannot open $f $!";
    $/=">";<IN>;

    my(%h,%t)=(); 
    my(%th,%tt)=();
    while(<IN>){chomp;  
        my($id,$seq)=(split /\n+/,$_); 
        my$len=length$seq; 
        $_=~/\_x(\d+)/;
        my$abd=$1;
        $id=~/\_(AM\w+)/;
        $h{$1}+=$len*$abd;
        $t{$1} +=$abd; 

        next unless exists $tmp{$id};
        $th{$1}+=$len*$abd;
        $tt{$1} +=$abd;
    }
    $/="\n";
    close IN;
    foreach my$k (keys %h){ 
        $com{$k}{$sam}=($h{$k}/$t{$k});
        $tcom{$k}{$sam}=($th{$k}/$tt{$k}) if exists $tt{$k};;
    }
}

my$title=join "\t",@sam; 
print "\t$title\n"; 
foreach my$k(keys %com){
    print"$k"; 
    foreach (@sam){
        print "\t$com{$k}{$_}" if exists $com{$k}{$_};
        print "\tNA" unless exists $com{$k}{$_};
    }
    print "\n"; 
} 


print STDERR "\t$title\n";
foreach my$k(keys %tcom){
    print STDERR "$k";
    foreach (@sam){
        print STDERR "\t$tcom{$k}{$_}" if exists $tcom{$k}{$_};
        print STDERR "\tNA" unless exists $tcom{$k}{$_};
    }
    print STDERR "\n";
}

