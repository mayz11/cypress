#!/usr/bin/perl

#This script performs PBS test to identify genes under positive selection within C. gigantea or northern C. duclouxiana. Cgi means C. gigantea. Cdun means northern C. duclouxiana. Cdus means sorthern C. duclouxiana.

use strict;
use warnings;

my($fileInFst1,$fileInFst2,$fileInFst3,$fileOut)=@ARGV;
die("usage: FstForCgi/Cdun-Cdus FstForCgi/Cdun-outgroup FstForCdus-outgroup out \n")unless($fileOut);

my %win1=ReadFst($fileInFst1);
my %win2=ReadFst($fileInFst2);
my %win3=ReadFst($fileInFst3);

my @pbs;
my %winTopbs;
foreach my $chr(keys %win1){
    my $t1=$win1{$chr};
    my $t2=$win2{$chr};
    my $t3=$win3{$chr};
    my $pbs;
    next if(!defined $t1 or !defined $t2 or !defined $t3);
    if($t1=~/Inf/ or $t2=~/Inf/ or $t3=~/Inf/){
      $pbs="($t1 + $t2 - $t3)/2";
    }
    else{
      $pbs=($t1+$t2-$t3)/2;
    }
    push(@pbs,$pbs) unless($pbs=~m/Inf/);
    $winTopbs{$chr}=$pbs unless($pbs=~m/Inf/);
}

@pbs=sort{$b<=>$a}@pbs;
my $total=@pbs;
my $ValueTop=$pbs[int($total*0.1)];

my $fileOutSig="$fileOut.sig";

open(Fo2,'>',$fileOut) or die("$!: $fileOut\n");
open(Fo1,'>',$fileOutSig) or die("$!: $fileOutSig\n");
foreach my $chr(keys %winTopbs){
    if($winTopbs{$chr} >= $ValueTop){
      print Fo1 "$chr\t$winTopbs{$chr}\n";
    }
    print Fo2 "$chr\t$winTopbs{$chr}\n";
}
close(Fo1);
close(Fo2);

sub ReadFst{
  my($fileInFst)=@_;
  my %win;
  open(F,'<',$fileInFst) or die("$!: $fileInFst\n");
  <F>;
  while(<F>){
    chomp;
    next if(/^\s*$/);
    next if(/^\s*#/);
    my($Chrom,$MeanFst)=split /\t/;
    my $fst=$MeanFst;
    my $t;
    if($fst<1){
      $t=-1*log(1-$fst);
    }
    else{
      $t=-1*log(1-0.9999);
    }
    $win{$Chrom}=$t;
  }
  close(F);
  return %win;
}
