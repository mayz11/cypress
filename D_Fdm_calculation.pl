#!/usr/bin/perl

#This script calculates D statistic and fdm statistic described in Martin et al. 2015, MBE and Malinsky, M. et al. 2015, Science. Pop1 is sorthern C. duclouxiana. Pop2 is northern C. duclouxiana. Pop3 is C. gigantea. Pop4 is outgroup.

use strict;
use warnings;

my($bcftools,$fileInVCFGZ,$inChr,$inStart,$inEnd,$inGroup1,$inGroup2,$inGroup3,$inGroup4)=@ARGV;
die("usage: PathToBcftools inVCFGZ SeqID StartPos EndPos SampleListForPOP1 SampleListForPOP2  SampleListForPOP3  SampleListForPOP4\n")unless($inGroup4);

my @vcf=`$bcftools filter -r '$inChr:$inStart-$inEnd' $fileInVCFGZ`;

my %nameToGroup;
my @name1=split /,/,$inGroup1;foreach my $name(@name1){$nameToGroup{$name}="P1";}
my @name2=split /,/,$inGroup2;foreach my $name(@name2){$nameToGroup{$name}="P2";}
my @name3=split /,/,$inGroup3;foreach my $name(@name3){$nameToGroup{$name}="P3";}
my @name4=split /,/,$inGroup4;foreach my $name(@name4){$nameToGroup{$name}="P4";}

my @name;
my %PosArray;
foreach (@vcf) {
  chomp;
  if (m|^#CHROM|) {
    my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@ind)=split /\t/;
    @name=@ind;
  }
  next if(/^\s*#/);
  my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@ind)=split /\t/;
  $PosArray{$POS}=[$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@ind];
}

my ($D,$FDM)=Calculate($inStart,$inEnd);
print "$inChr\t$D\t$FDM\n";

sub Calculate{
  my($PosStart,$PosEnd)=@_;
  my $numerator=0;
  my $denominator=0;
  my $FDnumerator=0;
  my $FDMdenominator=0;
  my $D="NA";
  my $FDM="NA";
  for(my $i=$PosStart;$i<=$PosEnd;$i++) {
    next if(!exists $PosArray{$i});
    my($ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@ind)=@{$PosArray{$i}};
    my %g;
    my %RefAltSupport;
    for (my $i=0;$i<@ind;$i++) {
      my $name=$name[$i];
      my $group=$nameToGroup{$name};
      if (!defined $group) {
        print "$name cannot assigned to group\n";
      }
      my($type,$other)=split/:/,$ind[$i];
      my @tp=split '\/',$type;
      foreach(@tp){
        $g{$group}{0}++ if($_ eq '0');
        $g{$group}{1}++ if($_ eq '1');
        $RefAltSupport{$group}{$_}{$name}=1;
      }
    }
    my $failGroupSupport=0;
    foreach my $group (keys %RefAltSupport) {
      if (exists $RefAltSupport{$group}{0} and keys %{$RefAltSupport{$group}{0}} < 2) {
        $failGroupSupport=1;
        last;
      }
      if (exists $RefAltSupport{$group}{1} and keys %{$RefAltSupport{$group}{1}} < 2) {
        $failGroupSupport=1;
        last;
      }
    }
    next if($failGroupSupport==1);

    foreach my $p ("P1","P2","P3","P4") {
      $g{$p}{0}=0 if(!defined $g{$p}{0});
      $g{$p}{1}=0 if(!defined $g{$p}{1});
    }
    next if ($g{P4}{0} == $g{P4}{1});

    my $DerivedAllele=$g{P4}{0}>$g{P4}{1}?'1':'0';

    my $Pi1=$g{P1}{$DerivedAllele}/($g{P1}{0}+$g{P1}{1});
    my $Pi2=$g{P2}{$DerivedAllele}/($g{P2}{0}+$g{P2}{1});
    my $Pi3=$g{P3}{$DerivedAllele}/($g{P3}{0}+$g{P3}{1});
    my $Pi4=$g{P4}{$DerivedAllele}/($g{P4}{0}+$g{P4}{1});
    my $PiD=$Pi2>$Pi3?$Pi2:$Pi3;
    my $PiDM=$Pi1>$Pi3?$Pi1:$Pi3;
    my $Cabbai=(1-$Pi1)*$Pi2*$Pi3*(1-$Pi4);
    my $Cbabai=$Pi1*(1-$Pi2)*$Pi3*(1-$Pi4);
    next if($Pi3==0);
    $numerator+=$Cabbai-$Cbabai;
    $denominator+=$Cabbai+$Cbabai;
    $FDnumerator+=$Cabbai-$Cbabai;
    if ($Pi2 >= $Pi1) {
      $FDMdenominator+=(1-$Pi1)*$PiD*$PiD*(1-$Pi4) - $Pi1*(1-$PiD)*$PiD*(1-$Pi4);
    } else {
      $FDMdenominator-=(1-$PiDM)*$Pi2*$PiDM*(1-$Pi4) - $PiDM*(1-$Pi2)*$PiDM*(1-$Pi4);
    }
  }
  if ($denominator!=0) {
    $D=$numerator/$denominator;
  }
  if ($FDMdenominator!=0) {
    $FDM=$FDnumerator/$FDMdenominator;
  }
  return ($D,$FDM);
}
