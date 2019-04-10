#!/usr/bin/perl

#This script performs HKA test to identify genes under positive selection within C. gigantea or northern C. duclouxiana. Cgi means C. gigantea. Cdun means northern C. duclouxiana. Cdus means sorthern C. duclouxiana. InFile2 contains sample ID information for each group.

use strict;
use warnings;

my($inFileVcfGz,$fileIn2,$inR,$fileOut,$tarName)=@ARGV;
die("usage: invcf sampleIdForEachGroup Path_To_R_Binary out cgi/cdun\n")unless($tarName);

my %h1;
open(F,'<',$fileIn2) or die("$!: $fileIn2\n");
while(<F>){
  chomp;
  next if(/^\s*$/);
  next if(/^\s*#/);
  my($group,$id)=split /\s/;
  $h1{$group}{$id}=1;
}
close(F);

my @outgroup=keys %{$h1{outgroup}};
my $outgroup=shift (@outgroup);

my $tarN=keys %{$h1{$tarName}};
my $tarHapN=2*$tarN;

my $conN=keys %{$h1{cdus}};
my $conHapN=2*$conN;

my %chromToNum;
my $allA=0;
my $allC=0;
my @name;

open(F,"gzip -dc $inFileVcfGz |") or die("$!: $inFileVcfGz\n");
while(<F>){
  chomp;
  next if(/^\s*$/);
  if (m|^#CHROM|) {
    my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@ind)=split /\t/;
    @name=@ind;
  }
  next if(/^\s*#/);
  my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@ind)=split /\t/;

  $chromToNum{$CHROM}{A}=0 if (!defined $chromToNum{$CHROM}{A});
  $chromToNum{$CHROM}{C}=0 if (!defined $chromToNum{$CHROM}{C});

  my %site;
  for (my $i=0;$i<@name;$i++){
    my @site= split (/;/,$ind[$i]);
    my $site=shift @site;
    $site{$name[$i]}=1 if ($site=~/0\/0/);
    $site{$name[$i]}=2 if ($site=~/0\/1/);
    $site{$name[$i]}=3 if ($site=~/1\/1/);
  }

  my ($tar0num,$tar1num,$tarR,$con0num,$con1num,$conR);
  $tar0num=0;
  $tar1num=0;
  $con0num=0;
  $con1num=0;
  foreach my $id(keys %{$h1{$tarName}}){
    if ($site{$id}==1){$tar0num+=2;}
    if ($site{$id}==2){$tar0num+=1;}
    if ($site{$id}==2){$tar1num+=1;}
    if ($site{$id}==3){$tar1num+=2;}
  }
  $tarR=$tar0num/$tarHapN;

  foreach my $id(keys %{$h1{cdus}}){
    if ($site{$id}==1){$con0num+=2;}
    if ($site{$id}==2){$con0num+=1;}
    if ($site{$id}==2){$con1num+=1;}
    if ($site{$id}==3){$con1num+=2;}
  }
  $conR=$con0num/$conHapN;

  if ($tarR<0.9 and $tarR>0.1){ $chromToNum{$CHROM}{A}++;$allA++;}
  if ($site{$outgroup}==1){
    if ($tarR<=0.1 and $conR>=0.5){$chromToNum{$CHROM}{C}++;$allC++;}
  }
  if ($site{$outgroup}==3){
    if ($tarR>=0.9 and $conR<=0.5){$chromToNum{$CHROM}{C}++;$allC++;}
  }
  if ($site{$outgroup}==2){
    if ($tarR>=0.9 and $conR<0.9 and $conR>0.1){$chromToNum{$CHROM}{C}++;$allC++;}
    if ($tarR<=0.1 and $conR<0.9 and $conR>0.1){$chromToNum{$CHROM}{C}++;$allC++;}
  }
}
close(F);

open(Fo,'>',$fileOut) or die("$!: $fileOut\n");
print Fo "ChromID\tAvalue\tCvalue\tallAvalue\tallCvalue\tchiCorSymbol\tchiCorPvalue\n";
foreach my $chrom(keys %chromToNum){
  my ($symChiCor,$pChiCor);
  my $chiCor=`$inR -e "data=matrix(c($chromToNum{$chrom}{A},$allA,$chromToNum{$chrom}{C},$allC),ncol=2);chisq.test(data,correct=T)"`;
  $chiCor=~m/p-value\s+(\S)\s+(\S+)/;
  $symChiCor=$1;
  $pChiCor=$2;
  print Fo "$chrom\t$chromToNum{$chrom}{A}\t$chromToNum{$chrom}{C}\t$allA\t$allC\t$symChiCor\t$pChiCor\n";
}
close(Fo);
