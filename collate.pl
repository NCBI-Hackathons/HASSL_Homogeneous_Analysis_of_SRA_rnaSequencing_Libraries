#!/usr/bin/env perl
use strict;
use warnings;

my $acc = shift;    # list of accessions in file
my @acc;            # array to save the order of Accessions
my $sra;            # hashref to hold info about each SRA
my $crsmheader;     # i.e. PF_BASES \t ... \t ... (taken from Picard output)
my $gotgenes = 0;
my @countsgenes;
my $genecounts;         # hashref to hold raw counts for each gene
open (my $ACC, "<", $acc);
while (<$ACC>) {
  chomp;
  my $acc = $_;
  push (@acc, $acc);
  
  # parse Picard CRSM output
  my $crsm = "$acc.GRCh38.p4.hisat.crsm";
  if (-e $crsm) {
    open (my $CRSM, "<", $crsm);
    my $hfound = 0;
    while (<$CRSM>) {
      chomp;
      if ($hfound) {
        $sra->{$acc}{crsm} = $_;
        last;
      }
      elsif (/^PF_BASES/) {
        $crsmheader = $_;
        $hfound = 1;   # qc header found
      }
    }
    close($CRSM);
  }

  # parse Hisat log for percent mapped  
  my $log =  "$acc.hisat.two.log";
  if (-e $log) {
    open (my $PM, "<", $log);
    while (<$PM>) {
      chomp;
      if (/^(.+)% overall alignment rate/) {
        $sra->{$acc}{pm} = $1/100;
      }
    }
    close($PM);
  }    

  # parse HTSeq counts
  my $counts = "$acc.GRCh38.p4.HTSeq.counts";
  if (-e $counts) {
    open (my $C, "<", $counts);
    while (<$C>) {
      chomp;
      next if /^__/;
      my ($gene, $cnt) = split /\t/,$_;
      if ($gotgenes == 0) {    # only need to get Gene names on the first SRA
        push (@countsgenes, $gene);
      }
      push @{$genecounts->{$gene}}, $cnt;  # this array of counts will be in the same order as @acc
    }
    if (@countsgenes) { $gotgenes++; }
  }
}
close($ACC);

unless ($sra) {
  print STDERR "did not find any *.crsm nor *.two.log files\n";
  exit;
}

# Output
open (my $QC, ">", "collate.qc.tsv");
print $QC join("\t", "Name", $crsmheader, "PCT_MAPPED"),"\n";
foreach (sort keys %$sra) {
  print $QC join("\t", $_, $sra->{$_}{crsm}, $sra->{$_}{pm}),"\n";
}
close ($QC);

open (my $COUNTS, ">", "collate.counts.tsv");
print $COUNTS join("\t", "Name", join ("\t", @acc)), "\n";
foreach my $g (@countsgenes) {
  print $COUNTS join("\t", $g, join("\t", @{$genecounts->{$g}} )  ), "\n";   # the array of counts will be in same order as @acc
}
close ($COUNTS);
