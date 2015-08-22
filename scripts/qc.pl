#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $o = {};
my @files = qw(maplogfile=s metricsfile=s);
my @qc = qw(pctmapped=s pctribo=s pctcoding=s pctutr=s pctintron=s
              pctintergenic=s pctmrna=s pctusable=s pctcorrect=s pctbias53=s);
my @global = qw(config=s help donotdelete debug sra=s);
GetOptions($o, @qc, @global, @files);

#print Dumper ($o),"\n";

die join(' ', @qc, @global, @files),$/ if $o->{help};
die "Need to define maplogfile (--maplogfile)\n" if !$o->{maplogfile};
die "Need to define metricsfile (--metricsfile)\n" if !$o->{metricsfile};
die "Need to define sra (--sra)\n" if !$o->{sra};

my $defaults = {pctmapped => {h => '', d => 0.65},
                  pctcoding => {h => 'PCT_CODING_BASES', d => 0},
                  pctutr => {h=>'PCT_UTR_BASES', d=>0},
                  pctintron => {h=>'PCT_INTRONIC_BASES',d=>0.15},
                  pctintergenic => {h=>'PCT_INTERGENIC_BASES',d=>1},
                  pctmrna => {h=>'PCT_MRNA_BASES', d=>0.65},
                  pctusable => {h=>'PCT_USABLE_BASES',d=>0.60},
                  pctcorrect => {h=>'PCT_CORRECT_STRAND_READS',d=>0},
                  pctbias53 => {h=>'MEDIAN_5PRIME_TO_3PRIME_BIAS',d=>0},
                  };
#print Dumper ($defaults), "\n";

foreach (@qc) {
  my ($qcname) = $_ =~ /^(.+?)=/;
  
  if (exists $o->{$qcname}) {
    $defaults->{$qcname}{d} = $o->{$qcname};
  }
}




#print Dumper ($defaults), "\n";
my %actual;
open (my $LOG, "<", $o->{maplogfile});
while (<$LOG>) {
  chomp;
  if (/^(.+)% overall alignment rate/) {
    $actual{pctmapped} = $1/100;
  }
}
close($LOG);

open (my $MET, "<", $o->{metricsfile});
while (<$MET>) {
  chomp;
  if (/^PF/) {
    my @h = split(/\t/, $_);
    my $metrics = <$MET>;
    chomp $metrics;
    my @f = split (/\t/, $metrics);
    for (my $i = 0; $i< @h; $i++) {
      foreach (keys %$defaults) {
        if ($defaults->{$_}{h} eq $h[$i]) {
          $actual{$_} = $f[$i];
          if ($actual{$_} eq '?'){
            $actual{$_}=1.0;
          }
        }
      }   
    }
    last
  }  
}
close ($MET);

#print Dumper (\%actual), "\n";

my $ok = 1;
foreach my $qcname (keys %actual) {
  my $actual = $actual{$qcname};
  my $default = $defaults->{$qcname}{d};
  
  if ($qcname eq 'pctintron' || $qcname eq 'pctintergenic') {
    if ($actual > $default) {
      $ok = 0;
      last;
    }
  }
  else {
    if ($actual < $default) {
      $ok = 0;
      last;
    }
  }
}
my $sra = $o->{sra};
unlink ("$sra.pass");
unlink ("$sra.fail");
($ok) ? `touch log/$sra.pass` : `touch log/$sra.fail`;


