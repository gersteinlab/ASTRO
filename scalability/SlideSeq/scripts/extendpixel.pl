#!/usr/bin/env perl
use strict;
use warnings;

my %IUPAC = (
    'A' => ['A'], 'C' => ['C'], 'G' => ['G'], 'T' => ['T'],
    'R' => ['A','G'], 'Y' => ['C','T'], 'S' => ['G','C'], 'W' => ['A','T'],
    'K' => ['G','T'], 'M' => ['A','C'],
    'B' => ['C','G','T'], 'D' => ['A','G','T'], 'H' => ['A','C','T'], 'V' => ['A','C','G'],
    'N' => ['A','C','G','T'],
);

#set default should use our! like:
#our $opt_i
use Getopt::Std;
use vars qw($opt_i $opt_o);
getopts('i:o:');

open IN, "$opt_i" or die;
open OUT, '>', "$opt_o" or die;
open OUTT, '>', "$opt_o.out" or die;

my $ii = 1;
print OUTT "$ii\n";
while (<IN>) {
    chomp;
    $ii += 1;

    my @parts = split /\s+/, $_;
    my $seq   = uc($parts[0]);
    my @rest  = @parts[1..$#parts];

    my @pos;
    for (my $i = 0; $i < length($seq); $i++) {
        my $ch = substr($seq, $i, 1);
        push @pos, $i if ($ch !~ /[ACGT]/);
    }
    
    print OUTT "$ii\n" if (@pos <= 1);
    if (@pos == 0) {
        print OUT join("\t", $seq, @rest), "\n";
        next;
    }
    
    next if (@pos > 1);
    

    my $i  = $pos[0];
    my $ch = substr($seq, $i, 1);
    my @choices = @{$IUPAC{$ch}};

    foreach my $b (@choices) {
        my $newseq = $seq;
        substr($newseq, $i, 1, $b);
        print OUT join("\t", $newseq, @rest), "\n";
    }
}
