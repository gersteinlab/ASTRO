#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile);
GetOptions(
    'inputfile|i=s'   => \$inputfile,
    'outputfile|o=s' => \$outputfile,
) or die "Usage: perl mergeReadname.pl -i <inputfile> -o <outputfile> \n";



open(my $fh, '<', $inputfile) or die "Cannot open $inputfile: $!\n";
open(my $out_fh, '>', $outputfile) or die "Cannot open temporary output file: $!\n";


my %outputhash;
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    my $valuetempit = ($fields[5] - ($fields[7] - $fields[5])) / $fields[6];
    my $infovector = splice(@fields, 8, 1);
    my $geneout = join("\t", @fields);

    if (exists($outputhash{$infovector})){
        if ($valuetempit > $outputhash{$infovector}->[0]) {
            $outputhash{$infovector}->[0] = $valuetempit;
            $outputhash{$infovector}->[1] = "$geneout";
        }elsif($valuetempit == $outputhash{$infovector}->[0]){
            $outputhash{$infovector}->[1] = $outputhash{$infovector}->[1] . "\t$geneout";
        }
    }else{
        $outputhash{$infovector}->[0] = $valuetempit;
        $outputhash{$infovector}->[1] = "$geneout";
    }
}

for my $var (keys(%outputhash)) {
    my $gene = $outputhash{$var}->[1];
    print $out_fh "$gene\n";
}

close($fh);
close($out_fh);

#unlink $interbed or warn "Could not delete $interbed: $!\n";


