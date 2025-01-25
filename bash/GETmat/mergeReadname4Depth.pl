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
    my $readname = splice(@fields, 10, 1);
    my ($Coord1, $Coord2, $gtf_name, $ASscore, $fastq_read_length, $overlap_length, $read_begin, $read_end, $gtf_start, $gtf_end) = split(/\t/, $line);
    my $read_length = $read_end-$read_begin;
    my $gtf_length = $gtf_end-$gtf_start+1;
    my $valuetempit = ($overlap_length - ($read_length - $overlap_length)) / $gtf_length;
    my $infovector = $readname;
    #my $geneout = "$Coord1\t$Coord2\t$gtf_name\t$ASscore\t$fastq_read_length\t$overlap_length\t$read_begin\t$read_end\t$gtf_start\t$gtf_end\t$infovector";
    my $geneout = "$Coord1\t$Coord2\t$gtf_name\t$ASscore\t$fastq_read_length\t$overlap_length\t$read_begin\t$read_end\t$gtf_start\t$gtf_end";
    if (exists($outputhash{$infovector})){
        if ($valuetempit > $outputhash{$infovector}->[0]) {
            $outputhash{$infovector}->[0] = $valuetempit;
            $outputhash{$infovector}->[1] = $geneout;
        }elsif($valuetempit == $outputhash{$infovector}->[0]){
            $outputhash{$infovector}->[1] = $outputhash{$infovector}->[1] . "\t$geneout";
        }
    }else{
        $outputhash{$infovector}->[0] = $valuetempit;
        $outputhash{$infovector}->[1] = $geneout;
    }
}

for my $var (keys(%outputhash)) {
    my $gene = $outputhash{$var}->[1];
    print $out_fh "$gene\n";
}

close($fh);
close($out_fh);

#unlink $interbed or warn "Could not delete $interbed: $!\n";


