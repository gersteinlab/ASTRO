#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my ($inputbed, $outputfile, $gtffile);
GetOptions(
    'inputbed|i=s'   => \$inputbed,
    'outputfile|o=s' => \$outputfile,
    'gtffile|g=s'    => \$gtffile
) or die "Usage: perl interBED2GeneFile.pl -i <inputbed> -o <outputfile> -g <gtffile> \n";
my $interbed = $inputbed . '.tempinter.bed';


my $cmd_intersect = "bedtools intersect -a \"$inputbed\" -b \"$gtffile\" -wao | "."awk 'BEGIN{OFS=\"\\t\"} { if (\$NF != 0) { print \$1, \$2, \$3, \$4, \$6, \$(NF-6), \$(NF-5), \$(NF-1), \$NF }}' > \"$interbed\"";
system($cmd_intersect) == 0 or die "Failed to run bedtools intersect: $!\n";


open(my $fh, '<', $interbed) or die "Cannot open $interbed: $!\n";
open(my $out_fh, '>', $outputfile) or die "Cannot open temporary output file: $!\n";



my %outputhash;
while (my $line = <$fh>) {
    chomp $line;
    my ($chrname, $read_begin, $read_end, $readname, $strand, $gtf_start, $gtf_end, $gtf_name, $overlap_length) = split(/\t/, $line);
    my $read_length = $read_end-$read_begin;
    my $gtf_length = $gtf_end-$gtf_start+1;
    my $valuetempit = ($overlap_length - ($read_length - $overlap_length)) / $gtf_length;
    my $infovector = $readname;
    my $geneout;
    my $geneout0;
    $read_begin = $read_begin + 1;
    if ($infovector =~ /(\d+)_(\d+)_(\d+):(.+)$/) {
        $infovector = $1;
        my $Coord1 = $2;
        my $Coord2 = $3;
        my @parts = split(/:/, $4);
        shift(@parts);
        my ($ASscore, $fastq_read_length) = (@parts, '.', '.');
        $geneout = "$Coord1\t$Coord2\t$gtf_name\t$ASscore\t$fastq_read_length\t$overlap_length\t$read_begin\t$read_end\t$gtf_start\t$gtf_end";
        $geneout0 = "$Coord1\t$Coord2\t$gtf_name\t$ASscore\t$fastq_read_length\t$overlap_length\t$read_begin\t$read_end\t$gtf_start\t$gtf_end\t$infovector";
    } else {
        die "wrong read name format here\n";
    }
    if (exists($outputhash{$infovector})){
        if ($valuetempit > $outputhash{$infovector}->[0]) {
            $outputhash{$infovector}->[0] = $valuetempit;
            $outputhash{$infovector}->[1] = $geneout0;
        }elsif($valuetempit == $outputhash{$infovector}->[0]){
            $outputhash{$infovector}->[1] = $outputhash{$infovector}->[1] . "\t$geneout";
        }
    }else{
        $outputhash{$infovector}->[0] = $valuetempit;
        $outputhash{$infovector}->[1] = $geneout0;
    }
}

for my $var (keys(%outputhash)) {
    my $gene = $outputhash{$var}->[1];
    print $out_fh "$gene\n";
}

close($fh);
close($out_fh);

#unlink $interbed or warn "Could not delete $interbed: $!\n";


