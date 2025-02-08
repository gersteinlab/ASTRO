use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_i $opt_o $opt_p);
getopts('i:o:p:');

my $in_bam = $opt_i;
my $out_bam = $opt_o;

my $regex = '.+\|:_:\|(.+)$';
my $command = "samtools markdup -u --barcode-rgx '$regex' -@ $opt_p $opt_i ${opt_o}.temp.sam";
system "$command";

open(my $bam, "<", "${opt_o}.temp.sam") or die "Cannot process input BAM file: ${opt_o}.temp.sam";
open(my $modified_bam, "|-", "samtools view -b -o $out_bam -") or die "Cannot write to output BAM file: $out_bam";
#open(my $modified_bam, ">", "$out_bam") or die "Cannot write to output BAM file: $out_bam";

my %dedup;
while (my $line = <$bam>) {
    unless ($line =~ /^\@/) {
        my @fields = split(/\t/, $line);
        if ($fields[1] & 1024) {
            my ($readname, $pos) = (split(/\|:_:\|/, $fields[0]))[0,1];
            $dedup{$readname} = 1;
        }
    }
}

close $bam;


sub to_base_62 {
    my $number = shift;
    my @digits = (0..9, 'A'..'Z', 'a'..'z');
    my $base_62 = "";
    while ($number) {
        my $remainder = $number % 62;
        $base_62 = $digits[$remainder] . $base_62;
        $number = int($number / 62);
    }
    return $base_62;
}

my $ii = 0;
close $bam;

my %rename;
open(my $bam2, "<", "${opt_o}.temp.sam") or die "Cannot process input BAM file: $!";
while (my $line = <$bam2>) {
    if ($line =~ /^\@/) {
        print $modified_bam "$line";
    } else {
        my @fields = split(/\t/, $line);
        my ($readname, $pos) = (split(/\|:_:\|/, $fields[0]))[0,1];
        if (exists($dedup{$readname})) {
            next;
        }
        my $thename;
        if (exists($rename{$readname})) {
            $thename = $rename{$readname};
        }else{
            $thename = $ii."_$pos";
            $rename{$readname} = $thename;
        }
        if ($fields[1] & 4) {
            next;
        }
        $ii += 1;
        $fields[0] = $thename;
        print $modified_bam join("\t", @fields);
    }
}


close $bam;
close $modified_bam;
unlink "${opt_o}.temp.sam";