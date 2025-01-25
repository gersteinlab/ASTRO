use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_i $opt_o $opt_p);
getopts('i:o:p:');

my $in_bam = $opt_i;
my $out_bam = $opt_o;

my $regex = '.+\|:_:\|(.+)$';
my $command = "samtools markdup -u -r --barcode-rgx '$regex' -@ $opt_p $opt_i ${opt_o}.temp.sam";
my $exit_status = system "$command";

if ($exit_status != 0) {
    my $exit_code = $exit_status >> 8;
    die "Error: Command failed with exit code $exit_code\n Command: $command\nVery possible due to memory restrictions\n";
}

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
        if ($fields[1] & 4) {
            next;
        }
        my $thename;
        if (exists($rename{$readname})) {
            $thename = $rename{$readname};
        }else{
            $thename = $ii."_$pos";
            $rename{$readname} = $thename;
            $ii += 1;
        }
        my $readlen =length($fields[9]);
        my $ASscore;
        if ($line =~ m/AS:i:([^\s]+)/) {
            $ASscore = $1;
        }else{
            print "$line\n";
            die 'wrong format of SAM about AS';
        }
        
        $fields[0] = "$thename:$ASscore:$readlen";
        print $modified_bam join("\t", @fields);
    }
}


close $bam;
close $modified_bam;
unlink "${opt_o}.temp.sam";