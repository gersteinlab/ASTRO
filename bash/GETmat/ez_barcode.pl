use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_u $opt_r $opt_b $opt_l);
getopts('i:r:o:u:b:l:');

my $position = 0;
my $length = 0;
if ($opt_b =~ m/^([0-9_]+)b$/) {
    my $catch = $1;
    if ($catch =~ m/\_/) {
        my ($pos, $len) = split(/_/, $catch);
        $position = $pos;
        $length = $len;
    }else{
        $length = $catch;
    }
}elsif ($opt_b =~ m/^b([0-9_]+)$/){ 
    my $catch = $1;
    if ($catch =~ m/\_/) {
        my ($pos, $len) = split(/_/, $catch);
        $position = -$pos;
        $length = $len;
    }else{
        $position = -$catch;
        $length = $catch;
    }
}else{
    die "wrong format $opt_b\n";
}

print "$position\n";
print "$length\n";


if (`file $opt_r` =~ /gzip compressed data/) {
    open IN, "gzip -dc $opt_r|" or die;
}else{
    open IN, "$opt_r" or die;
}

my %bar;
while (<IN>) {
	chomp;
	$_ =~ s/\r//g;
	my ($seq, $xx, $yy) = split(/\t/, $_);
	$bar{$seq} = "${xx}_$yy";
}
close IN;



if (`file $opt_i` =~ /gzip compressed data/) {
    open IN, "gzip -dc $opt_i|" or die;
}else{
    open IN, "$opt_i" or die;
}

if ($opt_i =~ m/(fq|fastq|faq|fq.gz|fastq.gz|faq.gz)$/) {
    $/='@';
}elsif ($opt_i =~ m/(fa|fasta|fas|fa.gz|fasta.gz|fas.gz)$/) {
    $/='>';
}
 
my $minlenofread = 0;
my $maxlenofread = 999999999;
if ($opt_l && $opt_l ne 'NA') {
   ($minlenofread, $maxlenofread) = (split(/_/, $opt_l))[0,1];
}


open OUT, '>', "$opt_o" or die;
open OUTT, '>', "$opt_u" or die;
my $line1 = <IN>;
while (<IN>) {
    $_ =~ s/\r//g;
    chomp($_);
	my ($name1, $seq1, $qual) = (split(/\n/, $_))[0,1];
	if (length($seq1) <= $minlenofread || length($seq1) >= $maxlenofread) {
	    next;
	}
	my $barcode = substr($seq1, $position, $length);
	if (exists($bar{$barcode})) {
	    $name1 =~ s/\s.+$//g;
		print OUT "$name1\t$bar{$barcode}\n";
	}else{
	    print OUTT "\@$_";
	}

}


close OUT;
close IN;
close OUTT;


