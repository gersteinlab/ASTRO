use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_r $opt_u);
getopts('i:r:o:u:');
#perl barcodedFq.pl -r index.out -u temp_3.fastq -o seesee.fq -i ../archie/test/raw_qc_linker1_R1.fastq.gz


my $indicate = `file $opt_r`;
if ($indicate =~ /gzip compressed data/) {
    $indicate = 1;
} else {
    $indicate = 0;
}
if ($indicate) {
    open REF, "gzip -dc $opt_r|" or die;
}else{
    open REF, "$opt_r" or die;
}

my %bar;
while (<REF>) {
        $_ =~ s/\r//g;
        $_ =~ s/\n+$//g;
        my ($name, $location) = split(/\t/, $_);
        if ($location ne '*') {
	        $name =~ s/\s.+$//g;
            $bar{$name} = $location;
	    }
        
}
close REF;



$indicate = `file $opt_u`;
if ($indicate =~ /gzip compressed data/) {
    $indicate = 1;
} else {
    $indicate = 0;
}

if ($indicate) {
    open UMI, "gzip -dc $opt_u|" or die;
}else{
    open UMI, "$opt_u" or die;
}
print "11\n";
my %ids;
$/='@';
while (<UMI>) {
        $_ =~ s/\n+$//g;
        $_ =~ s/\r//g;
        my ($name, $seq) = (split(/\n/, $_))[0,1];
        $name =~ s/\s.+$//g;
        if (exists($bar{$name}) && length($seq) >= 6) {
            $ids{$name} = "$bar{$name}:$seq";
        }
}
close UMI;
undef %bar;


$indicate = `file $opt_i`;
if ($indicate =~ /gzip compressed data/) {
    $indicate = 1;
} else {
    $indicate = 0;
}

if ($indicate) {
    open IN, "gzip -dc $opt_i|" or die;
}else{
    open IN, "$opt_i" or die;
}
 

open OUT, '>', "$opt_o" or die;
my $line1 = <IN>;
while (<IN>) {
    $_ =~ s/\r//g;
	my ($name1, $seq1, $qual) = (split(/\n/, $_))[0,1,3];
	if (length($seq1) == 0) {
	    next;
	}
	$name1 =~ s/\s.+$//g;
	if (exists($ids{$name1})) {
		my $location = $ids{$name1};
		print OUT  "\@$name1|:_:|$location\n$seq1\n+\n$qual\n";
	}

}


close OUT;
close IN;


