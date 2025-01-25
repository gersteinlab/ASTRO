use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_d $opt_f $opt_g $opt_e);
getopts('i:o:d:f:g:e');


my $ASthreshold; 
my $ratiothreshold_gtflen;
my $doit = 1;
$opt_f = $opt_f // '25:0.75';
if ($opt_f =~ m/([^:]+):([^:]+)/) {
	$ASthreshold = $1;
	$ratiothreshold_gtflen = $2;
}else{
	$doit = 0;
}




my %num2genehash;
my %gene2numhash;

if ($opt_g) {
	my $ii = 1;
	open DD, "$opt_g" or die;
	while (<DD>) {
		chomp;
		my $gene = (split(/\t/, $_))[8];
		if (!exists($gene2numhash{$gene})) {
			$gene2numhash{$gene} = $ii;
			$num2genehash{$ii} = $gene;
			$ii += 1;
		}
	}
	close DD;
}else{
	my $ii = 1;
	open DD, "$opt_i" or die;
	my %genes;
	while (<DD>) {
		chomp;
		my @fieldsin = split(/\t/, $_);
		while (@fieldsin) {
			my ($x, $y, $genename, $ASscore, $readlen, $overlaplen, $genelen, $read_length) = splice(@fieldsin, 0, 7);
			$genes{$genename} = 1;
		}
	}
	close DD;
	my @genes = sort {$a cmp $b} keys(%genes);
	for my $gene (@genes) {
		$gene2numhash{$gene} = $ii;
		$num2genehash{$ii} = $gene;
		$ii += 1;
	}
}




my @locations;
open DD, "$opt_d" or die;
while (<DD>) {
	chomp;
	my ($x, $y) = (split(/\t/, $_))[1,2];
    push(@locations, "${x}_$y");
}
close DD;


open IN, "$opt_i" or die;
open OUT, '>', "$opt_o" or die;


my %matrix;
my %genenames;
my $previous_connect;

while (<IN>) {
	chomp;
	my @fieldsin = split(/\t/, $_);
	my @genenames;
	my ($x, $y);
	while (@fieldsin) {
		($x, $y, my $genename, my $ASscore, my $readlen, my $overlaplen, my $genelen, my $read_length) = splice(@fieldsin, 0, 8);
		if ($doit && $ASscore ne '.' && $ASscore <= $ASthreshold && $ASscore <= $ratiothreshold_gtflen*$genelen) {
			next;
		}
		push @genenames, $gene2numhash{$genename};
	}
	
	my %seen;
	@genenames = grep { !$seen{$_}++ } @genenames;
	
	if (scalar(@genenames) == 1)  {
		$matrix{"${x}_${y}_$genenames[0]"} += 1;
		$genenames{$genenames[0]} += 1;
	}elsif (scalar(@genenames) > 1) {
		@genenames = sort {$a <=> $b} @genenames;
		my $genename;
		if ($opt_e) {
			$genename = $genenames[0];
		}else{
			$genename = join("-", @genenames);
		}
		$matrix{"${x}_${y}_$genename"} += 1;
		$genenames{$genename} += 1;
	}
}
close IN;



my @genenames1;
my @genenames2;

my @genenames = keys(%genenames);
for my $element (@genenames) {
    if ($element =~ /-/) {
        push @genenames2, $element;
    } else {
        push @genenames1, $element;
    }
}


@genenames1 = sort {$a <=> $b} @genenames1;
@genenames2 = sort {$a cmp $b} @genenames2;


for my $var (@genenames1) {
	print OUT "\t$num2genehash{$var}";
}
for my $var (@genenames2) {
	my @vars = split(/-/, $var);
	@vars = map { 
    	exists $num2genehash{$_} 
    	? $num2genehash{$_} 
    	: die "Error: No conversion found for $_\n" 
	} @vars;
	my $varout = join("-|:", @vars);
	print OUT "\t$varout";
}
print OUT "\n";
 

for my $var (@locations) {
	my $varout = $var;
	$varout =~ s/_/x/;
	print OUT "$varout";
	for my $var1 (@genenames1) {
		my $count;
		if (exists($matrix{"${var}_$var1"})) {
			$count = $matrix{"${var}_$var1"};
		}else{
			$count = 0;
		}
		print OUT "\t$count";
	}
	for my $var1 (@genenames2) {
		my $count;
		if (exists($matrix{"${var}_$var1"})) {
			$count = $matrix{"${var}_$var1"};
		}else{
			$count = 0;
		}
		print OUT "\t$count";
	}
	print OUT "\n";
}

close OUT;