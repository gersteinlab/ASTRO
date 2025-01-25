use strict;
use warnings;
use Getopt::Std;
use List::Util qw(max min);


use vars qw($opt_i $opt_o $opt_d $opt_f $opt_g $opt_e);
getopts('i:o:d:f:g:e');

sub merge_hashes {
    my ($hash1, $hash2) = @_;
    foreach my $key (keys %$hash2) {
        if (exists $hash1->{$key} && ref($hash1->{$key}) eq 'HASH' && ref($hash2->{$key}) eq 'HASH') {
            merge_hashes($hash1->{$key}, $hash2->{$key});
        } else {

            $hash1->{$key} = $hash2->{$key};
        }
    }
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
			my ($x, $y, $genename, $ASscore, $readlen, $overlaplen, $genelen) = splice(@fieldsin, 0, 7);
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
my %mirnames0;
my %mirnames;



while (<IN>) {
	chomp;
	my @fieldsin = split(/\t/, $_);
	my ($x, $y, @read_begins, @read_ends, @gtf_starts, @gtf_ends, @genenames);
	while (@fieldsin) {
		($x, $y, my $gtf_name, my $ASscore, my $fastq_read_length, my $overlap_length, my $read_begin, my $read_end, my $gtf_start, my $gtf_end) = splice(@fieldsin, 0, 10);
		my $genelen = $gtf_end - $gtf_start + 1;
		#print "$gtf_name  $ASscore\n";
		#print "$fastq_read_length $overlap_length\n";
		#print "$read_begin  $read_end\n";
		#print "$gtf_start  $gtf_end\n";

		if ($doit && $ASscore ne '.' && $ASscore <= $ASthreshold && $ASscore <= $ratiothreshold_gtflen*$genelen) {
			next;
		}
		
		
		push @read_begins, $read_begin;
		push @read_ends, $read_end;
		push @gtf_starts, $gtf_start;
		push @gtf_ends, $gtf_end;
		push @genenames, $gene2numhash{$gtf_name};
	}

	

	if (scalar(@genenames ) == 1 || $opt_e)  {
		if (scalar(@genenames ) > 1) {
			my @sorted_indices = sort { $genenames[$a] cmp $genenames[$b] } 0..$#genenames;
			@genenames = @genenames[@sorted_indices];
			@read_begins = @read_begins[@sorted_indices];
			@read_ends = @read_ends[@sorted_indices];
			@gtf_starts = @gtf_starts[@sorted_indices];
			@gtf_ends = @gtf_ends[@sorted_indices];
		}
		my $realbegin = max($read_begins[0], $gtf_starts[0]);
		my $realend = min($read_ends[0], $gtf_ends[0]);
		my %coordhash;
		$matrix{"${x}_${y}_${genenames[0]}"} += 1;
		for (my $var = $realbegin; $var <= $realend; $var++) {
			$matrix{"${x}_${y}_${genenames[0]}_${var}"} += 1;
			$coordhash{$var} = 1;
		}
		if (exists($mirnames0{${genenames[0]}})) {
			$mirnames0{${genenames[0]}} = { %{$mirnames0{${genenames[0]}}}, %coordhash};
		}else{
			$mirnames0{${genenames[0]}} = { %coordhash };
		}
	}elsif (scalar(@genenames ) > 1) {
		my @sorted_indices = sort { $genenames[$a] cmp $genenames[$b] } 0..$#genenames;
		@genenames = @genenames[@sorted_indices];
		@read_begins = @read_begins[@sorted_indices];
		@read_ends = @read_ends[@sorted_indices];
		@gtf_starts = @gtf_starts[@sorted_indices];
		@gtf_ends = @gtf_ends[@sorted_indices];

		my $genename;
		my %coordhash;
		my $gtf_name = join('-',@genenames);
		$matrix{"${x}_${y}_$genenames[0]"} += 1;
		for (my $ii = 0; $ii < scalar(@genenames); $ii++) {
			my $genename = $genenames[$ii];
			my $realbegin = max($read_begins[$ii], $gtf_starts[$ii]);
			my $realend = min($read_ends[$ii], $gtf_ends[$ii]);
			for (my $var = $realbegin; $var <= $realend; $var++) {
				$matrix{"${x}_${y}_${gtf_name}_${ii}_${var}"} += 1;
				$coordhash{$ii}{$var} = 1;
			}
			if (exists($mirnames{$gtf_name})) {
				for my $var (keys(%coordhash)) {
					if (exists($mirnames{$gtf_name}{$var})) {
						$mirnames{$gtf_name}{$var} = { %{$mirnames{$gtf_name}{$var}}, %{$coordhash{$var}} };
					}
					else{
						$mirnames{$gtf_name}{$var} = { %{$coordhash{$var}} };
					}
				}
			} else {
    			$mirnames{$gtf_name} = { %coordhash };
			}
		}
	}
}
	
close IN;


my @outmirnames;
my @mirnames = sort {$a cmp $b} keys(%mirnames0);
for my $var (@mirnames) {
	my @coords1 = keys(%{$mirnames0{$var}});
	@coords1 = sort {$a <=> $b} @coords1;
	push @outmirnames, $var;
	my $varout = $num2genehash{$var};
	print OUT "\t$varout";
	for my $vi (@coords1) {
		push @outmirnames, "${var}_${vi}";
		print OUT "\t${varout}_${vi}";
	}
}

@mirnames = sort {$a cmp $b} keys(%mirnames);
for my $var (@mirnames) {
	my @coords1 = keys(%{$mirnames{$var}});
	@coords1 = sort {$a <=> $b} @coords1;
	push @outmirnames, $var;
	my @vars = map { exists $num2genehash{$_} ? $num2genehash{$_} : die "Error: No conversion found for $_\n"} split(/-/, $var);
	my $varout = join("-|:", @vars);
	print OUT "\t$varout";
	for my $vi (@coords1) {
		my %coords2 = %{$mirnames{$var}{$vi}};
		my @coords2 = sort {$a <=> $b} keys(%coords2);
		#my $temp = join(" ", @coords2);
		#print scalar @coords2;
		#print " $vi\n$temp\n";
		for my $vi2 (@coords2) {
			push @outmirnames, "${var}_${vi}_${vi2}";
			#print " $vi2\n";
			print OUT "\t${varout}_${vi}_$vi2";
		}
	}
}
print OUT "\n";



for my $var (@locations) {
	my $varout = $var;
	$varout =~ s/_/x/;
	print OUT "$varout";	
	for my $var1 (@outmirnames) {
		my $count0;
		if (exists($matrix{"${var}_${var1}"})) {
			$count0 = $matrix{"${var}_${var1}"};
		}else{
			$count0 = '';
		}
		print OUT "\t$count0";
	}
	print OUT "\n";	
}
close OUT;