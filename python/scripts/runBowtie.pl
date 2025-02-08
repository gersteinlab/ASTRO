use strict;
use warnings;
use Getopt::Long;


my ($referencefile, $inputfa, $outputfile, $thread);
GetOptions(
    "r|referencefile=s" => \$referencefile,
    "i|inputfa=s"       => \$inputfa,
    "o|outputfile=s"    => \$outputfile,
    "t|thread=s"        => \$thread
);
$thread //= '8'; 


my @cmd = (
    'bowtie2', '-p', $thread, '-x', $referencefile, '-U', $inputfa, '-D', '20', '-R', '3', '-N', '1', '-L', '6',
    '-i', 'C,5', '--score-min', 'C,10', '--no-hd', '--no-sq', '-a', '--ma', '1', '--rdg', '0,1', '--rfg', '0,1',
    '--mp', '1,1', '--gbar', '1', '--local'
);


open my $process, '-|', @cmd or die "Could not start process: $!";


open my $output_file, '>', $outputfile or die "Could not open file '$outputfile': $!";
my ($major_read, $major_map, $major_value, $see1);

$see1 = 0;
while (my $line = <$process>) {

    chomp $line;
    my @parts = split(/\s+/, $line);
    my $current_read = $parts[0];
    my $current_map = $parts[2];
    my $current_value;
    if ($current_map ne '*') {
        unless ( ($current_value) = $parts[11] =~ /^AS:i:(\d+)$/) {
            print STDERR "Error in line: $line\n";
            die "Tag string does not match the required format AS:i:x";
        }
    }else{
        print $output_file "$current_read\t*\n";
        next;
    }

    if (!defined $major_read || $current_read ne $major_read) {
        print $output_file "$major_read\t$major_map\n" if defined $major_read;
        $see1 = 1;
        $major_read = $current_read;
        $major_value = $current_value;
        $major_map = $parts[2];
    }else{
        if ($see1) {
            if ($major_value <= $current_value) {
                $major_map = '*';
            }
            $see1 = 0;
        }
    }
}

if (defined $major_read) {
    print $output_file "$major_read\t$major_map\n";
}
close $output_file;
close $process;
