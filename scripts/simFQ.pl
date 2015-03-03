#!/usr/bin/perl

# Created Aug 14 - 2014 by Thomas Nordahl Petersen tnp@cbs.dtu.dk
# Technical University of Copenhagen, 2800 Lyngby, Denmark

my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Cwd;
use strict;

#
# parameters for pair-end reads
#
my $readLen=100;
my $insertSize=250;
my $std=$insertSize*0.1;

# Default base quality
my $base_qual='U';

# output fastq files
my $r1='F.fq.gz';
my $r2='R.fq.gz';

# parameters for a log normal distribution of depth also known as X
my $depth_mean=2;
my $depth_std=8;
my $depth=0;

my $seed=999;

# percentage is the fraction of entries of a fasta file that is selected
my $percentage=100;

my $tmpDir="simFQ_$$";

# simulate an error frequency. Negative number means no polymorphisms or errors are introduced
# a value og 1000 corresponds to Q3
my $error_frequency=-1;

#
my $skipN =1;
my @bad_chars=qw(N n);
my $skip_read_pairs=0;
my $Z = 4;
my $min = $insertSize - $Z * $std;
my $max = $insertSize + $Z * $std;
my $min_coverage=0.0;
#
# Process command line
#
getopts('hi:r:d:D:I:S:q:F:R:s:x:t:E:l:NZ:c:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-d number] [-D number] [-r number] [-I number] [-S number] [-q char] [-F file] [-R file] [-s number] [-x number] [-t string] [-E number] [-N] [-Z number] [-c fraction]\n");
  print ("Description:\n");
  print ("$0 - Simulate pair-end fastq files given a fasta file - fastq files are produced with a log normal distributed depth for each entry in fasta file\n");
  print ("\n");
  print ("Log normal distributed depth's are calculated; http://en.wikipedia.org/wiki/Log-normal_distribution#Generating_log-normally_distributed_random_variates\n\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input fasta file name [STDIN]\n");
  print ("  -d  : mean depth in output fastq [$depth_mean]\n");
  print ("  -D  : depth std fastq [$depth_std]\n");
  print ("  -r  : read length [$readLen]\n");
  print ("  -I  : ave insert size [$insertSize]\n");
  print ("  -S  : ave std on insert size [$std]\n");
  print ("  -Z  : Z-score - hard boundaries on insert size ([$min,$max]) - default Z-score [$Z]\n");
  print ("  -q  : base quality for all nucleotides [$base_qual]\n");
  print ("  -F  : output in temp dir Forward reads named [$r1]\n");
  print ("  -R  : output in temp dir reverse reads named [$r2]\n");
  print ("  -s  : seed random number [$seed]\n");
  print ("  -x  : percentage of entries to be randomly selected [$percentage]\n");
  print ("  -t  : temp dir [$tmpDir]\n");
  print ("  -E  : Error frequency - (1000 corresponds to Q30) no errors if negative [$error_frequency]\n");
  print ("  -N  : skip reads if [@bad_chars] are present [on]\n");
  print ("  -c  : minimum coverage in range from 0.0 - 1.0 [$min_coverage]\n");
  print ("  -l  : logfile  [$tmpDir/logfile]\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (not defined($Getopt::Std::opt_i)){
  # Read from standard input
  *INP = *STDIN;
} 
else{
  # Read from file
  if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open(INP,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
  else{
    open(INP,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
}
if (defined($Getopt::Std::opt_r)){
    $readLen=$Getopt::Std::opt_r;
}
if (defined($Getopt::Std::opt_I)){
    $insertSize=$Getopt::Std::opt_I;
}
if ($insertSize < $readLen){
    print STDERR "\n#\n# Error: insert size \($insertSize\) can not be shorter than read length \($readLen\)\n#\n";
    exit(1);
}
if (defined($Getopt::Std::opt_t)){
    $tmpDir=$Getopt::Std::opt_t;
}
if (! -d $tmpDir){
    system "mkdir -p $tmpDir";
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
else{
    open(LOG,">$tmpDir/logfile")
}

if (defined($Getopt::Std::opt_F)){
    $r1=$Getopt::Std::opt_F;
}
if ($r1 =~ /\.gz$/){
    open(FOUT,"| gzip > $tmpDir/$r1");
}
else{
    open(FOUT,"| cat > $tmpDir/$r1");
}
if (defined($Getopt::Std::opt_R)){
    $r2=$Getopt::Std::opt_R;
}

if ($r2 =~/\.gz$/){
    open(ROUT,"| gzip > $tmpDir/$r2");
}
else{
    open(ROUT,"| cat > $tmpDir/$r2");
}
if (defined($Getopt::Std::opt_d)){
    $depth_mean=$Getopt::Std::opt_d;
}
if (defined($Getopt::Std::opt_D)){
    $depth_std=$Getopt::Std::opt_D;
}
if (defined($Getopt::Std::opt_S)){
    $std=$Getopt::Std::opt_S;
}
if (defined($Getopt::Std::opt_q)){
    $base_qual=$Getopt::Std::opt_q;
}
if (defined($Getopt::Std::opt_s)){
    $seed=$Getopt::Std::opt_s;
}
if (defined($Getopt::Std::opt_E)){
    $error_frequency=$Getopt::Std::opt_E;
}
if (defined($Getopt::Std::opt_x)){
    $percentage=$Getopt::Std::opt_x;
}
if (defined($Getopt::Std::opt_Z)){
    $Z=$Getopt::Std::opt_Z;
}
if (defined($Getopt::Std::opt_N)){
    $skipN=0;
}
if (defined($Getopt::Std::opt_c)){
    $min_coverage=$Getopt::Std::opt_c;
}
$min = $insertSize - $Z * $std;
$max = $insertSize + $Z * $std;
###############################################################################
# Main
#
###############################################################################
srand($seed);
my $thisDir=cwd();
my $datestring = localtime();
print LOG "## Local date and time $datestring - start program\n";
print LOG "# command: $command\n";
print LOG "# min allowed insert size: $min\n";
print LOG "# max allowed insert size: $max\n";
print LOG "# working dir: $thisDir\n";

my @bases=();
$bases[0]='A';
$bases[1]='T';
$bases[2]='C';
$bases[3]='G';

my %rec;
my $name;
my $qual='';
my $seq;
my $mutated_sequence;
for (my $i=1;$i<=$readLen;$i++){
    $qual .= $base_qual;
}
my $fastaEntryCounter=0;
my $x;
my $skipped=0;
my $selected=0;
my $val=0;
my $nucleotides_selected=0;
my $total_polymorphisms=0;
my $total_reads=0;

#
# log normal distributed depth for an organism
# http://en.wikipedia.org/wiki/Log-normal_distribution#Generating_log-normally_distributed_random_variates
# depth=exp(my + sigma*Z);
my $Z;
my $meanD=$depth_mean;
#my $varD=2*$depth_std * 2*$depth_std;
my $varD=$depth_std * $depth_std;
my $my=log($meanD*$meanD/sqrt($varD+$meanD*$meanD));
my $sigma=sqrt( log(1+($varD/($meanD*$meanD))) );

while (defined($_=<INP>)){
    chomp;
    if (m/^>(\S+)/){

	#
	# handle the previous sequence that is stored in $seq
	#
	if ($fastaEntryCounter>=1){

	    my $len=length($seq);
	    my $ratio = $len / $max;
	    my $ratio_flag=0;
	    if ($ratio < 1){
		$fastaEntryCounter++;
		$seq='';
		$name=$1;		
		next;
	    }
	    
	    
	    my $min_possible_depth = 2 * $readLen / $len;
	    my $depth_flag=1;
	    while ($depth_flag){
		$Z=gaussian_rand(0,1);
		$depth=exp($my+$sigma*$Z);
		if ($depth >= $min_possible_depth){
		    $depth_flag=0;
		}
	    }

	    #
	    # Select sequence for printout if random number is <= percent of entries to be selected
	    #
	    $x=int(rand(101));

	    #
	    # Should the sequence be selected ?
	    #
	    if ($x <= $percentage){
		$val=printout($seq,$name,$depth);
		if (! $val){
		    $selected++;
		}
	    }
	}
	$fastaEntryCounter++;
	$seq='';
	$name=$1;

	next;
    }

    $seq .= $_;
}
#
# print out for last read sequence
#
my $len=length($seq);
my $ratio = $len / $max;
my $ratio_flag=0;
if ($ratio >= 1){
    my $min_possible_depth = 2 * $readLen / $len;
    my $depth_flag=1;
    while ($depth_flag){
	$Z=gaussian_rand(0,1);
	$depth=exp($my+$sigma*$Z);
	if ($depth >= $min_possible_depth){
	    $depth_flag=0;
	}
    }
    
    #
    # Select sequence for printout if random number is <= percent of entries to be selected
    #
    $x=int(rand(101));
    
    #
    # Should the sequence be selected ?
    #
    if ($x <= $percentage){
	$val=printout($seq,$name,$depth);
	if (! $val){
	    $selected++;
	}
    }
}
$fastaEntryCounter++;



print LOG "# Selected $selected entries out of $fastaEntryCounter possible fasta entries\n";
$skipped = $fastaEntryCounter-$selected;
print LOG "# Skipped $skipped entries out of $fastaEntryCounter possible fasta entries\n";
print LOG "# Read pairs skipped due to presence of bad chars \(@bad_chars\): $skip_read_pairs\n";
print LOG "# Total number of reads: $total_reads\n";

my $total_bases = $total_reads * $readLen;
print LOG "# Total number of nucleotides: $total_bases\n";
print LOG "# Introduced sequence errors: $total_polymorphisms\n";

my $perc_polymorphisms = $total_polymorphisms*100/$total_bases;
printf LOG "# Percentage introduced sequencing errors: %.4f\n",$perc_polymorphisms;

$datestring = localtime();
print LOG "## Local date and time $datestring - end program\n";

sub printout{
    my ($sequence,$entry,$depth)=@_;
    my $i=0;
    my $counter=1;
    my $from;
    #
    # store time af beginning
    #
    my $datestring_start = localtime();

    #
    # Length of sequence
    #
    my $len=length($seq);


    #
    # sequencing depth is continuesly stored
    #
    my $current_depth=0;

    #
    # Number of reads simulated from this sequence
    #
    my $readsInEntry=0;

    #
    # Number of polymorphisms introduced
    #
    my $polymorphisms=0;

    #
    # whenever 2 reads are made, the depth is increased by the number add_depth
    #
    my $add_depth = 2*$readLen/$len;

    # Max number of re-tries in case insert size becomes lower than 2*read length
    my $maxTries=10000;
    my $tries=0;

    #
    # Coverage parameters
    #
    my %rec=();
    my $coverage_count=0;
    my $coverage=0;
    my $continue=1;

#    PE: while ( ($current_depth <= ($depth - $add_depth)) && ($coverage < $min_coverage)){
  PE: while ($continue){

	if ($current_depth == 0){
	    $tries++;
	    if ($tries >= $maxTries){
		print LOG "# Could not finish for entry $entry after $tries attempts\n";
		return(1);
	    }
	}

	
	#
	# chose a gausian distributed insert size - only >= 0 values
	#
	my $insert=abs(int(gaussian_rand($insertSize,$std)));
	my $rightmost_position  = $len - $insert;

	#
	# chose a random position where forward read must start
	#
	my $sub_sequence_from = int(rand($rightmost_position+1));
	my $check = $len - ($sub_sequence_from + $insert);

	while (($rightmost_position < 0) || ($check < 0) || ($insert < $readLen) || ($insert < $min) || ($insert > $max)){
#	    print STDERR "I= $insert\trightmost=  $rightmost_position\n";
	    $insert=abs(int(gaussian_rand($insertSize,$std)));
	    $rightmost_position  = $len - $insert;
	    $sub_sequence_from = int(rand($rightmost_position+1));
	    $check = $len - ($sub_sequence_from + $insert);
	}


	my $forward_seq = substr($sequence, $sub_sequence_from, $readLen);

	#
	# first position in reverse strand before its reverse complimented
	#
	my $sub_sequence_end = $sub_sequence_from + $insert - $readLen;

	#
	# reverse_seq will be reverse complimented after potential errors are introduced
	#
	my $reverse_seq = substr($sequence, $sub_sequence_end, $readLen);

	#
	# skip reads if an N is present
	#
	if ($skipN){
	    foreach my $char (@bad_chars){
		if (index($forward_seq,$char) != -1){
		    $skip_read_pairs++;
		    next PE;
		}
		if (index($reverse_seq,$char) != -1){
		    $skip_read_pairs++;
		    next PE;
		}
	    }
	}
	if ($error_frequency > 0){
	    my $variants = 0;
	    ($forward_seq, $variants) = introduce_errors($forward_seq,$error_frequency);
	    #print LOG "$variants\t$forward_seq\n";
	    $total_polymorphisms += $variants;
	    $polymorphisms += $variants;

	    ($reverse_seq, $variants) = introduce_errors($reverse_seq,$error_frequency);
	    $total_polymorphisms += $variants;
	    $polymorphisms += $variants;
	}

	my $Rrev=reverse($reverse_seq);
	$Rrev =~ tr/ACGTacgt/TGCAtgca/;
#	my $lenF=length($forward_seq);
#	my $lenR=length($Rrev);
#	if ($lenF != $lenR){
#	    print STDERR "Error - sequence length for forward and reverse differ\n";
#	    print STDERR "\@$entry:F$sub_sequence_from:L$readLen:I$insert:$counter 1\n$forward_seq\n+\n$qual\n";
#	    print STDERR "\@$entry:F$sub_sequence_from:L$readLen:I$insert:$counter 2\n$Rrev\n+\n$qual\n";	    
#	    exit(1);
#	}

	#
	# add 1 to the zero-based position such that it fits with normal counting in bwa or blast
	#
	$from = $sub_sequence_from +1;
	print FOUT "\@$entry:F$from:L$readLen:I$insert:$counter 1\n$forward_seq\n+\n$qual\n";
	print ROUT "\@$entry:F$from:L$readLen:I$insert:$counter 2\n$Rrev\n+\n$qual\n";

	#
	# Forward strand
	#
	my $forward_to = $from + $readLen;
	for (my $i=$from;$i<=$forward_to;$i++){
	    if (! exists($rec{$i})){
		$rec{$i}=1;
		$coverage_count++;
	    }
	}
	#
	# reverse strand
	#
	my $rev_from = $from + $insert - $readLen;
	my $rev_to = $rev_from + $readLen;
	for (my $i=$rev_from; $i<= $rev_to;$i++){
	    if (! exists($rec{$i})){
		$rec{$i}=$i;
		$coverage_count++;
	    }
	}
	$current_depth += $add_depth;
	$readsInEntry += 2;
	$total_reads += 2;
	$counter++;
	$coverage = $coverage_count/$len;

	if (($current_depth > ($depth - $add_depth)) && ($coverage >= $min_coverage)){
	    $continue=0;
	}
    }
    
    if ($readsInEntry >0){
	printf LOG "$entry\tLen= %d\tDepth= %.4f\tReads= %d\tPolymorphisms= %d",$len,$current_depth,$readsInEntry,$polymorphisms;
	printf LOG "\tCov= %.4f\n",$coverage;
	return(0);
    }
    else{
	print LOG "# Failed to make any pair-end reads for entry $entry\tLen=$len\ttarget_depth=$depth\n";
	return(1);
    }
}

sub gaussian_rand {
    my ($mean,$std)=@_;

    # uniformly distributed random numbers
    my ($u1, $u2);

    # variance, then a weight
    my $w;

    # gaussian-distributed numbers
    my ($g1, $g2);
    do { $u1 = 2 * rand() - 1; $u2 = 2 * rand() - 1; $w = $u1*$u1 + $u2*$u2; } while ( $w >= 1 ); $w = sqrt( (-2 * log($w)) / $w ); $g2 = $u1 * $w; $g1 = $u2 * $w; # return both if wanted, else just one
return wantarray ? ($g1, $g2) : $g1*$std+$mean; }


sub introduce_errors{
    my ($seq, $freq) = @_;
    my $len = length($seq);
    my $polymorphisms=0;
    my $polymorphism;

    
    for (my $i=0;$i<$len;$i++){
	
	if (int(rand($freq)) == 0){
	    my $continue=1;
	    $polymorphisms++;
	    while ($continue){
		$polymorphism=$bases[int(rand(4))];
		if ($polymorphism ne substr($seq,$i,1)){
		    $continue=0;
		}
	    }
	    substr($seq,$i,1,$polymorphism);
	}
    }
    return($seq, $polymorphisms);
}
