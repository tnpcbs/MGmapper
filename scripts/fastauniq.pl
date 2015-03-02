#!/usr/bin/perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Digest::MD5 qw(md5_hex);
use strict;
use Cwd;
my $columns=60;
# Default parameters
*LOG=*STDERR;
my $verbose=0;
#
# Process command line
#
getopts('hi:o:vl:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-l logfile] [-v]\n");
  print ("Description:\n");
  print ("$0 - Rejects identical fasta entries based on md5sum on the sequence only\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name [STDIN]\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
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
#
# If not file name is given, use standard output
#
if (not defined($Getopt::Std::opt_o)){
  # Output goes to std output
  *OUT = *STDOUT;
} else {
  # Open file to write to
  open(OUT, "| cat >$Getopt::Std::opt_o") || die ("can't open file $Getopt::Std::opt_o: $!");
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
print LOG "## Local date and time $datestring - start program\n" if ($verbose);
print LOG "# $command\n" if ($verbose);
my $thisDir=cwd();
print LOG "# working dir: $thisDir\n" if ($verbose);
my $seq='';
my $first=1;
my $line;
my $newline;
my $md5;
my %rec=();
my $entries=0;
my $redundant=0;
print LOG "# Kept entry\tRejected 100% identical sequence\n" if ($verbose);
my %entryNames=();
my @duplicateEntryNames=();

while (defined($_=<INP>)){
    chomp;
    if (m/^>/){
	if ($first){
	    $first=0;
	}
	else{
	    my @tmp=split(/\s+/,$line);
	    my $name=$tmp[0];

	    #
	    # If duplicate entry names exists, then only write out one of them
	    #
	    if (! exists($entryNames{$name})){
		$entryNames{$name}=1;
		$md5 = md5_hex($seq);
		if (! exists($rec{$md5})){
		    $entries++;
		    $rec{$md5}=$line;
		    
		    print OUT "\>$line\n";
		    $seq =~ s/\s+//g;
		    while (my $chunk = substr($seq, 0, $columns, "")) {
			print OUT "$chunk\n";
		    }
		}
		else{
		    $redundant++;
		    my @w=split(/\s+/,$rec{$md5});
		    printf LOG "%s\t%s\n",$w[0],$line if ($verbose);
		}
	    }
	    else{
		push(@duplicateEntryNames,$name);
	    }
	    $seq='';
	}

	#
	# replace special characters in name of the fasta sequence
	#
	my @w=split(/\s+/);
	my $id = $w[0];
	my $replace=0;
#	if ($id =~ /[\'\(\)]/){
	if ($id =~ /\'/){
#	    print STDERR "replacing chars in $id\n";
#	    $id =~ s/\)/_/g;
#	    $id =~ s/\(/_/g;
	    $id =~ s/\'/_/g;
	    $replace=1;
	}
	$newline = "$id";
	if ($replace){
	    my $name=substr($w[0],1);
	    $newline .= " original name='$name'";
	}
	my $j=1;
	while (defined($w[$j])){
	    $newline .= " $w[$j]";
	    $j++;
	}
	$line=$newline;
	$line = substr($line,1);
	next;
    }
    $seq .= $_;
}
$md5 = md5_hex($seq);

my @tmp=split(/\s+/,$line);
my $name=$tmp[0];

if (! exists($entryNames{$name})){
    $entryNames{$name}=1;

    if (! exists($rec{$md5})){
	$entries++;
	$rec{$md5}=$line;
	
	print OUT "\>$line\n";
	$seq =~ s/\s+//g;
	while (my $chunk = substr($seq, 0, $columns, "")) {
	    print OUT "$chunk\n";
	}
    }
    else{
	$redundant++;
	my @w=split(/\s+/,$rec{$md5});
	printf LOG "%s\t%s\n",$w[0],$line if ($verbose);
    }
}
else{
    push(@duplicateEntryNames,$name);
}

my $duplicateEntries=$#duplicateEntryNames + 1;
if ($duplicateEntries > -1){
    foreach my $id (@duplicateEntryNames){
	print LOG "# Duplicate entries observed for: $id\n" if ($verbose);
    }
}
print LOG "# ok entries: $entries\trejected entries: $redundant + $duplicateEntries duplicate sequence names\n" if ($verbose);
$datestring = localtime();
print LOG "## Local date and time $datestring - done\n" if ($verbose);
