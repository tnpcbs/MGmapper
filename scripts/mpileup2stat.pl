#!/usr/bin/perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use strict;

# Default parameters
#$tmpdir=$$;
*LOG=*STDERR;
#
# Process command line
#
getopts('hi:o:vl:a:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-a indexFile] [-l logfile] [-v]\n");
  print ("Description:\n");
  print ("$0 - Read and write files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input mpileup file name [STDIN]\n");
  print ("  -a  : bwa index file\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (! defined($Getopt::Std::opt_a)){
    print "Option -a has not been defined\n";
    exit;
}
else{
    open(IDX,"<$Getopt::Std::opt_a") || die ("File not found: $Getopt::Std::opt_a");
}

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
  open(OUT, ">$Getopt::Std::opt_o") || die ("can't open file $Getopt::Std::opt_o: $!");
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
###############################################################################
# Main
#
###############################################################################
if (defined($Getopt::Std::opt_v)){
    print LOG ("# $command\n");
}
my %rec=();
#
# Read info in Index file
#
$_=<IDX>;
while (! eof (IDX)){
    $_=<IDX>;
    chomp;
    my $name;
    my $desc;
    my $size;
    if (m/^\d+ (\S+) (.+)/){
	$name=$1;
	$desc=$2;
	$rec{$name}{name}=$name;
	$rec{$name}{desc}=$desc;
    }
    elsif (m/^\d+ (\S+)/){
	$name=$1;
	$rec{$name}{name}=$name;
	$rec{$name}{desc}='NotDefined';
    }
    else{
	print LOG "Error parsing file: $Getopt::Std::opt_a\n";
    }
    $_=<IDX>;
    chomp;
    if (m/^\d+ (\d+) \d+/){
	$size=$1;
    }
    else{
	print LOG  "Error parsing File: $Getopt::Std::opt_a\n";
    }
    $rec{$name}{size}=$size;
}
close(IDX);

#
# Read the mpileup file
#
my $nucTotal=0;
my $readTotal=0;
while (! eof (INP)){
    $_=<INP>;
    chomp;
    if (m/^(\S+)\s+(\d+)\s+\w+\s+(\d+)\s+(\S+)\s+\S+\s+(\d+)\s+\S+\s+\S+/){
	my $name = $1;
	my $pos = $2;
	my $line_ALT = $4;
	my $count = $3 + $5;

	if ($count >= 1){
	    $rec{$name}{hits}++;
	}

        my @tmp = split(/\^!/,$line_ALT);
	my $countRead = scalar(@tmp) -1;

	$rec{$name}{reads} += $countRead;
	$readTotal += $countRead;

	$rec{$name}{nucleotides} += $count;
	$nucTotal += $count;
    }
    elsif (m/^(\S+)\s+(\d+)\s+\w+\s+(\d+)\s+(\S+)\s+\S+/){
	my $name = $1;
	my $pos = $2;
	my $count = $3;
	my $line_ALT = $4;
	if ($count >= 1){
	    $rec{$name}{hits}++;
	}
	$rec{$name}{nucleotides} += $count;
	$nucTotal += $count;

        my @tmp = split(/\^/,$line_ALT);
	my $countRead = scalar(@tmp) -1;
	$rec{$name}{reads} += $countRead;
	$readTotal += $countRead;
    }
}

foreach my $key (keys %rec){
    if (exists ($rec{$key}{hits})){
	my $depth = $rec{$key}{nucleotides}/$rec{$key}{size};
	my $coverage = $rec{$key}{hits}/$rec{$key}{size};
	my $percentNuc = $rec{$key}{nucleotides}*100/$nucTotal;
	my $reads = $rec{$key}{reads};
	my $percentReads = $rec{$key}{reads}*100/$readTotal;
	printf OUT "%-25s\t%7.4f\t%7.4f\t%d\t%7.4f\t$rec{$key}{hits}\t$rec{$key}{size}\t%7.4f\t%d\t$rec{$key}{desc}\n",$rec{$key}{name},$percentNuc,$depth, $rec{$key}{nucleotides},$coverage,$percentReads,$reads;
    }
}
