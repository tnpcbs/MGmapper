#!/usr/bin/perl
use Getopt::Std;
use strict;
#
# Configuration
#

my $maxreadentries = 250; # In 1000


#
# Process command line
#
sub Usage {
  print ("Usage: $0 [-h] [-s readsize] -f <name> -r <name> -a <name> -b <name>\n");
  print ("Description:\n");
  print ("$0 - find reads in commen given 2 input fastq files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -s  : read size (memory use and performance option, the bigger the better) default 250 which uses 1 GB RAM\n");
  print ("  -f  : input forward fastq\n");
  print ("  -r  : input reverse fastq\n");
  print ("  -a  : output forward fastq\n");
  print ("  -b  : output reverse fastq\n");
  print ("  -A  : output singletons from forward fastq, optional\n");
  print ("  -B  : output singletons from reverse fastq, optional\n");
  print ("\n");
  exit;
} # Usage
getopts('hf:r:a:b:s:A:B:')||Usage();
#
# Usage
#
&Usage($Getopt::Std::opt_h) if defined($Getopt::Std::opt_h);
&Usage() unless defined($Getopt::Std::opt_f) and defined($Getopt::Std::opt_r);
&Usage() unless defined($Getopt::Std::opt_a) and defined($Getopt::Std::opt_b);
&Usage() if defined($Getopt::Std::opt_s) and $Getopt::Std::opt_s !~ m/^\d+$/;

# calulate read size (performance)
# The larger, the faster, the more memory. 250 is around 1 GB
$maxreadentries = $Getopt::Std::opt_s if defined $Getopt::Std::opt_s;
$maxreadentries *= 1000;

#
# Open forward fastq file
#
if (($Getopt::Std::opt_f=~/\.gz$/) or ($Getopt::Std::opt_f=~/\.Z$/)){
   open(F,'-|', 'gunzip', '-c', $Getopt::Std::opt_f) or die ("can't open file $Getopt::Std::opt_f: $!");
} else {
   open(F,'<', $Getopt::Std::opt_f) or die ("can't open file $Getopt::Std::opt_f: $!");
}

#
# Open reverse fastq file
#
if (($Getopt::Std::opt_r=~/\.gz$/) or ($Getopt::Std::opt_r=~/\.Z$/)){
   open(R,'-|', 'gunzip', '-c', $Getopt::Std::opt_r) or die ("can't open file $Getopt::Std::opt_r: $!");
} else {
   open(R,'<', $Getopt::Std::opt_r) or die ("can't open file $Getopt::Std::opt_r: $!");
}


# No logging
#if (defined($Getopt::Std::opt_l)){
#    open(LOG,">$Getopt::Std::opt_l");
#}


# Output files
if (($Getopt::Std::opt_a=~/\.gz$/) or ($Getopt::Std::opt_a=~/\.Z$/)){
   open(OUTF,'|-', "gzip -c > $Getopt::Std::opt_a") or die ("can't write file $Getopt::Std::opt_a: $!");
} else {
   open(OUTF, '>', $Getopt::Std::opt_a) or die ("can't write file $Getopt::Std::opt_a: $!");
}
if (($Getopt::Std::opt_b=~/\.gz$/) or ($Getopt::Std::opt_b=~/\.Z$/)){
   open(OUTR,'|-', "gzip -c > $Getopt::Std::opt_b") or die ("can't write file $Getopt::Std::opt_b: $!");
} else {
   open(OUTR, '>', $Getopt::Std::opt_b) or die ("can't write file $Getopt::Std::opt_b: $!");
}
# optional singleton output
unless (defined $Getopt::Std::opt_A) {
} elsif (($Getopt::Std::opt_A=~/\.gz$/) or ($Getopt::Std::opt_A=~/\.Z$/)){
   open(OUTSF,'|-', "gzip -c > $Getopt::Std::opt_A") or die ("can't write file $Getopt::Std::opt_A: $!");
} else {
   open(OUTSF, '>', $Getopt::Std::opt_A) or die ("can't write file $Getopt::Std::opt_A: $!");
}
unless (defined $Getopt::Std::opt_B) {
} elsif (($Getopt::Std::opt_B=~/\.gz$/) or ($Getopt::Std::opt_B=~/\.Z$/)){
   open(OUTSR,'|-', "gzip -c > $Getopt::Std::opt_B") or die ("can't write file $Getopt::Std::opt_B: $!");
} else {
   open(OUTSR, '>', $Getopt::Std::opt_B) or die ("can't write file $Getopt::Std::opt_B: $!");
}

###############################################################################
# Main
#
###############################################################################

# Initialise variables, make perl aware of the sizes
my $done = 0;
=pod
my $buf_f = 'x' x $buffersize;
my $buf_r = 'x' x $buffersize;
my @f;
$f[$arraysize] = 1;
my @r;
$r[$arraysize] = 1;
my @fkey;
$fkey[$arraysize] = 1;
my @rkey;
$rkey[$arraysize] = 1;
my %hash = map($_ => $_, 0..$arraysize);
$buf_f = '';
$buf_r = '';
@f = ();
@r = ();
@fkey = ();
@rkey = ();
=cut
my ($buf_f, $buf_r) = ('', '');
my (@f, @r, @rkey, @fkey, %hash);
until ($done) {
   # Read a big block from both filehandles
   my $countf = $maxreadentries;
   while (--$countf and defined(my $line = <F>)) {
      next unless $line =~ m/^@(\S+)[\/|#|\t ]/;
      push(@fkey, $1);
      $line .= <F> . <F> . <F>;
      push(@f, $line);
   }
   my $countr = $maxreadentries;
   while (--$countr and defined(my $line = <R>)) {
      next unless $line =~ m/^@(\S+)[\/|#|\t ]/;
      push(@rkey, $1);
      $line .= <R> . <R> . <R>;
      push(@r, $line);
   }
   $done = 1 if $countr and $countf;
   # Find commons quick
   my ($indexf, $indexr) = (-1, -1);
   my (@fspos, $rspos);
   %hash = ();
   for (my $i = $#rkey; $i >= 0; $i--) {
      $hash{$rkey[$i]} = $i;
   }
   for (my $i = 0; $i <= $#fkey; $i++) {
      next unless exists $hash{$fkey[$i]};
      $buf_f .= $f[$indexf = $i];
      $buf_r .= $r[$indexr = $hash{$fkey[$i]}];
   }
   # Find singletons
   if (defined $Getopt::Std::opt_A) {
      my $target = $done ? $#fkey : $indexf-1;
      my $buffer = '';
      for (my $i = 0; $i <= $target; $i++) {
         $buffer .= $f[$i] unless exists $hash{$fkey[$i]};
      }
      print OUTSF $buffer;
   }
   if (defined $Getopt::Std::opt_B) {
      %hash = ();
      for (my $i = $#fkey; $i >= 0; $i--) {
         $hash{$fkey[$i]} = $i;
      }
      my $target = $done ? $#rkey : $indexr-1;
      my $buffer = '';
      for (my $i = 0; $i <= $target; $i++) {
         $buffer .= $f[$i] unless exists $hash{$rkey[$i]};
      }
      print OUTSR $buffer;
   }
   # Clean arrays
   if ($indexf >= $#fkey) {
      @f = ();
      @fkey = ();
   } else {
      splice(@f, 0, $indexf+1);
      splice(@fkey, 0, $indexf+1);
   }
   if ($indexr >= $#rkey) {
      @r = ();
      @rkey = ();
   } else {
      splice(@r, 0, $indexr+1);
      splice(@rkey, 0, $indexr+1);
   }
   # Output
   print OUTF $buf_f;
   $buf_f = '';
   print OUTR $buf_r;
   $buf_r = '';
}

close(F);
close(R);
close(OUTF);
close(OUTR);
close(OUTSF) if defined $Getopt::Std::opt_A;
close(OUTSR) if defined $Getopt::Std::opt_B;

