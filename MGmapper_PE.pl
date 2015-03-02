#!/usr/bin/perl
BEGIN{
    $ENV{LANG}='C';
    *LOG=*STDERR;
}

umask 0022;
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Cwd;
use strict;
use lib "$ENV{MGmapper}/modules";
use Parallel::ChildManager;

#
# program calls
#

my $prog_xlsWrite = "$ENV{MGmapper}/scripts/MGcsv2xls.py";
my $prog_pileup2fasta = "$ENV{MGmapper}/scripts/pileup2fasta.pl";
my $prog_mpileup2stat = "$ENV{MGmapper}/scripts/mpileup2stat.pl";
my $prog_readsInCommon = "$ENV{MGmapper}/scripts/readsInCommon.pl";

#######################################################################
#
# THESE PROGRAMS ARE NEEDED FOR MGmapper TO RUN - MODIFY PATH IF NEEDED
#
#######################################################################
my $prog_bwa = " /tools/ngs/bwa-0.7.10/bwa";
my $prog_pigz = 'pigz';
my $prog_samtools = "/tools/bin/samtools";
my $prog_bamtools = "/tools/bin/bamtools";
my $prog_cutadapt = "/tools/bin/cutadapt";
#######################################################################
#
# NOTHING BELOW NEEDS TO BE CHANGED
#
#######################################################################
my $origDir = cwd();
#
# Process command line
#
getopts('hi:j:I:J:d:rkKq:m:pF:C:xc:t:Te:EISbVNoVn:ws:a:X:vP:Q:M:D:H:B:R')||Usage();
my $parameterFile = "$ENV{MGmapper}/databases.txt";
my $parameterFile_adapter = "$ENV{MGmapper}/adapter.txt";
my $datestring;

if (defined($Getopt::Std::opt_P)){
    $parameterFile = $Getopt::Std::opt_P;
}
if (! -e "$parameterFile"){
    print LOG "can't open file $parameterFile\n";
    print LOG "Done!\n";
    die;
}
open(PARAMETER,"<$parameterFile");
if (defined($Getopt::Std::opt_Q)){
    $parameterFile_adapter = $Getopt::Std::opt_Q;
}
my $databaseCounter=0;
my %databaseNames=();
my @dbNames=();
my %idx=();
my %databaseCount=();
my $meanInsertSize;
my %countEntries=();

while (defined($_=<PARAMETER>)){
    if (/^\#/){
	next;
    }
    chomp;
    my @w=split(/\s+/);
    my $file = $w[0];
    my $id = $w[1];
    $databaseNames{$databaseCounter} = $id;
    $databaseCount{$id} = $databaseCounter;
    push(@dbNames,$id);
    
    my $i=1;
    my $splitFile = $file . ".$i";
    my $flag=1;
    #
    # Big fasta files can be split into smaller chunks and called with a suffix being 1 2 3 ... N i.e. Bacteria.1 Bacteria.2  
    # Each fasta file must of-course be indexed with both bwa index and samtools faidx
    #
    if (-e $splitFile){
	while ($flag){
	    my $split = "$file" . ".$i";
	    if (-e $split){
		$idx{index}{$id}[$i] = $split;
		$idx{ann}{$id}[$i] = "$split" . '.ann';
	    }
	    else{
		$flag=0;
	    }
	    $i++;
	}
    }
    else{
	$idx{index}{$id}[$i] = $file;
	$idx{ann}{$id}[$i] = "$file" . '.ann';	
    }
    my $dbNum = $databaseCount{$id};
    my $dbName = $databaseNames{$dbNum};
    $idx{dbNum}{$dbName}=$dbNum;
    $idx{dbName}{$dbNum}=$dbName;

    my $i=2;
    while (defined($w[$i])){
	$idx{desc}{$dbNum} .= "$w[$i] ";
	$i++;
    }
    # Lookup file for renaming a fasta entry and showing location of index files
    # for each individual fasta entry.
    $idx{lookup}{$dbName} = $file . '.lookup';

    #
    # A logfile showing redundant entries that were removed from the database
    #
    my $fastauniq = $file . '.fastauniq.log';
    if (-e $fastauniq){
	$idx{fastauniq}{$dbName} = $fastauniq;
    }
    $databaseCounter++;
}
$databaseCounter -= 1;

my %insertSizeDistrib=();
my $insertSizeAverage;
my $insertSizeStd=0;
my $i=0;
my @rmMapNotList=();

#
# Default parameters
#

my $wwwMode=0;
my $fileIn;
my $inFileF;
my $inFileR;
my $inFileI;
my $inFileJ;
my $listF;
my $listR;
my $listIn;
my $workDir = "mapper_$$";
my $cleanFile;
my $redo=0;
my $commonName;
my $commonNameF;
my $commonNameR;
my $qualityCut = 30;
my $minLen = 30;
my $cores = 1;
my $matchRatio = 0.8;
my $matchRatioFlag = 1;
my $absoluteMatchCount = $minLen;
my $absoluteMatchCountFlag = 0;
my $mapNotPhiXReads=0;
my $sampleName = 'sample';
my $SNP_threshold=4;
my $skip_chain_fraction = 0.5;
my $XA_hits =5;
my $QUALITY_BASE = 33;
#
# Total number of reads that are ok either from the start or after cutadapt
#
my $okReads=0;
my $mapNotReads=0;
my @databases_fullmode=();
my @databases_chainmode=();
my @databases_all=();
my %abundance=();
my %fastaOut=();
my $lastMapNot;
my @rmFiles;
my $verbose=1;
my $cleanup=1;
my $cleanupCutadapt=1;
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-P parameter file] [-i fastqFile] [-j fastqFile [-I fastqList] [-J fastqList] [-d outDir] [-q number] [-m number] [-F num serarate w comma] [-C numbers separate w comma] [-c cores] [-t number] [-T] [-e number] [-E] [-s number] [-r] [-k] [-K] [-S] [-x] [-o] [-p] [-n name] [-D number] [-H number] [-M number] [-R]\n");
  print ("Description:\n");
  print ("$0 - Map Pair-end reads against reference sequence databases\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -P  : parameter file with reference to database [$parameterFile]\n");
  print ("  -Q  : parameter file with adapters removed by cutadapt [$parameterFile_adapter]\n");
  print ("  -c  : cores to use [$cores]\n");
  print ("  Input files:\n");
  print ("  -i  : input forward fastq file\n");
  print ("  -j  : input reverse fastq file\n");
  print ("  -I  : input forward fastq list\n");
  print ("  -J  : input reverse fastq list\n");
  print ("  Parameters to filter hits:\n");
  print ("  -t  : true hits have a Match/readLen >= fraction [$matchRatio]\n");
  print ("  -T  : option -t is active [on]\n");
  print ("  -e  : minimum number of Matches+Mismatches for a valid read [$absoluteMatchCount] - only active if -E is defined\n");
  print ("  -E  : option -e is active [off]\n");
  print ("  -D  : drop chains shorter than FLOAT fraction of the longest overlapping chain [$skip_chain_fraction]\n");
  print ("  -H  : if there are <INT hits with score >80% of the max score, output all in XA [$XA_hits]\n");
  print ("  Cutadapt parameters:\n");
  print ("  -S  : skip cutadapt [off]\n");
  print ("  -m  : minimum read length in cutadapt [$minLen]\n");
  print ("  -q  : Quality q cutoff in cutadapt [$qualityCut]\n");
  print ("  -B  : QUALITY_BASE, quality values are encoded as ascii(quality + QUALITY_BASE) [$QUALITY_BASE]\n");
  print ("  Mapping mode:\n");
  print ("  -F  : map reads in Full mode against these databases - comma separated numbers\n");
  print ("  -C  : map reads in Chain mode against these databases - comma separated numbers - order is important\n");
  print ("  -p  : properly paired reads only in mapTo [on]\n");
  print ("  Output files\/directory\n");
  print ("  -d  : output directory [$workDir]\n");
  print ("  -o  : output unmapped fastq reads as 1 interleaved file [off] or default 2 fastq files with F and R reads\n");
  print ("  -n  : Sample name [$sampleName]\n");
  print ("  Assembeled fasta files and parameters:\n");
  print ("  -a  : make assembled fasta files for reads mapping to these databases - comma separated numbers - order dont matter ex.'1,2'\n");
  print ("  -s  : Min number of occurences to accept a variation in contig fasta file [$SNP_threshold]\n");
  print ("  -R  : make Read count matrices [off]\n");
  print ("  -k  : keep all files [off]\n");
  print ("  -K  : keep all cutadapt files [off]\n");
  print ("  -x  : dont run - just show settings and exit [off]\n");
  print ("---Databases---\n");
  for (my $i=1;$i<=$databaseCounter;$i++){
      print ("   $i\t$databaseNames{$i}\t$idx{desc}{$i}\n");
  }
 exit;
} # Usage
#
# A working directory for all output files
#

if (defined($Getopt::Std::opt_M)){
    $meanInsertSize = $Getopt::Std::opt_M;
}
if (defined($Getopt::Std::opt_D)){
    $skip_chain_fraction = $Getopt::Std::opt_D;
}
if (defined($Getopt::Std::opt_H)){
    $XA_hits = $Getopt::Std::opt_H;
}
if (defined($Getopt::Std::opt_d)){
    $workDir = $Getopt::Std::opt_d;
}

if ( (! -d $workDir) && (! defined($Getopt::Std::opt_x)) ){
    system ("mkdir -p $workDir; chmod 0755 $workDir");
}

#if (defined($Getopt::Std::opt_v)){
#    $verbose=0;
#}

if ( ($verbose) && (! defined($Getopt::Std::opt_x))){
    open(LOG,">$workDir/MGmapper.log");
}
if (defined($Getopt::Std::opt_v)){
    print LOG "Option -v has no functionality anymore - it is always on and will produce a logfile called MGmapper.log in output directory\n";
}
if (defined($Getopt::Std::opt_w)){
    $wwwMode=1;
}
if (defined($Getopt::Std::opt_s)){
    $SNP_threshold=$Getopt::Std::opt_s;
}

if (defined($Getopt::Std::opt_r)){
    $redo=1;
}
if (defined($Getopt::Std::opt_n)){
    $sampleName=$Getopt::Std::opt_n;
}

if (defined($Getopt::Std::opt_c)){
    $cores = $Getopt::Std::opt_c;
}
my $cm = new ChildManager($cores);

if (defined($Getopt::Std::opt_i)){
    $inFileF = $Getopt::Std::opt_i;
    if (! -e $inFileF){
	print LOG "File not found: $inFileF\n";
	print STDERR "File not found: $inFileF\nDone!\n";
	exit;
    }
}

if (defined($Getopt::Std::opt_j)){
    $inFileR = $Getopt::Std::opt_j;
    if (! -e $inFileR){
	print LOG "File not found: $inFileR\n";
	print STDERR "File not found: $inFileR\nDone!\n";
	die;
    }
    if ($inFileF eq $inFileR){
	print LOG "File names for input forward and reverse fastq names are identical\n";
	die;
    }
}
#
# Quality score q in cutadapt 
#
if (defined($Getopt::Std::opt_q)){
    $qualityCut = $Getopt::Std::opt_q;
}
#
# minimum length of read after cutadapt
#
if (defined($Getopt::Std::opt_m)){
    $minLen = $Getopt::Std::opt_m;
}
#
# option to activate option -t, default is on
#
if (defined($Getopt::Std::opt_T)){
    $matchRatioFlag=0;
}
if (defined($Getopt::Std::opt_e)){
    $absoluteMatchCount=$Getopt::Std::opt_e;
}

#
# activate option -e, default is off
#
if (defined($Getopt::Std::opt_E)){
    $absoluteMatchCountFlag=1;
}
#
# Databases to map against phiX
#
my %runningMode=();
if (defined($Getopt::Std::opt_F)){
    @databases_fullmode=();
    my @choices = split(/\,/,$Getopt::Std::opt_F);
    my $i=0;
    while (defined $choices[$i]){
	if ($choices[$i] eq '0'){
	    @databases_fullmode=();
	    $i++;
	    next;
	}
	if (! exists $databaseNames{$choices[$i]}){	    
	    print "Error with options to -F\n";
	    print STDERR "Error with options to -F\nDone!";
	    exit;
	}
	$runningMode{$databaseNames{$choices[$i]}}='Fullmode';
	push(@databases_fullmode,$databaseNames{$choices[$i]});
	my $j=0;
	while (defined($databases_all[$j])){
	    if ($databases_all[$j] eq $databaseNames{$choices[$i]}){
		print LOG "Database $databaseNames{$choices[$i]} has been selected more than 1 time - i exit now\n";
		print STDERR "Done!\n";
		exit;
	    }
	    $j++;
	}
	push(@databases_all,$databaseNames{$choices[$i]});
	$i++;
    }
}

#
# Databases to map against in a chain
#
if (defined($Getopt::Std::opt_C)){
    @databases_chainmode=();
    my @choices = split(/\,/,$Getopt::Std::opt_C);
    my $i=0;
    while (defined $choices[$i]){
	if ($choices[$i] eq '0'){
	    @databases_chainmode=();
	    $i++;
	    next;
	}
	if (! exists $databaseNames{$choices[$i]}){
	    print LOG "Error with options to -C\n";
	    print STDERR "Error with options to -C\nDone!\n";
	    exit;
	}
	$runningMode{$databaseNames{$choices[$i]}}='Chainmode';
	push(@databases_chainmode,$databaseNames{$choices[$i]});
	my $j=0;
	while (defined($databases_all[$j])){
	    if ($databases_all[$j] eq $databaseNames{$choices[$i]}){
		print LOG "Database $databaseNames{$choices[$i]} has been selected more than 1 time - i exit now\n";
		print STDERR "Done!\n";
		exit;
	    }
	    $j++;
	}
	push(@databases_all,$databaseNames{$choices[$i]});
	$i++;
    }
}

if (defined($Getopt::Std::opt_t)){
    $matchRatio=$Getopt::Std::opt_t;
}
#if (defined($Getopt::Std::opt_M)){
#    $maxInsertSize=$Getopt::Std::opt_M;
#}
my $noclean="$workDir/noclean";
my $currentDir = cwd();

if (defined($Getopt::Std::opt_k)){
    $cleanup=0;
}
if (defined($Getopt::Std::opt_K)){
    $cleanupCutadapt=0;
}
if (defined($Getopt::Std::opt_B)){
    $QUALITY_BASE=$Getopt::Std::opt_B;
}

print LOG ("# command: $command\n\n");
print LOG ("# Current directory: $currentDir\n");
if (defined($Getopt::Std::opt_i)){
    print LOG "# -i: file with reads: $inFileF\n";
}

if (defined($Getopt::Std::opt_I)){
    $listF = $Getopt::Std::opt_I;
    if (! -e $listF){
	print "File not found: $listF\n";
	print STDERR "File not found: $listF\nDone!\n";
	exit;
    }
    open(LISTF,"<$listF");
    print LOG "# -I: list with reads: $listF\n";
}

if (defined($Getopt::Std::opt_j)){
    print LOG "# -j: file with mate reads: $inFileR\n";
}
if (defined($Getopt::Std::opt_J)){
    $listR = $Getopt::Std::opt_J;
    if (! -e $listR){
	print  LOG "File not found: $listR\n";
	print STDERR "File not found: $listR\nDone!\n";
	exit;
    }
    open(LISTR,"<$listR");
    print LOG "# -J: list with mate reads: $listR\n";
}

if ($#databases_fullmode >= 0){
    print LOG "# -F: Full mode: map reads against these databases: @databases_fullmode\n";
}
if ($#databases_chainmode >= 0){
    print LOG "# -C: Chain mode: map reads against these databases: @databases_chainmode\n";
}
if ($verbose){
    print LOG "# -c: Number of cores: $cores\n";
    print LOG "# -t: Matches\/readLength: $matchRatio\n";
    print LOG "# -v: logfile: $workDir/MGmapper.log\n";
}
if (defined($Getopt::Std::opt_S)){
    print LOG "# -S: Skip running cutadapt\n";
}
print LOG "# -d: Output directory: $workDir\n\n";

if (defined($Getopt::Std::opt_a)){
    my @w=split(/\,/,$Getopt::Std::opt_a);
    my $i=0;
    while (defined($w[$i])){
	if (! exists($databaseNames{$w[$i]})){
	    print STDERR "Assembled fasta files requested for an unexisting database referenced with: '$w[$i]'\n";
	    die;
	}
	$fastaOut{$databaseNames{$w[$i]}}=1;
	print LOG "fasta: $w[$i] $databaseNames{$w[$i]}\n" if ($verbose);
	$i++;
    }
}
if (defined($Getopt::Std::opt_x)){
    foreach my $key (keys %databaseNames){
	print "$key\t$databaseNames{$key}\n";
    }
    exit;
}
###############################################################################
#
# Main program start
#
###############################################################################
my $mappingSummary = "$workDir/MGmapper.summary";
my $fileStatSummary = "$workDir/filestat.tab";
my $thisDir=cwd();
$datestring = localtime();
if ($wwwMode){
    open(FILESTAT,">$fileStatSummary");
    open(TIME,">$workDir/time.tab");
    print TIME "#Job start\tJob end\n";
    print TIME "$datestring\t";
}

open(STAT,">$mappingSummary");
print STAT "## Local date and time $datestring - start main program\n";
print STAT "# command: $command\n";
print STAT "# working dir: $thisDir\n\n";

my $cmd;
if (! defined($Getopt::Std::opt_S)){
#    $prog_cutadapt .= " -f fastq -q $qualityCut -m $minLen -b ATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTATCTCGTATGCC -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGAGCATCTCGTATGC -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA";
    $prog_cutadapt .= " --quality-base=$QUALITY_BASE -f fastq -q $qualityCut -m $minLen";
    if (! -e "$parameterFile_adapter"){
	print LOG "can't open file $parameterFile_adapter\n";
	print LOG "Done!\n";
	die;
    }
    open(ADAPTER,"<$parameterFile_adapter") || die ("can't open file $parameterFile_adapter: $!");
    print LOG "Adapters removed by cutadapt:\n" if ($verbose);
    my $i=0;
    my $adapter;
    while (defined($_=<ADAPTER>)){
	if (/^\#/){
	    next;
	}
	chomp;
	if (m/(\S+)/){
	    $adapter = $1;
	}
	else{
	    next;
	}
	my $len=length($adapter);
	if ($len >1){
	    $i++;
	    print LOG "$i\t$adapter\n" if ($verbose);
	    $prog_cutadapt .= " -b $adapter";
	}
    }
    print LOG "\n" if ($verbose);
    close(ADAPTER);
}
my %readCount=();
my $readCountOpenFlag=0;
my $readCountFile = "$workDir/readCount.txt";


my $readsInCommonF = "$workDir/readsInCommon.F.fq";
my $readsInCommonR = "$workDir/readsInCommon.R.fq";
my $readsInCommonFgz = "$workDir/readsInCommon.F.fq.gz";
my $readsInCommonRgz = "$workDir/readsInCommon.R.fq.gz";

my $sumReadsBeforeF=0;
my $sumReadsBeforeR=0;
my $sumReadsBefore=0;

my $sumReadsAfterF=0;
my $sumReadsAfterR=0;
my $sumReadsAfter=0;

my $fastqFiles=0;
my $dummy='dummy';

if (! defined($Getopt::Std::opt_S)){
    if ($wwwMode){
	print FILESTAT "#Counter\tFileName\tReads in\tReads after cutadapt\n";
    }
}
else{
    if ($wwwMode){
	print FILESTAT "#Counter\tFileName\tReads in\n";
    }
    print STAT "Counter\tFileName\tReads\n";
}

my $runReadsInCommonFlag=1;
my $cleanFileFlag=0;
if (-e "$workDir/cleaned.nophiX.bam"){
    $cleanFileFlag=1;
    if (! exists($readCount{$cleanFile})){
	$okReads = bamCountReads($cleanFile);
    }
}

if ( ((-e $readsInCommonF) && (-e $readsInCommonR)) || ((-e $readsInCommonFgz) && (-e $readsInCommonRgz)) || ($cleanFileFlag)){
    $runReadsInCommonFlag=0;
    if (-e $readsInCommonFgz){
	$readsInCommonF=$readsInCommonFgz;
    }
    if (-e $readsInCommonRgz){
	$readsInCommonR=$readsInCommonRgz;
    }
}

if (-e $readCountFile){
    open(COUNT,"<$readCountFile");
    while (defined($_=<COUNT>)){
	chomp;
	my @w=split(/\t+/);
	my $file = $w[0];
	my $reads = $w[1];
	$readCount{$file}=$reads;
	print LOG "# Reading from $readCountFile: $_\n" if ($verbose);
    }
    close(COUNT);
}

my %readsBefore=();
my %readsAfter=();
my %cutadaptFile=();

my $sumReadsBefore=0;
my $sumReadsAfter=0;

my @fileHolder=();
my @fileHolderF=();
my @fileHolderR=();
my $cutF = "$workDir/cutadapt.all.F.fq";
my $cutR = "$workDir/cutadapt.all.R.fq";

if ( defined($Getopt::Std::opt_i) && defined($Getopt::Std::opt_j) ){
    push(@fileHolder,$inFileF);
    push(@fileHolder,$inFileR);
    push(@fileHolderF,$inFileF);
    push(@fileHolderR,$inFileR);
}
elsif ( defined($Getopt::Std::opt_I) && defined($Getopt::Std::opt_J) ){    
    while (defined($_=<LISTF>)){
	if (m/^\#/){
	    next;
	}
	chomp;
	my $len = length($_);
	if ($len < 1){next}
	my $file = $_;
	if (! -e $file){
	    print LOG "File not found; $file\n";
	    print STDERR "File not found; $file\nDone!\n";
	    die;
	}
	push(@fileHolder, $file);
	push(@fileHolderF, $file);
    }
    close(LISTF);

    while (defined($_=<LISTR>)){
	if (m/^\#/){
	    next;
	}
	chomp;
	my $len = length($_);
	if ($len < 1){next}
	my $file = $_;
	if (! -e $file){
	    print LOG "File not found; $file\n";
	    print STDERR "File not found; $file\nDone!\n";
	    die;
	}
	push(@fileHolder, $file);
	push(@fileHolderR, $file);
    }
    close(LISTR);
}

foreach my $file (@fileHolder){
    if (exists($readCount{$file})){
	$readsBefore{$file}=$readCount{$file};
    }
    else{
	open(COUNT,">>$readCountFile");
	$readsBefore{$file}=countReads($file);
	print COUNT "$file\t$readsBefore{$file}\n";
	close(COUNT);
    }
    $sumReadsBefore += $readsBefore{$file};
}


if (! defined($Getopt::Std::opt_S)){
    my $first=1;
    my @commands=();

    foreach my $file (@fileHolder){
	if ($first){
	    $datestring = localtime();	    
	    print STAT "## Local date and time $datestring - cutadapt start\n";
	    $first=0;
	}

	my @w=split(/\s+/,$command);
	my $command = cutadapt_SE($file);
	my @w=split(/\s+/,$command);
	$cutadaptFile{$file}=$w[-1];

	if (-e "${cutadaptFile{$file}}.gz"){
	    $cutadaptFile{$file}= "${cutadaptFile{$file}}.gz";
	}

	if ($command ne ''){
	    if ((! -e $cutadaptFile{$file}) || (-z $cutadaptFile{$file})){
		push(@commands,$command);
	    }
	}
    }
    
    foreach my $command (@commands){
	print LOG "$command\n" if ($verbose);
	$cm->start("$command");		
    }
    #
    # Wait for all jobs to finish
    #
    $cm->wait_all_children;
    
    foreach my $file (@fileHolder){
	$fastqFiles++;
	if (exists($readCount{$cutadaptFile{$file}})){
	    $readsAfter{$file}=$readCount{$cutadaptFile{$file}};
	}
	else{
	    open(COUNT,">>$readCountFile");
	    $readsAfter{$file}=countReads($cutadaptFile{$file});
	    print COUNT "$cutadaptFile{$file}\t$readsAfter{$file}\n";
	    close(COUNT);
	}
	
	$sumReadsAfter += $readsAfter{$file};
	
	if ($wwwMode){
	    my @FS=split(/\//,$file);
	    print FILESTAT "$fastqFiles\t$FS[-1]\t$readsBefore{$file}\t$readsAfter{$file}\n";	    
	}
	print STAT "$fastqFiles\t$file\t$readsBefore{$file}\t$readsAfter{$file}\n";
    }
}
else{
    foreach my $file (@fileHolder){
	$fastqFiles++;
	if ($wwwMode){
	    my @FS=split(/\//,$file);
	    print FILESTAT "$fastqFiles\t$FS[-1]\t$readsBefore{$file}\n";
	}
	print STAT "$fastqFiles\t$file\t$readsBefore{$file}\n";
    }
}


if (($#fileHolderF >= 0) && ($#fileHolderR >=0)){
    #
    # Forward reads in several files
    #
    if ($runReadsInCommonFlag){
	$datestring = localtime();
	print STAT "## Local date and time $datestring - start joining files before readsInCommon\n";
	my $cmd;

	if (-e "${cutF}.gz"){
	    $cutF="${cutF}.gz";
	}
	if (-e "${cutR}.gz"){
	    $cutR="${cutR}.gz";
	}

	if (! -e $cutF){
	    foreach my $file (@fileHolderF){
		if (-e $cutadaptFile{$file}){
		    if (($cutadaptFile{$file} =~/\.gz$/) || ($cutadaptFile{$file}=~/\.Z$/)){
			$cmd="gunzip -c $cutadaptFile{$file} >> $cutF"; 
		    }
		    else{
			$cmd="cat $cutadaptFile{$file} >> $cutF";
		    }
		}
		elsif ( ($file =~/\.gz$/) || ($file=~/\.Z$/) ){
		    $cmd = "gunzip -c $file >> $cutF";
		}
		else{
		    $cmd = "cat $file >> $cutF";
		}
		print LOG "# Doing: $cmd\n" if ($verbose);
		system("$cmd");
	    }
	}
	
	if (! -e $cutR){
	    foreach my $file (@fileHolderR){
		if (-e $cutadaptFile{$file}){
		    if (($cutadaptFile{$file} =~/\.gz$/) || ($cutadaptFile{$file}=~/\.Z$/)){
			$cmd="gunzip -c $cutadaptFile{$file} >> $cutR"; 
		    }
		    else{
			$cmd="cat $cutadaptFile{$file} >> $cutR"; 
		    }
		}
		elsif (($file =~/\.gz$/) || ($file=~/\.Z$/)){
		    $cmd = "gunzip -c $file >> $cutR";
		}
		else{
		    $cmd = "cat $file >> $cutR";
		}
		print LOG "# Doing: $cmd\n" if ($verbose);
		system("$cmd");
	    }
	}
   
	if ($cleanupCutadapt){
	    push(@rmFiles,$cutF);
	    push(@rmFiles,$cutR);
	}


	$datestring = localtime();
	print STAT "## Local date and time $datestring - start readsInCommon\n";
	$cmd = "$prog_readsInCommon -f $cutF -r $cutR -a $readsInCommonF -b $readsInCommonR";
	print LOG "# Doing: $cmd\n" if ($verbose);
	system("$cmd");
	
	push(@rmFiles,$readsInCommonF);
	push(@rmFiles,$readsInCommonR);
	if (-z $readsInCommonF){
	    print LOG "No reads were found in pairs with $prog_readsInCommon\n" if ($verbose);
	    print STDERR "No reads were found in pairs with $prog_readsInCommon\n";
	    print STDERR "\nDone!\n";
	    die;
	}
    }
    else{
	print STAT "Files exist - Skipping readsIncommon\n";
    }
    my $count;
    if (exists($readCount{$readsInCommonF})){
	$count= $readCount{$readsInCommonF};
	$okReads= 2 * $count;
    }
    else{
	if (-e "$readsInCommonF"){
	    $count = countReads($readsInCommonF);
	    open(COUNT,">>$readCountFile");
	    print COUNT "$readsInCommonF\t$count\n";
	    close(COUNT);	
	    $okReads= 2 * $count;
	}
	else{
	    if (exists($readCount{$cleanFile})){
		$okReads = $readCount{$cleanFile};
	    }
	    else{
		print LOG "Problem occured counting reads from readsInCommon or the cleaned bam file- I die now\n";
		print LOG "Done!\n";
		die;
	    }
	}
    }
    print STAT "\nTotal\tReads read\treads after cutadapt\treads after readsInCommon\n";
    print STAT "Total\t$sumReadsBefore\t$sumReadsAfter\t$okReads\n\n";	    
    if ($wwwMode){
	open(FILESTAT2,">$workDir/filestat2.tab");
	print FILESTAT2 "#Files\tReads read\tReads after cutadapt\tReads in common\n";
	print FILESTAT2 "$fastqFiles\t$sumReadsBefore\t$sumReadsAfter\t$okReads\n";	    
	close(FILESTAT2);
    }
}
else{
    print LOG "Empty fileholders - I die now\n";
    print LOG "Done!\n";
    die;
}
if (-e $cleanFile){
    print LOG "Calling mapping with parameters: $cleanFile,$dummy\n" if ($verbose);
    $lastMapNot=mapping($cleanFile,$dummy);
}
else{
    print LOG "Calling mapping with parameters: $readsInCommonF,$readsInCommonR\n" if ($verbose);
    $lastMapNot=mapping($readsInCommonF,$readsInCommonR);
}
############################################################################
#
# Make overall abundance statistics
#
############################################################################
open(ABUNDANCE,">$workDir/abundance.databases.txt");
open(DBCount,">$workDir/dbCount.tab");
print DBCount "#Database\tNo of sequence\tNo of nucleotides\n";

if ($wwwMode){
    open(FILESTAT3,">$workDir/filestat3.tab");
    print FILESTAT3 "#Mapping mode\tDatabase\tPercent mapped\tNo reads\tReads available\n";
    open(FILESTAT4,">$workDir/filestat4.tab");
    print FILESTAT4 "#Mapping mode\tDatabase\tPercent mapped\tNo reads\tReads available\n";
}
my $unmapped = $mapNotPhiXReads;

printf ABUNDANCE "Fullmode\tnotPhiX\t100.00\t%d\t%d\n",$mapNotPhiXReads,$okReads;
#if ($wwwMode){
    #printf FILESTAT3 "Fullmode\tnotPhiX\t100.00\t%d\t%d\n",$mapNotPhiXReads,$okReads;
    #printf FILESTAT4 "Fullmode\tnotPhiX\t100.00\t%d\t%d\n",$mapNotPhiXReads,$okReads;
#}
print LOG "mapNotPhiXReads = $mapNotPhiXReads used as reference to max possible reads that can map to a database\n" if ($verbose);

#foreach my $key (keys %abundance){
#
foreach my $key (@databases_all){
    my @tmp=split(/\./,$key);
    if (($tmp[1] ne 'uniq') && ($tmp[2] ne 'uniq')){
	my $perc = $abundance{$key}{reads}*100/$mapNotPhiXReads;
	if ($runningMode{$key} eq 'Fullmode'){
	    printf ABUNDANCE "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$key}{reads},$mapNotPhiXReads;
	    print DBCount "$key\t$countEntries{$key}{sequences}\t$countEntries{$key}{nuc}\n";
	    if ($wwwMode){
		printf FILESTAT3 "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$key}{reads},$mapNotPhiXReads;
	    }
	}
	else{
	    printf ABUNDANCE "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$key}{reads},$abundance{$key}{readsAvailable};
	    print DBCount "$key\t$countEntries{$key}{sequences}\t$countEntries{$key}{nuc}\n";
	    if ($wwwMode){
		printf FILESTAT3 "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$key}{reads},$abundance{$key}{readsAvailable};
	    }
	}
	if ($runningMode{$key} eq 'Chainmode'){
	    $unmapped -= $abundance{$key}{reads};
	}
    }
}

my $perc = $unmapped*100/$mapNotPhiXReads;
my $str="Unmapped";
printf ABUNDANCE "Chainmode\t%s\t%7.3f\t%d\t%d\n",$str,$perc,$unmapped,$mapNotPhiXReads;
close(ABUNDANCE);
if ($wwwMode){
    printf FILESTAT3 "Chainmode\t%s\t%7.3f\t%d\t%d\n",$str,$perc,$unmapped,$mapNotPhiXReads;
    close(FILESTAT3);
}
close(DBCount);
#
# Make overall statistics for uniq reads only
#
open(ABUNDANCE,">$workDir/abundance.databases.uniq.txt");
printf ABUNDANCE "Fullmode\tnotPhiX\t100.00\t%d\t%d\n",$mapNotPhiXReads,$okReads;

my $unmapped = $mapNotPhiXReads;
foreach my $key (@databases_all){
    my $db = $key . '.uniq';
    my $perc = $abundance{$db}{reads}*100/$mapNotPhiXReads;
    if ($runningMode{$key} eq 'Fullmode'){
	printf ABUNDANCE "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$db}{reads},$mapNotPhiXReads;
	if ($wwwMode){
	    printf FILESTAT4 "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$db}{reads},$mapNotPhiXReads;
	}
    }
    else{
	printf ABUNDANCE "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$db}{reads},$abundance{$db}{readsAvailable};
	if ($wwwMode){
	    printf FILESTAT4 "%s\t%s\t%7.3f\t%d\t%d\n",$runningMode{$key},$key,$perc,$abundance{$db}{reads},$abundance{$db}{readsAvailable};
	}
    }
    
    if ($runningMode{$key} eq 'Chainmode'){
	$unmapped -= $abundance{$db}{reads};
    }
}

my $perc = $unmapped*100/$mapNotPhiXReads;
my $str="Unmapped";
printf ABUNDANCE "Chainmode\t%s\t%7.3f\t%d\t%d\n",$str,$perc,$unmapped,$mapNotPhiXReads;
if ($wwwMode){
    printf FILESTAT4 "Chainmode\t%s\t%7.3f\t%d\t%d\n",$str,$perc,$unmapped,$mapNotPhiXReads;
    close(FILESTAT4);
}
close(ABUNDANCE);

open(DISTRIB,"| sort -k 1 -n > $workDir/insertSizeDistrib.txt");
my $insertSizeSum=0;
my $insertSizeCount=0;
foreach my $insertSize (keys %insertSizeDistrib){
    print DISTRIB "$insertSize\t$insertSizeDistrib{$insertSize}\n";
    $insertSizeSum += $insertSizeDistrib{$insertSize}*$insertSize;
    $insertSizeCount += $insertSizeDistrib{$insertSize};
}
close(DISTRIB);
if ($insertSizeCount > 0){
    $insertSizeAverage = $insertSizeSum/$insertSizeCount;
    my $nn=0;
    $insertSizeCount=0;
    foreach my $insertSize (keys %insertSizeDistrib){
	my $n =  $insertSizeAverage - $insertSize;
	$nn += $insertSizeDistrib{$insertSize} * ($n *$n);
	$insertSizeCount += $insertSizeDistrib{$insertSize};
    }
    if ($insertSizeCount > 1){
	$insertSizeStd= sqrt($nn/($insertSizeCount - 1));
    }
    printf LOG "Average insert size: %.2f\n",$insertSizeAverage if ($verbose);
}
    
#
# Join nucleotide statistics file with read abundance statistics file
#
$datestring = localtime();
print STAT "## Local date and time $datestring - Start making excel files\n";

my $xlsWriteUniq = "/tools/bin/python $prog_xlsWrite $workDir/$sampleName.uniq.xlsx $workDir/abundance.databases.uniq.txt";
my $xlsWrite = "/tools/bin/python $prog_xlsWrite $workDir/$sampleName.xlsx $workDir/abundance.databases.txt";
if ($insertSizeCount >0){
    $xlsWrite .= " $workDir/insertSizeDistrib.txt";
}
foreach my $db (@databases_all){
    my $stat = "$workDir/stat.$db.txt";
    if ((! -z $stat) && (-e $stat)){
	$xlsWrite .= " $stat";
    }

    my $statUniq = "$workDir/stat.$db.uniq.txt";
    if ((! -z $statUniq) && (-e $statUniq)){
	$xlsWriteUniq .= " $statUniq";
    }
}
if ($verbose){
    print LOG "$xlsWrite\n";
    print LOG "$xlsWriteUniq\n";
}
if (-e "$workDir/$sampleName.xlsx"){
    system("rm $workDir/$sampleName.xlsx");
}
if (-e "$workDir/$sampleName.uniq.xlsx"){
    system("rm $workDir/$sampleName.uniq.xlsx");
}
system("$xlsWrite");
system("$xlsWriteUniq");

#
# 1 fastq file with unmapped reads
#
if ($insertSizeCount > 0){
    printf STAT "Ave insert size: %.2f\n",$insertSizeAverage;
    printf STAT "Std insert size: %.2f\n\n",$insertSizeStd;
    if ($wwwMode){
	printf "\nAve insert size: %.2f\n",$insertSizeAverage;
	printf "Std insert size: %.2f\n\n",$insertSizeStd;
    }
}

if ($#databases_chainmode >=0){
    if (($wwwMode) && ($insertSizeCount > 0)){
        open(INSERTSIZE,">$workDir/InsertSize.tab");
	print INSERTSIZE "#Ave insert size\tStd insert size\n";
	printf INSERTSIZE"%.2f\t%.2f\n",$insertSizeAverage,$insertSizeStd;
        close(INSERTSIZE);
    }

    if (defined($Getopt::Std::opt_o)){
	#                                                                                                                                                                                
	# 1 interleaved fastq file with unmapped reads                                                                                                                                   
	#
	my $count;
	my $readsUnmapped = "$workDir/unmapped.fq";
	my $readsUnmappedgz = "$workDir/unmapped.fq.gz";
	if (! -e $readsUnmappedgz){
	    bamTofastq($lastMapNot,$readsUnmapped);
	    $count = `wc -l $readsUnmapped`;
	    $count /= 4;
	    chomp($count);
	    print STAT ("File with unmapped reads : $readsUnmapped.gz\n");
	    if ($count == 0){
		push(@rmFiles,$readsUnmapped);
	    }
	}
	else{
	    $count = `gunzip -c $readsUnmappedgz | wc -l`;
	    $count /= 4;
	    chomp($count);	
	}
	print STAT ("Number of unmapped reads : $count\n");
    }
    else{
	#                                                                                                                                                                                
	# 2 files with pair and mate reads respectively                                                                                                                                  
	#                                                                                                                                                                                
	my $unmappedF = "$workDir/unmapped.F.fq";
	my $unmappedR = "$workDir/unmapped.R.fq";
	
	my $unmappedFgz = "$workDir/unmapped.F.fq.gz";
	my $unmappedRgz = "$workDir/unmapped.R.fq.gz";
	
	my $count;
	if (! -e $unmappedFgz){
	    my $cmd = "$prog_samtools view -f64 $lastMapNot -b | $prog_bamtools convert -format fastq | cat > $unmappedF";
	    print LOG "# doing: $cmd\n" if ($verbose);
	    
	    system("$cmd");
	    $count = `wc -l $unmappedF`;
	    $count /= 4;
	    chomp($count);
	    print STAT ("File with unmapped reads in : $unmappedFgz\n");
	    print STAT ("Number of unmapped reads : $count\n");
	    if ($count == 0){
		push(@rmFiles,$unmappedF);
	    }
	}
	else{
	    $count = `gunzip -c $unmappedFgz | wc -l`;
	    $count /= 4;
	    chomp($count);
	    print STAT ("File with unmapped reads in : $unmappedFgz\n");
	    print STAT ("Number of unmapped reads : $count\n");
	}
	
	if (! -e $unmappedRgz){
	    $cmd = "$prog_samtools view -f128 $lastMapNot -b | $prog_bamtools convert -format fastq | cat > $unmappedR";
	    print LOG "# doing: $cmd\n" if ($verbose);
	    
	    system("$cmd");
	    $count = `wc -l $unmappedR`;
	    $count /= 4;
	    chomp($count);
	    print STAT ("File with unmapped reads in : $unmappedRgz\n");
	    print STAT ("Number of unmapped reads : $count\n");
	    if ($count == 0){
		push(@rmFiles,$unmappedR);
	    }
	}
	else{
	    $count = `gunzip -c $unmappedRgz | wc -l`;
	    $count /= 4;
	    chomp($count);
	    print STAT ("File with unmapped reads in : $unmappedRgz\n");
	    print STAT ("Number of unmapped reads : $count\n");
	}
    }
}
if ($wwwMode){
    my $maxShowLines=5;
    my $rc="$workDir/Chainmode.tab";
    my $rf="$workDir/Fullmode.tab";
    my $rcUniq="$workDir/Chainmode.uniq.tab";
    my $rfUniq="$workDir/Fullmode.uniq.tab";
    open(RC,">$rc");
    open(RF,">$rf");
    open(RCUNIQ,">$rcUniq");
    open(RFUNIQ,">$rfUniq");
    my %writeToFile=();
    $writeToFile{Chainmode}=0;
    $writeToFile{Fullmode}=0;
    $writeToFile{ChainmodeUniq}=0;
    $writeToFile{FullmodeUniq}=0;
    my $statFile;
    foreach my $key (@databases_all){	
	print LOG "Making sbippet files for $key: $rc and $rf: wwwmode=$wwwMode\n";
	$statFile = "$workDir/stat.$key.txt";
	if (($runningMode{$key} eq 'Chainmode') && (-e $statFile) && (! -z $statFile)){
	    print RC "#$key\n";
	    print RC "#Ref Seq\t% Nucleotides\tDepth\tSum of Nucleotides\tCoverage\tUniq positions\tSize\t%Reads\tNo Reads\tDescription\n";
	    open(A,"<$statFile");
	    my $i=0;
	    while (defined($_=<A>) && ($i<$maxShowLines)){
		print RC "$_";
		$i++;
		$writeToFile{Chainmode}=1;
	    }
	    close(A);
	    print RC "\n"; 
	}
	$statFile = "$workDir/stat.$key.uniq.txt";
	if (($runningMode{$key} eq 'Chainmode') && (-e $statFile) && (! -z $statFile)){
	    print RCUNIQ "#$key\n";
	    print RCUNIQ "#Ref Seq\t% Nucleotides\tDepth\tSum of Nucleotides\tCoverage\tUniq positions\tSize\t%Reads\tNo Reads\tDescription\n";
	    open(A,"<$statFile");
	    my $i=0;
	    while (defined($_=<A>) && ($i<$maxShowLines)){
		print RCUNIQ "$_";
		$i++;
		$writeToFile{ChainmodeUniq}=1;
	    }
	    close(A);
	    print RCUNIQ "\n"; 
	}
	$statFile = "$workDir/stat.$key.txt";
	if (($runningMode{$key} eq 'Fullmode') && (-e $statFile) && (! -z $statFile)){
	    print RF "#$key\n";
	    print RF "#Ref Seq\t% Nucleotides\tDepth\tSum of Nucleotides\tCoverage\tUniq positions\tSize\t%Reads\tNo Reads\tDescription\n";
	    open(A,"<$statFile");
	    my $i=0;
	    while (defined($_=<A>) && ($i<$maxShowLines)){
		print RF "$_";
		$i++;
		$writeToFile{Fullmode}=1;
	    }
	    close(A);
	    print RF "\n"; 
	}
	$statFile = "$workDir/stat.$key.uniq.txt";
	if (($runningMode{$key} eq 'Fullmode') && (-e $statFile) && (! -z $statFile)){
	    print RFUNIQ "#$key\n";
	    print RFUNIQ "#Ref Seq\t% Nucleotides\tDepth\tSum of Nucleotides\tCoverage\tUniq positions\tSize\t%Reads\tNo Reads\tDescription\n";
	    open(A,"<$statFile");
	    my $i=0;
	    while (defined($_=<A>) && ($i<$maxShowLines)){
		print RFUNIQ "$_";
		$i++;
		$writeToFile{FullmodeUniq}=1;
	    }
	    close(A);
	    print RFUNIQ "\n"; 
	}
    }
    close(RC);
    close(RF);
    close(RCUNIQ);
    close(RFUNIQ);
    if ($writeToFile{Chainmode} == 0){
	push(@rmFiles,$rc);
    }
    if ($writeToFile{ChainmodeUniq} == 0){
	push(@rmFiles,$rcUniq);
    }
    if ($writeToFile{Fullmode} == 0){
	push(@rmFiles,$rf);
    }
    if ($writeToFile{FullmodeUniq} == 0){
	push(@rmFiles,$rfUniq);
    }
}
if ($cleanup){
    my $cmd;
    if (! -e "$noclean"){
	foreach my $id (@rmFiles){
	    if (-e $id){
		$cmd = "rm $id";
		print LOG "Doing: $cmd\n" if ($verbose);

		system("$cmd");
	    }
	}
    }
}
#
# gzip all fastq files
#
my $cmd="$prog_pigz -p $cores $workDir/*.fq $workDir/*.cutadapt";
print LOG "Doing: $cmd\n" if ($verbose);

$datestring = localtime();
print STAT "## Local date and time $datestring - gzip files\n";

system("$cmd");

if ($wwwMode){
    close(FILESTAT);
    close(LOG);
    print STDERR "\nDone!\n";
    close(STAT);
    exit;
}

#
# Make a header filethat explains columns in $database.summary files
#
my $readme="$workDir/Readme";
open(README,">$readme");
print README "Explanation of columns in 'stat.databaseName.txt' files:\n";
print README "Column 1: Ref sequence name\n";
print README "Column 2: % of nucleotides that map to Ref sequence, compared to number of nucleotides that mapped to one paticular database\n";
print README "Column 3: Depth or X,  i.e. total number of nucleotides mapped to Ref sequence divided with length of Ref sequence\n";
print README "Column 4: Total number of nucleotides mapped to Ref sequence\n";
print README "Column 5: Coverage, fraction of Ref sequence that is covered\n";
print README "Column 6: Number of unique positions in Ref sequence that is covered\n";
print README "Column 7: Size of Ref sequence\n";
print README "Column 8: % of reads that map to Ref sequence, compared to number of reads that mapped to one paticular database\n";
print README "Column 9: Total number of reads mapped to Ref sequence\n";
print README "Column 10: Description of Ref sequence\n";
print README "\n";
print README "A special file is 'abundance.databases.txt' which shows read aboundance for all databases\n";
print README "Percentages are calculated based on reads available after removing those that map to PhiX \(notPhiX\). notPhiX count is set\n";
print README "Unmapped reads is calculated only on runs in Chainmode.\n";
print README "\n";
if (defined($Getopt::Std::opt_R)){
    print README "A readcount matrix file: readCount.database.matrix\n";
    print README "A matrix showing ALL entries that are present in a database, also entries where not reads mapped to that entry\n";
    print README "Column 1: Name of reference sequence\n";
    print README "Column 2: Number of reads that mapped to that reference sequence\n";
    print README "Column 3: Length of reference sequence\n";
    print README "Column 4: Description of reference sequence\n";
    print README "\n";
}
print README "Files named *uniq* only include uniquely mapped reads i.e. sam tags AS > XS meaning there is only one best hit. Only a criteria for one of the paired reads\n";
print README "AS is the alignment-score and AX is the second best alignment score for one paticular read mapping to a reference sequence database\n";
print README "\n";
print README "All files are combined in two excel workbooks named sample.xlsx and sample.uniq.xlsx which include headers and figures";

close(README);
$datestring = localtime();
print STAT "## Local date and time $datestring - done\n";
close(STAT);
exit(0);
###############################################################################
#
# Main program end
#
###############################################################################
sub single_bwa{
    my ($bam, $list, $db) = @_;
    
    my $remapDir="${db}_remapped";
    my $newDir="$workDir/$remapDir";
    $datestring = localtime();

    if (! -d $newDir){
	system("mkdir -p $newDir");
    }
    my $logfile = "$workDir/remapping_${db}.log";
    if (-e $logfile){
	open(LOGFILE,">>$logfile");
	print LOGFILE "## Local date and time $datestring - appending to logfile - db= $db\n";
    }
    else{
	open(LOGFILE,">$logfile");
	print LOGFILE "## Local date and time $datestring - start making scripts to remap single entries -  db= $db\n";
    }
    #
    # $bam is a sorted mapTo bam file
    #

    # Read the lookup file with following format: 3 columns being 'newName' 'newDir' 'targetName'
    #
    # ResFinder_1   singleFasta/1   aac(2')-Ia_1_L06156
    # ResFinder_2   singleFasta/1   aac(2')-Ib_1_U41471
    #
    # ...
    #
    # ResFinder_839 singleFasta/9   blaOXA-327_1_KF203101
    my $lookupPath='';
    my %lookup=();
    my @pileupToDepthCommands=();
    my @pileupToFastaCommands=();
    my @bamToPileupCommands=();

    if (-e "$idx{lookup}{$db}"){
	my $pos = rindex($idx{lookup}{$db},'/');
	$lookupPath=substr($idx{lookup}{$db},0,$pos);
    }
    else{
	print LOGFILE "Could not find database lookup file $idx{lookup}{$db}\n";
	close(LOGFILE);
	return(1);
    }
    my $singleIndexPath='';    
    my $targetName;

    open(LOOKUP,"cat $idx{lookup}{$db} |");
    while (defined($_=<LOOKUP>)){
	chomp;
	my @w=split(/\s+/);
	$targetName=$w[2];
	$lookup{newName}{$targetName}=$w[0];
	$lookup{index}{$targetName}="$lookupPath/$w[1]/$w[0]";
	$singleIndexPath="$lookupPath/$w[1]";
    }
    close(LOOKUP);

    # Get list of target names and map all reads from bam-file against each target name.
    # Each target name exists as bwa indexed files in directory for the database $db
    
    
    #
    # The input bam file is sorted and Interleaved bam file or fastq files where pairs are ordered in F.fq and R.fa
    # So do something to the bam file !!! puuuh
    #
    my $newBam= "$newDir/$db.SortByName.bam";
    my $newBamName= "$newDir/$db.SortByName";
    my %ALscores=();
    my $targetBam;

    if (! -e $newBam){
	my $cmd="$prog_samtools sort -n $bam $newBamName";
	print LOG ("# $cmd\n") if ($verbose);
	system("$cmd");
    }

    open(INBAM,"$prog_samtools view $newBam |");
    print LOG "# in single_bwa: reading alignment scores in $newBam\n" if ($verbose);
    #
    # read alignment scores from the bam-file and store them in hash
    #
    # Also determine the insertsize which will be used when remapping against reference sequences, where
    # the name is specified in $list
    #
    my $pairs=0;
    my $insert_size_sum=0;
    while (defined($_=<INBAM>)){
	my @w=split(/\s+/);
	my $read=$w[0];
	my $flag=$w[1];
	my $insert_size=abs($w[8]);

	if ($flag & 64){
	    #
	    # first read in the pair
	    #
	    $pairs++;
	    $insert_size_sum += $insert_size;
	}
	    
	my $as = substr($w[13],5);
	$ALscores{$read}{count}=0;
	    
	#
	# if flag <128 the the read is the first read in the pair otherwise second in pair
	#
	if ($flag & 64){ 
	    $ALscores{$read}{First}=$as;
	}
	else{
	    $ALscores{$read}{Second}=$as;
	}	    
    }
    close(INBAM);

    my $is_mean = $insert_size_sum/$pairs;
    my $is_std = $is_mean*0.1;
    my $is_max = $is_mean + 4*$is_std;
    my $is_min = $is_mean - 4*$is_std;
    
    print LOGFILE "# Using following parameters to call bwa with insert size information: mean, std=mean*0.1, max and min 4 std from mean\n";
    print LOGFILE "# mean: $is_mean\n";
    print LOGFILE "# std: $is_std\n";
    print LOGFILE "# max: $is_max\n";
    print LOGFILE "# min: $is_min\n";


    my @bwaMem=();
    open(LIST,"<$list");

    while (defined($_=<LIST>)){
	chomp;
	my $targetName = $_;
	my $index = $lookup{index}{$targetName};
	my $targetBam= "$newDir/$lookup{newName}{$targetName}.bam";

	#
	# Make 
	#
	if (! -e $targetBam){
	    $cmd = "$prog_bamtools convert -format fastq -in  $newBam | $prog_bwa mem -t $cores -I $is_mean,$is_std,$is_max,$is_min -M -p $index -| $prog_samtools view -h -F256 -f2 -Sb -o $targetBam -";
	    print LOG "# Doing: $cmd\n" if ($verbose);
	    $datestring = localtime();
	    print LOGFILE "## Local date and time $datestring - start making $targetBam  db= $db\n";

	    system($cmd);
	    $datestring = localtime();
	    print LOGFILE "## Local date and time $datestring - done making $targetBam  db= $db\n";
	}
    }
    close(LIST);

    open(LIST,"<$list");
    while (defined($_=<LIST>)){
	chomp;
	my $targetName = $_;
	my $index = $lookup{index}{$targetName};
	my $bamIn= "$newDir/$lookup{newName}{$targetName}.bam";
	my $bamOut= "$newDir/$lookup{newName}{$targetName}.1.bam";

	if (! -e $bamIn){
	    next;
	}
	open(INBAM,"$prog_samtools view -h $bamIn |");
	open(OUTBAM,"| $prog_samtools view -h -Sb - | cat > $bamOut");
	
	my $reads=0;
	if (! -e $bamOut){
	    my $header=1;
	    while (defined($_=<INBAM>)){
		if ((m/^\@/) && ($header)){
		    print OUTBAM "$_";
		    next;
		}
		else{
		    $header=0;
		}
		my $read_first=0;
		my $read_second=0;
		
		my @w=split(/\s+/);
		my $read=$w[0];
		my $flag=$w[1];
		my $as= substr($w[13],5);


		my $read_line=$_;
		if ($flag & 64){
		    if ($as >= $ALscores{$read}{First}){
			$read_first=1;
		    }
		}
		else{
		    if ($as >= $ALscores{$read}{Second}){
			$read_second=1;
		    }
		}
		
		#
		# Reading the paired read
		#
		$_=<INBAM>;	    
		@w=split(/\s+/);
		$read=$w[0];
		$flag=$w[1];
		$as = substr($w[13],5);
		my $mate_line=$_;

		if ($flag & 64){
		    if ($as >= $ALscores{$read}{First}){
			$read_first=1;
		    }
		}
		else{
		    if ($as >= $ALscores{$read}{Second}){
			$read_second=1;
		    }
		}
		
		if (($read_first == 1) && ($read_second==1)){
		    print OUTBAM "$read_line";
		    print OUTBAM "$mate_line";
		    $reads +=2;
		    $ALscores{$read}{count}++;
		}
	    }
	    close(INBAM);
	    close(OUTBAM);
	}

	if ($reads != 0){
	    my $prefix = "$bamOut" . '.sort';
	    my $sortBam= "${prefix}.bam";
	    if (! -e $sortBam){
		print LOG "# Calling sortBam - parameters: $targetBam $prefix\n" if ($verbose);
		sortBam($bamOut, $prefix);
	    }
	    
	    my $dummy=0;
	    my $pileup="$newDir/$lookup{newName}{$targetName}.pileup";
	    my $summary="$newDir/$lookup{newName}{$targetName}.summary";
	    if (! -e $pileup){
		$cmd=bamToPileupCmd($sortBam,$dummy,$dummy,$pileup,$index);
		print LOG "Calling bamToPileupCmd: $sortBam, $dummy, $dummy, $pileup, $index\n" if ($verbose);
		if ($cmd ne ''){
		    push(@bamToPileupCommands,$cmd);
		}
	    }
	    

	    my $ann= $index . '.ann';	    
	    if (! -e $summary){
		print LOG "Calling pileupToDepthCmd: $pileup, $dummy, $dummy, $summary, $ann\n" if ($verbose);
		$cmd = pileupToDepthCmd($pileup,$dummy,$dummy,$summary,$ann);	    
		if ($cmd ne ''){
		    push(@pileupToDepthCommands,$cmd);
		}
	    }

	    my $fsa = "$newDir/$lookup{newName}{$targetName}.fsa";
	    $ann = "$lookup{index}{$targetName}" . '.ann';
	    my $refSeqLen = `tail -1 $ann| cut -f2 -d\' \'`;
	    chomp($refSeqLen);

	    if (! -e $fsa){
		$cmd = pileupToFastaCmd($pileup,$fsa,$targetName,$refSeqLen);
		push(@pileupToFastaCommands,$cmd);
	    }
	}
    }
    close(LIST);
    

    $datestring = localtime();
    if ($#bamToPileupCommands >=0){
	print LOGFILE "## Local date and time $datestring - start mapping against single entries for db= $db , bamToPileup\n";
    }
    else{
	print LOGFILE "## Local date and time $datestring - bamToPileupArray is empty\n";
    }
    foreach my $command (@bamToPileupCommands){
	print LOGFILE "'$command'\n" if ($verbose);
	$cm->start("$command");
    }
    #
    # Wait for all jobs to finish
    #
    $cm->wait_all_children;
    

    $datestring = localtime();
    if ($#pileupToDepthCommands >=0){
	print LOGFILE "## Local date and time $datestring - start mapping against single entries for db= $db , pileupToDepth\n";
    }
    else{
	print LOGFILE "## Local date and time $datestring - pileupToDepthArray  is empty\n";
    }
    foreach my $command (@pileupToDepthCommands){
	print LOGFILE "'$command'\n" if ($verbose);
	$cm->start("$command");
    }
    #
    # Wait for all jobs to finish
    #
    $cm->wait_all_children;

#    if (exists($fastaOut{$db})){
    $datestring = localtime();
    if ($#pileupToFastaCommands >=0){
	print STAT "## Local date and time $datestring - start mapping against single entries for db= $db , pileupToFasta\n";
    }
    else{
	print LOGFILE "## Local date and time $datestring - pileupToFastaArray  is empty\n";
    }
    my $makeFastaFiles=0;
    foreach my $command (@pileupToFastaCommands){
	print LOGFILE "'$command'\n" if ($verbose);
	$cm->start("$command");
	$makeFastaFiles=1;
    }
    #
    # Wait for all jobs to finish
    #
    $cm->wait_all_children;
    
    if ($makeFastaFiles){
	my $fasta = "$workDir/stat.$db.uniq_remapped.fasta";
	system("cat $newDir/*.fsa > $fasta");
	
	my $stat="$workDir/stat.$db.uniq_remapped.txt";
	open(SUMMARY,"cat $newDir/*.summary |");
	my $nuc_sum=0;
	my $reads_sum=0;
	while (defined($_=<SUMMARY>)){
	    my @w=split(/\s+/);
	    $nuc_sum+=$w[3];
	    $reads_sum+=$w[8];	
	}
	close(SUMMARY);
	open(SUMMARY,"cat $newDir/*.summary | sort -k 4 -nr |");
	open(STATISTICS,">$stat");
	my @nameList=();
	while (defined($_=<SUMMARY>)){
	    my @w=split(/\t+/);
	    my @tmp=split(/\s+/,$w[9]);
	    my $name=$tmp[0];
	    my $j=1;
	    my $desc='';
	    while (defined($tmp[$j])){
		$desc .= "$tmp[$j] ";
		$j++;
	    }
	    my $perc_nuc=$w[3]*100/$nuc_sum;
	    my $perc_reads=$w[8]*100/$reads_sum;
	    my $depth=$w[2];
	    my $nucleotides=$w[3];
	    my $coverage=$w[4];
	    my $uniq_positions=$w[5];
	    my $seq_len=$w[6];
	    my $reads=$w[8];
	    printf STATISTICS "%-25s\t%7.4f\t%7.4f\t%d\t%7.4f\t%d\t%d\t%7.4f\t%d\t%s\n",$name,$perc_nuc,$depth,$nucleotides,$coverage,$uniq_positions,$seq_len,$perc_nuc,$reads,$desc;
	}
	close(SUMMARY);
	close(STATISTICS);
	
	if (-e  $idx{fastauniq}{$db}){
	    if (! -z $stat){
		my $redundant="$workDir/redundant.$db.uniq_remapped.txt";
		my @list =();
		open(TMP,"<$stat");
		while (defined($_=<TMP>)){
		    chomp;
		    my @w=split(/\s+/);
		    push(@list,$w[0]);
		}
		close(TMP);
		foreach my $id (@list){
		    my $cmd = "grep '^$id' $idx{fastauniq}{$db} >> $redundant";
		    print LOGFILE "# Doing: $cmd\n" if ($verbose);
		    system "$cmd";
		}
	    }
	    print LOGFILE "# fastauniq file with redundant entries removed from database $db\n" if ($verbose);
	}
    }
#    my $cmd = "rm -r $newDir";
#    print LOG "# Doing $cmd\n" if ($verbose);
#   system "$cmd";

    #
    # Now identify all reads that have not yet been used
    #
    my $unused="$workDir/unused.$db.bam";
    open(BAMOUT,"| samtools view -h -Sb -o $unused - ");
    open(BAM,"$prog_samtools view -h $newBam |");
    my $header=1;
    while (defined($_=<BAM>)){
	if ((m/^SQ/) && ($header)){
	    print BAMOUT "$_";
	    next;
	}
	$header=0;
	my @w=split(/\s+/);
	my $read=$w[0];

	if ($ALscores{$read}{count} == 0){
	    print BAMOUT "$_";
	}
    }
    close(BAM);
    close(BAMOUT);
   return(0);
}

sub cutadapt_SE{
    my ($fastq) = @_;
    #
    # Identify input file for forward reads 
    #
    my $pos1 = rindex($fastq,'/') +1;
    my $name = substr $fastq , $pos1;

#    my $pos2 = index($name,'.');    
#    my $commonName = substr $name , 0, $pos2;
    

    #
    # remove adapter, trimm for forward reads
    #
    my $cutAdaptFile = "$workDir/$name" . '.cutadapt';
    my $gz = $cutAdaptFile . '.gz';
    if (-e $gz){
	$cutAdaptFile = $gz;
    }
    
    if ($cleanupCutadapt){
	push(@rmFiles,$cutAdaptFile);
    }
    my $cmd = "$prog_cutadapt $fastq > $cutAdaptFile";

    return("$cmd");
}

sub mergeSummaryFiles{
    my ($in,$out,$count) = @_;

    my $cmd;
    if ($count == 1){
        $cmd="mv $in $out";
	print LOG "# In mergeSummeryFiles: $cmd\n" if ($verbose);
        system("$cmd");
	return;
    }
    else{
        #
        # % nucleotides and reads need to be corrected
        #
	my $nucSum=0;
	my $readSum=0;
        $cmd="cat $in | cut -f4,9 |";
	print LOG "# mergeSummaryFiles: $cmd\n" if ($verbose);
        open(TMP,"$cmd");
        while (defined($_=<TMP>)){
            chomp;
	    my @w=split(' ');
            $nucSum += $w[0];
            $readSum += $w[1];
        }
        close(TMP);


	open(TMP,"<$in");
        open(TMPOUT,"| sort -k 4 -nr > $out");
        while (defined($_=<TMP>)){
            chomp;
            my @w=split(/\t+/);
	    my $name = $w[0];
	    my $nucleotides=$w[3];
	    my $readCount = $w[8];
            my $nucPerc=$nucleotides*100/$nucSum;
            my $depth = $nucleotides/$w[6];
            my $coverage = $w[5]/$w[6];
	    my $readPerc = $readCount*100/$readSum;
	    my $desc=$w[9];
            printf TMPOUT "%s\t%7.3f\t%7.3f\t%d\t%6.4f\t%d\t%d\t%7.4f\t%d\t%s\n",$name,$nucPerc,$depth,$nucleotides,$coverage,$w[5],$w[6],$readPerc,$readCount,$desc;
        }
        close(TMPOUT);
        close(TMP);
        $cmd = "rm $in";
	print LOG "# Doing: $cmd\n" if ($verbose);
	system("$cmd");
    }
    return;
}


sub mapping{
    my ($readF, $readR) = @_;
    #
    # bwa against phiX
    #
    $cleanFile = "$workDir/cleaned.nophiX.bam";

    my $sortBamFile;
    my $summaryFile;

    if (($redo==1) || (! -e $cleanFile) || (-z $cleanFile)){
	print LOG "Calling cleanData with parameters: $readF , $readR , $cleanFile\n" if ($verbose);
	$datestring = localtime();
	print STAT "## Local date and time $datestring - start mapping against phiX\n";
	cleanData($readF, $readR, $cleanFile);
    }
    else{
	print LOG "Skipping 'cleanData' as file exists: $cleanFile\n" if ($verbose);
    }
    
#    push(@rmFiles,$cleanFile);

    my $reads = bamCountReads($cleanFile);
    print STAT "Reads not mapped to db= phiX:\t$reads\n\n";

    $mapNotPhiXReads = $reads;
    $mapNotReads = $mapNotPhiXReads;
    if (! defined($Getopt::Std::opt_S)){    
	if ($cleanup){
	    if ( ((-e $readF) && (-e $readR)) && (! -e $noclean)){
		if ($readF ne $cleanFile){
		    push(@rmFiles,$readF);
		    push(@rmFiles,$readR);
		}
	    }
	}
    }
    #
    # map reads in Full mode against databases
    #
    my $size_databases_fullmode = $#databases_fullmode;
    if ($size_databases_fullmode >=0){
	print STAT "\nFull mode mapping starts here:\n";
	print STAT "------------------------------\n\n";
    }
    
    foreach my $db (@databases_fullmode){
	$datestring = localtime();
	print STAT "## Local date and time $datestring - start mapping against $db\n";

	#
	# print to summary file
	#
	my @w = split(/\//, $cleanFile);

	my $lastBamFile=$w[-1];
	if ($w[-1] =~ m/(\S+)\.bam\.(\d+)/){
	    $lastBamFile = "$1" . '.bam';
	}
	print STAT "Unmapped reads from U= $lastBamFile\n";
	print STAT "Mapping reads from U against db= $db:\n";

	#
	# Start mapping
	#
	my $ra_mapToList;
	my $ra_mapToListUniq;
	print LOG "calling mapToDatabases with parameters: '$cleanFile' '$db'\n" if ($verbose);
	my ($mapNot, $mapTo, $ra_mapToList, $mapToUniq, $ra_mapToListUniq) = mapToDatabases($cleanFile,$db);

	if ($verbose){
	    print LOG "after mapToDatabases fullmode\n";
	    print LOG "mapNot: '$mapNot'\n";
	    print LOG "mapTo: '$mapTo'\n";
	}
	my $ii=0;
	while (defined($ra_mapToList->[$ii])){
	    print LOG "mapToList: $ii $ra_mapToList->[$ii]\n" if ($verbose);
	    $ii++;
	}
	print LOG "mapToUniq: $mapToUniq\n" if ($verbose);
	$ii=0;
	while (defined($ra_mapToListUniq->[$ii])){
	    print LOG "mapToListUniq: $ii $ra_mapToListUniq->[$ii]\n" if ($verbose);
	    $ii++;
	}

	#
	# Count reads in mapTo and mapNot
	#
	my $reads = bamCountReads($mapTo);
	$abundance{$db}{reads}=$reads;

	my $readsUniq = bamCountReads($mapToUniq);
	my $uniqdb=$db . ".uniq";
	$abundance{$uniqdb}{reads}=$readsUniq;
	print LOG "uniq reads in $db: $readsUniq\n";
	print LOG "uniq reads $uniqdb in abundance hash: $readsUniq\n";

	printf STAT ("Map to db:\t%d\n",$reads);
	$mapNotReads = ($mapNotPhiXReads - $reads);

	if ($mapNot ne 'dummy'){
	    printf STAT ("Dont map to db:\t%d\t\(%s\)\n\n",$mapNotReads,$mapNot);
	}
	else{
	    print STAT ("Dont map to db:\t$mapNotReads\t\(bam not available\)\n\n");
	}


	#######################################################################################
	#
	# Make readCount matrices
	#
	#######################################################################################
	if (defined($Getopt::Std::opt_R)){
	    my $readCountMatrix ="$workDir/readCount.$db.matrix";
	    my $readCountMatrixUnsorted  ="$workDir/readCount.$db.matrix.Unsorted";
	    if ((! -e "$readCountMatrix") && ($reads > 0)){
		my $i=0;
		while (defined($ra_mapToList->[$i])){
		    my $j=$i+1;
		    readCountAnalysis($ra_mapToList->[$i], $j, $db, $readCountMatrixUnsorted);    
		    $i++;
		}
		#
		# sort the readCount matrix according to fastaEntry name
		#
		
		$cmd = "sort -k 1 $readCountMatrixUnsorted > $readCountMatrix";
		print LOG "# Doing: $cmd\n" if ($verbose);
		
		system($cmd);
		system("rm $readCountMatrixUnsorted");
	    }
	    
	    my $readCountMatrixUniq ="$workDir/readCount.$db.uniq.matrix";
	    my $readCountMatrixUnsortedUniq  ="$workDir/readCount.$db.uniq.matrix.Unsorted";
	    
	    if ((! -e "$readCountMatrixUniq") && ($readsUniq > 0)){
		my $i=0;
		while (defined($ra_mapToListUniq->[$i])){
		    my $j=$i+1;
		    readCountAnalysis($ra_mapToListUniq->[$i], $j, $db, $readCountMatrixUnsortedUniq);    
		    $i++;
		}
		#
		# sort the readCount matrix according to fastaEntry name
		#
		
		$cmd = "sort -k 1 $readCountMatrixUnsortedUniq > $readCountMatrixUniq";
		print LOG "# Doing: $cmd\n" if ($verbose);
		
		system($cmd);
		system("rm $readCountMatrixUnsortedUniq");
	    }
	}
	#####################################################################################
	#
	# calc coverage, depth ... using mpileup
	#
	#####################################################################################

	my $summaryFileMain ="$workDir/stat.$db" . ".txt";
	my $fsaMain ="$workDir/contig.$db.fasta";
	if ((! -e $summaryFileMain) || (! -e $fsaMain)){
	    my @summaryFiles=();
	    my @Fsa=();
	    my $i=0;
	    while (defined($ra_mapToList->[$i])){
		my $j=$i+1;
		my $summaryFile ="$workDir/stat.$db" . ".txt.$j";
		push(@summaryFiles,$summaryFile);

		my $pileup ="$workDir/pileup.$db.$j";
		bamToPileup($ra_mapToList->[$i], $j, $db, $pileup);
		push(@rmFiles,$pileup);

		if ((! -z $pileup) && (! -e $summaryFile)){
		    pileupToDepth($pileup, $j, $db, $summaryFile);
		}
		$i++;
	    }
	    #
	    # make assembled fasta files if requested - option -a
	    #
	    if (exists($fastaOut{$db})){
		my $cmd='';
		if ($#Fsa == 0){
		    $cmd = "mv $Fsa[0] $fsaMain";
		}
		elsif ($#Fsa >0){
		    my $i=0;
		    $cmd = 'cat ';
		    while (defined($Fsa[$i])){
			$cmd .= "$Fsa[$i] ";
			$i++;
		    }
		    $cmd .= " > $fsaMain";
		}
		print LOG "Doing: $cmd\n" if ($verbose);
		system("$cmd");

		push(@rmFiles,@Fsa);
	    }
	    #
	    # cat all summaryfiles for database db
	    #
	    my $cmd = "cat ";
            my $tmpfile = $summaryFileMain . '.tmp';
            my $count=0;

	    foreach my $id (@summaryFiles){
		$count++;
		chomp($id);
		$cmd .= "$id ";
	    }
	    $cmd .= "| sort -k 4 -nr > $tmpfile";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system($cmd);
	    mergeSummaryFiles($tmpfile,$summaryFileMain,$count);	        
            if (-z $summaryFileMain){
                print LOG "$summaryFileMain is empty - I remove it\n";
                system("rm $summaryFileMain");
            }

	    if ($cleanup){
		foreach my $id (@summaryFiles){
		    if (-e $id){
			system("rm $id");
		    }
		}
	    }
	}

	#
	# make summary file for uniq mapped reads
	# 
	my $summaryFileMainUniq ="$workDir/stat.$db.uniq" . ".txt";

	if (! -e $summaryFileMainUniq){
	    my @summaryFiles=();
	    my @Fsa=();
	    my $i=0;
	    while (defined($ra_mapToListUniq->[$i])){
		my $j=$i+1;
		my $summaryFile ="$workDir/stat.$db.uniq.txt.$j";
		push(@summaryFiles,$summaryFile);
		my $pileup ="$workDir/pileup.$db.uniq.$j";
		bamToPileup($ra_mapToListUniq->[$i], $j, $db, $pileup);
		push(@rmFiles,$pileup);

		if ((! -z $pileup) && (! -e $summaryFile)){
		    pileupToDepth($pileup, $j, $db, $summaryFile);
		}
		$i++;
	    }
	    #
	    # cat all summaryfiles for database db
	    #
            my $tmpfile = $summaryFileMainUniq . '.tmp';
            my $count=0;
	    my $cmd = "cat ";
	    foreach my $id (@summaryFiles){
		$count++;
		chomp($id);
		$cmd .= "$id ";
	    }
	    $cmd .= "| sort -k 4 -nr > $tmpfile";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system($cmd);
	    mergeSummaryFiles($tmpfile,$summaryFileMainUniq,$count);	        
            if (-z $summaryFileMainUniq){
                print LOG "$summaryFileMainUniq is empty - I remove it\n";
                system("rm $summaryFileMainUniq");
            }

	    if ($cleanup){
		foreach my $id (@summaryFiles){
		    if (-e $id){
			system("rm $id");
		    }
		}
	    }
	}
	
	$abundance{$db}{summary}=$summaryFileMain;
	$abundance{$uniqdb}{summary}=$summaryFileMainUniq;

	$lastMapNot = $mapNot;
	if ($cleanup){
	    push(@rmFiles,$mapNot);
	}
	
	if ($#rmMapNotList > 0){
	    $cmd="rm $rmMapNotList[0]";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system("$cmd");
	    shift(@rmMapNotList);
	}
    }

    ###########################################
    #
    # map reads in Chain mode against databases
    #
    ###########################################
    my $startFile = $cleanFile;
    my $size_databases_chainmode = $#databases_chainmode;
    if ($size_databases_chainmode >=0){
	print STAT "\nChain mode mapping starts here:\n";
	print STAT "-------------------------------\n\n";
	$mapNotReads = $mapNotPhiXReads;
    }
    foreach my $db (@databases_chainmode){
	$datestring = localtime();
	print STAT "## Local date and time $datestring - start mapping against $db\n";
	#
	# print to summary file
	#
	my @w = split(/\//, $startFile);
	my $lastBamFile=$w[-1];
	if ($w[-1] =~ m/(\S+)\.bam\.(\d+)/){
	    $lastBamFile = "$1" . '.bam';
	}
	print STAT "Unmapped reads from U= $lastBamFile\n";
	print STAT "Mapping reads from U against db= $db:\n";

	#
	# Start mapping
	#
	my $ra_mapToList;
	my $ra_mapToListUniq;
	print LOG "calling mapToDatabases with parameters: '$startFile' '$db'\n" if ($verbose);

	my ($mapNot, $mapTo, $ra_mapToList, $mapToUniq, $ra_mapToListUniq) = mapToDatabases($startFile,$db);

	if ($verbose){
	    print LOG "after mapToDatabases chainmode\n";
	    print LOG "mapNot: '$mapNot'\n";
	    print LOG "mapTo: '$mapTo'\n";
	}
	my $ii=0;
	while (defined($ra_mapToList->[$ii])){
	    print LOG "mapToList: $ii $ra_mapToList->[$ii]\n" if ($verbose);
	    $ii++;
	}
	print LOG "mapToUniq: $mapToUniq\n" if ($verbose);
	$ii=0;
	while (defined($ra_mapToListUniq->[$ii])){
	    print LOG "mapToListUniq: $ii $ra_mapToListUniq->[$ii]\n" if ($verbose);
	    $ii++;
	}
	#
	# Next start file are those reads that did not map
	#
	$startFile=$mapNot;

	#
	# Count inferred insertSize  
	#
	print LOG "Storing insert size distribution for file: $mapTo\n" if ($verbose);

	$cmd = "$prog_samtools view $mapTo | cut -f9 |";
	open(IN,"$cmd");
	while (defined($_=<IN>)){
	    chomp;
	    if ($_>=0){
		$insertSizeDistrib{$_}++;
	    }
	}
	close(IN);
	#
	# Count reads in mapTo and mapNot
	#
	my $reads = bamCountReads($mapTo);
	my $readsUniq = bamCountReads($mapToUniq);
	my $uniqdb= $db . '.uniq';

	$abundance{$db}{readsAvailable}=$mapNotReads;
	$abundance{$uniqdb}{readsAvailable}=$mapNotReads;
	$mapNotReads -= $reads;

	printf STAT ("Map to db:\t%d\n",$reads);
	$abundance{$db}{reads}=$reads;
	$abundance{$uniqdb}{reads}=$readsUniq;

	print LOG "uniq reads in $db: $readsUniq\n" if ($verbose);

	if ($mapNot ne 'dummy'){
	    printf STAT ("Dont map to db:\t%d\t\(%s\)\n\n",$mapNotReads,$mapNot);
	}
	else{
	    print STAT ("Dont map to db:\t$mapNotReads\t\(bam not available\)\n\n");
	}

	#######################################################################################
	#
	# Make readCount matrices
	#
	#######################################################################################
	if (defined($Getopt::Std::opt_R)){
	    my $readCountMatrix ="$workDir/readCount.$db.matrix";
	    my $readCountMatrixUnsorted  ="readCount.$db.matrix.Unsorted";
	    if ((! -e "$readCountMatrix") && ($reads > 0)){
		my $i=0;
		while (defined($ra_mapToList->[$i])){
		    my $j=$i+1;
		    readCountAnalysis($ra_mapToList->[$i], $j, $db, $readCountMatrixUnsorted);	    
		    $i++;
		}
		#
		# sort the readCount matrix according to fastaEntry name
		#
		$cmd = "sort -k 1 $readCountMatrixUnsorted > $readCountMatrix";
		print LOG "# Doing: $cmd\n" if ($verbose);
		
		system($cmd);
		system("rm $readCountMatrixUnsorted");
	    }
	    
	    
	    my $readCountMatrixUniq ="$workDir/readCount.$db.uniq.matrix";
	    my $readCountMatrixUnsortedUniq  ="readCount.$db.uniq.matrix.Unsorted";
	    
	    if ((! -e "$readCountMatrixUniq") && ($readsUniq > 0)){
		my $i=0;
		while (defined($ra_mapToListUniq->[$i])){
		    my $j=$i+1;
		    readCountAnalysis($ra_mapToListUniq->[$i], $j, $db, $readCountMatrixUnsortedUniq);	    
		    $i++;
		}
		#
		# sort the readCount matrix according to fastaEntry name
		#
		$cmd = "sort -k 1 $readCountMatrixUnsortedUniq > $readCountMatrixUniq";
		print LOG "# Doing: $cmd\n" if ($verbose);
		
		system($cmd);
		system("rm $readCountMatrixUnsortedUniq");
	    }
	}
	#####################################################################################
	#
	# calc coverage, depth ... using mpileup
	#
	#####################################################################################

	my $summaryFileMain ="$workDir/stat.$db.txt";
	my $fsaMain ="$workDir/contig.$db.fasta";
	if ((! -e $summaryFileMain) || (! -e $fsaMain)){
	    my @summaryFiles=();
	    my @Fsa=();
	    my $i=0;
	    while (defined($ra_mapToList->[$i])){
		my $j=$i+1;
		my $summaryFile ="$workDir/stat.$db.txt.$j";
		push(@summaryFiles,$summaryFile);

		my $pileup ="$workDir/pileup.$db.$j";
		bamToPileup($ra_mapToList->[$i], $j, $db, $pileup);
		push(@rmFiles,$pileup);
		
		if (! -z $pileup){
		    if (! -e $summaryFile){
			pileupToDepth($pileup, $j, $db, $summaryFile);
		    }
		    #
		    # make assembled fasta files if requested - option -a
		    #
		    if (exists($fastaOut{$db})){		    
			my $fsa ="$workDir/contig.$db.$j.fasta";
			push(@Fsa,$fsa);
			pileupToFasta($pileup,$fsa);
		    }
		}
		$i++;
	    }

	    #
	    # make assembled fasta files if requested - option -a
	    #
	    if (exists($fastaOut{$db})){		    
		my $cmd='';
		if ($#Fsa == 0){
		    $cmd = "mv $Fsa[0] $fsaMain";
		}
		elsif ($#Fsa >0){
		    my $i=0;
		    $cmd = "cat";
		    while (defined($Fsa[$i])){
			$cmd .= " $Fsa[$i]";
			$i++;
		    }
		    $cmd .= " > $fsaMain";
		}
		print LOG "Doing: $cmd\n" if ($verbose);
		system("$cmd");
		
		push(@rmFiles,@Fsa);
	    }
	    #
	    # cat all summaryfiles for database db
	    #
	    my $tmpfile = $summaryFileMain . '.tmp';
            my $count=0;
	    my $cmd = "cat ";
	    foreach my $id (@summaryFiles){
		$count++;
		chomp($id);
		$cmd .= "$id ";
	    }
	    $cmd .= "| sort -k 4 -nr > $tmpfile";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system($cmd);
	    mergeSummaryFiles($tmpfile,$summaryFileMain,$count);
            if (-z $summaryFileMain){
                print LOG "$summaryFileMain is empty - I remove it\n";
                system("rm $summaryFileMain");
            }
	    if ($cleanup){
		foreach my $id (@summaryFiles){
		    if (-e $id){
			system("rm $id");
		    }
		}
	    }
	}

	#
	# make summary file for uniq mapped reads
	# 
	my $summaryFileMainUniq ="$workDir/stat.$db.uniq.txt";

	if ((! -e $summaryFileMainUniq) || (-z $summaryFileMainUniq)){
	    my @summaryFiles=();
	    my $i=0;
	    my $cmd;

            while (defined($ra_mapToListUniq->[$i])){
                my $j=$i+1;
                my $summaryFile ="$workDir/stat.$db.uniq" . ".txt.$j";
                push(@summaryFiles,$summaryFile);
                my $pileup ="$workDir/pileup.$db.uniq.$j";
                bamToPileup($ra_mapToListUniq->[$i], $j, $db, $pileup);
                push(@rmFiles,$pileup);

		if (! -z $pileup){
                    pileupToDepth($pileup, $j, $db, $summaryFile);
                }
                $i++;
	    }
	    #
	    # cat all summaryfiles for database db
	    #
	    my $tmpfile = $summaryFileMainUniq . '.tmp';
	    my $count=0;
	    my $cmd = "cat ";
	    foreach my $id (@summaryFiles){
		$count++;
		chomp($id);
		$cmd .= "$id ";
	    }
	    $cmd .= "| sort -k 4 -nr > $tmpfile";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system($cmd);
	    mergeSummaryFiles($tmpfile,$summaryFileMainUniq,$count);
            if (-z $summaryFileMainUniq){
                print LOG "$summaryFileMainUniq is empty - I remove it\n";
                system("rm $summaryFileMainUniq");
            }
	    if ($cleanup){
		foreach my $id (@summaryFiles){
		    if (-e $id){
			system("rm $id");
		    }
		}
	    }
	}



	$abundance{$db}{summary}=$summaryFileMain;
	$abundance{$uniqdb}{summary}=$summaryFileMainUniq;

	$lastMapNot = $mapNot;
	if ($cleanup){
	    push(@rmFiles,$mapNot);
	}
	
	if ($#rmMapNotList > 0){
	    $cmd="rm $rmMapNotList[0]";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system("$cmd");
	    shift(@rmMapNotList);
	}
    }
    return($lastMapNot);
}


sub mapToDatabases{
    my ($startFile, $db) = @_;

    my $mapToFinal = "$workDir/mapTo.$db.bam";
    my $mapToFinalUniq = "$workDir/mapTo.$db.uniq.bam";

    my $lnMapToFinal = "mapTo.$db.bam";
    my $lnMapToFinalUniq = "mapTo.$db.uniq.bam";
    my $lnMapNotFinal = "mapNot.$db.bam";
    my $mapNotFinal = "$workDir/mapNot.$db.bam";
    my $newMapNot = $startFile;
    my $mapTo;
    my $mapToSort;

    my $mapToUniq;
    my $mapToSortUniq;
    my $sorted;
    my @mapToList=();
    my @mapToListUniq=();

    my $i=1;
    my $cmd;

    if ((-e $mapToFinal) && (-e $mapToFinalUniq)){
	while (defined($idx{index}{$db}[$i])){
	    $mapToSort = "$workDir/mapTo.$db.$i.sort.bam";
	    $mapToSortUniq = "$workDir/mapTo.$db.uniq.$i.sort.bam";
	    print LOG "Adding to mapToList; $mapToSort\n";
	    push(@mapToList,$mapToSort);
	    push(@mapToListUniq,$mapToSortUniq);

	    my $prelimMapNot = "$workDir/mapNot.$db.$i.bam";
	    if (-e $prelimMapNot){
		$mapNotFinal = $prelimMapNot;
	    }
	    else{
		$mapNotFinal = 'dummy';
	    }
	    $i++;
	}
	print LOG "File already exists: $mapToFinal\n" if ($verbose);

	return($mapNotFinal,$mapToFinal, \@mapToList, $mapToFinalUniq, \@mapToListUniq);
    }

    while (defined($idx{index}{$db}[$i])){
	my $j=$i-1;
	my $previousMapNot="$workDir/mapNot.$db.$j.bam";
	my $map = "$workDir/map.$db.$i.bam";
	my $mapNot = "$workDir/mapNot.$db.$i.bam";
	$mapTo = "$workDir/mapTo.$db.$i.bam";
	$mapToUniq = "$workDir/mapTo.$db.uniq.$i.bam";
	$mapToSortUniq = "$workDir/mapTo.$db.uniq.$i.sort.bam";
	$mapToSort = "$workDir/mapTo.$db.$i.sort.bam";

        if ((-e $mapToSort) && (-e $mapToSortUniq)){
            push(@mapToList,$mapToSort);
	    print LOG "Adding to mapToList; $mapToSort\n";
            push(@mapToListUniq,$mapToSortUniq);
            if (-e $mapNot){
                $newMapNot=$mapNot;
            }
            else{
                $newMapNot='dummy';
            }
            $i++;
            next;
        }

        if ($newMapNot eq 'dummy'){
            print LOG "Can not recover from this state as last mapNot is set to 'dummy'\n";
            print LOG "db=$db i=$i\n";
            print LOG "I will exit now\n";
	    print STDERR "\nDone!\n";
            exit;
        }
	my $mapToSortPrefix = "$workDir/mapTo.$db.$i.sort";
	my $mapToSortPrefixUniq = "$workDir/mapTo.$db.uniq.$i.sort";

	#
	# map against databases
	#
	my $index = $idx{index}{$db}[$i];
	my $annFile = $idx{ann}{$db}[$i];
	
	#
	# extract reads from $startFile and run bwa against database $index
	#    
	print LOG "calling mapToIndex with parameters: '$newMapNot' '$index' '$map' for i=$i\n" if ($verbose);

	mapToIndex($newMapNot,$index,$map);
	my $annotationFile = $index . '.ann';
        open(ANN,"head -1 $annotationFile |");
        $_=<ANN>;
        chomp;
        my @w=split(/\s+/);
        close(ANN);
        $countEntries{$db}{sequences} += $w[1];
        $countEntries{$db}{nuc} += $w[0];


	if (! defined($Getopt::Std::opt_p)){
	    if ($verbose){
		print LOG "# Calling subroutine analyzeTrueHit with these parameters:\n";
		print LOG "# $map\n# $mapTo\n# $mapNot\n# $matchRatio\n# $annFile\n"; 
	    }
	    analyzeTrueHit($map, $mapTo, $mapToUniq, $mapNot, $matchRatio);
	}
	else{
	    mapPaired($map, $mapTo);
	    mapNotPaired($map, $mapNot);
	}

	if (! -e $mapToSort){
	    sortBam($mapTo, $mapToSortPrefix);
	    $cmd = "rm $mapTo";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system($cmd);
	}
	if (! -e $mapToSortUniq){
	    sortBam($mapToUniq, $mapToSortPrefixUniq);
	    $cmd = "rm $mapToUniq";
	    print LOG "# Doing: $cmd\n" if ($verbose);

	    system($cmd);
	}
		
	if ($cleanup){
	    if (-e $previousMapNot){
		$cmd="rm $previousMapNot";
		print LOG "# Doing: $cmd\n" if ($verbose);
		system($cmd);
	    }
	    if (-e $map){
		$cmd="rm $map";
		print LOG "# Doing: $cmd\n" if ($verbose);
		system($cmd);
	    }
	}
	push(@mapToList,$mapToSort);
	print LOG "Adding to mapToList; $mapToSort\n";

	push(@mapToListUniq,$mapToSortUniq);
	print LOG "Adding to mapToListUniq; $mapToSortUniq\n";
	$newMapNot=$mapNot;
	$i++;
    }
    
    if (($#mapToList == 0) && (! -e $mapToFinal)){
	my @tmp=split(/\//,$mapToList[0]);

	chdir($workDir);
	$cmd="ln -s $tmp[-1] $lnMapToFinal";
	print LOG "# Doing: $cmd\n" if ($verbose);

	system($cmd);
	chdir($origDir);
    }
    else{
	my $merge = "$workDir/merge.$db.sam";
	mergeBam(\@mapToList,$merge,$mapToFinal);
    }

    if (($#mapToListUniq == 0) && (! -e $mapToFinalUniq)){
	my @tmp=split(/\//,$mapToListUniq[0]);

	chdir($workDir);
	$cmd="ln -s $tmp[-1] $lnMapToFinalUniq";
	print LOG "# Doing: $cmd\n" if ($verbose);

	system($cmd);
	chdir($origDir);
    }
    else{
	my $merge = "$workDir/merge.$db.uniq.sam";
	mergeBam(\@mapToListUniq,$merge,$mapToFinalUniq);
    }

    $mapNotFinal = $newMapNot;
    if (($mapNotFinal ne "$workDir/cleaned.nophiX.bam") && ($mapNotFinal ne 'dummy')){
	push(@rmMapNotList,$mapNotFinal);
    }
    return($mapNotFinal,$mapToFinal, \@mapToList, $mapToFinalUniq, \@mapToListUniq);
}

sub mergeBam{
    my ($ra, $sam, $bam)=@_;
    if (-e $sam){
	system("rm $sam");
    }
    my $i=0;
    my $cmd;
    while (defined($ra->[$i])){
	$cmd = "$prog_samtools view -H $ra->[$i] >> $sam";
	print LOG "# In mergeBam: $cmd\n" if ($verbose);

	system("$cmd");
	$i++;
    }
    $i=0;
    while (defined($ra->[$i])){
	$cmd = "$prog_samtools view $ra->[$i] | cat >> $sam";
	print LOG "# In mergeBam: $cmd\n" if ($verbose);

	system("$cmd");
	$i++;
    }
    
    $cmd = "$prog_samtools view $sam -Sb | cat > $bam";
    print LOG "# In mergeBam: $cmd\n" if ($verbose);

    system("$cmd");
    
    $cmd="rm $sam";
    print LOG "# In mergeBam: $cmd\n" if ($verbose);

    system("$cmd");
    return;
}

sub mapPaired{
    my ($inFile, $outFile) = @_;

    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
        return;
    }

    my $cmd = "$prog_samtools view -f1 -F12 -b $inFile > $outFile";

    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    return;
}

sub mapNotPaired{
    my ($inFile, $outFile) = @_;

    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
        return;
    }

    my $cmd = "$prog_samtools view -f12 -b $inFile > $outFile";

    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    return;
}


sub mapProperlyPaired{
    my ($inFile, $outFile) = @_;

    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
        return;
    }

    my $cmd = "$prog_samtools view -f2 -b $inFile | cat > $outFile";

    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    return;
}

sub mapNotProperlyPaired{
    my ($inFile, $outFile) = @_;

    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
        return;
    }

    my $cmd = "$prog_samtools view -F2 -b $inFile | cat > $outFile";

    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    return;
}

sub analyzeTrueHit{
    my ($map, $mapTo, $mapToUniq, $mapNot, $frac) = @_;

    if (($redo==0) && (-e $mapTo) && (-e $mapToUniq) && (-e $mapNot)){
        return;
    }
    my $cmd;
    my $excludeFromPreNotMap=0;
    my $line;
    my $reject=0;

    #
    # No criteria on reads other that they are properly paired
    #
    if (! defined($Getopt::Std::opt_p)){
	if (($absoluteMatchCountFlag == 0) && ($matchRatioFlag == 0)){
	    mapProperlyPaired($map, $mapTo);
	    mapNotProperlyPaired($map, $mapNot);
	    return;
	}
	else{
	    $cmd = "$prog_samtools view -h -f2 $map |";
	}
    }
    else{
	if (($absoluteMatchCountFlag == 0) && ($matchRatioFlag == 0)){
	    mapPaired($map, $mapTo);
	    mapNotPaired($map, $mapNot);
	    return;
	}
	else{
	    $cmd = "$prog_samtools view -h -f1 $map |";
	}
    }

    if ($absoluteMatchCountFlag){
	print LOG "Paired reads are only accepted if both have a sum of Matches+Mismatches >= $absoluteMatchCount\n";
    }
    if ($matchRatioFlag){
	print LOG "Paired reads are only accepted if both have a sum of Matches+Mismatches divided by length of read >= $matchRatio\n";
    }
    my %rec=();

    my $samNotOk="$mapTo.notOk.sam";
    open(NO,">$samNotOk");

    my $samOk="$mapTo.Ok.sam";
    my $samOkUniq="$mapTo.Ok.Uniq.sam";
    open(UNIQ,">$samOkUniq");
    open(YES,">$samOk");

    if ($verbose){
	print LOG "Doing: $cmd\n";
	print LOG "# Doing $cmd and writing to $samOk and $samNotOk\n";
    }
    open(FILE,"$cmd");
    my $line_read;
    my $line_mate;
    my $reads=0;
    my $headerLine=1;
    my @matchRatios=();
    my @matchCounts=();
    my @as=();
    my @xs=();
    my @uniqHit=();
    my @rejectCause=();
    for (my $i=0;$i<=3;$i++){
	$rejectCause[$i]=0;
    }
    while (defined($_=<FILE>)){

	if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
	    # this is a header line
	    print YES "$_";
	    print UNIQ "$_";
	    next;
	}
	$headerLine=0;
	$reads++;
	$line_read=$_;
        chomp($line_read);
        my @w_read=split(/\s+/,$line_read);

        my $len=length($w_read[9]);
        my $cigar = $w_read[5];
#
	$xs[0] = substr($w_read[14],5);
	$as[0] = substr($w_read[13],5);
	#
        # Now count number of matches in cigar string
        #
	my $matchSum=0;
        my @tmp = $cigar =~ m/(\d+)M/g;
        foreach my $s (@tmp) {
            $matchSum += $s;
        }

	$matchRatios[0] = $matchSum/$len;
	$matchCounts[0] = $matchSum;
	#
	# calc the same values as above for the mate
	#
	$line_mate=<FILE>;
        chomp($line_mate);
	$reads++;
        my @w_mate=split(/\s+/,$line_mate);

        $len=length($w_mate[9]);
        $cigar = $w_mate[5];
	$xs[1] = substr($w_mate[14],5);	
	$as[1] = substr($w_mate[13],5);

	$matchSum=0;
        @tmp = $cigar =~ m/(\d+)M/g;
        foreach my $s (@tmp) {
            $matchSum += $s;
        }

	$matchRatios[1] = $matchSum/$len;
	$matchCounts[1] = $matchSum;

	#
	# Init flags as true
	#
	my @matchRatioFlags=();
	$matchRatioFlags[0]=1;
	$matchRatioFlags[1]=1;

	my @matchCountFlags=();
	$matchCountFlags[0]=1;
	$matchCountFlags[1]=1;

	my $uniqHit=();;
	$uniqHit[0]=0;
	$uniqHit[1]=0;

	for (my $i=0;$i<=1;$i++){
	    if ($matchRatioFlag){
		if ($matchRatios[$i] < $frac){
		    $matchRatioFlags[$i] = 0;
		    $rejectCause[$i]++;
		}
	    }
	    if ($absoluteMatchCountFlag){
		if ($matchCounts[$i] < $absoluteMatchCount){
		    my $j=2+$i;
		    $matchCountFlags[$i] = 0;
		    $rejectCause[$j]++;
		}
	    }
	    if ($as[$i] > $xs[$i]){
		$uniqHit[$i]=1;
	    }
	}
	
	my $readsOk = 1;
	for (my $i=0;$i<=1;$i++){
	    $readsOk *= $matchRatioFlags[$i];
	    $readsOk *= $matchCountFlags[$i];
	}

	if ($readsOk){
	    print YES "$line_read\n";
	    print YES "$line_mate\n";
	    if ( ($uniqHit[0]== 1) || ($uniqHit[1] == 1) ){
		print UNIQ "$line_read\n";
		print UNIQ "$line_mate\n";
	    }
	}
	else{
	    print NO "$line_read\n";
	    print NO "$line_mate\n";
	    $reject += 2;
	}
    }
    close(FILE);
    if ($verbose){
	print LOG "REJECT: db= $mapTo\n";
	print LOG "REJECT: Removed $reject out of $reads reads\n";
	print LOG "REJECT: Encoding of reads marked to be rejected\n";
	print LOG "REJECT: 0: matchRatioFlags for read, 1: matchRatioFlags for mate, 2: absoluteMatchFlags for read, 3: absoluteMatchFlag for mates\n";
	for (my $i=0;$i<=3;$i++){
	    print LOG "REJECT:\t$i\t$rejectCause[$i] reads marked to be rejected\n";
	}
	print LOG "REJECT: A read pair will be rejected if either read or mate are maked to be rejected\n";
    }
    close(YES);
    close(UNIQ);
    close(NO);

    $cmd="$prog_samtools view $samOk -Sb | cat > $mapTo";
    print LOG "Doing: $cmd\n" if ($verbose);

    system("$cmd");

    $cmd="$prog_samtools view $samOkUniq -Sb | cat > $mapToUniq";
    print LOG "Doing: $cmd\n" if ($verbose);

    system("$cmd");

    my $cmdRead;
    if (! defined($Getopt::Std::opt_p)){
	#
	# Make mapNot file i.e. reads that did not map properly and reads that were rejected dur to filtering criteria
	#
	$cmdRead = "$prog_samtools view -h -F2 $map |";
    }
    else{
	#
	# Make mapNot file i.e. reads that were not paired and reads that were rejected dur to filtering criteria
	#
	$cmdRead = "$prog_samtools view -h -F1 $map |";
    }
    my $cmdWrite = "| $prog_samtools view -Sb - | cat > $mapNot";

    if ($verbose){
	print LOG "# Reading via $cmdRead\n";
	print LOG "# Writing via $cmdWrite\n";
    }
    open(MAPNOT,"$cmdWrite");
    open(A,"$cmdRead");
    while (defined($_=<A>)){
	print MAPNOT ("$_");
    }
    close(A);
    
    if ($verbose){
	print LOG "# Reading from $samNotOk\n";
	print LOG "# Writing via $cmdWrite\n";
    }
    
    #
    # now include reads that were removed dut to filtering criteria
    #
    open(A,"<$samNotOk");
    while (defined($_=<A>)){
	print MAPNOT ("$_");
    }
    close(A);
    close(MAPNOT);

    my @rmNow=();
    push(@rmNow,$samOk);
    push(@rmNow,$samOkUniq);
    push(@rmNow,$map);
    push(@rmNow,$samNotOk);
    if ($cleanup){
	foreach my $id (@rmNow){
	    if ((-e $id) && ($id ne "$workDir/cleaned.nophiX.bam")){
		$cmd = "rm $id";
		print LOG "# Doing: $cmd\n" if ($verbose);

		system("$cmd");
	    }
	}
    }
    return;
}


sub cleanData{
    my ($fastqF, $fastqR, $outFile) = @_;

    my $dbName = $idx{index}{phiX174}[1];
    my $cmd = "$prog_bwa mem  -t $cores $dbName $fastqF $fastqR | $prog_samtools view -f1 -f12 -Sb - | cat > $outFile";
    
    print LOG "# Doing: $cmd\n\n" if ($verbose);

    system("$cmd");
    return;
}

sub cutadapt{
    my ($inFile, $outFile) = @_;

    my $cmd = "$prog_cutadapt $inFile > $outFile";
    

#    system("$cmd");
    return($cmd);
}

sub sortBam{
    my ($inFile, $outFile) = @_;

    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
	print LOG "# Skipping as file exists: $outFile\n" if ($verbose);
	return($outFile);
    }

    my $cmd = "$prog_samtools sort $inFile $outFile";
    
    print LOG "# Doing: $cmd\n\n" if ($verbose);
    system("$cmd");
    
    return;
}
sub bamTofastq{
    my ($bam,$fq) = @_;

    $cmd = "$prog_bamtools convert -in $bam -format fastq | cat > $fq";
    
    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    if (-z $fq){
	print LOG "converted $bam to a fastq file which is empty: $fq\n";
	print STDERR "Dome!\n";
	exit;
    }
    return;
}

sub mapToIndex{
    my ($inFile, $index, $outFile) = @_;
    my $cmd;

    if ( ($redo==0) && (-e $outFile) && (! -z $outFile) ){
	print LOG "# Skipping mapToIndex as file already exists: $outFile\n" if ($verbose);
	return;
    }
    

    $cmd = "$prog_bamtools convert -format fastq -in $inFile | $prog_bwa mem -h $XA_hits -D $skip_chain_fraction -t $cores -M -p $index";
    if (defined($Getopt::Std::opt_M)){
	$cmd .= " -I $Getopt::Std::opt_M ";
    }
    $cmd .= " - | $prog_samtools view -F 256 -Sb - | cat > $outFile";
    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    return;
}


sub mapToBam{ 
    my ($inFile, $outFile) = @_;
    
    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
	return;
    }
    
    my $cmd = "$prog_samtools view -F 4 -b $inFile | cat > $outFile";
    
    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");    
    return;
}

sub mapNotBam{ 
    my ($inFile, $outFile) = @_;
    
    if ( ($redo == 0) && (-e $outFile) && (! -z $outFile) ){
	return;
    }
    
    my $cmd = "$prog_samtools view -f 4 -b $inFile | cat > $outFile";
    
    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");    
    return;
}

sub bamToPileup{
    my ($bam, $i, $db, $outFile, $indexFile) = @_;
    my $index;
    my $cmd;
    if (-e $indexFile){
	$index=$indexFile;
    }
    else{
	$index = $idx{index}{$db}[$i];
    }
    
    if ((! -e $outFile) || (-z $outFile)){
	$cmd = "$prog_samtools mpileup -f $index $bam | cat > $outFile";	
	print LOG "# Doing: $cmd\n" if ($verbose);
	system("$cmd");
    }

    return;
}
sub bamToPileupCmd{
    my ($bam, $i, $db, $outFile, $indexFile) = @_;
    my $index;
    if (-e $indexFile){
	$index=$indexFile;
    }
    else{
	$index = $idx{index}{$db}[$i];
    }
    my $cmd = '';

    # the -A option in samtools mpileup has the effect that also not properly paired reads are used
    # that can be usefull if bamtools convert actually make fastq files where reads are reversed if stated in bam file
    # and when its reversed it can no longer be properly paired as that only mean reads like this: ----->   <------
    if ((! -e $outFile) || (-z $outFile)){
#	$cmd = "$prog_samtools mpileup -A -f $index $inFile | cat > $outFile";
	$cmd = "$prog_samtools mpileup -f $index $bam | cat > $outFile";
    }

    return($cmd);
}

sub pileupToDepth{
    my ($inFile, $i, $db, $summaryFile, $annFile) = @_;
    my $ann;
    if (-e $annFile){
	$ann=$annFile;
    }
    else{
	$ann = $idx{ann}{$db}[$i];
    }

    my $cmd = "$prog_mpileup2stat -i $inFile -a $ann  | cat > $summaryFile";
    
    if (! -e $summaryFile){
	print LOG "# Doing: $cmd\n" if ($verbose);	
	system("$cmd");    
    }
    return;
}

sub pileupToDepthCmd{
    my ($inFile, $i, $db, $summaryFile, $annFile) = @_;
    my $ann;
    if (-e $annFile){
	$ann=$annFile;
    }
    else{
	$ann = $idx{ann}{$db}[$i];
    }

    my $cmd = "$prog_mpileup2stat -i $inFile -a $ann  | cat > $summaryFile";
    
    return($cmd);
}

sub pileupToFasta{
    my ($inFile, $fasta, $sequenceName, $refseqLen) = @_;

    my $cmd='';
    if (! -e $fasta){
	$cmd = "$prog_pileup2fasta -i $inFile -c $fasta -b $SNP_threshold";
	if (defined($sequenceName)){
	    $cmd .= " -n '$sequenceName'";
	}
	if (defined($refseqLen)){
	    $cmd .= " -L $refseqLen";
	}

	print LOG "# Doing: $cmd\n" if ($verbose);
	system("$cmd");	
    }
    return;
}

sub pileupToFastaCmd{
    my ($inFile, $fasta, $sequenceName, $refseqLen) = @_;

    $cmd = "$prog_pileup2fasta -i $inFile -c $fasta -b $SNP_threshold";
    if (defined($sequenceName)){
	$cmd .= " -n '$sequenceName'";
    }
    if (defined($refseqLen)){
	$cmd .= " -L $refseqLen";
    }
    return($cmd);
}

sub readCountAnalysis{
    my ($inFile, $i, $db, $matrix) = @_;
    my $annFile = $idx{ann}{$db}[$i];
    my %rec=();


    open(ANN,"<$annFile");
    if (! eof ANN){
        # Read header line
	$_=<ANN>;
    }
    my $fastaName;
    my $desc;
    my @w=();
    while (defined($_=<ANN>)){
	if (m/^(\d+) (\S+) (.*)/){
	    $fastaName=$2;
	    $desc=$3;
	}

        if (! exists($rec{$fastaName}{count})){
            $rec{$fastaName}{count}=0;
            $rec{$fastaName}{desc}=$desc;
        }

        # read nest line which include length of that fasta entry
	$_=<ANN>;
        @w=split(/\s+/);
        $rec{$fastaName}{len}=$w[1];
    }
    close(ANN);

    my $cmd="$prog_samtools view $inFile | cut -f3 | sort | uniq -c |";
    print LOG "# Doing: $cmd\n" if ($verbose);

    open(READCOUNT,"$cmd ");
    while (defined($_=<READCOUNT>)){
        my @w=split(' ',$_);
        my $count=$w[0];
        my $fastaName=$w[1];
        if (! exists($rec{$fastaName}{count})){
            print LOG "Error In readCountAnalysis: $inFile,$i,$db\n";
            print LOG "A new fasta entry name is found '$fastaName' which was not in the file: $annFile\n";
            print LOG "$_";
        }
        else{
            $rec{$fastaName}{count}=$count;
        }
    }
    close(READCOUNT);
    if ($i == 1){
	print LOG "# writing to $matrix\n" if ($verbose);
        open(MATRIX,"| cat > $matrix");
    }
    else{
	print LOG "# adding to $matrix\n" if ($verbose);
        open(MATRIX,"| cat >> $matrix");
    }

    foreach my $id (keys %rec){
        print MATRIX "$id\t$rec{$id}{count}\t$rec{$id}{len}\t$rec{$id}{desc}\n";
    }
    close(MATRIX);
    return;
}

sub countReads{
    my ($fastq) = @_;

    my $cmd;
    chomp($fastq);
    if ($fastq =~ '.gz$'){
#	$cmd = "gunzip -c $fastq | grep -c \'^+$\'";
	$cmd = "gunzip -c $fastq | wc -l";
    }
    else{
#	$cmd = "grep -c \'^+$\' $fastq";
	$cmd = "wc -l $fastq";
    }
    print LOG "# Doing: $cmd\n" if ($verbose);

    my $reads = `$cmd`;
    chomp($reads);

    $reads /= 4;

    return($reads);
}

sub bamCountReads{
    my ($bam) = @_;

    my $count;
    chomp($bam);
    if (exists($readCount{$bam})){
        $count=$readCount{$bam};
    }
    else{
	print LOG "# Doing: $prog_samtools flagstat $bam | head -1 | cut -f1 -d' '\n" if ($verbose);

        $count=`$prog_samtools flagstat $bam | head -1 | cut -f1 -d' '`;
        chomp($count);
        open(TMP,">>$readCountFile");
        print TMP "$bam\t$count\n";
        close(TMP);
    }
    return($count);
}
