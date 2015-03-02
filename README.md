# MGmapper
A metagenomics mapping package
# First the environment vairable needs to be set
setenv MGmapper /home/projects9/tnp/programs/MGmapper.v1.4

# In directory $MGmapper/test_run the MGmapper programs have been tested. Two simulated fastq files are generated using data from the database:
# $MGmapper/db/Bacteria_test/
# The two fastq files have been used to test the pair-end version (MGmapper_PE.pl), where output were directed to outDir_PE. Try to run the command below,
# where the output directory have been changed to testing_PE and check the file MGmapper.log for and error messages.

cd $MGmapper/test_run
$MGmapper/MGmapper_PE.pl -i simulated.F.fq.gz -j simulated.R.fq.gz -d testing_PE -C 1

# For the single-end version, edit the location of the files in: $MGmapper/test_run/fastq.list
$MGmapper/MGmapper_PE.pl -f fastq.list -d testing_SE -C 1

# compare results (MGmapper.summary and MGmapper.log) to the output generated in $MGmapper/test_run/outDir_SE and output you just generated in testing_SE
