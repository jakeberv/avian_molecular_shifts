#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# perl script to use mitofinder with a set assembled sequences
#   - see mitofinder call below, next to comment that reads "Use mitofinder 
#     to conduct the annotation"
#   - will organize the large gene regions (protein coding and rRNAs) into
#     a tab-delimited format that includes the sequences
#   - sequences can be extracted from the appropriate columns, aligned, and
#     concatenated
############################################################################

############################################################################
# Set the global variables
############################################################################
my($progname) = $0;

# Set for vertebrate mitochondrial code with 4 processors and 10 GB memory
# Change for other run modes
my($settings) = "-o 2 -p 4 -m 10";

my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($nter);
my($qter);

my($yter);
my($zter);

my($tempch1);
my($tempch2);
my($tempvar);
my($tempstr);
my @temparray;

my($verbose) = 0;

if ( @ARGV != 1 ) {
	print STDERR "Usage:\n  \$ $progname <ctlfile>\n";
	print STDERR "  ctlfile   = tab delimited list of mitogenomes to annotate\n";
	print STDERR "              (read the comments in this code for format details)\n";
	print STDERR "exiting...\n";
	exit;
}

my($ctlfile) = $ARGV[0];

############################################################################
# Read the ctlfile - the format has information about one mitogenome per 
#                    line in this order:
# 	0. Outfile name (e.g., Atlapetes_gutturalis_UWMB93636)
# 	1. Source note - any notes regarding source
# 	2. Source - N = NCBI; otherwise the path to the file
# 	3. Read file or accession number
# 	4. Reference mitogenome (e.g., Passer_montanus_JX486030.gb)
#   5. Source of reference mitogenome - N = NCBI; F = File
############################################################################
my @inlist;
my($com);
open (my $INF, "$ctlfile") or die "Could not open file $ctlfile for input.\n";
@inlist = <$INF>; # Read the tab-delimited file of mitogenomes
close($INF) or die "Could not close file $ctlfile\n";
my($inlistnum) = $#inlist + 1;

print "\n";
print "########################################\n";
print "Will analyze $inlistnum datasets:\n";
for ($iter=0; $iter<$inlistnum; $iter++) {
	chomp($inlist[$iter]);
	(@temparray) = split(/\s+/, $inlist[$iter]);
	print "  -- $temparray[0] - $temparray[3]\n"
}
print "########################################\n";
print "\n";
print "\n";

############################################################################
# Read the inputfile
############################################################################
my @mitofinderoutput;
my($found);

open (my $OUTF, ">$ctlfile.OUTPUT.txt") or die "Could not open file $ctlfile.OUTPUT.txt for output.\n";

for ($iter=0; $iter<$inlistnum; $iter++) {
	chomp($inlist[$iter]);
	(@temparray) = split(/\s+/, $inlist[$iter]);
	
	if ( uc($temparray[1]) eq "S" ) { my($assembler) = "--megahit"; }
	$tempvar = localtime();
	print "\n";
	print "****************************************\n";
	print " Annotating data for:\n";
	print " 	$temparray[0]\n";
	print " Using fasta file:\n";
	print " 	$temparray[3]\n";
	print " Note regarding source:\n";
	print " 	$temparray[1]\n";
	if ( uc($temparray[2]) eq "N" ) { print " 	-- downloaded from SRA\n"; }
	else { print " 	-- local file\n"; }  # Default is local file
	print " Annotate using reference:\n";
	print " 	$temparray[4]\n";
	print " Analysis start:\n";
	print " 	$tempvar\n";
	print "****************************************\n";
	print "\n";
	
	if ( uc($temparray[2]) eq "N" ) { # download from NCBI
		system("efetch -db nuccore -id $temparray[3] -format fasta > temp.fasta");
	}
	else {
		system("cp $temparray[2]/$temparray[3] temp.fasta")
	}
	
	if ( uc($temparray[5]) eq "N" ) { # download genbank file from NCBI
		system("efetch -db nuccore -id $temparray[4] -format gb > temp.gb");
		$temparray[4] = "temp.gb";
	}
	
	# Use mitofinder to conduct the annotation
	system("mitofinder -j $temparray[0] -a temp.fasta -r $temparray[4] $settings");
	unlink("temp.fasta");
	
	$tempvar = localtime();
	print "\n";
	print " Assembly complete at time:\n";
	print " 	$tempvar\n";
	print "****************************************\n";
	print "\n";
	
	# Output the mitofinder data in the order:
	#	0.	Notes
	#	1.  source (SRA or file)
	#	2.  Accession (file if from file)
	#	3.  Species_plus
	#	4.	Number_of_genes_found
	#	5.  rRNA seq 1
	#	6.  rRNA seq 2
	#	7.  ND1_seq
	#	8.  ND2_seq
	#	9.  COX1_seq
	#	10. COX2_seq
	#	11. ATP8_seq
	#	12. ATP6_seq
	#	13. COX3_seq
	#	14. ND3_seq
	#	15. ND4L_seq
	#	16. ND4_seq
	#	17. ND5_seq
	#	18. CYTB_seq
	#	19. ND6_seq
	# Elements 5 through 19 are actually triple columns. The first column is 0 for "not 
	# found" and 1 for "found", the second is the sequence length, and the third is the
	# sequence. The sequence is listed as "not_found" if the gene was not found
	if ( uc($temparray[2]) eq "S" ) { print $OUTF "SRA\t"; }
	else { print $OUTF "file\t"; }
	print $OUTF "$temparray[1]\t";
	print $OUTF "$temparray[3]\t";
	print $OUTF "$temparray[0]\t";
	
	# Copy the final fasta file to a temporary file:
	$tempstr = "sed \"s/ND4L/nd4l/g\" $temparray[0]" . "/" . "$temparray[0]" . "*_Final_Results/" . "$temparray[0]" . "_final_genes_NT.fasta > $temparray[0]" . ".TEMP.fasta";
	system("$tempstr");
	$tempstr = "$temparray[0]" . ".TEMP.fasta";
	open (my $MTF, "$tempstr") or die "Could not open file $tempstr for input.\n";
	@mitofinderoutput = <$MTF>; # Read the mitofinder output
	close($MTF) or die "Could not close file $tempstr\n";
	$zter = $#mitofinderoutput + 1;
	
	# Count the number of sequences in the temporary fasta  file
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/^>/) { $found++; }
	}
	print $OUTF "$found\t";
	
	# Extract small subunit rRNA (rrnS)
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/rrnS/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract large subunit rRNA (rrnL)
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/rrnL/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND1
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ND1/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND2
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ND2/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }

	# Extract COX1
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/COX1/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract COX2
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/COX2/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }

	# Extract ATP8
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ATP8/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ATP6
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ATP6/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract COX3
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/COX3/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND3
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ND3/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND4L
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/nd4l/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND4
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ND4/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND5
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ND5/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract CYTB
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/CYTB/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\t";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\t"; }
	
	# Extract ND6
	$found = 0;
	for ($jter=0; $jter<$zter; $jter++) {
		if ( $mitofinderoutput[$jter] =~ m/ND6/) {
			chomp($mitofinderoutput[$jter+1]);
			$qter = length($mitofinderoutput[$jter+1]);
			print $OUTF "1\t$qter\t$mitofinderoutput[$jter+1]\n";
			$found = 1;
			$jter = $zter;
		}
	}
	if ( $found == 0 ) { print $OUTF "0\t0\tnot_found\n"; }
	
}

close($OUTF) or die "Could not close file $ctlfile.OUTPUT.txt\n";

