#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);
use Carp;
use FindBin;
use File::Basename;

my $param_string = <<__PARAMS__;

################# Evidence Modeler ##############################
#
#  parameters:
#
#  Required:
# 
#  --partitions               partitions file 
#
#  Base File Names:
#
#  --genome            | -G   genome sequence in fasta format
#  --gene_predictions       | -g   gene predictions gff3 file
#  --protein_alignments | -p  protein alignments gff3 file
#  --transcript_alignments   | -e    est alignments gff3 file
#  --repeats          | -r    gff3 file with repeats masked from genome file
#  --terminalExons    | -t   supplementary file of additional terminal exons to consider (from PASA long-orfs)
#  --output_file_name | -O   name of output file to be written to in the directory containing the inputs
#
#  ** Files with full and fixed paths ** :
#
#  --weights          | -w    weights for evidence types file
#
#  Optional arguments:
#
#  --stop_codons             list of stop codons (default: TAA,TGA,TAG)
#                               *for Tetrahymena, set --stop_codons TGA
#  --min_intron_length       minimum length for an intron (default 20 bp)
#
# flags:
#
#  --forwardStrandOnly   runs only on the forward strand
#  --reverseStrandOnly   runs only on the reverse strand
#  -S                    verbose flag
#  --debug               debug mode, writes lots of extra files.
#
#  Misc:
#
#  --search_long_introns  <int>  when set, reexamines long introns (can find nested genes, but also can result in FPs) (default: 0 (off))
#
#
#  --re_search_intergenic <int>  when set, reexamines intergenic regions of minimum length (can add FPs) (default: 0  (off))
#  --terminal_intergenic_re_search <int>   reexamines genomic regions outside of the span of all predicted genes (default: 10000)
#
#
#################################################################

__PARAMS__

    ;


#  --stitch_ends             file listing source types to apply end stitching 
#                                 into existing exons (ie. 'genewise,alignAssembly')
#  --extend_to_terminal      file listing source types to extend terminal alignment segment into a terminal exon (ie. 'genewise')



## Input parameters
my $partitions_file;
my $genomicSeqFile;
my $genePredictionsFile;
my $estAlignmentsFile;
my $proteinAlignmentsFile;
my $repeatsFile;
my $weightsFile;
my $terminalExonsFile;
my $stitch_ends;

my $SEE;
my $DEBUG;
my $extend_to_terminal; # extend termini of genewise alignments to find a start or stop codon, forming a terminal exon.

# stop codon issues
my $DEFAULT_STOP_CODONS = "TAA,TGA,TAG";
my $stop_codons_arg;
my @STOP_CODONS;
my $FORWARD_STRAND_ONLY_FLAG;
my $REVERSE_STRAND_ONLY_FLAG;
my $min_intron_length;
my $output_file_name;
my $help_flag;

my $MIN_LONG_INTRON_LENGTH;
my $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH;
my $MIN_INTERGENIC_SIZE_ON_RE_SEARCH;


&GetOptions ("partitions=s" => \$partitions_file,
             "output_file_name|O=s" => \$output_file_name,
             "genome|G=s"=>\$genomicSeqFile,
             "gene_predictions=s"=>\$genePredictionsFile,
             "transcript_alignments|e=s"=>\$estAlignmentsFile,
             "protein_alignments|p=s" => \$proteinAlignmentsFile,
             "repeats|r=s"=>\$repeatsFile,
             "weights|w=s"=>\$weightsFile,
             "forwardStrandOnly"=>\$FORWARD_STRAND_ONLY_FLAG,
             "reverseStrandOnly"=>\$REVERSE_STRAND_ONLY_FLAG,
             "terminalExonsFile|t=s"=>\$terminalExonsFile,
             "stop_codons=s" => \$stop_codons_arg,
             "help|h" => \$help_flag,
             "stitch_ends=s" => \$stitch_ends,
             "extend_to_terminal=s" => \$extend_to_terminal,
             "min_intron_length=i" => \$min_intron_length,

             "search_long_introns=i" => \$MIN_LONG_INTRON_LENGTH,
             "re_search_intergenic=i" => \$MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH,
             "terminal_intergenic_re_search=i" => \$MIN_INTERGENIC_SIZE_ON_RE_SEARCH,
             
             );

unless ($output_file_name && $partitions_file && $genomicSeqFile && $genePredictionsFile && $weightsFile) {
    die $param_string;
}
if ($help_flag) { 
    die $param_string;
}

my $evm_bin_dir = $FindBin::Bin . "/../.";
my $evm_cmd = "$evm_bin_dir/evidence_modeler.pl";

## ensure basenames
$genomicSeqFile = basename($genomicSeqFile);
$genePredictionsFile = basename($genePredictionsFile);
$proteinAlignmentsFile = basename($proteinAlignmentsFile) if $proteinAlignmentsFile;
$estAlignmentsFile = basename($estAlignmentsFile) if $estAlignmentsFile;
$repeatsFile = basename($repeatsFile) if $repeatsFile;
$terminalExonsFile = basename($terminalExonsFile) if $terminalExonsFile;

## ensure full paths for those that requie it:
unless ($weightsFile =~ /^\//) {
    die "Error, weights file $weightsFile must be specified by the full path\n";
}
if ($stitch_ends && $stitch_ends !~ /^\//) {
    die "Error, stitch_ends file $stitch_ends must be specified by the full path\n";
}
if ($extend_to_terminal && $extend_to_terminal !~ /^\//) {
    die "Error, extend_to_terminal file $extend_to_terminal must be specified by the full path\n";
}


## add options to the command:
$evm_cmd .= " -G $genomicSeqFile -g $genePredictionsFile -w $weightsFile"; #required

$evm_cmd .= " -e $estAlignmentsFile " if $estAlignmentsFile;
$evm_cmd .= " -p $proteinAlignmentsFile " if $proteinAlignmentsFile;
$evm_cmd .= " -r $repeatsFile " if $repeatsFile;
$evm_cmd .= " -t $terminalExonsFile " if $terminalExonsFile;
$evm_cmd .= " --stitch_ends $stitch_ends " if $stitch_ends;
$evm_cmd .= " --extend_to_terminal $extend_to_terminal " if $extend_to_terminal;
$evm_cmd .= " --min_intron_length $min_intron_length " if $min_intron_length;
$evm_cmd .= " --stop_codons $stop_codons_arg " if $stop_codons_arg;

# add flags
$evm_cmd .= " --forwardStrandOnly " if $FORWARD_STRAND_ONLY_FLAG;
$evm_cmd .= " --reverseStrandOnly " if $REVERSE_STRAND_ONLY_FLAG;

# add misc opts
$evm_cmd .= " --search_long_introns $MIN_LONG_INTRON_LENGTH " if $MIN_LONG_INTRON_LENGTH;
$evm_cmd .= " --re_search_intergenic $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH " if $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH;
$evm_cmd .= " --terminal_intergenic_re_search $MIN_INTERGENIC_SIZE_ON_RE_SEARCH " if $MIN_INTERGENIC_SIZE_ON_RE_SEARCH;


open (my $fh, $partitions_file) or die "Error, cannot open $partitions_file";
while (<$fh>) {
    chomp;
    my ($accession, $acc_dir, $is_partitioned, $partition_dir) = split (/\t/);
    
    my $data_dir = $acc_dir;
    
    my $core_cmd = $evm_cmd;


    if ($is_partitioned eq "Y") {
        $data_dir = $partition_dir;
    }
    
    unless (-d $data_dir) {
        die "Error, cannot locate data directory $data_dir";
    }
    
    $core_cmd .= " --exec_dir $data_dir";

    # add output direction
    $core_cmd .= " > $data_dir/$output_file_name 2> $data_dir/$output_file_name.log";
    
    #my $cmd = "cd $data_dir; if [ \$? -eq 0 ]; then $evm_cmd; E=\$?; else E=9; fi; exit \$E;";
    
    print "$core_cmd\n";
    
    
}

exit(0);



