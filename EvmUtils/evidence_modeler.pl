#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);
use Carp;



my $param_string = <<__PARAMS__;

################# Evidence Modeler ##############################
#
#  parameters:
#
#  Required:
# 
#  --genome              genome sequence in fasta format
#  --weights              weights for evidence types file
#  --gene_predictions     gene predictions gff3 file
#
#  Optional but recommended:
#  --protein_alignments   protein alignments gff3 file
#  --transcript_alignments       transcript alignments gff3 file
#
#  Optional and miscellaneous
#
#  --repeats               gff3 file with repeats masked from genome file
#
#  
#  --terminalExons         supplementary file of additional terminal exons to consider (from PASA long-orfs)
#
#  --stop_codons             list of stop codons (default: TAA,TGA,TAG)
#                               *for Tetrahymena, set --stop_codons TGA
#  --min_intron_length       minimum length for an intron (default 20 bp)
#  --exec_dir                directory that EVM cd's to before running.
#
# flags:
#
#  --forwardStrandOnly   runs only on the forward strand
#  --reverseStrandOnly   runs only on the reverse strand
#
#  -S                    verbose flag
#  --debug               debug mode, writes lots of extra files.
#  --report_ELM          report the eliminated EVM preds too.
#
#  --search_long_introns  <int>  when set, reexamines long introns (can find nested genes, but also can result in FPs) (default: 0 (off))
#
#
#  --re_search_intergenic <int>  when set, reexamines intergenic regions of minimum length (can add FPs) (default: 0  (off))
#  --terminal_intergenic_re_search <int>   reexamines genomic regions outside of the span of all predicted genes (default: 10000)
#
#################################################################################################################################



__PARAMS__

    ;

#  --stitch_ends             file listing source types to apply end stitching 
#                                 into existing exons (ie. 'genewise,alignAssembly')
#  --extend_to_terminal      file listing source types to extend terminal alignment segment into a terminal exon (ie. 'genewise')
#
#  --INTERGENIC_SCORE_ADJUST_FACTOR    value <= 1 applied to the calculated intergenic scores  (default 1)


## hidden opts:
my $limit_range_lend;
my $limit_range_rend;


my $usage = <<__EOUSAGE__;

    ## Placeholder for additional usage info.


__EOUSAGE__


    ;



## Allowable evidence types
my %ALLOWABLE_EVIDENCE_CLASSES = ( 'PROTEIN' => 1,  # ie. AAT_nap or genewise data
                                   'TRANSCRIPT' => 1, # ie. AAT_gap2, sim4, or pasa data
                                   'ABINITIO_PREDICTION' => 1, # ie. genscan abinitio gene prediction
                                   'OTHER_PREDICTION' => 1, # ie. ensembl known genes (homology-based prediction). Count like protein alignments, but include well defined initial/terminal exons too!
                                   ); 
my %EV_TYPE_TO_EV_CLASS;



## Flags
my $SCORE_FORWARD_STRAND_FLAG = 1;
my $SCORE_REVERSE_STRAND_FLAG = 1;

## Input parameters
my $genomicSeqFile;
my $genePredictionsFile;
my $transcriptAlignmentsFile;
my $proteinAlignmentsFile;
my $repeatsFile;
my $weightsFile;
my $terminalExonsFile;
my $prediction_progs;
my $stitch_ends;
my $avg_intergenic_length;

my $FORWARD_STRAND_ONLY_FLAG;
my $REVERSE_STRAND_ONLY_FLAG;
my $SEE;
my $DEBUG;
my $extend_to_terminal; # extend termini of genewise alignments to find a start or stop codon, forming a terminal exon.
my $min_intron_length = 20; # smallest valid intron length I'm aware of (paramecium)
my $MIN_CODING_LENGTH = 300; # smallest coding gene length (sum exons) to be reported. (applied to nested and 'intergenic' genes only!).

my $MIN_CODING_NONCODING_SCORE_RATIO = 0.75; # predictions with score ratios less than this are eliminated. (use --report_ELM to still report them).

# stop codon issues
my $DEFAULT_STOP_CODONS = "TAA,TGA,TAG";
my $stop_codons_arg;
my @STOP_CODONS;
my $exec_dir = undef;

my $RECURSION_COUNT = 0;
my $RECURSE_FLAG = 0;

my $help_flag;
my $report_ELM_flag = 0;  #if set, the ELIMINATED EVM preds are reported too (those that score better as noncoding or are ultra short).

my $INTERGENIC_SCORE_ADJUST_FACTOR = 1;  # can change the percent of the actual intergenic score that's used in the DP scoring function.

my $MIN_LONG_INTRON_LENGTH = 0; # minimum intron length to go searching for a gene within an intron
my $MIN_INTERGENIC_SIZE_ON_RE_SEARCH = 10000;
my $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH = 0; ## off


my $MAX_NUM_PREV_EXONS_COMPARE = 500; # heuristic: arbitrary value, but hopefully long enough to give the same result as complete DP.


&GetOptions ("genome|G=s"=>\$genomicSeqFile,
             "gene_predictions|g=s"=>\$genePredictionsFile,
             "transcript_alignments|e=s"=>\$transcriptAlignmentsFile,
             "protein_alignments|p=s" => \$proteinAlignmentsFile,
             "repeats|r=s"=>\$repeatsFile,
             "weights|w=s"=>\$weightsFile,
             "forwardStrandOnly"=>\$FORWARD_STRAND_ONLY_FLAG,
             "reverseStrandOnly"=>\$REVERSE_STRAND_ONLY_FLAG,
             "terminalExonsFile|t=s"=>\$terminalExonsFile,
             "S"=>\$SEE,
             "debug" => \$DEBUG,
             "stop_codons=s" => \$stop_codons_arg,
             "help|h" => \$help_flag,
             "stitch_ends=s" => \$stitch_ends,
             "extend_to_terminal=s" => \$extend_to_terminal,
             "min_intron_length=i" => \$min_intron_length,
             "exec_dir=s" => \$exec_dir,
             "report_ELM" => \$report_ELM_flag,
             "INTERGENIC_SCORE_ADJUST_FACTOR=f" => \$INTERGENIC_SCORE_ADJUST_FACTOR,

             "search_long_introns=i" => \$MIN_LONG_INTRON_LENGTH,
             "re_search_intergenic=i" => \$MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH,
             "terminal_intergenic_re_search=i" => \$MIN_INTERGENIC_SIZE_ON_RE_SEARCH,
             
             # hidden opts
             "limit_range_lend=i" => \$limit_range_lend,
             "limit_range_rend=i" => \$limit_range_rend,

             'trellis_search_limit=i' => \$MAX_NUM_PREV_EXONS_COMPARE,
             
             
             );

if ($help_flag) {
    die "$param_string\n\n$usage\n";
}

## check required options.
unless ($genomicSeqFile && $genePredictionsFile && $weightsFile) {
    die $param_string;
}

if ($DEBUG) {
    open (STDERR, ">&STDOUT");
}



if ($MIN_LONG_INTRON_LENGTH) {
    $RECURSE_FLAG = 1;
}

## track prediction progs under consideration.

my %PREDICTION_PROGS; #populated on reading weights (all gene predictions).
my %PREDICTION_PROGS_CONTRIBUTE_INTERGENIC;  #only those prediction programs that score intergenic regions. OTHER_PREDICTION types do not count.  Only the ab initio preds should be included here.


if ($exec_dir) {
    chdir $exec_dir or die "Error, cannot cd to exec dir: $exec_dir ";
    print STDERR "-successfully changed directories to $exec_dir\n";
}


# process end-stitching info
my %END_STITCHING_PROG_NAMES;
if ($stitch_ends) {
    open (my $fh, $stitch_ends) or die "Error, cannot open file $stitch_ends";
    while (<$fh>) {
        chomp;
        my @x = split (/\s+/);
        my $ev_type = $x[1];
        if ($ev_type =~ /\w/) {
            $END_STITCHING_PROG_NAMES{$ev_type} = 1;
        }
    }
    close $fh;
    
    unless (%END_STITCHING_PROG_NAMES) {
        confess "Error, supposed to have stitch types, but none parsed from file $stitch_ends";
    }
    
    print STDERR "stitching ev_types: " . join (",", keys %END_STITCHING_PROG_NAMES) . "\n";
    
}

# process extend-to-termini info
my %EXTEND_TO_TERMINI_TYPES;
if ($extend_to_terminal) {
    open (my $fh, $extend_to_terminal) or die "Error, cannot open file $extend_to_terminal";
    while (<$fh>) {
        chomp;
        my @x = split (/\s+/);
        my ($ev_type) = $x[1];
        if ($ev_type) {
            $EXTEND_TO_TERMINI_TYPES{$ev_type}=1;
        }
    }
    close $fh;

    unless (%EXTEND_TO_TERMINI_TYPES) {
        confess "Error, supposed to have extend_to_terminal types, but none parsed from file $extend_to_terminal";
    }
    
    print STDERR "stitching ev_types: " . join (",", keys %EXTEND_TO_TERMINI_TYPES) . "\n";
    
}

# process stop codons:
if ($stop_codons_arg) {
    $stop_codons_arg =~ s/\s//g; #rid whitespace
    my @stops = split (/,/, $stop_codons_arg);
    foreach my $stop (@stops) {
        if (length($stop) != 3 || $DEFAULT_STOP_CODONS !~ /$stop/) {
            die "Error, stop codon ($stop) is not acceptable";
        }
    }
    @STOP_CODONS = @stops;
} else {
    @STOP_CODONS = split (/,/, $DEFAULT_STOP_CODONS);
}
unless (@STOP_CODONS) {
    die "Error, no stop codons in use.\n";
}


if ($FORWARD_STRAND_ONLY_FLAG && $REVERSE_STRAND_ONLY_FLAG) {
    die "Error, don't specify both forward strand only and reverse strand only. They are mutually exclusive.";
} elsif ($FORWARD_STRAND_ONLY_FLAG) {
    $SCORE_REVERSE_STRAND_FLAG = 0;
} elsif ($REVERSE_STRAND_ONLY_FLAG) {
    $SCORE_FORWARD_STRAND_FLAG = 0;
}


## Global variables:
my $genomicSeq;
my $genomic_seq_length;
my @GENOME_FEATURES; #array of genome positions, features set by constants below, default zero.

############################
#constants                 #
my $START = 1;             #
my $DONOR = 2;             #
my $ACCEPTOR = 3;          #
my $STOP = 4;              #
############################

my $MASK_EXON_MIN_PERCENT = 40; # minimum amount of an exon to overlap masked region to be removed from candidacy
my $INTRON_MEDIAN_FACTOR = 2;  ## protein-alignements with length at most this factor x median intron length decrement the coding region.

my %EVIDENCE_WEIGHTING; #populated from evidence weight file.

# coding score
my @CODING_SCORES; #each nucleotide incremented if covered by a match evidence.

## containers and lookups for exons
our @EXONS; #holds exon objects
my %EXONS_VIA_COORDS;


## other configurable parameters:
my $START_STOP_RANGE = 500;
my $CHAIN_TERMINI_WINDOW = 250;
my $START_STOP_MIN_EVIDENCE_THRESHOLD = 0; # at least this many pieces of evidence supporting start or stop


## other useful global vars:
my %ACCEPTABLE_EXON_LINKAGES;
my %PHASED_CONNECTIONS;
my %INTERGENIC_CONNECTIONS;
my %GFF_PHASE_CONVERSIONS;  #here, we use phases [123456], in gff we have [012] and orient [+-]

## introns
my %INTRONS_TO_SCORE;  #holds intron end5_end3 and score of intron.

my @FORWARD_PRED_INTRON_VEC;  # holds individual base pair scores for introns.
my @REVERSE_PRED_INTRON_VEC;


## repeat masking: no repeat masked region should contribute to score. Entries derived from repeat annotations or N's in sequence.
my @MASK;

# intergenic regions
my @INTERGENIC_SCORES;


# tracking start/end peaks supported by evidence:
my @START_PEAKS;
my @END_PEAKS;
my @begins;
my @ends;


# misc
my $mock_exon_counter = 0;
my @INTRON_FORWARD_CONSENSUS_SCORES; # for calculating offsets when rescore introns as intergenic 
my @INTRON_REVERSE_CONSENSUS_SCORES;
my %intergenic_cache; # for debugging purposes, extract the scores calc'd for intergenics.

my %PREDICTED_INTRONS;
my %INTRONS_TO_EVIDENCE;  # retains list of [ [acc, ev_type], ...] keyed on intron key.

my $SUM_GENEPRED_WEIGHTS = 0; # maximum protein coding by genefinders alone.

my $MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON = 30; #used for estimating intron counts within est and protein alignment chains, NOT a restriction on perfect intron lengths.  Alignment gaps smaller than this length are assumed to be simple alignment gaps and not introns when evaluating alignment chains to decrement intergenic scores.

my $AUGMENT_INTERGENIC_FROM_START_STOP_PEAKS_FLAG = 1;
my $TRIM_PARTIAL_GENES_TO_INITIAL_EXONS_FLAG = 1; # for 5' partials, alternative exons are chosen to provide an initial exon


my $PERCENT_EXON_START_EXPLORE = 30; # used with --extend_to_terminal option, exploring the first 30% of an initial segment for a start codon if it cannot be extend upstream to a start codon.
my $PERCENT_EXON_STOP_EXPLORE = 200; # as above, but for terminal segments extending to a stop codon.


my @LONG_INTRONS;  ## store them for processing later.

## private vars
my $_TRELLIS_ROUND = 0;

my $RESET_INTERGENIC_REGIONS_FLAG = 0; # reset intergenic regions at first filtering of low quality predictions. Undoes the evidence-based adjustments to intergenic scores after the first trellis build.


main: {
     
    &readWeights();
    foreach my $pred_type (keys %PREDICTION_PROGS) {
        $SUM_GENEPRED_WEIGHTS += $EVIDENCE_WEIGHTING{$pred_type};
    }
    
    ## read in the genomic sequence:
    print STDERR "\n-reading genomic sequence.";
    &read_genomic_sequence();
        
    ## populate any required data structures.
    &initialization();
    
    &repeatMask('+');  # both forward and reverse strand repeats are masked.

    
    
    if ($SCORE_FORWARD_STRAND_FLAG) {
        ## do forward strand
        &process_features('+');
    }
    
    
     ## cache features.
    my @cached_forward_exons = @EXONS;
    
    
    ## prepare for reverse strand analysis
    &initialization(); #to clear globals
    
    if ($SCORE_REVERSE_STRAND_FLAG) {
        ## do reverse strand analysis
                
        ## reverse complement the sequence
        $genomicSeq = &reverse_complement($genomicSeq);
        
        &repeatMask('-');
        &process_features('-');  ## @EXONS is populated here, now.
        
        ## bring everything back to the forward strand:

        $genomicSeq = &reverse_complement($genomicSeq);
        &transpose_exons_back_to_forward_strand(); ## operates on the @EXONS list
        
        &repeatMask('+');  # both forward and reverse strand repeats are masked. 
    
    }
    
    ## add the cached exons from the forward strand run back to the list.
    if (@cached_forward_exons) {
        push (@EXONS, @cached_forward_exons);  
    }
    @cached_forward_exons = (); #clear
    

    ## populate the intron score vectors (used during the filtering of preds w/ low support).
    &populate_forward_reverse_pred_intron_vectors();
    
    if ($DEBUG) {
        &write_introns_and_exons_to_file();
        &write_start_and_stop_peaks_to_file();
        &write_intron_vectors_to_file();
    }
        
    ## mask repeat-containing exons, removing them as candidates.
    if ($repeatsFile) {
        if ($DEBUG) {
            &dump_repeat_mask();
        }
        &mask_exons();
    }
    
    &populate_intergenic_regions(); # mask applied to regions.
        
    # &decrement_intergenic_for_evidence_spans(); # only applies to the first trellis run.
    
    if ($AUGMENT_INTERGENIC_FROM_START_STOP_PEAKS_FLAG) {
        &augment_intergenic_from_start_stop_peaks();
    }
    
    ## Do DP to find the highest scoring path thru all exons, introns, and intergenic regions.
    if (defined $limit_range_lend && defined $limit_range_rend) {
        &generate_consensus_gene_predictions ($limit_range_lend, $limit_range_rend, "STANDARD");
    }
    else {
        ## whole sequence processed, as default
        &generate_consensus_gene_predictions(1, $genomic_seq_length, "STANDARD"); # set bounds to full sequence 
    }
    
    ## Run now just within the long introns.  Also, recreate the INTERGENIC scores w/o decrementing for evidence spans.
    if (($RECURSE_FLAG) && @LONG_INTRONS) {
        
        # recreate the intergenic vector:
        &populate_intergenic_regions(); # mask applied to regions.
        
        foreach my $long_intron (@LONG_INTRONS) {
            my ($intron_lend, $intron_rend) = @$long_intron;
            print STDERR "-running in long intron: $intron_lend-$intron_rend\n";
            &generate_consensus_gene_predictions($intron_lend, $intron_rend, "INTRON");
        }
    }
    
    exit(0);
}

####
sub write_start_and_stop_peaks_to_file {
    &write_peaks_to_file("start_peaks", \@START_PEAKS);
    &write_peaks_to_file("end_peaks", \@END_PEAKS);
    
    return;
}


####
sub generate_consensus_gene_predictions {
    my ($range_lend, $range_rend, $mode) = @_;
    
    print "## SEARCHING ($RECURSION_COUNT) $range_lend-$range_rend (mode: $mode)\n" if $DEBUG;
    
    $RECURSION_COUNT++;

    my $region_length = $range_rend - $range_lend + 1; 
    if ($region_length < $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH) { 
        $RECURSION_COUNT--;
        return;
    }
        
    print "** generate_consensus_gene_predictions ($range_lend, $range_rend) (mode: $mode)\n" if $DEBUG;

    my @exons_within_range;
    
    foreach my $exon (@EXONS) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();

        if ($lend >= $range_lend && $rend <= $range_rend) {
            push (@exons_within_range, $exon);
        }
    }

    unless (@exons_within_range) {
        print STDERR "\t-no exons in range $range_lend - $range_rend\n";
        $RECURSION_COUNT--;
        return;
    }
    
    print STDERR "\n-building DP trellis";
    my $top_scoring_exon = &build_trellis($range_lend, $range_rend, \@exons_within_range);
    
    my @EVM_predictions = &traverse_path($top_scoring_exon);
    print STDERR "-found " . scalar (@EVM_predictions) . " predictions from path traversal\n";
    
    unless (@EVM_predictions) {
        ## nothing to report
        print STDERR "\t-no predictions reported within region $range_lend - $range_rend\n";
        $RECURSION_COUNT--;
        return;
    }
    
    @EVM_predictions = &filter_predictions_low_support(\@EVM_predictions, $mode);
    
    
    unless (@EVM_predictions) {
        ## must have been filtered
        print STDERR "\t-no predictions after filtering those with low support\n";
        $RECURSION_COUNT--;
        return;
    }
    

     if ($TRIM_PARTIAL_GENES_TO_INITIAL_EXONS_FLAG) {
         &convert_5prime_partials_to_complete_genes_where_possible(@EVM_predictions);
     }

    ## report the final predictions:
        
    print STDERR "\n-reporting highest scoring path thru the trellis.\n";
    
    my ($predictions_span_lend, $predictions_span_rend) = &get_range_covered_by_predictions(\@EVM_predictions);
    
	my $preds_remain_after_filter = 0;
    foreach my $prediction (@EVM_predictions) {
		unless ($prediction->is_eliminated()) {
			$preds_remain_after_filter++;
		}
	}
	
    print STDERR "-found $preds_remain_after_filter predictions persist after filtering for low support\n";
	print "!! Predictions spanning range $predictions_span_lend - $predictions_span_rend [R$RECURSION_COUNT]\n" if ($preds_remain_after_filter || $report_ELM_flag);
	foreach my $pred (@EVM_predictions) {
		if ($pred->is_eliminated() && ! $report_ELM_flag) { next; } # not reporting the eliminated ones.
		
		$pred->{mode} = $mode;
		print $pred->toString() . "\n";
	}
    
    if ($mode ne 'INTRON') {
        ## examine the tail of predictions.
        ## this is only important when some very high-scoring partial gene internally located in the contig
        ## prevents trellis extension.
        
        ## left tail
        my ($left_range_lend, $left_range_rend) = ($range_lend, $predictions_span_lend -1);
        my $region_len = $left_range_rend - $left_range_lend;
        if ($region_len >= $MIN_INTERGENIC_SIZE_ON_RE_SEARCH) {
            &generate_consensus_gene_predictions($left_range_lend, $left_range_rend, $mode);
        }
        
        ## right tail
        my ($right_range_lend, $right_range_rend) = ($predictions_span_rend + 1, $range_rend);
        $region_len = $right_range_rend - $right_range_lend;
        if ($region_len >= $MIN_INTERGENIC_SIZE_ON_RE_SEARCH) {
            &generate_consensus_gene_predictions($right_range_lend, $right_range_rend, $mode);
        }
    }
    

	## Identify Long introns and store them for later:
	if ( ($RECURSE_FLAG) && $mode ne 'INTRON') {
		## pursue long introns:
		my @long_introns = &get_long_introns(\@EVM_predictions);
		if (@long_introns) {
			print STDERR "storing long introns: " . Dumper (\@long_introns);
			push (@LONG_INTRONS, @long_introns); # store for later.
		}
	}
	
    
    if ($mode ne 'INTRON' && $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH) { 
        	
        ## pursue the intergenic regions again, just in case.  
        my @intergenic_regions = &get_intergenic_regions(\@EVM_predictions);
        if ($DEBUG) {
            print "Got Intergenic regions: " . Dumper(\@intergenic_regions);
        }
        
        foreach my $intergenic_region (@intergenic_regions) {
            my ($intergenic_lend, $intergenic_rend) = @$intergenic_region;
            print STDERR "+++++++++++ searching intergenic region[$RECURSION_COUNT]: $intergenic_lend-$intergenic_rend\n";
            &generate_consensus_gene_predictions($intergenic_lend, $intergenic_rend, $mode);
        }
        
        
        unless (defined $predictions_span_lend && defined $predictions_span_rend) { return; }
        
        ## Recurse on unpredicted regions:
        
        my ($left_side_lend, $left_side_rend) = ($range_lend, $predictions_span_lend - 1);
        
        my ($right_side_lend, $right_side_rend) = ($predictions_span_rend + 1, $range_rend);
        
        my $new_range_lend_length = $left_side_rend - $left_side_lend + 1;
        
        if ($new_range_lend_length >= $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH) {
            print STDERR "-search left intergenic[$RECURSION_COUNT]: $left_side_lend-$left_side_rend\n";
            &generate_consensus_gene_predictions($left_side_lend, $left_side_rend, $mode);
        }
        
        my $new_range_rend_length = $right_side_rend - $right_side_lend + 1;
        
        if ($new_range_rend_length >= $MIN_GENE_LENGTH_SIZE_ON_RE_SEARCH) {
            print STDERR "-search right intergenic[$RECURSION_COUNT]: $right_side_lend-$right_side_rend\n"; 
            &generate_consensus_gene_predictions($right_side_lend, $right_side_rend, $mode);
        }

    }

        
    $RECURSION_COUNT--;
    return;
}    


####
sub get_range_covered_by_predictions {
    my ($predictions_aref) = @_;
    
    my $min_lend = undef;
    my $max_rend = undef;
    foreach my $prediction (@$predictions_aref) {
        my ($lend, $rend) = sort {$a<=>$b} $prediction->get_span();
        if (!defined $min_lend) {
            $min_lend = $lend;
            $max_rend = $rend;
        }
        else {
            if ($lend < $min_lend) {
                $min_lend = $lend;
            }
            if ($rend > $max_rend) {
                $max_rend = $rend;
            }
        }
    }

    unless (defined $min_lend && defined $max_rend) {
        die "Error, didn't parse min_lend && max_rend from predictions.\n";
    }
    
    return ($min_lend, $max_rend);

}

####
sub write_peaks_to_file {
    my ($filename, $peaks_aref) = @_;

    open (my $fh, ">$filename") or die "Error, cannot write to file $filename";
    @$peaks_aref = sort {$a->{position}<=>$b->{position}} @$peaks_aref;
    
    foreach my $peak (@$peaks_aref) {
        my ($position, $score, $strand) = ($peak->{position}, $peak->{score}, $peak->{strand});
        if ($strand eq '-') {
            $score *= -1; # position score towards bottom strand in graph
        }
        print $fh "$position\t$score\n";
    }
    close $fh;

    return;
}


####
sub process_features {
    
    my $genomic_strand = shift;
    print STDERR "\n\n** Processing features ($genomic_strand)\n";

    ## get splice site info:
    print STDERR "\n-finding all potential splice sites, looking for GT/GC donors and AG acceptors.";
    &populate_splice_sites();
    
    ## get the starts and stops:
    print STDERR "\n-finding all potential starts and stops. (strand: $genomic_strand)";
    &populate_starts_and_stops();
    
    ## load predictions:
    print STDERR "\n-loading gene prediction evidence. (strand: $genomic_strand)";
    &load_prediction_data($genomic_strand);
    	
    if ($terminalExonsFile) {
        print STDERR "\n-loading the terminal exon supplement (strand: $genomic_strand)";
        &supplement_terminal_exons($genomic_strand);
        # all pasa-based terminal exons are treated as high confidence starting and ending points.
        &include_pasa_terminal_exons_in_gene_boundary_estimation($terminalExonsFile, $genomic_strand, \@begins, \@ends);
    }
    
    ## load the search evidence
    print STDERR "\n-loading the search evidence (protein and EST search results, strand: $genomic_strand)";
    
    if ($transcriptAlignmentsFile) {
        ## IMPORTANT!!! Load transcript alignments before protein alignments, since start/stop applications are done during protein phase.
        &load_evidence_data($genomic_strand, $transcriptAlignmentsFile);
    }
    
    if ($proteinAlignmentsFile) {
        &load_evidence_data($genomic_strand, $proteinAlignmentsFile);
        &decrement_coding_using_protein_alignment_introns($genomic_strand);
    }

    &analyze_gene_boundaries($genomic_strand, \@begins, \@ends);
    
    print STDERR "\n-scoring exons (strand: $genomic_strand)";
    &score_exons($genomic_strand); 
    

    return;

}


####
sub dump_prediction {
    my @coordsets = @_;
    
    foreach my $coordset (@coordsets) {
        my ($end5, $end3) = @$coordset;
        
        my $potential_start = substr ($genomicSeq, $end5-1, 3);
        my $acceptor = substr ($genomicSeq, $end5-2-1, 2);
        my $donor = substr ($genomicSeq, $end3, 2);
        my $potential_stop = substr ($genomicSeq, $end3-2, 3);
        print STDERR "($potential_start)\t $acceptor\t$end5-$end3\t$donor\t$potential_stop\n";
    }

    return;
}


####
sub dump_exons {
    @EXONS = sort {$a->{end5} <=> $b->{end5}} @EXONS;
    
    foreach my $exon (@EXONS) {
        my $base_score = $exon->{base_score};
        my $length = $exon->length();
        my $score_per_base = sprintf ("%.2f", $base_score / $length);

        print $exon->toString() . " base score: $base_score, score_per_base: $score_per_base\n";
    }

    return;
}


####
sub read_genomic_sequence {
    open (FILE, $genomicSeqFile) or die "Error, cannot open $genomicSeqFile\n\n";
    while (<FILE>) {
        if (/^>/) { next;} #ingore header:
        s/\s+//g;
        $genomicSeq .= $_;
    }
    close FILE;
    $genomicSeq =~ s/\s//g;
    $genomicSeq = uc $genomicSeq;
    
    $genomic_seq_length = length ($genomicSeq);
    
    return;
}


####
sub populate_splice_sites() {
    
    # positions are stored in GENOME_FEATURES starting at position 1 (instead of zero).

    ## get the donors
    my @gt_positions = &find_all_positions('GT');
    my @gc_positions = &find_all_positions('GC');
    foreach my $pos (@gt_positions, @gc_positions) {
        $GENOME_FEATURES[$pos+1] = $DONOR;
    }
    
    ## get the acceptors
    my @ag_positions = &find_all_positions('AG');
    foreach my $pos (@ag_positions) {
        $GENOME_FEATURES[$pos+1] = $ACCEPTOR;
    }
    
    return;

}


####
sub populate_starts_and_stops {
    
    ## get the starts
    my @start_pos = &find_all_positions("ATG");
    foreach my $start (@start_pos) {
        $GENOME_FEATURES[$start+1] = $START;
    }
    
    ## get the stops
    my @stop_pos;
    
    foreach my $stop (@STOP_CODONS) {
        push (@stop_pos, &find_all_positions($stop));
    }
    
    foreach my $stop (@stop_pos) {
        $GENOME_FEATURES[$stop+1] = $STOP;
    }

    return;
}


####
sub find_all_positions {
    my $string = shift;
    
    my @positions;
    
    my $start = 0;
    while ($start != -1) {
        $start = index($genomicSeq, $string, $start);
        if ($start != -1) {
            push (@positions, $start);
            $start+= 1; #begin search right after found position.
        }
    }
    
    return (@positions);
}



####
sub load_prediction_data {
    
    my $genomic_strand = shift;
    
    my %model_to_exon_coord_list;
    my %model_to_evidence_weight;
    my %model_to_ev_class;
    my %model_to_ev_type;

    open (FILE, $genePredictionsFile) or die "Error, cannot open file $genePredictionsFile";
    while (<FILE>) {
        my $line = $_;
        chomp;
        unless (/\w/) { next;}
        if (/^\#/) { next;}
        
        my @x = split (/\t/);
        my ($feat_type, $predType, $id, $lend, $rend, $orient, $phase) = ($x[2], $x[1], $x[8],$x[3], $x[4], $x[6], $x[7]);
        unless ($feat_type eq 'CDS') { next;}
        
        my $ev_class = $EV_TYPE_TO_EV_CLASS{$predType};
        
        unless (&in_range_of_genomic_sequence($lend, $rend) ) {
            print STDERR "-WARNING, IGNORING line: $line\tBECAUSE NOT IN RANGE OF GENOMIC SEQUENCE (1-$genomic_seq_length)\n";
            next;
        }
                

        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

        # only consider those predictions specified as prediction progs:
        unless (exists $PREDICTION_PROGS{$predType}) {
            print STDERR "-WARNING, skipping predType: $predType, since not specified in weights file.\n";
            next;
        }
        
        
        my $evidence_weight = $EVIDENCE_WEIGHTING{$predType};
        
        if (! defined ($evidence_weight)) {
            confess "Error, no weight defined for prediction type: $predType ";  # only processing those entries defined in the weights file.
        }
        
        unless ($orient eq $genomic_strand) {
            next;
        }
        
        if ($genomic_strand eq '-') {
            ($end5, $end3) = &revcomp_coordinates($end5, $end3);
        }
        
		$id =~ /Parent=([^;\s]+)/ or confess "Error, cannot parse parent feature from $id";
		$id = $1;
		
        my $model_id =  $predType . "_" . $id;
        $model_to_evidence_weight{$model_id} = $evidence_weight;
        $model_to_ev_class{$model_id} = $ev_class;
        $model_to_ev_type{$model_id} = $predType;

        my $list_ref = $model_to_exon_coord_list{$model_id};
        unless (ref $list_ref) {
            $list_ref = $model_to_exon_coord_list{$model_id} = [];
        }
        
        push (@$list_ref, [$end5, $end3]);
                
    }
    close FILE;
    

    ## Go thru and assign exons
    foreach my $model_id (keys %model_to_exon_coord_list) {
        
        print "\n\nInputting $model_id\n" if $SEE;

        my $evidence_weight = $model_to_evidence_weight{$model_id};
        unless (defined $evidence_weight) {
            confess "Error, no weight assigned to $model_id\n";
        }
        
        my $ev_class = $model_to_ev_class{$model_id};
        my $ev_type = $model_to_ev_type{$model_id};
        
        my @coordsets = @{$model_to_exon_coord_list{$model_id}};
        
        @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;
        
        my $found_single = 0;
        
        ## process single-exon gene
        if (scalar(@coordsets) == 1) {
            my ($end5, $end3) = @{$coordsets[0]};
        
            if ($GENOME_FEATURES[$end5] == $START && $GENOME_FEATURES[ $end3-2 ] == $STOP) {
                
                &add_exon ($model_id, $end5, $end3, "single", 1, $ev_type);
                &add_match_coverage($end5, $end3, $evidence_weight, $model_id, $ev_type);
                $found_single = 1;
            } 
        }
        
        
        if (! $found_single) {
            
            ## process multi-exon genes and single-partials.
            
            my $valid_prediction = 1;
            my $cds_length = 0;
            my @classified_coordsets;
            
            ## check each exon
            
          
            for (my $i=0; $i <= $#coordsets; $i++) {
                my ($end5, $end3) = @{$coordsets[$i]};
                
                my $exon_length = abs ($end3 - $end5) + 1;

                &add_match_coverage($end5, $end3, $evidence_weight, $model_id, $ev_type);
                
                my $type = undef;
                
                print "InputPrediction end5: $end5, end3: $end3\t" if $SEE;
                my $has_start = $GENOME_FEATURES[$end5] == $START;
                my $has_stop =  $GENOME_FEATURES[$end3-2] == $STOP;
                my $has_acceptor = $GENOME_FEATURES[$end5-2] == $ACCEPTOR;
                my $has_donor = $GENOME_FEATURES[$end3+1] == $DONOR;

                print "(i:$i, ex_len($exon_length)) start($has_start), acceptor($has_acceptor), donor($has_donor), stop($has_stop)\t" if $SEE;
                
                
                if ($i==0 && $GENOME_FEATURES[$end5] == $START && $GENOME_FEATURES[$end3+1] == $DONOR) {
                    ## don't allow splitting of start codons:
                    if ($exon_length >= 3) {
                        $type = "initial";
                    }
                } 
                
                elsif ($i == $#coordsets && $GENOME_FEATURES[$end5-2] == $ACCEPTOR && $GENOME_FEATURES[$end3-2] == $STOP) {
                    # don't split stop codons.
                    if ($exon_length >= 3) {
                        $type = "terminal";
                    }
                }
                
                elsif ($i != 0 && $i != $#coordsets && $GENOME_FEATURES[$end5-2] == $ACCEPTOR && $GENOME_FEATURES[$end3+1] == $DONOR) {
                    $type = "internal";
                }
                
                print "type: $type\n" if $SEE;
                
                if (! $type) {
                    print STDERR "Couldn't classify $end5, $end3\n";
                    ### Should try to recover partials later
                    $valid_prediction = 0;
                }
                
                my $start_frame = $cds_length % 3 + 1;
                $cds_length += abs ($end3 - $end5) + 1;
                push (@classified_coordsets, [$end5, $end3, $type, $start_frame]);
				
            }
            print "\n" if $SEE;
            if ($valid_prediction) {
                print STDERR "Valid prediction for $model_id\n";
                my @coords;
                foreach my $classified_coordset (@classified_coordsets) {
                    my ($end5, $end3, $type, $start_frame) = @$classified_coordset;
                    &add_exon($model_id, $end5, $end3, $type, $start_frame, $ev_type);
                    push (@coords, [$end5, $end3]);
                }
                # extract intron info based on exon coordinates:
                &add_introns($model_id, \@coordsets, $genomic_strand, $evidence_weight, $ev_type);
                
            } else {
                ## Recover what we can from partials, or proper exons within otherwise corrupt gene predictions.
                warn "Sorry, prediction $model_id " . #Dumper (\@coordsets) . 
                    "fails validation.\n";
                &dump_prediction(@coordsets);
                &try_recover_partial_prediction(\@coordsets, $model_id, $evidence_weight, $genomic_strand, $ev_type);
            } 
        }
    }
    
    return;
}


####
sub try_recover_partial_prediction {
    my ($coordsets_aref, $model_id, $evidence_weight, $genomic_strand, $ev_type) = @_;
    
    my @candidate_initials;
    my @candidate_internals;
    my @candidate_terminals;
    
    ## Simply see what exon category each coordset could fit:
    
    foreach my $coordset (@$coordsets_aref) {
        my ($end5, $end3) = @$coordset;
        my $exon_length = abs ($end3-$end5) + 1;
        
        if ($GENOME_FEATURES[$end5] == $START && $GENOME_FEATURES[$end3+1] == $DONOR) {
            # don't split start codons
            if ($exon_length >= 3) {
                push (@candidate_initials, $coordset);
            }
        } 
    
        if ($GENOME_FEATURES[$end5-2] == $ACCEPTOR && $GENOME_FEATURES[$end3-2] == $STOP) {
            # don't split stop codons:
            if ($exon_length >= 3) {
                push (@candidate_terminals, $coordset);
            }
        }
        
        if ($GENOME_FEATURES[$end5-2] == $ACCEPTOR && $GENOME_FEATURES[$end3+1] == $DONOR) {
            push (@candidate_internals, $coordset);
        }
    }

    ## load candidate initials:
    foreach my $initial_exon_coordset (@candidate_initials) {
        my ($end5, $end3) = @$initial_exon_coordset;
        &add_exon($model_id, $end5, $end3, "initial", 1, $ev_type);
        print STDERR "-recovered $model_id, initial, $end5, $end3\n";
    }

    ## load candidate terminals:
    foreach my $terminal_exon_coordset (@candidate_terminals) {
        my ($end5, $end3) = @$terminal_exon_coordset;
        my $length = $end3 - $end5 + 1;
        my $num_prev_codon_chars = $length % 3;
        my $phase;
        if ($num_prev_codon_chars == 0) {
            $phase = 1;
        } 
        elsif ($num_prev_codon_chars == 1) {
            $phase = 3;
        } 
        elsif ($num_prev_codon_chars == 2) {
            $phase = 2;
        }

        &add_exon($model_id, $end5, $end3, "terminal", $phase, $ev_type);
        print STDERR "-recovered $model_id, terminal, $end5, $end3\n";
    }

    ## load internals
    foreach my $internal_exon_coordset (@candidate_internals) {
        my ($end5, $end3) = @$internal_exon_coordset;
        # treat like evidence exon, add it any way we can
        &add_internal_evidence_exon($model_id, $end5, $end3, $ev_type);
        print STDERR "-recovered $model_id, internal, $end5, $end3\n";
    }
    
    ## add in the introns:
    &add_introns($model_id, $coordsets_aref, $genomic_strand, $evidence_weight, $ev_type);

    return;
}


####
sub add_introns {
    my ($accession, $coords_list_aref, $genomic_strand, $weight, $intron_type) = @_;
    
    ## intron type must be one of the following: prediction, protein, est
    my $intron_ev_class = $EV_TYPE_TO_EV_CLASS{$intron_type};
    
    unless ($ALLOWABLE_EVIDENCE_CLASSES{$intron_ev_class}) {
        confess "Error, don't recognize intron_type: $intron_ev_class from type $intron_type\n";
    }
    
    ## prediction flag indicates that the intron is from a prediction.
    
    my @coords_list = sort {$a->[0]<=>$b->[0]} @$coords_list_aref;
    
    for (my $i=0; $i < $#coords_list; $i++) {
        my $first_end3 = $coords_list[$i]->[1];
        my $next_end5 = $coords_list[$i+1]->[0];
        
        if ($next_end5 < $first_end3) {
            warn "ERROR adding intron for $accession: next end5 $next_end5 < this end3 $first_end3 " . #Dumper (\@coords_list) . 
                " mode $genomic_strand\n";
            next;
        }
        
        if ($DEBUG) {
            my $intron_bound_end5 = $first_end3 + 1;
            my $intron_bound_end3 = $next_end5 - 1;
            if ($genomic_strand eq '-') {
                ## revcomp them:
                ($intron_bound_end5, $intron_bound_end3) = &revcomp_coordinates($intron_bound_end5, $intron_bound_end3);
            }
            print "ADDING_INTRON\t$intron_bound_end5-$intron_bound_end3\t$accession\tweight:$weight\n";
        }
        

        my $potential_donor = $first_end3+1;
        my $potential_acceptor = $next_end5-2;
        
        my $intron_length = abs ($potential_acceptor - $potential_donor) + 1;
        if ($intron_length < $min_intron_length) {
            warn "Sorry, intron length ($intron_length) is less than minimum required ($min_intron_length)\n";
            next;
        }
        
        if ($GENOME_FEATURES[$potential_donor]==$DONOR && $GENOME_FEATURES[$potential_acceptor]==$ACCEPTOR) {
            # got an intron.
            
            
            # score the intron array:
            my $intron_score = 0;
            for (my $i = $potential_donor; $i <= $potential_acceptor+1; $i++) {
                unless ($MASK[$i]) {
                    $intron_score += $weight;
                }
            }
            
            my ($intron_end5, $intron_end3) = ($potential_donor, $potential_acceptor);
            
            # Introns coordinate pairs always stored based on forward reference coordinates.
            if ($genomic_strand eq '-') {
                ($intron_end5, $intron_end3) = &revcomp_coordinates($intron_end5, $intron_end3);
            }
            
            # store entry for coordinate pair
            my $intron_key = "$intron_end5" . "_" . "$intron_end3";
            $INTRONS_TO_SCORE{$intron_key} += $intron_score;
        
            ## add evidence to the intron key:
            my $evidence_list_aref = $INTRONS_TO_EVIDENCE{$intron_key};
            unless (ref $evidence_list_aref) {
                $evidence_list_aref = $INTRONS_TO_EVIDENCE{$intron_key} = [];
            }
            push (@$evidence_list_aref, [ $accession, $intron_type ] );
            
            if ($intron_ev_class eq "ABINITIO_PREDICTION") {
                $PREDICTED_INTRONS{$intron_key} += $intron_score;
            }
        }
    }
    
    return;
}


####
sub load_evidence_data {
    
    my $genomic_strand = shift;
    my $evidence_filename = shift;

    ## ToDo: decouple start/end analysis from evidence-based exon instantiation.
    

    my @evidence_chains = &parse_evidence_chains($genomic_strand, $evidence_filename);
    &instantiate_evidence_based_exons(\@evidence_chains, \@begins, \@ends, $genomic_strand, $evidence_filename);
    
    return;
}


####
sub analyze_gene_boundaries {
    my ($genomic_strand, $begins_aref, $ends_aref) = @_;

    print "Analyzing gene boundaries (strand: $genomic_strand)\n" if $DEBUG;

    my @startPeaks = &analyze_peaks($begins_aref);
    my @endPeaks = &analyze_peaks($ends_aref);
    if ($DEBUG) {
        &dump_begins_and_ends($genomic_strand, $begins_aref, $ends_aref);
    }
    
    ## add peak info to global vars and include strand:
    foreach my $start_peak (@startPeaks) {
        my ($position, $score) = ($start_peak->{position}, $start_peak->{score});
        if ($genomic_strand eq '-') {
            ($position) = &revcomp_coordinates($position);
        }
        push (@START_PEAKS, { position => $position,
                              score => $score,
                              strand => $genomic_strand } );
    }
    
    foreach my $end_peak (@endPeaks) {
        my ($position, $score) = ($end_peak->{position}, $end_peak->{score});
        if ($genomic_strand eq '-') {
            ($position) = &revcomp_coordinates($position);
        }
        push (@END_PEAKS, { position => $position,
                            score => $score,
                            strand => $genomic_strand } );
    }
    
    return;
    
}


####
sub dump_begins_and_ends {
    my ($genomic_strand, $begins_aref, $ends_aref) = @_;
    
    open (my $fh, ">begins.$genomic_strand.coords") or die $!;
    for (my $i = 1; $i <= $#$begins_aref; $i++) {
        if (my $val = $begins_aref->[$i]) {
            print $fh "$i\t$val\n";
        }
    }
    close $fh;

    open ($fh, ">ends.$genomic_strand.coords") or die $!;
    for (my $i = 1; $i <= $#$ends_aref; $i++) {
        if (my $val = $ends_aref->[$i]) {
            my ($rev_coord) = &revcomp_coordinates($i);
            print $fh "$rev_coord\t$val\n";
        }
    }
    close $fh;

    return;
}


####
sub include_pasa_terminal_exons_in_gene_boundary_estimation {
    my ($term_exons_file, $genomic_strand, $begins_aref, $ends_aref) = @_;

    open (my $fh, $term_exons_file) or die "Error, cannot open file $term_exons_file";
    while (<$fh>) {
        chomp;
        my $line = $_;
        unless (/\w/) { next;}
        if (/^\#/) { next;}
        my @x = split (/\t/);
        my ($exonType, $ev_type, $acc, $lend, $rend, $orient, $phase) = ($x[2], $x[1], $x[8],$x[3], $x[4], $x[6], $x[7]);
        unless ($orient eq $genomic_strand) { next;}
        my $evidence_weight = $EVIDENCE_WEIGHTING{$ev_type};
        unless (defined ($evidence_weight)) {
            next; 
        }
        
        unless (&in_range_of_genomic_sequence($lend, $rend) ) {
            print STDERR "-WARNING, IGNORING line: $line\tBECAUSE NOT IN RANGE OF GENOMIC SEQUENCE (1-$genomic_seq_length)\n";
            next;
        }


        if ($orient eq '-') {
            ($lend, $rend) = sort {$a<=>$b} &revcomp_coordinates($lend, $rend);
        }
            
        if ($exonType eq 'initial') {
            $begins_aref->[$lend] += $evidence_weight;
        }
        elsif ($exonType eq 'terminal') {
            $ends_aref->[$rend] += $evidence_weight;
        }
        elsif ($exonType eq 'single') {
            $begins_aref->[$lend] += $evidence_weight;
            $ends_aref->[$rend] += $evidence_weight;
        }
        else {
            confess "Error, cannot parse exon type $exonType";
        }
    }
    close $fh;
    
    return;
    
}


####
sub instantiate_evidence_based_exons {
    my ($evidence_chains_aref, $begins_aref, $ends_aref, $genomic_strand, $evidence_filename) = @_;
    
    foreach my $chain (@$evidence_chains_aref) {
		        
        my $accession = $chain->{accession};
        
        if (my $target = $chain->{target}) {
            ## include target name in the accession as stored in the various data structures below.
            # this is very useful in the evm output.
            $accession .= "/Target=$target";
        }
        
        my $ev_type = $chain->{ev_type};
        my $ev_class = $chain->{ev_class};
        my $evidence_weight = $EVIDENCE_WEIGHTING{$ev_type};
        unless (defined $evidence_weight) {
            confess "Error, no weight for $ev_type\n";
        }
        
        ## increment support for beginning and end of gene
        if ($ev_class eq "PROTEIN") {
            $begins_aref->[$chain->{lend}] += $evidence_weight;
            $ends_aref->[$chain->{rend}] += $evidence_weight;
        }
        
        my @links = sort {$a->[0]<=>$b->[0]} @{$chain->{links}};
        
        my $link_count = 0;
        my $num_links = scalar @links;
        foreach my $link (@links) {
            my ($end5, $end3) = @$link;
                      
            $link_count++;
                       

            ## Check for good exon boundaries and add exon as needed.
            my $got_donor = 0;
            my $got_acceptor = 0;
            
            ## Not adding terminal segments as internal exons even if they end at splice boundaries.  These are sometimes artifacts.

            if ($GENOME_FEATURES[$end5-2] == $ACCEPTOR && $link_count != 1 && $link_count != $num_links) {
                $got_acceptor = 1;
            }
            
            
            if ($GENOME_FEATURES[$end3+1] == $DONOR && $link_count != 1 && $link_count != $num_links) {
                $got_donor = 1;
            }
            
            if ($got_donor && $got_acceptor) {
                my $found_exon_orf = &add_internal_evidence_exon ($accession, $end5, $end3, $ev_type);
                if ($found_exon_orf && $ev_class eq "TRANSCRIPT") {
                    ## allow ESTs to contribute to coding score if they provide internal exons with orfs.
                    &add_match_coverage($end5, $end3, $evidence_weight, $accession, $ev_type); #### Set weight to evidence!!!
                }

            }
            
            if ($ev_class eq "PROTEIN") {
                ## add range to coverage incrementer
                ## proteins count regardless of proper splice boundaries.
                &add_match_coverage($end5, $end3, $evidence_weight, $accession, $ev_type); #### Set weight to evidence!!!
            }
            
        }
        &add_introns($accession, \@links, $genomic_strand, $evidence_weight, $ev_type);
        
        if ($END_STITCHING_PROG_NAMES{$ev_type}) {
            &try_stitching_evidence_ends_into_existing_exons($chain);
        }
        if ($EXTEND_TO_TERMINI_TYPES{$ev_type}) {
            &try_extending_termini_to_terminal_exons($chain, $begins_aref, $ends_aref);
        }
    }
}


####
sub try_stitching_evidence_ends_into_existing_exons {
    my ($alignment_chain) = @_;
    
    my $accession = $alignment_chain->{accession};
    my $ev_type = $alignment_chain->{ev_type};
    
    my @links = sort {$a->[0]<=>$b->[0]} @{$alignment_chain->{links}};
    
    unless (scalar (@links) > 1) {
        # for now, only try stitching spliced alignments
        return;
    }

    my $ev_info_href = { accession => $alignment_chain->{accession},
                         ev_type => $alignment_chain->{ev_type},
                     };
    
    my $first_link = shift @links;
    {
        my ($end5, $end3) = @$first_link;
        &try_stitch_5prime_exon($end5, $end3, $ev_info_href);
    }

    my $last_link = pop @links;
    {
        my ($end5, $end3) = @$last_link;
        &try_stitch_3prime_exon($end5, $end3, $ev_info_href);
    }
    
    
}


####
sub try_stitch_5prime_exon {
    my ($end5, $end3, $ev_info_href) = @_;

    unless ($end5 && $end3 && ref $ev_info_href) {
        confess "Error, params needed: end5, end3, ev_info_href\n";
    }
        
    my $accession = $ev_info_href->{accession};
    my $ev_type = $ev_info_href->{ev_type};
    
    
    if ($DEBUG) {
        print "Trying to stitch 5prime: $end5,$end3\t$accession;$ev_type\n";
    }
    
    ## require that end3 correspond to a donor:
    unless ($GENOME_FEATURES[$end3+1] == $DONOR) {
        print "\tsorry, stitched $end3 + 1 != donor\n" if $DEBUG;
        return; # not interested
    }
    
    my @overlapping_exons = &find_overlapping_exons($end5, $end5); # just want end5 to anchor
   
    foreach my $exon (@overlapping_exons) {
        my $exon_type = $exon->getExonType();
        unless ($exon_type eq "initial" || $exon_type eq "internal") {
            next;
        }
        
        my ($exon_end5, $exon_end3) = $exon->get_coords();
        if ($exon_end5 < $end5) {
            ## good candidate for stitching
            my ($new_end5, $new_end3) = ($exon_end5, $end3);
            foreach my $phase (&determine_good_phases($new_end5, $new_end3)) {
                $mock_exon_counter++;    
                my $exon_name = "mock_exon_5prime_stitch_${new_end5}-${new_end3}_${accession}_${ev_type}";
                if ($exon_type eq "initial") {
                    if ($phase == 1) {
                        &try_add_stitched_exon($exon_name, $new_end5, $new_end3, "initial", 1, $ev_info_href);
                    }
                }
                else {
                    # internal exon
                    &try_add_stitched_exon($exon_name, $new_end5, $new_end3, "internal", $phase, $ev_info_href);
                }
            }
        }
    }
}


####
sub try_stitch_3prime_exon {
    my ($end5, $end3, $ev_info_href) = @_;
    
    
    ## require that end5 correspond to an acceptor:
    unless ($GENOME_FEATURES[$end5-2] == $ACCEPTOR) {
        return; # not interested
    }

    my $accession = $ev_info_href->{accession};
    my $ev_type = $ev_info_href->{ev_type};
    
    print "trying to stitch 3prime exon ($end5, $end3)\t$accession;$ev_type\n" if $DEBUG;
    
    my @overlapping_exons = &find_overlapping_exons($end3, $end3); # just want end3 to anchor
       
    foreach my $exon (@overlapping_exons) {
        my $exon_type = $exon->getExonType();
        unless ($exon_type eq "internal" || $exon_type eq "terminal") {
            next;
        }
        
        my ($exon_end5, $exon_end3) = $exon->get_coords();
        if ($exon_end3 > $end3) { 
            ## good candidate for stitching
            print "3prime stitching: found candidate ($exon_type: $exon_end5,$exon_end3)\n" if $DEBUG;
            my ($new_end5, $new_end3) = ($end5, $exon_end3);
            
            my $phase_check_end3 = $new_end3;
            if ($exon_type eq "terminal") {
                ## knock off the stop codon, otherwise it'll get trashed as stop-containing as you'd expect!
                $phase_check_end3 -= 3;
                if ($phase_check_end3 < $new_end5) {
                    next; #doh! too short
                }
            }
            
            foreach my $phase (&determine_good_phases($new_end5, $phase_check_end3)) {
                $mock_exon_counter++;    
                my $exon_name = "mock_exon_3prime_stitch_${new_end5}-${new_end3}_${accession}_${ev_type}";
                print "OKphase($phase) for $exon_name\n" if $DEBUG;
                if ($exon_type eq "terminal") {
                    # make sure stop codon remains in frame:
                    my $new_exon_length = $new_end3 - $new_end5 + 1;
                    if ( ($new_exon_length + $phase - 1) % 3 == 0) { # phase is [1.2.3], make adjustment for integral codons
                        &try_add_stitched_exon($exon_name, $new_end5, $new_end3, "terminal", $phase, $ev_info_href);
                    }
                    else {
                        print "sorry, phase($phase) and exon length ($new_exon_length) are not terminal-friendly\n" if $DEBUG;
                    }
                }
                else {
                    # internal exon
                    &try_add_stitched_exon($exon_name, $new_end5, $new_end3, "internal", $phase, $ev_info_href);
                }
            }
        }
    }
}



####
sub try_add_stitched_exon {
    my ($exon_name, $new_end5, $new_end3, $exon_type, $phase, $ev_info_href) = @_;
    
    my ($accession, $ev_type) = ($ev_info_href->{accession}, $ev_info_href->{ev_type});
    
    print "method try_add_stitched_exon($exon_name, $new_end5, $new_end3, $exon_type, $phase)\n" if $DEBUG;
    
    my $exon = &get_exon_via_coords($new_end5, $new_end3, $exon_type, $phase);
    if ($exon) {
        ## don't add it if it already exists as a candidate
        if ($DEBUG) {
            print "\talready got this exon ($new_end5, $new_end3, $exon_type, $phase:\n";
            print "\t" . $exon->toString() . "\n";
        }
        ## add to evidence:
        $exon->appendEvidence($exon_name, $ev_type);
    }
    else {
        &add_exon($exon_name, $new_end5, $new_end3, $exon_type, $phase, $ev_type);
    }
}


####
sub add_internal_evidence_exon {
    my ($accession, $end5, $end3, $ev_type) = @_;
    
    ## determine which, if any frames allow addition of exon; no intervening stops allowed.
    
    my @good_frames = &determine_good_phases($end5, $end3);
    
    if (@good_frames) {
        
        foreach my $frame (@good_frames) {
            print "Adding evidence ($accession) based exon $end5-$end3, frame: $frame\n" if $SEE;
            &add_exon($accession, $end5, $end3, "internal", $frame, $ev_type);
        }

        return (1); # found at least one good frame
    }

    else {
        ## no good frames. 
        return (0);
    }

}


####
sub determine_good_phases {
    my ($end5, $end3) = @_;
    
    ## checks all frames (phases) starting with phase 1 at bp = end5
    
    
    my %phase_ok = ( 1 => 1,
                     2 => 1,
                     3 => 1 );
    
    for (my $i = $end5; $i <= $end3-2; $i++) {
        if ($GENOME_FEATURES[$i] == $STOP) {
            
            # calc phase of exon
            my $delta = $i - $end5;
            my $mod = $delta % 3;
            
            my $phase;
            
            if ($mod == 0) {
                # xxx xxx STOP
                $phase = 1;
            }
            
            elsif ($mod == 1) {
                # x xxx xxx STOP
                $phase = 3;
            }

            elsif ($mod == 2) {
                # xx xxx xxx STOP
                $phase = 2;
            }
            
            $phase_ok{$phase} = 0;
            # print "\tgot stop in phase ($phase) coord ($i)\n" if $DEBUG;
        }
    }
    
    my @good_phases;
    
    foreach my $phase (sort {$a<=>$b} keys %phase_ok) {
        if ($phase_ok{$phase}) {
            push (@good_phases, $phase);
            
        }
    }
    
    if ($DEBUG) {
        if (@good_phases) {
            print "good phases: " . join (",", @good_phases) . "\n";
        }
        else {
            print "\tsorry, no good phases between $end5 and $end3\n";
        }
    }
    
    return (@good_phases);
}


####
sub parse_evidence_chains {
    
    my $genomic_strand = shift;  #  allowed [+|-|?]   if (?), applied_orient is original orient.  if (+|-), strand is restrictive to evidence gathering, and (-) are transposed to the (+) strand.
    my $evidence_filename = shift;

    ## structure evidence chain:
    ##    accession => string,   (combo of ev_type, accession, and chainID)
    ##    ev_type => string,
    ##    lend => int
    ##    rend => int
    ##    links => [ [end5, end3], [end5, end3], ...]
    ##    gaps => [ [gap_end5, gap_end3], ...]
    ##    applied_orient => +|-  
    
    my %acc_to_chain;
    
    open (FILE, $evidence_filename) or die "Error, cannot open file $evidence_filename";
    while (<FILE>) {
        my $line = $_;
        chomp;
        unless (/\w/) { next;}
        if (/^\#/) { next;}
        my @x = split (/\t/);
        
        my ($ev_type, $lend, $rend, $attributes, $orient) = ($x[1], $x[3], $x[4], $x[8], $x[6]);

        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend); # just in case...
        
        unless (defined ($EVIDENCE_WEIGHTING{$ev_type})) {
            print STDERR "WARNING: not considering ev_type: $ev_type since not included in weights file\n";
            next; # not processing anything not indicated in the weights file.
        }
        
        unless (&in_range_of_genomic_sequence($lend, $rend) ) {
            print STDERR "-WARNING, IGNORING line: $line\tBECAUSE NOT IN RANGE OF GENOMIC SEQUENCE (1-$genomic_seq_length)\n";
            next;
        }
        
        $attributes =~ /ID=([^; ]+)/ or confess "error, no chainID in attributes $attributes of line $_\n";
        my $chainID = $1;
        
        
        my $accession = "";
        if ($attributes =~ /(Target|Query)=([^; ]+)/) {
            $accession = $2;
        }
        
        my $is_child = 0;
        if ($attributes =~ /Parent=([^; ]+)/) {
            $is_child = 1;
            $chainID = $1;
        }
        
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        unless ($genomic_strand eq '?' || $orient eq $genomic_strand) { next; }
        
        my $ev_class = $EV_TYPE_TO_EV_CLASS{$ev_type} or confess "Error, cannot determine ev_class from ev_type: $ev_type";
        
        if ($genomic_strand eq '-') {
            ($end5, $end3) = &revcomp_coordinates($end5, $end3);
            $orient = '+'; # reset
        }
        
        #my $key = $ev_type . "/" . $accession . "/" . $chainID;  # some users aren't including Target|Query in top 'protein_match' line.
        my $key = "ev_type:$ev_type" . "/" . "ID=$chainID";

        my $href = $acc_to_chain{$key};
        unless (ref $href) {
            $href = $acc_to_chain{$key} = { accession => $key,
                                            target => $accession, # we'll put it here for now.
                                            ev_type => $ev_type,
                                            ev_class => $ev_class,
                                            lend => undef,
                                            rend => undef,
                                            links => [],  ## end5,end3 values
                                            gaps => [],  ## gap coords
                                            applied_orient => $orient,
                                            
                                            ## allow for match/match_part representation
                                            children_links => [],
                                        };
        }
        
        if ($accession) {
            ## ensure the top chain record has it recorded under target
            $href->{target} = $accession;
        }
        

        if ($is_child) {
            push (@{$href->{children_links}}, [$end5, $end3]);
        }
        else {
            push (@{$href->{links}}, [$end5, $end3]);
        }
    }
    close FILE;
    
    
    ## order the links, set lend and rend
    my @evidence_chains;
    foreach my $key (keys %acc_to_chain) {
        my $href = $acc_to_chain{$key};
        
        ## give children links priority
        if (@{$href->{children_links}}) {
            $href->{links} = $href->{children_links};
        }
        
        my $links_aref = $href->{links};
        @$links_aref = sort {$a->[0]<=>$b->[0]} @$links_aref;
        
        my @coords;
        foreach my $link (@$links_aref) {
            push (@coords, @$link);
        }
        my $lend = min(@coords);
        my $rend = max(@coords);

        
        $href->{lend} = $lend;
        $href->{rend} = $rend;

        ## proces gaps:
        my $gaps_aref = $href->{gaps};
        for (my $i = 1; $i <= $#$links_aref; $i++) {
            my ($prev_end5, $prev_end3) = sort {$a<=>$b} @{$links_aref->[$i-1]};
            my ($curr_end5, $curr_end3) = sort {$a<=>$b} @{$links_aref->[$i]};
            my $gap_end5 = $prev_end3+1;
            my $gap_end3 = $curr_end5 - 1;
            if ($gap_end3 > $gap_end5) {
                push (@$gaps_aref, [$gap_end5, $gap_end3]);
            }
            else {
                print STDERR "alignment gap has misordered coords: $gap_end5 - $gap_end3\n" . Dumper ($links_aref);
            }
        }
        
        push (@evidence_chains, $href);
    }
    

    return (@evidence_chains);
}

####
sub fragment_evidence_chains_using_max_gap {
    my ($max_gap_length, $evidence_chains_aref) = @_;

    my @ret_evidence_chains;
    foreach my $chain (@$evidence_chains_aref) {
        my ($accession, $ev_type, $lend, $rend, $links, $applied_orient, $gaps) = ($chain->{accession},
                                                                                   $chain->{ev_type},
                                                                                   $chain->{lend},
                                                                                   $chain->{rend},
                                                                                   $chain->{links},
                                                                                   $chain->{applied_orient},
                                                                                   $chain->{gaps} );
        ## see if requires fragmentation:
        my $requires_fragmentation = 0;
        foreach my $gap (@$gaps) {
            if ($gap >= $max_gap_length) {
                $requires_fragmentation = 1;
                last;
            }
        }
        unless ($requires_fragmentation) {
            push (@ret_evidence_chains, $chain);
            next;
        }
        
        ## fragment it:
        my $first_link = shift @$links;
        unless (ref $first_link) {
            confess "error, no link in " . Dumper ($chain);
        }
        
        my @link_clusters = ( [$first_link] );
        
        foreach my $link (@$links) {
            my ($prev_lend, $prev_rend) = (@{$link_clusters[$#link_clusters]->[-1]}); #last element of the last array ref in clusters
            my ($curr_lend, $curr_rend) = (@$link);
                       
            if ($curr_lend - $prev_rend -1 >= $max_gap_length) {
                push (@link_clusters, [$link]); # new cluster
            }
            else {
                push (@{$link_clusters[$#link_clusters]}, $link); # add to existing cluster
            }
        }
        
        ## create new evidence chains for each:
        foreach my $cluster (@link_clusters) {
            my @coords;
            foreach my $link (@$cluster) {
                push (@coords, @$link);
            }
            my $lend = min (@coords);
            my $rend = max (@coords);
            
            my $chain = { accession => $accession,
                          ev_type => $ev_type,
                          lend => $lend,
                          rend => $rend,
                          links => $cluster,
                          applied_orient => $applied_orient,
                          gaps => undef, # should recalc if needed later on...
            };
            push (@ret_evidence_chains, $chain);
        }


    }
    
    return (@ret_evidence_chains);
}


#### 

sub analyze_peaks {
    my $vector_aref = shift;
    
    my @foundPeaks;
    
    my $bestPeakScoreSoFar = 0;
    my $bestPeakPositionSoFar = 0;

    my $currentPeakScore = 0;
    
    for (my $i=1; $i <= length($genomicSeq); $i++) { 
        
        my $leadingEdge = $i;
        my $trailingEdge = $i - $CHAIN_TERMINI_WINDOW;
        
        
        ## Check to see if we're outside the window.
        if ($leadingEdge - $bestPeakPositionSoFar > $CHAIN_TERMINI_WINDOW) {
            if ($bestPeakScoreSoFar > $SUM_GENEPRED_WEIGHTS) {
                push (@foundPeaks, 
                      { position => $bestPeakPositionSoFar,
                        score => $bestPeakScoreSoFar } 
                      );
                if ($DEBUG) {
                    print "($i) ADDING peak for pos: $bestPeakPositionSoFar, score: $bestPeakScoreSoFar, SUMPRED: $SUM_GENEPRED_WEIGHTS\n";
                }
            }
            # reset
            $bestPeakScoreSoFar = 0;
            $bestPeakPositionSoFar = $i;
        }
        
        my $count = $vector_aref->[$leadingEdge];
        
        $currentPeakScore += $vector_aref->[$leadingEdge];
        if ($trailingEdge > 0) {
            $currentPeakScore -= $vector_aref->[$trailingEdge];
        }
        
        
        if ($currentPeakScore > $bestPeakScoreSoFar) {
            
            ## update current statistics
            $bestPeakScoreSoFar = $currentPeakScore;
            $bestPeakPositionSoFar = $leadingEdge;
        }
        
    }
    
    return (@foundPeaks);
}


####
sub find_overlapping_introns {
    my ($genomic_strand, $range_lend, $range_rend) = @_;
 
    ## If strand is (-), the coordinates returned are revcomped to refer to the minus strand flipped to the forward orientation.

    my @intron_coords;
    
    foreach my $intron_key (keys %INTRONS_TO_SCORE) {
        my ($intron_end5, $intron_end3) = split (/_/, $intron_key);
        
        my $orient = ($intron_end5 < $intron_end3) ? '+' : '-';
        
        if ($genomic_strand ne $orient) {
            next; #not interested
        }
        
        ## need reference coordinates, not actual
        if ($orient eq '-') {
            ($intron_end5, $intron_end3) = sort {$a<=>$b} &revcomp_coordinates($intron_end5, $intron_end3);
        }
        # check for overlap
        if ($intron_end5 < $range_rend && $intron_end3 > $range_lend) {
            # got overlap
            push (@intron_coords, [$intron_end5, $intron_end3]);
        }
    }
    
    return (@intron_coords);
}


####
sub extend_exon_downstream_to_stop {
    my $exon = shift;
    
    my ($end5, $end3) = $exon->get_coords();
    my $endFrame = $exon->getEndFrame();
    
    return (&extend_downstream_to_stop($end3, $endFrame));
    
}

####
sub extend_downstream_to_stop {
    my ($end3, $endFrame) = @_;
    
    my $adjust = $endFrame - 1;
    $end3 -= $adjust;
    
    my $new_end3 = undef;
    while ( ($end3 < length($genomicSeq)) && (!defined($new_end3)) ) {
        if ($GENOME_FEATURES[$end3] == $STOP) {
            $new_end3 = $end3+2; #make end of stop codon.
            last;
        }
        $end3+=3;
    }
    
    return ($new_end3);
}



####
sub extend_exon_upstream_to_start {
    my $exon = shift;
    
    my ($end5, $end3) = $exon->get_coords();
    my $startFrame = $exon->getStartFrame();
    
    return (&extend_upstream_to_start_codon($end5, $startFrame));
    
}


####
sub extend_upstream_to_start_codon {
    my $end5 = shift;
    my $startFrame = shift;
    my $first_start_encountered_flag = shift; #if set, the first start found is returned.  Otherwise, keeps climbing
    
    if ($startFrame == 2) {
        $end5 += 2;
    } elsif ($startFrame == 3) {
        $end5 += 1;
    }
    
    my $upstream_start = undef;
    while ($end5 >= 0 && $GENOME_FEATURES[$end5] != $STOP) {
        if ($GENOME_FEATURES[$end5] == $START) {
            $upstream_start = $end5;
            if ($first_start_encountered_flag) { last; } # done searching.
        }
        $end5 -= 3;
    }
    
    return ($upstream_start);
}

####
sub extend_exon_internally_to_first_start {
    my $exon = shift;
    my ($end5, $end3) = $exon->get_coords();
    
    my $frame = $exon->getStartFrame();
    if ($frame == 2) {
        $end5 += 2;
    }
    elsif ($frame == 3) {
        $end5 += 1;
    }
    
    
    my $firstStart = undef;
    while ($end5 <= $end3-2) {
        if ($GENOME_FEATURES[$end5] == $START) {
            $firstStart = $end5;
            last;
        }
        $end5 += 3;
    }
    
    return ($firstStart);
    
}


####
sub find_overlapping_exons {
    my ($lendRange, $rendRange) = @_;
    
    ($lendRange, $rendRange) = sort {$a<=>$b} ($lendRange, $rendRange); #just to be sure.


    my @overlapping_exons;
    foreach my $exon (@EXONS) {
        my ($end5, $end3) = sort {$a<=>$b} $exon->get_coords();
        
        #print "overlapping exon analysis ($end5-$end3) vs. ($lendRange-$rendRange)\n" if $DEBUG;
        
        if ($end5 < $rendRange && $end3 > $lendRange) { #overlap
            print "\t** got overlapping exon ($end5,$end3)\n" if $DEBUG;
            push (@overlapping_exons, $exon);
        }
    }
    
    if ($DEBUG) {
        if (@overlapping_exons) {
            print "Got the following overlapping exons for range: ($lendRange, $rendRange):\n";
            foreach my $overlapping_exon (@overlapping_exons) {
                print "exon: " . $overlapping_exon->toString() . "\n";
            }
        }
        else {
            print "sorry, no overlapping exons for coords ($lendRange, $rendRange)\n";
        }
    }
        
    return (@overlapping_exons);
}



####
sub add_exon {
    my ($evidence_acc, $end5, $end3, $type, $start_phase, $ev_type) = @_;
    
    unless ($ev_type) {
        confess "Error, must specify ev_type";
    }
    
    my $ev_class = $EV_TYPE_TO_EV_CLASS{$ev_type} or confess "Error, cannot determine ev_class from ev_type $ev_type";

    print "\n\n\nmethod add_exon($evidence_acc, $end5, $end3, $type, phase: $start_phase, $ev_type, $ev_class)\n" if $DEBUG || $SEE;
        
    ## verify frame is OK. (probably redundant check, but just in case)
    my ($check_end5, $check_end3) = ($end5, $end3);
    if ($type eq "terminal" || $type eq "single") { 
        # trim expected stop codon
        $check_end3 -= 3;
    }
    # adjust check_end5 so begins at phase 1
    my @good_phases = &determine_good_phases($check_end5, $check_end3);
    print "Incoming phase: $start_phase, Good phases are: @good_phases\n" if $DEBUG || $SEE;
    my $good_phase_flag = 0; 
    foreach my $phase (@good_phases) {
        if ($phase == $start_phase) {
            $good_phase_flag = 1;
            last;
        }
    }
    
    unless ($good_phase_flag) {
        print STDERR "add_exon() Sorry, $evidence_acc $end5-$end3-$type-$start_phase is invalid.\n";
        return;
    }
    
    my $exon = &get_exon_via_coords($end5, $end3, $type, $start_phase);
    
    
    if ($evidence_acc =~ /Mock/i && $exon) {
        confess "Error, already got exon, $end5, $end3 for $evidence_acc\n" . Dumper ($exon);
    }

    unless ($exon) {
        
        ## instantiate exon object
        $exon = Exon->new($end5, $end3);
        push (@EXONS, $exon);
		my $exon_coord_token = &coord_token($end5, $end3, $type, $start_phase);
        $EXONS_VIA_COORDS{$exon_coord_token} = $exon;
        
        ## set exon attributes
        $exon->setExonType($type);
        
        $exon->setLeftSeqBoundary( substr($genomicSeq, $end5-1, 2));
        $exon->setRightSeqBoundary( substr($genomicSeq, $end3-2, 2));
        
    }
    
    
    ## set evidence:
    $exon->appendEvidence($evidence_acc, $ev_type);
    
    ## set the frames:
    $exon->setStartFrame($start_phase);
    
}

sub get_exon_via_coords {
	my ($end5, $end3, $type, $start_phase)  = @_;
	my $exon_coord_token = &coord_token($end5, $end3, $type, $start_phase);
    
	my $exon = $EXONS_VIA_COORDS{$exon_coord_token};
	return ($exon);
}


#### 
sub coord_token {
    return (join ("_", @_));
}



#### 
sub score_exons {
    my ($genomic_strand) = @_;

    if ($DEBUG) {
    
        ## write the coding vector to a file.
        open (my $logfh, ">coding_vector.$genomic_strand.dat") or die $!;

        for (my $i = 1; $i <= $#CODING_SCORES; $i++) {
            my $coord = $i;
            my $coding_value = $CODING_SCORES[$coord];
            if ($genomic_strand eq '-') {
                ($coord) = &revcomp_coordinates($coord);
                $coding_value = -1 * $coding_value; 
            }
            print $logfh "$coord\t$coding_value\n";
        }
        
        close $logfh;
    }
    

    foreach my $exon (@EXONS) {
        my ($end5, $end3) = $exon->get_coords();
        my $type = $exon->getExonType();
        my $length = abs ($end3 - $end5) + 1;
        
        my $coding_score = 0;
        
        ## add evidence contribution to score
        for (my $i=$end5; $i <= $end3; $i++) {
            $coding_score += &max(0, $CODING_SCORES[$i]);  # not adding negative coding scores; prevent negative coding score from defeating the predicted coding potential.
        }
        
        ## add evidence-specific scores derived from homology-prediction and transcript alignments supporting exact exon structures.
        my @evidence = $exon->get_evidence();
        foreach my $ev (@evidence) {
            my ($accession, $ev_type) = @$ev;
            my $ev_class = $EV_TYPE_TO_EV_CLASS{$ev_type} || confess "Error, cannot determine class from type $ev_type";
            if ($ev_class =~ /^(OTHER_PREDICTION|TRANSCRIPT)$/) { # protein and prediction types contribute to coding vector scored above.
                my $ev_weight = $EVIDENCE_WEIGHTING{$ev_type};
                my $coding_contribution = 0;
                ## walk the coordinates and ignore masked values
                for (my $i = $end5; $i <= $end3; $i++) {
                    unless ($MASK[$i]) {
                        $coding_contribution += $ev_weight;
                    }
                }
            
                # add to coding score
                $coding_score += $coding_contribution;
                if ($DEBUG) {
                    print "EXON-SPECIFIC CODING CONTRIBUTION: $accession, $ev_type, $ev_class contributes $coding_contribution to exon: $end5-$end3\n";
                }
            }
        }
        
        $exon->{base_score} = $exon->{sum_score} = $coding_score;
        
        #store score contributions
        $exon->{coding_score} = $coding_score;
    }
}


####
sub build_trellis {
    my ($range_lend, $range_rend, $exons_within_range_aref) = @_;
    
    print STDERR "\n\n-building trellis in range $range_lend - $range_rend\n";
    
    local @EXONS = @$exons_within_range_aref;
    
    ## sort by 5' end
    @EXONS = sort {$a->{end5}<=>$b->{end5}} @EXONS;
    
    my $exon_count = scalar (@EXONS);
    print STDERR "-exon list size: $exon_count\n";
    
    ## init exon vals
    foreach my $exon (@EXONS) {
        $exon->{sum_score} = $exon->{base_score};
        $exon->{link} = 0;
    }
    
    ## add boundaries
    my $left_bound = Exon->new($range_lend, $range_lend);
    $left_bound->setExonType("bound");
    $left_bound->setStartFrame(1);
    unshift (@EXONS, $left_bound);
    
    my $right_bound = Exon->new($range_rend, $range_rend);
    $right_bound->setExonType("bound");
    $right_bound->setStartFrame(1);
    push (@EXONS, $right_bound);
    
    my $num_exons = scalar (@EXONS);
    
    my $highest_score = 0;
    my $highest_scoring_exon = 0;
    
    if ($DEBUG) {
        unless ($_TRELLIS_ROUND) {
            # want to append to the trellis for each round.
            if (-e "trellis") { unlink "trellis";};
            if (-e "final_path") { unlink "final_path";}
        }
        
        $_TRELLIS_ROUND++;
        open (TRELLIS, ">>trellis") or die $!;
        print TRELLIS "\n\n\n# Trellis round: $_TRELLIS_ROUND\n";
    }
    
    print STDERR "\n\nbuilding trellis\n";
    ## build trellis from left to right
    for (my $i=1; $i <= $#EXONS; $i++) {
        
        my $percentageFinished = sprintf ("%.1f    ", $i/$num_exons * 100);
        print STDERR "\r$percentageFinished";
        
        my $exon_i = $EXONS[$i];
        
        my $exon_i_base_score = $exon_i->{base_score};
        my $exon_i_sum_score = $exon_i->{sum_score};
        
        my $exon_i_type = $exon_i->getExonType();
        my $exon_i_name = $exon_i_type . $exon_i->get_orientation() 
            . "~" . $exon_i->getStartFrame() . "~" . $exon_i->getEndFrame() . " " 
            . join (", " , $exon_i->get_coords()) . " ($exon_i_sum_score) ";
        
        my $compare_count = 0;
        my $found_compatible_comparison_flag = 0;
        
        my ($exon_i_lend, $exon_i_rend) = sort {$a<=>$b} $exon_i->get_coords();


        for (my $j = $i-1; $j >=0 
             && ($compare_count < $MAX_NUM_PREV_EXONS_COMPARE || !$found_compatible_comparison_flag)
             ; $j--) {
                

            $compare_count++;        
            my $exon_j = $EXONS[$j];
            
            my ($exon_j_lend, $exon_j_rend) = sort {$a<=>$b} $exon_j->get_coords();

            my $exon_j_sum_score = $exon_j->{sum_score};
            
            my $exon_j_type = $exon_j->getExonType();
            
            my $exon_j_start_frame = $exon_j->getStartFrame() or confess "Error, missing start frame for exon " . Dumper ($exon_j);
            
            my $exon_j_name = $exon_j_type . $exon_j->get_orientation() 
                . "~$exon_j_start_frame" . "~" . $exon_j->getEndFrame() . " "  
                .  join (", ", $exon_j->get_coords()) 
                . " ($exon_j_sum_score) ";
            
            
            my $score = $exon_i_base_score + $exon_j_sum_score;
            my $compare_text = "$exon_i_name to $exon_j_name, ";
            
            my $join_score = -1;
            
            ## check for boundary condition:
            if ($exon_i_type eq "bound" || $exon_j_type eq "bound") {
                $join_score = &score_boundary_condition($exon_j, $exon_i);
            }
            else {
                $join_score = &are_compatible_exons($exon_j, $exon_i);
            }
        
            if ($join_score >= 0) { 
                
                $found_compatible_comparison_flag = 1;

                $compare_text .= "compatible ($join_score) ";
                
                $score += $join_score;
                
                $compare_text .= " TotalScore ($score) ";
                
                if ($score > $exon_i_sum_score) {
                    ## better link:
                    $exon_i->{link} = $exon_j;
                    $exon_i_sum_score = $exon_i->{sum_score} = $score;
                    
                    $compare_text .= " [BESTSOFAR]";
                }
                else {
                    
                    $compare_text .= " , but best so far = $exon_i_sum_score when linked to: " . &get_exon_name($exon_i->{link});
                }
                
            }
            
            $compare_text .= "\n";
            print $compare_text if $SEE;
            if ($DEBUG) {
                print TRELLIS $compare_text;
            }
        }
        
        if ($exon_i_sum_score >= $highest_score) {
            $highest_score = $exon_i_sum_score;
            $highest_scoring_exon = $exon_i;
        }
        
    }
    
    print STDERR "\n";
    if ($DEBUG) {
        close TRELLIS;
        open (INTERGENIC, ">intergenic_feat_scores") or die $!;
        foreach my $intergenic (keys %intergenic_cache) {
            my $score = $intergenic_cache{$intergenic};
            print INTERGENIC "$intergenic\t$score\n";
        }
        close INTERGENIC;
    }
    
    return ($highest_scoring_exon);
    
}

####
sub get_exon_name {
    my ($exon) = @_;
    
    if (ref $exon) {

        my $exon_name = $exon->getExonType() . $exon->get_orientation() 
            . $exon->getStartFrame() . "~" . $exon->getEndFrame() . " "  
            .  join (", ", $exon->get_coords());
        
        return($exon_name);
    }

    else {
        return("boundary");
    }
}


        

#### 
sub initialization {
    
    ## clear data structures
    @EXONS = ();
    %EXONS_VIA_COORDS = ();
    @GENOME_FEATURES = ();
    @CODING_SCORES = ();
    @begins = ();
    @ends = ();
    
    $#GENOME_FEATURES = $genomic_seq_length; #preallocate size
    $#CODING_SCORES = $genomic_seq_length; #preallocate size
    
    for (my $i=0; $i <= $genomic_seq_length + 1; $i++) {
        $GENOME_FEATURES[$i] = 0; #init
        $CODING_SCORES[$i] = 0;
        $begins[$i] = 0;
        $ends[$i] = 0;
    }
    
    unless (%ACCEPTABLE_EXON_LINKAGES) { ## this needs initing only once
        ## format (typeA, typeB, phasedFlag)
        my @acceptableExonLinkages = ( 
                                       
                                       ## Forward strand
                                       ["initial+", "terminal+", 1],
                                       ["initial+", "internal+", 1],
                                       
                                       
                                       ["internal+", "internal+", 1],
                                       ["internal+", "terminal+", 1],
                                       
                                       ["terminal+", "initial+", 0],
                                       ["terminal+", "single+", 0], 
                                       
                                       ["single+", "single+", 0],
                                       ["single+", "initial+", 0],
                                       
                                       ## reverse strand
                                       ["terminal-", "initial-", 1],
                                       ["internal-", "initial-", 1],
                                       
                                       ["internal-", "internal-", 1],
                                       ["terminal-", "internal-", 1],
                                       
                                       ["initial-", "terminal-", 0],
                                       ["single-", "terminal-", 0],
                                       
                                       ["single-", "single-", 0],
                                       ["initial-", "single-", 0],
                                       
                                       ## transitions forward to reverse
                                       ["single+", "single-", 0],
                                       ["single+", "terminal-", 0],
                                       ["terminal+", "terminal-", 0],
                                       ["terminal+", "single-", 0],
                                       
                                       
                                       ## transitions reverse to forward
                                       ["single-", "single+", 0],
                                       ["single-", "initial+", 0],
                                       ["initial-", "initial+", 0],
                                       ["initial-", "single+", 0],
                                       
                                       ## include frame info for internal exon connects:
                                       ## forward strand
                                       [1,2],
                                       [2,3],
                                       [3,1],
                                       ## reverse strand
                                       [4,5],
                                       [5,6],
                                       [6,4]
                                       
                                       
                                       );
        
        
        foreach my $pair (@acceptableExonLinkages) {
            my ($typeA, $typeB, $phased) = @$pair;
            $ACCEPTABLE_EXON_LINKAGES{$typeA}->{$typeB} = 1;
            
            if ($phased) {
                $PHASED_CONNECTIONS{$typeA}->{$typeB} = 1;
            }
            
        }
        
        
        my @intergenic_connects = (
                                   
                                   # forward strand
                                   ["terminal+", "initial+"],
                                   ["terminal+", "single+"], 
                                   
                                   ["single+", "single+"],
                                   ["single+", "initial+"],
                                   
                                   # reverse strand
                                   ["initial-", "terminal-"],
                                   ["single-", "terminal-"],
                                   
                                   ["single-", "single-"],
                                   ["initial-", "single-"],
                                   
                                   # transitions forward to reverse:
                                   ["single+", "single-"],
                                   ["single+", "terminal-"],
                                   ["terminal+", "terminal-"],
                                   ["terminal+", "single-"],
                                   
                                   ## transitions reverse to forward
                                   ["single-", "single+"],
                                   ["single-", "initial+"],
                                   ["initial-", "initial+"],
                                   ["initial-", "single+"],
                                   
                                   );
        
        foreach my $intergenic_pair (@intergenic_connects) {
            my ($typeA, $typeB) = @$intergenic_pair;
            $INTERGENIC_CONNECTIONS{$typeA}->{$typeB} = 1;
        }
        
        
    }
    
    unless (%GFF_PHASE_CONVERSIONS) {
        ## only need do this once
        %GFF_PHASE_CONVERSIONS = ( 
                                   ## forward strand phases
                                   '0+' => 1,
                                   '1+' => 2,
                                   '2+' => 3,

                                   ## reverse strand phases
                                   '0-' => 4,
                                   '1-' => 5,
                                   '2-' => 6,
                                   );
    }
    
}

####
sub add_match_coverage {
    my ($end5, $end3, $weight, $accession, $ev_type) = @_;
    
    unless ($ev_type) {
        confess "Error, must specify ev_type in param";
    }
    my $ev_class = $EV_TYPE_TO_EV_CLASS{$ev_type};
    unless ($ev_class) {
        confess "Error, cannot determine ev class from type: $ev_type";
    }
    
    ## allowing only PROTEIN and ABINITIO_PREDICTION types to contribute to match coverage.  OTHER_PREDICTION and TRANSCRIPT exons will contribute to exon-specific scores only.
    unless ($ev_class =~ /^(PROTEIN|ABINITIO_PREDICTION)$/) {
        return;
    }
    
    
    
    if ($end5 > $end3) {
        confess "Error, trying to add match coverage from $end5-$end3 but end5 > end3 and we're only operating on forward strand!";
    }

    if ($end5 > $genomic_seq_length || $end3 > $genomic_seq_length
        ||
        $end5 < 1 || $end3 < 1) {
        confess "ERROR, coordinates of feature are not within range of genomic sequence (1-$genomic_seq_length) ";
    }
    
    print "COVERAGE: $end5 - $end3 ($weight) $accession\n" if $DEBUG;
    
    for (my $i= $end5; $i <= $end3; $i++) {
        unless ($MASK[$i]) {
            my $current_coding_score = $CODING_SCORES[$i];
            if ($weight < 0) {
                # don't employ negative coding scores.
                $CODING_SCORES[$i] = max($current_coding_score + $weight, 0);
            }
            else {
                $CODING_SCORES[$i] = $current_coding_score + $weight;
            }
        }
    }
}



####
sub get_converted_GFF_phase {
    my ($phase, $orient) = @_;
    
    my $ret_phase = $GFF_PHASE_CONVERSIONS{ "$phase$orient" } or confess "error, no phase conversion for (phase: $phase, orient: $orient) ";
    
    return ($ret_phase);
}

####
sub are_compatible_exons {
    
    ## Exon A must come before Exon B
    
    # returns -1 on 'No'
    # returns intron score on Yes.  If not an intron, returns zero.
    
    my ($exonA, $exonB) = @_;
    

    my $SEE_incompatible = 0; #$DEBUG;

    if ($SEE_incompatible) {
        print "TESTING_COMPATIBILITY between " . $exonA->toString() . " and " . $exonB->toString() . "\n";
    }
    
    
    my $exonAtype = $exonA->getExonType() . $exonA->get_orientation();
    my $exonBtype = $exonB->getExonType() . $exonB->get_orientation();
    
    my $intron_score = 0;
    
    ## are types compatible?
    unless ($ACCEPTABLE_EXON_LINKAGES{$exonAtype}->{$exonBtype}) {
        
        print "incompatibile: unacceptable linkage of $exonAtype to $exonBtype\n" if $SEE_incompatible;
        
        return (-1); #false
    }
    
    my $endFrame = $exonA->getEndFrame();
    my $startFrame = $exonB->getStartFrame();
    
    my ($exonA_end5, $exonA_end3) = $exonA->get_coords();
    my ($exonA_lend, $exonA_rend) = sort {$a<=>$b} ($exonA_end5, $exonA_end3);
    
    my ($exonB_end5, $exonB_end3) = $exonB->get_coords();
    my ($exonB_lend, $exonB_rend) = sort {$a<=>$b} ($exonB_end5, $exonB_end3);
    
    
    ## No overlap allowed:
    if ($exonA_lend <= $exonB_rend && $exonA_rend >= $exonB_lend) { 
        print "incompatible: exons overlap\n" if $SEE_incompatible;
        return (-1); #false
    }
    
    ## check maximum intron length
    if ($PHASED_CONNECTIONS{$exonAtype}->{$exonBtype}) {
        
        
        my ($intron_end5, $intron_end3) = ($exonA_rend+1, $exonB_lend-2);
        
        ## Check to see if stop codon is created across junction:
        my ($exon_before, $exon_after) = ($exonA, $exonB);
        if ($exonAtype =~ /\-$/) { #reverse strand
            ($exon_before, $exon_after) = ($exon_after, $exon_before);
            ($intron_end5, $intron_end3) = ($exonB_lend-1, $exonA_rend+2);
        }
        
        my $intron_key = "$intron_end5" . "_" . "$intron_end3";
        $intron_score = $INTRONS_TO_SCORE{$intron_key};
        unless (defined $intron_score) {
            ## invalid intron
            print "\tINVALID_INTRON: $intron_key\n" if $SEE_incompatible;
            return (-1);
        }
        
        ## check that the exon out/in frames are compatible.
        unless ($ACCEPTABLE_EXON_LINKAGES{$exon_before->getEndFrame()}->{$exon_after->getStartFrame()}) {
            print "incompatible: frames incompatible: " . $exon_before->getEndFrame() . " -> " . $exon_after->getStartFrame() . "\n" if $SEE_incompatible;
            return (-1); #false
        }
        
        my $seqJunction = $exon_before->{rightSeqBoundary} . $exon_after->{leftSeqBoundary};
        
        #print STDERR "SEQ_JUNCTION: $seqJunction\n";
        
        my $endFrame = $exon_before->getEndFrame() % 3;
        if ($endFrame == 0) { $endFrame = 3;}
        
        my $potential_stop = "XXX";
        if ($endFrame == 1) {
            $potential_stop = substr($seqJunction, 1, 3);
        } elsif ($endFrame == 2) {
            $potential_stop = substr($seqJunction, 0, 3);
        }
        
        if (&is_stop_codon($potential_stop)) {
            #print STDERR "Created a stop across the junction : $seqJunction\n";
            print "incompatible: stop codon created across junction ($potential_stop)\n" if $SEE_incompatible;
            return (-1); #created a stop codon across the junction.
        }
        
    }
    
    
    else {
        ## Unphased connection.  Must be intergenic.
        unless ($INTERGENIC_CONNECTIONS{$exonAtype}->{$exonBtype}) {
            confess "Error, supposedly connectable but cannot find ($exonAtype,$exonBtype) categorized specifically as intergenic connections.\n";
        }
        
        unless ($exonA_rend < $exonB_lend) {
            confess "Error, why is exonA_rend ($exonA_rend) ! < exonB_lend ($exonB_lend) ? ";
        }
        my $intergenic_score = &calc_intergenic_score ($exonA_rend+1, $exonB_lend-1);
        
        
        return ($intergenic_score);
    }
    
    
    ## If got here, passed all the above tests.
    return ($intron_score); #true
    
}


####
sub traverse_path {
    my $top_exon = shift;
    
    my $link = $top_exon;
    
    my @entries;
    
    while (ref $link) {
        push (@entries, $link);
        $link = $link->{link};
    }
    
    ## join segments into complete predictions
    my @predictions;
    my $current_prediction = [];

    @entries = reverse @entries;
    
    my $final_path_fh;
    if ($DEBUG) {
        open ($final_path_fh, ">>final_path") or confess "error, cannot open file final_path for writing";
        print $final_path_fh "\n# path chosen for trellis: $_TRELLIS_ROUND\n";
    }
    
    foreach my $entry (@entries) {
        
        if ($DEBUG) {
            my $exon_name = $entry->getExonType() . $entry->get_orientation() . "~" . $entry->getStartFrame() . " " . join (", " , $entry->get_coords());
            print $final_path_fh "$exon_name\n";
        }


        my $type = $entry->getExonType() . $entry->get_orientation();
        unless ($type =~ /bound/) {
            #print $entry->toString()  . "\n";
            push (@$current_prediction, $entry);
        }
		
        if ($type =~ /terminal\+|single|initial\-/) {
            #separate gene entries.
            #print "\n";
            push (@predictions, $current_prediction);
            $current_prediction = [];
        }
        
    }
    if (@$current_prediction) {
        # don't forget the last one
        push (@predictions, $current_prediction);
    }
    
    close $final_path_fh if $DEBUG;

    ###############################
    # Create Prediction Objects ###
    ###############################
    
    my @gene_prediction_objects;
    foreach my $prediction (@predictions) {
        my @segments = @$prediction;
        my $prediction_object = new EVM_prediction (@segments);
        push (@gene_prediction_objects, $prediction_object);
    }
    
    return (@gene_prediction_objects);
    
}


####
sub revcomp_coordinates {
    my @coords = @_;
    
    my @ret_coords;
    foreach my $coord (@coords) {
        my $rev_coord = $genomic_seq_length - $coord + 1;
        push (@ret_coords, $rev_coord);
    }
    
    return (@ret_coords);
}

sub reverse_complement { 
    my($s) = @_;
    my ($rc);
    $rc = reverse ($s);
    $rc =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    
    return($rc);
}



####
sub transpose_exons_back_to_forward_strand {
    
    my %forward_to_reverse_frame = ( 1 => 4,
                                     2 => 5,
                                     3 => 6 );
    
    
    foreach my $exon (@EXONS) {
    
        my ($end5, $end3) = $exon->get_coords();
        my $startFrame = $exon->getStartFrame();
        my $endFrame = $exon->getEndFrame();
        
        ($exon->{end5}) = &revcomp_coordinates($end5);
        ($exon->{end3}) = &revcomp_coordinates($end3);
        
        $exon->{startFrame} = $forward_to_reverse_frame{$startFrame};
        $exon->{endFrame} = $forward_to_reverse_frame{$endFrame};
        $exon->{orientation} = '-';
    }
    
}

####
sub readWeights {
    open (WEIGHTS, $weightsFile) or die "Error, cannot open $weightsFile";
    while (<WEIGHTS>) {
        if (/^\#/) { next; } #comment line
        unless (/\w/) { next;}
        chomp;
        my ($ev_class, $ev_type, $weight) = split (/\s+/);
        
        unless ($ev_class && $ev_type && $weight =~ /\d/) {
            confess "cannot properly parse line of weights file: $_";
        }
        
        unless ($ALLOWABLE_EVIDENCE_CLASSES{$ev_class}) {
            confess "Sorry, [$ev_class] evidence class is not a recognized class type.\n";
        }
        
        unless ($weight =~ /\d/) {
            confess "Error, $ev_type ev-type in weights file $weightsFile is given non-numeric value $weight.\n";
        }
        
        $EV_TYPE_TO_EV_CLASS{$ev_type} = $ev_class;
        
        # store weight value.
        $EVIDENCE_WEIGHTING{$ev_type} = $weight;
        
        # examine prediction type.
        if ($ev_class =~ /_PREDICTION$/) {
            $PREDICTION_PROGS{$ev_type} = 1;
            if ($ev_class eq 'ABINITIO_PREDICTION') {
                $PREDICTION_PROGS_CONTRIBUTE_INTERGENIC{$ev_type} = 1; # others just count towards introns and exons, not intergenic.
            }
        }
                    
    }
}

####
sub write_introns_and_exons_to_file {
    open (EXONS, ">exon_list.out") or die $!;
    select EXONS;
    &dump_exons();
    close EXONS;
    select STDOUT;
    
    open (INTRONS, ">introns.out");
    my @introns;
    foreach my $key (keys %INTRONS_TO_SCORE) {
        my ($end5, $end3) = split (/_/, $key);
        my $score = $INTRONS_TO_SCORE{$key};
        my $evidence_list_aref = $INTRONS_TO_EVIDENCE{$key};
        my $evidence_text = "";
        foreach my $ev (@$evidence_list_aref) {
            my ($acc, $ev_type) = @$ev;
            $evidence_text .= "{$acc; $ev_type},";
        }
        chop $evidence_text; # remove trailing comma
        
        my $orient = ($end5 < $end3) ? "+" : '-';
        $end3 = ($orient eq '+') ? ++$end3 : --$end3; # adjust for dinuc splice offset in introns_to_score hash.
        
        push (@introns, [$end5, $end3, $score, $evidence_text]);
    }
    @introns = sort {$a->[0]<=>$b->[0]} @introns;
    foreach my $intron (@introns) {
        my ($end5, $end3, $score, $evidence_text) = @$intron;
        my $intron_length = abs($end3-$end5) + 1;
        my $per_bp_score = sprintf ("%.2f", $score / $intron_length);

        print INTRONS "$end5-$end3 (score: $score) (score_per_base: $per_bp_score) $evidence_text\n";
    }
    close INTRONS;
    
    return;

}


####
sub populate_forward_reverse_pred_intron_vectors {
    
    $#FORWARD_PRED_INTRON_VEC = $genomic_seq_length; #prealloc
    $#REVERSE_PRED_INTRON_VEC = $genomic_seq_length; #prealloc

    ## init
    for (my $i=0; $i <= $genomic_seq_length; $i++) {
        $FORWARD_PRED_INTRON_VEC[$i] = 0;
        $REVERSE_PRED_INTRON_VEC[$i] = 0;
    }
    
    ## propagate discrete intron scores into basepair values:
    foreach my $intron (keys %PREDICTED_INTRONS) {
        my ($end5, $end3) = &intron_key_to_intron_span($intron);
        
        my $orient = ($end5 < $end3) ? "+" : "-";
        
        my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
        
        my $adj_intron_length = &adjust_feature_length_for_mask($lend, $rend);
        
        unless ($adj_intron_length > 0) {
            next;
        }

        my $intron_score = $PREDICTED_INTRONS{$intron};
        my $score_per_bp = $intron_score / $adj_intron_length;
        
        # contribute scores to the corresponding vector:
        my $scoring_vec = ($orient eq '+') ? \@FORWARD_PRED_INTRON_VEC : \@REVERSE_PRED_INTRON_VEC;
        
        for (my $i = $lend; $i <= $rend; $i++) {
            unless ($MASK[$i]) {
                $scoring_vec->[$i] += $score_per_bp;
            }
        }
    }
    
    return;

}



####
sub write_intron_vectors_to_file {
    
    print STDERR "\n-writing intron vector to file.\n";

    ## create a forward and reverse strand vector:
    my @forward_intron_vec;
    my @reverse_intron_vec;
    $#forward_intron_vec = $genomic_seq_length; #prealloc
    $#reverse_intron_vec = $genomic_seq_length; #prealloc

    for (my $i=0; $i <= $genomic_seq_length; $i++) {
        $forward_intron_vec[$i] = 0;
        $reverse_intron_vec[$i] = 0;
    }

    ## propagate discrete intron scores into basepair values:
    foreach my $intron (keys %INTRONS_TO_SCORE) {
        my ($end5, $end3) = split (/_/, $intron);
        
        my $orient = ($end5 < $end3) ? '+' : '-';
        $end3 = ($orient eq '+') ? ++$end3 : --$end3; # adjust for dincleotide start pos in intron_to_score hash.
        
        my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);

        my $intron_length = $rend - $lend + 1;
        
        my $intron_score = $INTRONS_TO_SCORE{$intron};
        my $score_per_bp = $intron_score / $intron_length;

        my $scoring_vec = ($orient eq '+') ? \@forward_intron_vec : \@reverse_intron_vec;
        
        for (my $i = $lend; $i <= $rend; $i++) {
            $scoring_vec->[$i] += $score_per_bp;
        }
    }

    ## convert reverse intron vector to negative:
    foreach my $value (@reverse_intron_vec) {
        $value = -1 * $value;
    }
    

    ## write vectors to file:
    &write_vector_to_file("introns_decomposed_to_vec.+.dat", \@forward_intron_vec);
    &write_vector_to_file("introns_decomposed_to_vec.-.dat", \@reverse_intron_vec);
    
    return;
    
}

####
sub write_vector_to_file {
    my ($filename, $vector_aref) = @_;
    
    open (my $fh, ">$filename") or die "Error, cannot write to file $filename";
    for (my $i =0; $i <= $#$vector_aref; $i++) {
        print $fh "$i\t" . $vector_aref->[$i] . "\n";
    }

    close $fh;
    return;
}



####
sub repeatMask {
    my $genomic_strand = shift;

    ## init the array
    @MASK = ();
    for (my $i=0; $i <= $genomic_seq_length + 1; $i++) {
        $MASK[$i] = 0;
    }
    
    ## mask N's in the genome sequence:
    my @chars = split (//, uc $genomicSeq);
    for (my $i = 0; $i <= $#chars; $i++) {
        if ($chars[$i] eq 'N') {
            $MASK[$i] = 1;
        }
    }
    
    unless ($repeatsFile) {
        return;
    }
    
    # a repeat is always masked, regardless of orientation.
    # if reading the minus strand, must revcomp all the coords
    
    open (REPEAT, $repeatsFile) or die "Error, cannot open repeats file ($repeatsFile)";
    while (<REPEAT>) {
        my $line = $_;
        chomp;
        unless (/\w/) { next;}
        if (/^\#/) { next;}
        my @x = split (/\t/);
        
        my ($lend, $rend, $orient) = ($x[3], $x[4], $x[6]);
        
        
        unless (&in_range_of_genomic_sequence($lend, $rend) ) {
            print STDERR "-WARNING, IGNORING line: $line\tBECAUSE NOT IN RANGE OF GENOMIC SEQUENCE (1-$genomic_seq_length)\n";
            next;
        }

        if ($genomic_strand eq '-') {
            ($lend, $rend) = sort {$a<=>$b} &revcomp_coordinates($lend, $rend);
        }
        for (my $i = $lend; $i <= $rend; $i++) {
            $MASK[$i] = 1;
        }
    }
    close REPEAT;
}

sub get_gene_predictions {
    my %pred_data;
    open (FILE, $genePredictionsFile) or die $!;
    while (<FILE>) {
        my $line = $_;
        chomp;
        unless (/\w/) { next;}
        if (/^\#/) { next;}
        
        my @x = split (/\t/);
        
        my ($predType, $feat_type, $id, $lend, $rend, $orient) = ($x[1], $x[2], $x[8], $x[3], $x[4], $x[6]);
        
        unless ($feat_type eq 'CDS') { next;}
        
        unless (exists $PREDICTION_PROGS{$predType}) {
            print STDERR "-WARNING: IGNORING prediction type: $predType, since not specified in weights file.\n";
            next;
        }
        
        unless (defined ($EVIDENCE_WEIGHTING{$predType})) {
            confess "no weight for pred-type $predType"; # not examining anything not defined in the weights file.
        }
        
        
        unless (&in_range_of_genomic_sequence($lend, $rend) ) {
            print STDERR "-WARNING, IGNORING line: $line\tBECAUSE NOT IN RANGE OF GENOMIC SEQUENCE (1-$genomic_seq_length)\n";
            next;
        }


        ## if only doing one strand, then only count genes on that strand, and 
        ##   intergenic in between.
        if ($FORWARD_STRAND_ONLY_FLAG && $orient eq '-') { next; }
        if ($REVERSE_STRAND_ONLY_FLAG && $orient eq '+') { next; }
        
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
		
        my $coord_list_ref = $pred_data{$predType}->{$id};
        unless (ref $coord_list_ref) {
            $coord_list_ref = $pred_data{$predType}->{$id} = [];
        }
        push (@$coord_list_ref, [$end5, $end3]);
    }
    close FILE;
    return (%pred_data);
}



####
sub populate_intergenic_regions {
    
    # init the intergenic array
    for (my $i=0; $i <= $genomic_seq_length; $i++) {
        $INTERGENIC_SCORES[$i]=0;
    }
    
    # read in the coordinate list for each gene
    my %pred_data = &get_gene_predictions();
        
    # find the intergenic regions for each prediction type as the 
    # distance between neighboring genes.
    foreach my $predType (keys %PREDICTION_PROGS_CONTRIBUTE_INTERGENIC) {
        
        print "Processing Intergenic: $predType\n" if $SEE;
        my $evidence_weight = $EVIDENCE_WEIGHTING{$predType};
        unless (defined $evidence_weight) {
            confess "Error, weight not defined for type: $predType\n";
        }
        my @coordsets;
        # get the span of each prediction
        foreach my $pred_id (keys %{$pred_data{$predType}}) {
            
            #print STDERR "$pred_id: " . Dumper ($pred_data{$predType}->{$pred_id}) if $DEBUG;

            my @coords;
            foreach my $coordset ( @{$pred_data{$predType}->{$pred_id}}) {
                push (@coords, @$coordset);
            }
            
            @coords = sort {$a<=>$b} @coords;
            unless (@coords) {
                confess "Error, why no coords for $pred_id?\n";
            }
            my $lend = shift @coords;
            my $rend = pop @coords;
            push (@coordsets, [$lend, $rend]);
        }
        
        # add the boundaries including genome begin and end.
        unshift (@coordsets, [0,0]);
        push (@coordsets, [$genomic_seq_length, $genomic_seq_length]);
        @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;
        for (my $i =0; $i < $#coordsets; $i++) {
            my $lend_intergenic = $coordsets[$i]->[1];
            my $rend_intergenic = $coordsets[$i+1]->[0];
            if ($lend_intergenic > $rend_intergenic) {
                warn "Error with prediction: $predType" . #Dumper (\@coordsets) . 
                    " lend_intergenic: $lend_intergenic, rend_intergenic: $rend_intergenic ";
                next;
            }
            print "$predType\tintergenic: $lend_intergenic - $rend_intergenic\n" if $SEE;
            # add the intergenic weighted score
            for (my $j=$lend_intergenic+1; $j <= $rend_intergenic-1; $j++) {
                unless ($MASK[$j]) {
                    #print STDERR "$j intergenic\n";
                    $INTERGENIC_SCORES[$j] += $evidence_weight;
                }
            }
        }
        
        
    }
    
    if ($DEBUG) {
        open (my $fh, ">intergenic.bps") or die $!;
        for (my $i = 0; $i <= $#INTERGENIC_SCORES; $i++) {
            print $fh "$i\t$INTERGENIC_SCORES[$i]\n";
        }
        close $fh;
    }
}


####
sub score_boundary_condition {
    my ($exonA, $exonB) = @_;
    ## exonA comes before exonB
    
    my $exonA_type = $exonA->getExonType();
    my $exonB_type = $exonB->getExonType();
    
    my $exonA_orient = $exonA->get_orientation();
    my $exonB_orient = $exonB->get_orientation();
    
    my ($exonA_end5, $exonA_end3) = $exonA->get_coords();
    my ($exonB_end5, $exonB_end3) = $exonB->get_coords();
    
    if ($exonA_type eq "bound" && $exonB_type eq "bound") {
        # rediculous, don't allow it.
        return (-1);
    }
    
    ## analyze left-bound condition first:
    
    ## !!!! Should really score terminal introns and intergenic separately
    ##       and only allow for known terminal introns.
    #   for now, allowing every feature to extend to the sequence termini
    #   scoring for intergenic.
    #   If a better path exists (which should for all non terminal features), it will score higher than this and be taken.
    
    if ($exonA_type eq "bound") {
        
        

        my ($lend, $rend) = sort {$a<=>$b} ($exonB_end5, $exonB_end3);
        return (&calc_intergenic_score($exonA_end3, $lend-1)); # bound has same coords for end5 and end3, using end3 below
    }
    
    elsif ($exonB_type eq "bound") {
        my ($lend, $rend) = sort {$a<=>$b} ($exonA_end5, $exonA_end3);
        return (&calc_intergenic_score ($rend+1, $exonB_end5));
    }
    
    else {
        confess "one of the exons was supposed to be the bound type...";
    }
    
}

####
sub calc_intergenic_score {
    my ($lend, $rend) = @_;
    
    if ($lend > $rend + 1) { #allow for identical coords
        confess "Error, calc_intergenic_score($lend, $rend)  where $lend > $rend ! ";
    }
    
    my $intergenic_score = 0;
    
    for (my $i = $lend; $i <= $rend; $i++) {
        unless (defined $INTERGENIC_SCORES[$i]) {
            print STDERR "Error, no intergenic score defined for position $i\n";
        }
        $intergenic_score += $INTERGENIC_SCORES[$i]; #should already be masked.
    }
    
    ## apply intergenic adjustment factor:
    $intergenic_score *= $INTERGENIC_SCORE_ADJUST_FACTOR;
    
    if ($DEBUG) {
        $intergenic_cache{"$lend-$rend"} = $intergenic_score;
    }
    
    return ($intergenic_score);
}


####
sub supplement_terminal_exons {
    my $genomic_strand = shift;
    
    my %reversePhase = ( 4 => 1,
                         5 => 2,
                         6 => 3 ); #must flip it around if reverse strand being processed.
    
    open (my $fh, $terminalExonsFile) or die $!;
    while (<$fh>) {
        chomp;
        my $line = $_;
        unless (/\w/) { next;}
        if (/^\#/) { next;}
        my @x = split (/\t/);
        my ($exonType, $ev_type, $acc, $lend, $rend, $orient, $phase) = ($x[2], $x[1], $x[8],$x[3], $x[4], $x[6], $x[7]);
        
        
        unless ($orient eq $genomic_strand) { next;}
        $phase = &get_converted_GFF_phase($phase, $orient);
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
                
        my $evidence_weight = $EVIDENCE_WEIGHTING{$ev_type};
        unless (defined ($evidence_weight)) {
            next;
        }

        my $ev_class = $EV_TYPE_TO_EV_CLASS{$ev_type};        


        my $exon_length = abs ($end3 - $end5) + 1;
        unless ($exon_length >= 3) {
            # no splitting of start/stop codons:
            next;
        }
        
        if ($genomic_strand eq '-') {
            ($end5, $end3) = &revcomp_coordinates($end5, $end3);
            $phase = $reversePhase{$phase} or confess "Error, no reverse phase for phase($phase) line: $line\n";
        }
        
        &add_exon($acc, $end5, $end3, $exonType, $phase, $ev_type);
        
        ## this is the only way that EST-based terminal exons can contribute to the coding score:
        if ($evidence_weight) { 
            &add_match_coverage($end5, $end3, $evidence_weight, "$ev_type,$acc", $ev_type);
        }
    }
    close $fh;
}

####
sub is_stop_codon {
    my ($potential_stop_codon) = @_;
    foreach my $stop_codon (@STOP_CODONS) {
        if ($potential_stop_codon eq $stop_codon) {
            return (1); #yes
        }
    }
    
    return (0); # no
}

####
sub mask_exons {
    print STDERR "mask_exons()\n" if $SEE||$DEBUG;
    
    unless (@EXONS) {
        return;
    }
    
    my @masked_exons;

    my @unmasked_exons; # hold those that lack any masking:
    foreach my $exon (@EXONS) {
        my ($end5, $end3) = sort {$a<=>$b} $exon->get_coords();
        my $exon_length = $end3 - $end5 + 1;
        my $mask_count = 0;
        for (my $i=$end5; $i <= $end3; $i++) {
            if ($MASK[$i]) {
                $mask_count++;
            }
        }
        
        my $mask_exon_percent = $mask_count / $exon_length * 100;
        if ($mask_exon_percent < $MASK_EXON_MIN_PERCENT) {
            push (@unmasked_exons, $exon);
        } # otherwise avoided.
        else {
            $exon->set_masked_percent($mask_exon_percent);
            push (@masked_exons, $exon);
        }
    }

    
    my $total_exons = scalar (@EXONS);
    my $total_after_masking = scalar (@unmasked_exons);
    my $percent_masked = sprintf ("%.2f", ($total_exons - $total_after_masking) / $total_exons * 100);
    print STDERR "\nmask_exons: $percent_masked % of exons are masked.\n";
    
    # replace exons with unmasked ones:
    @EXONS = @unmasked_exons;

    if ($DEBUG) {
        open (my $masked_exons_fh, ">masked_exons.out") or die $!;
        @masked_exons = sort {$a->{end5}<=>$b->{end5}} @masked_exons;
        foreach my $exon (@masked_exons) {
            print $masked_exons_fh $exon->toString . "\t" . $exon->get_masked_percent() . "\n";
        }
    }
}



####
sub filter_predictions_low_support {
    my ($prediction_objects_aref, $mode) = @_;
    
    my @prediction_objects = @$prediction_objects_aref;
    
    ## reset the intergenic regions on first application of this method.  After resetting them, they remain constant for future recursive runs.
    unless ($RESET_INTERGENIC_REGIONS_FLAG) {
        &populate_intergenic_regions();
        $RESET_INTERGENIC_REGIONS_FLAG = 1;
    }
    
    ## Check each prediction to see if its noncoding score exceeds its prediction score:
    
    my @filtered_preds;

    foreach my $prediction (@prediction_objects) {
                
        ## compare coding to recomputed noncoding score
        
        my ($prediction_lend, $prediction_rend) = $prediction->get_span();
        my $prediction_score = $prediction->get_score();
        
        my $noncoding_score = &calc_intergenic_score ($prediction_lend, $prediction_rend);
        
        ## including predicted introns scored as intergenic.
        ## get a noncoding contribution from all introns, plus and minus strands
        my $noncoding_intron_addition = 0;
        for (my $i = $prediction_lend; $i <= $prediction_rend; $i++) {
            $noncoding_intron_addition += $FORWARD_PRED_INTRON_VEC[$i] + $REVERSE_PRED_INTRON_VEC[$i];
        }
        
        ## offset by intron consensus on the prediction-oriented strand.
        my $prediction_orient = $prediction->get_orientation();
        
        my @exons = $prediction->get_exons();
        my $num_exons = scalar (@exons);
        
        my @introns = $prediction->get_intron_coords();
        
        if ($num_exons > 1 && ! @introns) {
            confess "Error, no introns stored but have $num_exons exons. ";
        }
        
        my $offset = 0;
        foreach my $intron (@introns) {
            my ($intron_lend, $intron_rend) = @$intron;
            my $intron_key = &get_intron_key($intron_lend, $intron_rend, $prediction_orient);
            ## get the list of ab initio predictions yielding this intron.
            my $predicted_intron_score_contribution = &get_predicted_intron_score_contribution($intron_key);
            ($intron_lend, $intron_rend) = sort {$a<=>$b} &intron_key_to_intron_span($intron_key); # now set to exon adjacent bases.
            my $intron_len = &adjust_feature_length_for_mask($intron_lend, $intron_rend);
            
            unless ($intron_len > 0) { next; }

            my $existing_score_contribution_per_base = $predicted_intron_score_contribution / $intron_len;
			
			## adjust for intergenic score 
			$offset += &calc_intergenic_score($intron_lend, $intron_rend);
			
            for (my $i=$intron_lend; $i <= $intron_rend; $i++) {
                if ($MASK[$i]) { next; }
                
                if ($prediction_orient eq "+") {
                    $offset += $FORWARD_PRED_INTRON_VEC[$i] - $existing_score_contribution_per_base;
                }
                else {
                    $offset += $REVERSE_PRED_INTRON_VEC[$i] - $existing_score_contribution_per_base;
                }
            }
        }
        
        if ($offset < 0) { 
            confess "Error, offset to augment prediction score for predicted introns is negative: $offset\n";
        }
        
        my $raw_noncoding_score = $noncoding_score + $noncoding_intron_addition;
                
        $prediction->{raw_noncoding} = $raw_noncoding_score;
        $prediction->{offset_noncoding} = $offset; 
        
        my $noncoding_equivalent = $raw_noncoding_score - $offset;
        
        if ($noncoding_equivalent <= 0) {
			$noncoding_equivalent = 0.0001 * $prediction_score; # make it a small number to avoid div-by-zero errors.
        }
        
        $prediction->set_noncoding_equivalent($noncoding_equivalent);
        
        my $coding_noncoding_score_ratio = 0;
        unless ($prediction_score == 0 && $noncoding_equivalent == 0) {
            # avoid div/zero error
            $coding_noncoding_score_ratio = sprintf ("%.2f", $prediction_score / $noncoding_equivalent);
        }
        $prediction->{score_ratio} = $coding_noncoding_score_ratio;
        
        ## first, check the coding length:
        my $coding_length = $prediction->get_coding_length();
        
        if ($DEBUG) {
            print "** Examining prediction for removal:\nMode: $mode\n" . $prediction->toString()
                . "noncoding_intron_addition: $noncoding_intron_addition\n"
                . "noncoding: $noncoding_score\n"
                . "raw_noncoding{sum noncoding + intron_addition}: $raw_noncoding_score\n"
                . "offset from same strand predicted introns: $offset\n"
                . "noncoding equivalent = $noncoding_equivalent\n"
                . "prediction_score: $prediction_score + offset( $offset )\n";
        }
        
        my $min_coding_length = ($mode eq 'STANDARD') ? 150 : $MIN_CODING_LENGTH; #  by default, don't report anything less than 50 aa (150 nt) as a gene.
		

        if ( ($coding_noncoding_score_ratio < $MIN_CODING_NONCODING_SCORE_RATIO) || ($coding_length < $min_coding_length) ) {
            print "\t** eliminating prediction.\n" if $DEBUG;
            $prediction->eliminate();
        }
        
        push (@filtered_preds, $prediction); 
    }
    
    return (@filtered_preds);
    
}

####
sub get_intron_key {
    my ($intron_lend, $intron_rend, $orient) = @_;

    ($intron_lend, $intron_rend) = sort {$a<=>$b} ($intron_lend, $intron_rend); # just be sure they're sorted.

    my ($intron_end5, $intron_end3) = ($orient eq '+') ? ($intron_lend, $intron_rend-1) : ($intron_rend, $intron_lend+1);
    
    my $intron_key = "$intron_end5" . "_" . "$intron_end3";
    
    unless (exists $INTRONS_TO_SCORE{$intron_key}) {
        confess "Error, intron key: $intron_key does not exist";
    }

    return ($intron_key);
}
 

####
sub intron_key_to_intron_span {
    my ($intron_key) = @_;

    my ($intron_end5, $intron_end3) = split (/_/, $intron_key);
    my $orient = ($intron_end5 < $intron_end3) ? "+" : "-";

    if ($orient eq '+') {
        return ($intron_end5, $intron_end3 - 1);
    }
    else {
        return ($intron_end5, $intron_end3 + 1);
    }
}

####
sub get_predicted_intron_score_contribution {
    my ($intron_key) = @_;
   
    my ($intron_lend, $intron_rend) = &intron_key_to_intron_span($intron_key);
    my $intron_len = &adjust_feature_length_for_mask($intron_lend, $intron_rend);
    
    unless ($intron_len > 0) { 
        return(0);
    }
        
    my @evidence_list = @{$INTRONS_TO_EVIDENCE{$intron_key}};
    my @pred_types;
    
    foreach my $ev_pair (@evidence_list) {
        my ($acc, $ev_type) = @$ev_pair;
        if ($PREDICTION_PROGS_CONTRIBUTE_INTERGENIC{$ev_type}) {
            push (@pred_types, $ev_type);
        }
    }
    unless (@pred_types) {
        return (0);
    }
    if (@pred_types) {
        my $existing_contribution = 0;
        foreach my $pred_type (@pred_types) {
            my $weight = $EVIDENCE_WEIGHTING{$pred_type};
            $existing_contribution += $weight;
        }
        return ($existing_contribution * $intron_len);
    }
    
    confess "should have returned something already, never reaching here.";

}

####
sub adjust_feature_length_for_mask {
    my ($lend, $rend) = @_;
    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend); # just to be sure.

    my $len = 0;
    for (my $i = $lend; $i <= $rend; $i++) {
        unless ($MASK[$i]) {
            $len++;
        }
    }
    return ($len);
}

####
sub convert_5prime_partials_to_complete_genes_where_possible {
    my @EVM_predictions = @_;

    foreach my $prediction (@EVM_predictions) {
        
        if ($prediction->is_5prime_partial()) {
          
            print STDERR "GOT 5prime partial gene.\n";
            
            my $orientation = $prediction->get_orientation();
            
            my @exons = sort {$a->{end5}<=>$b->{end5}} $prediction->get_exons();
            
            ## only doing this for multi-exon genes for now.
            if (scalar @exons == 1) { next; }
            
            ## order from gene beginning to gene end
            if ($orientation eq '-') {
                @exons = reverse @exons;
            }

            # walk each exon and see if there's an initial exon that can replace an internal exon
            my @replacement_exons;
            
          INIT_EXON_SEARCH:
            while (my $exon = shift @exons) {
                print "-examining exon for partial->complete conversion: " . $exon->toString() . "\n" if $DEBUG;

                unless ((my $exon_type = $exon->getExonType()) eq 'internal') { 
                    print STDERR "exon type not internal: $exon_type\n"; 
                    last; 
                }
                my ($end5, $end3) = $exon->get_coords();
                my @overlapping_exons = &find_overlapping_exons($end5, $end3);
                my $best_initial_candidate_exon = undef;
                foreach my $overlapping_exon (@overlapping_exons) {
                    print STDERR "\t-found overlapping exon: " . $overlapping_exon->toString() . "\n" if $DEBUG;

                    if ($overlapping_exon->getExonType() eq 'initial' 
                        && $overlapping_exon->get_orientation() eq $orientation
                        && $overlapping_exon->{end3} == $end3
                        && $overlapping_exon->getEndFrame() == $exon->getEndFrame()
                        ) {
                        
                        if ( (!defined($best_initial_candidate_exon)) 
                             ||  $overlapping_exon->get_exon_score() > $best_initial_candidate_exon->get_exon_score()) {
                            $best_initial_candidate_exon = $overlapping_exon;
                        }
                    }
                }

                if ($best_initial_candidate_exon) {
                    ## set replacement exons
                    print STDERR "-found replacement initial exon " . $best_initial_candidate_exon->toString() . "\n" if $DEBUG;
                    @replacement_exons = ($best_initial_candidate_exon);
                    if (@exons) {
                        push (@replacement_exons, @exons);
                    }
                    last INIT_EXON_SEARCH;
                }
                else {
                    ## do not continue search for now.
                    last INIT_EXON_SEARCH;
                    
                    ## Note: breaking the while loop either way!
                }
            }
            if (@replacement_exons) {
                $prediction->replace_exons(@replacement_exons);
            }
        }
    }

    return;
}


sub dump_repeat_mask {
    open (my $fh, ">mask_coords.out") or die $!;
    for (my $i = 0; $i <= $#MASK; $i++) {
        print $fh "$i\t$MASK[$i]\n";
    }
    close $fh;
}


####
sub try_extending_termini_to_terminal_exons {
    my ($chain, $begins_aref, $ends_aref) = @_;
    
    if ($DEBUG) {
        print "-trying to extend termini to exons: " . $chain->{accession} . " " . $chain->{lend} . "-" . $chain->{rend} . "\n";
    }
    
    my $evidence_weight = $EVIDENCE_WEIGHTING{$chain->{ev_type}};
    unless (defined $evidence_weight) { confess "Error, no evidence wieght found for chain: " . Dumper ($chain); }
    
    my $accession = $chain->{accession};
    my $ev_type = $chain->{ev_type};
    
    my $ev_info_href = { accession => $accession,
                         ev_type => $ev_type, 
                     };
    
    my $chain_links_aref = $chain->{links};
        
    my @coordsets = @$chain_links_aref;
    my $num_coordsets = scalar (@coordsets);
    
    ## only trying this for multi-segment alignment chains
    ## where the terminal segments are at consensus splice junctions
    
    if ($num_coordsets < 2) { return; } # no introns!

    ## sort by coordinates so first segment is first and last is last.
    @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;

    my $first_segment = shift @coordsets;
    my $last_segment = pop @coordsets;
    
    {  ## try making an initial exon
                
        my $added_initial_exon_flag = 0;
        my $found_existing_initial_exon_flag = 0;
            
        my ($lend, $rend) = @$first_segment;
        ## must have a splice boundary:
        if ($GENOME_FEATURES[ $rend + 1 ] == $DONOR) {
           
            ## check to make sure that we don't already have a good candidate initial exon with this donor:
            $found_existing_initial_exon_flag = 0;
            my @overlapping_exons = &find_overlapping_exons($lend, $rend);
            foreach my $exon (@overlapping_exons) {
                my $exon_type = $exon->getExonType();
                my ($exon_end5, $exon_end3) = $exon->get_coords();
                if ($exon_type eq 'initial' && $exon_end3 == $rend) {
                    $found_existing_initial_exon_flag = 1;
                    last;
                }
            }
            
            unless ($found_existing_initial_exon_flag) {
                $added_initial_exon_flag = &try_extend_initial_segment_to_initial_exon($first_segment, $ev_info_href);
                
                ## if one couldn't be extended, try searching internally for the first in-frame start codon
                #unless ($added_initial_exon_flag) {
                #    &try_intend_initial_segment_to_initial_exon(@$first_segment); # intend is opposite of extend :)
                #}

                # no longer need the above since we're doing a second search from the end3 to upstream for the start. 
                ## this is superior since some termini might have intervening stops...

            }
        }
    
        if ($added_initial_exon_flag || $found_existing_initial_exon_flag) {
            ## add to the begins collection.
            $begins_aref->[$lend] += $evidence_weight;
        }
    }
    
    
    {  ## try making a terminal exon
        
        my ($lend, $rend) = @$last_segment;
        if ($GENOME_FEATURES[ $lend - 2 ] == $ACCEPTOR) {
            print "Try to create terminal exon from $accession, $lend-$rend\n" if $DEBUG;
            my $got_terminal_exon_flag = &try_extend_terminal_segment_to_terminal_exon($last_segment, $ev_info_href);
            if ($got_terminal_exon_flag) {
                print "\tfound a terminal exon.  Adding weight $evidence_weight to end series at coord: $rend\n" if $DEBUG;
                $ends_aref->[$rend] += $evidence_weight;
            }
        }
    }
    
    return;
}

####
sub try_extend_initial_segment_to_initial_exon {
    my ($coords_aref, $ev_info_href) = @_;

    my ($end5, $end3) = @$coords_aref;
    
    my $accession = $ev_info_href->{accession};
    my $ev_type = $ev_info_href->{ev_type};
    
    my $added_exon_flag = 0;

    ## if cannot extend to start, try finding a more internal start:
    
    my $exon_found_flag = 0;

  MOCK_INIT_EXON:
    foreach my $nucleate_pos ($end5, $end3 - 3) {
        
        my $first_start_flag = ($nucleate_pos == $end5) ? "FIRST_START_PLEASE" : undef; ## find first start only if climbing from end5
        
        foreach my $start_frame (1..3) {
            
            my $candidate_start_pos = &extend_upstream_to_start_codon($nucleate_pos, $start_frame, $first_start_flag);
            
            if ($candidate_start_pos) {
               
                my $exon_name = "mock_initial_seg_extend_${candidate_start_pos}-${end3}_${accession}_${ev_type}";

                ## check to see that the entire orf is open:
                my @good_phases = &determine_good_phases($candidate_start_pos, $end3);
                if (grep { $_ == $start_frame } @good_phases) {
                    
                    ## make sure this doesn't already exist:
                    my $existing_exon = &get_exon_via_coords($candidate_start_pos, $end3, "initial", 1);
                    if ($existing_exon) { 
                        $existing_exon->appendEvidence($exon_name, $ev_type);
                    } 
                    else {
                        ## add it here:
                        
                        $mock_exon_counter++;
                        &add_exon($exon_name, $candidate_start_pos, $end3, "initial", 1, $ev_type);
                        $added_exon_flag = 1;
                    }
                    $exon_found_flag = 1;
                }
            }
        }
        
        if ($exon_found_flag) {
            last MOCK_INIT_EXON; ## if found an initial exon during extension from end5, then stop.
        }
    }
    
    return ($added_exon_flag);
    
}


sub try_intend_initial_segment_to_initial_exon {
    ## intend opposite of extend.  Walking inward looking for a start:
    my ($end5, $end3) = @_;

    my $adj_length = int( ($PERCENT_EXON_START_EXPLORE / 100) * ($end3 - $end5 + 1) );
    
    my $stop_search_coord = $end5 + $adj_length;
    
    foreach my $frame (1..3) {
        
        my $coord = $end5;
        $coord += $frame - 1;
        
        my $first_start;
        while ($coord < $stop_search_coord) {
            if ($GENOME_FEATURES[$coord] == $START) {
                $first_start = $coord;
                last;
            }
            $coord += 3; #codon hop
        }

        if ($first_start) {
            ## check for good phases.
            my @good_phases = grep {$_ == 1} &determine_good_phases($first_start, $end3);
            if (@good_phases) {
                
                $mock_exon_counter++;
                my $accession = "mock_initial_seg_intend_$mock_exon_counter";
                &add_exon($accession, $first_start, $end3, "initial", 1);
            }
        }
    }
    

    return;
}


####
sub try_extend_terminal_segment_to_terminal_exon {
    my ($coords_aref, $ev_info_href) = @_;
    
    my ($end5, $end3) = @$coords_aref;

    my $accession = $ev_info_href->{accession};
    my $ev_type = $ev_info_href->{ev_type};
    
    print "-trying to extend terminal segment $end5-$end3, $accession, $ev_type to terminal exon\n" if $DEBUG;
    
    my $segment_length = $end3 - $end5 + 1;
    
    my $max_exon_length = int ($segment_length + ( ($PERCENT_EXON_STOP_EXPLORE/100) * $segment_length));
    
    my $found_exon_flag = 0;

  MOCK_TERM_SEARCH:
    foreach my $nucleate_pos ($end3, $end5) {  ## start with end3, if don't find anything, try again from end5
        
        foreach my $end_frame (1..3) {

            my $new_end3 = &extend_downstream_to_stop($nucleate_pos, $end_frame);

            unless ($new_end3) {
                print STDERR "no stop codon reached from position: $nucleate_pos\n";
                next;
            }
            print "\textending to stop from $nucleate_pos, endframe: $end_frame, found stop at $new_end3\n" if $DEBUG;
            
            ## given stop position and end5 value plus end frame, compute the start frame.
            my $exon_length = ($new_end3 - $end5) + 1;
            
            if ($exon_length > $max_exon_length) {
                print "\t\tsorry, exon length: $exon_length exceeds max_exon_length: $max_exon_length\n" if $DEBUG;
                next; 
            }
            
            my $start_frame = $exon_length % 3;
            
            if ($start_frame == 0) {
                $start_frame = 1;
            }
            elsif ($start_frame == 1) {
                $start_frame = 3;
            }
            
            my $exon_name = "mock_term_seg_extend_${end5}-${new_end3}_${accession}_${ev_type}";

            if ($new_end3) {
                ## check for open orf
                my @good_phases = &determine_good_phases($end5, $new_end3 - 3); #trim stop codon from search!
                print "-extending $end5-$end3, start_frame=$start_frame, end_frame=$end_frame, new_end3=$new_end3, good_phases: @good_phases\n" if $DEBUG;
                
                if (grep { $_ == $start_frame } @good_phases) {
                    ## make sure this doesn't already exist:
                    my $existing_exon = &get_exon_via_coords($end5, $new_end3, "terminal", $start_frame);
                    if ($existing_exon) { 
                        print "\t-have this exon already: " . $existing_exon->toString() . "\n" if $DEBUG;
                        $existing_exon->appendEvidence($exon_name, $ev_type);
                    } 
                    else {
                        ## add it here:
                        print "\t-adding a new terminal exon.\n" if $DEBUG;
                        $mock_exon_counter++;
                        
                        &add_exon($exon_name, $end5, $new_end3, "terminal", $start_frame, $ev_type);
                    }
                    $found_exon_flag = 1;
                }
            }
            
        }
        if ($found_exon_flag) {
            last MOCK_TERM_SEARCH;
        }
    }
    return ($found_exon_flag);
    
}



####
sub get_intergenic_regions {
    my ($predictions_aref) = @_;
    
    print STDERR "-get_intergenic_regions()\n";


    if (scalar (@$predictions_aref) < 2) {
        # no intergenic regions to report.
        print STDERR "\tno intergenic regions to report.\n";
        return ();
    }
    
    my @preds = @$predictions_aref;
    
    @preds = sort {$a->{prediction_lend} <=> $b->{prediction_lend}} @preds;

    my $prev_pred = shift @preds;
   
    my @intergenic_regions;

    

    while (@preds) {
        my $next_pred = shift @preds;

        my ($prev_lend, $prev_rend) = $prev_pred->get_span();
        my ($next_lend, $next_rend) = $next_pred->get_span();
        
        print STDERR "INTERGENIC SPACE between ($prev_lend-$prev_rend) and ($next_lend-$next_rend)\n";
        

        my ($intergenic_lend, $intergenic_rend) = ($prev_rend + 1, $next_lend - 1);
        
        push (@intergenic_regions, [$intergenic_lend, $intergenic_rend]);
    
        $prev_pred = $next_pred;
    }

    print STDERR "\t-intergenic regions: " . Dumper (\@intergenic_regions);
    
    return (@intergenic_regions);
}


####
sub get_long_introns {
    my ($preds_aref) = @_;

    my @long_introns;
    
    foreach my $pred (@$preds_aref) {
        my @introns = $pred->get_intron_coords();
        foreach my $intron (@introns) {
            my ($intron_lend, $intron_rend) = sort {$a<=>$b} (@$intron);
            my $intron_length = $intron_rend - $intron_lend + 1;
            
            print STDERR "INTRON LENGTH: $intron_length from ($intron_lend - $intron_rend)\n";
            
            if ($intron_length >= $MIN_LONG_INTRON_LENGTH) {
                push (@long_introns, [$intron_lend, $intron_rend]);
            }
        }
    }

    return (@long_introns);
    
}




## some useful math functions:

sub avg {
    my @nums = @_;
    my $total = $#nums + 1;
    my $sum = 0;
    foreach my $num (@nums) {
        $sum += $num;
    }
    my $avg = $sum/$total;
    return ($avg);
}                 


sub stDev {
    # standard deviation calculation
    my @nums = @_;
    my $avg = avg(@nums);
    my $total = $#nums + 1;
    
    ## sum up the sqr of diff from avg
    my $sum_avg_diffs_sqr = 0;
    foreach my $num (@nums) {
        my $diff = $num - $avg;
        my $sqr = $diff**2;
        $sum_avg_diffs_sqr += $sqr;
    }
    my $stdev = sqrt ( (1/($total-1)) * $sum_avg_diffs_sqr);
    return ($stdev);
}                                     

sub median {
    my @nums = @_;

    unless (@nums) {
        return (0);
    }
    
    @nums = sort {$a<=>$b} @nums;
    
    my $num_eles = scalar (@nums);
        
    my $midpt = int (scalar (@nums) / 2);
    if ($num_eles % 2 == 0) {
        # even
        return ( avg($nums[$midpt-1], $nums[$midpt]) );
    }
    else {
        return ($nums[$midpt]);
    }
}


####
sub decrement_intergenic_for_evidence_spans {
    my @evidence_filenames;
    if ($transcriptAlignmentsFile) {
        push (@evidence_filenames, $transcriptAlignmentsFile);
    }
    if ($proteinAlignmentsFile) {
        push (@evidence_filenames, $proteinAlignmentsFile);
    }

    unless (@evidence_filenames) {
        return;
    }

    print STDERR "-decrementing intergenic scores for evidence spans.\n";
    
    my @evidence_chains;
    foreach my $evidence_filename (@evidence_filenames) {
        my @ev_chains = &parse_evidence_chains('?', $evidence_filename);
        push (@evidence_chains, @ev_chains);
    }
    
    ## determine gap sizes:
    my @gap_lengths = &determine_gap_lengths_from_evidence_chains(@evidence_chains);
    my $median_gap_length = median(@gap_lengths);
    my $max_gap_length = $INTRON_MEDIAN_FACTOR * $median_gap_length;

    @evidence_chains = &fragment_evidence_chains_using_max_gap($max_gap_length, \@evidence_chains);

    my $dec_int_fh;
    if ($DEBUG) {
        open ($dec_int_fh, ">decremented_intergenic_from_evspans.log") or die $!;
    }
    
    foreach my $evidence_chain (@evidence_chains) {
        my ($ev_type, $lend, $rend) = ($evidence_chain->{ev_type}, $evidence_chain->{lend}, $evidence_chain->{rend});
        if ($lend > $rend) { 
            print STDERR "Error, got evidence chain with lend > rend: ($lend > $rend) " . Dumper ($evidence_chain);
            next;
        }
        my $weight = $EVIDENCE_WEIGHTING{$ev_type};
        if (! defined ($weight)) {
            confess "Error, no weight for ev_type. " . Dumper($evidence_chain);
        }
        
        print $dec_int_fh "$lend-$rend\t$ev_type\t$weight\n" if $DEBUG;
   
        ## adjust intergenic scores in this region:
        for (my $i = $lend; $i <= $rend; $i++) {
            my $intergenic_score = $INTERGENIC_SCORES[$i];
            $INTERGENIC_SCORES[$i] = max(0, $intergenic_score - $weight);
        }
        
        
    }
    
    close $dec_int_fh if $DEBUG;
    
    return;
}

####
sub determine_gap_lengths_from_evidence_chains {
    my @evidence_chains = @_;
    
    my @gap_lengths;
    foreach my $evidence_chain (@evidence_chains) {
        my $gaps_aref = $evidence_chain->{gaps};
        foreach my $gap (@$gaps_aref) {
            my ($gap_lend, $gap_rend) = @$gap;
            my $gap_length = $gap_rend - $gap_lend + 1;
            push (@gap_lengths, $gap_length) if $gap_length >= $MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON;
        }
    }
    
    return (@gap_lengths);
}

####
sub decrement_coding_using_protein_alignment_introns {
    my ($genomic_strand) = @_;
    

    ## To Do:
    ## Should do this for all alignment gaps, not just introns, and should do it based on whole genome region stats independent of strand
    ######
    my @evidence_chains = &parse_evidence_chains('?', $proteinAlignmentsFile);
    ## determine gap sizes:
    my @gap_lengths = &determine_gap_lengths_from_evidence_chains(@evidence_chains);
    my $median_gap_length = median(@gap_lengths);
    my $max_gap_length = $INTRON_MEDIAN_FACTOR * $median_gap_length;
    
    ## get chains again, this time strand-specifically:
    @evidence_chains = &parse_evidence_chains($genomic_strand, $proteinAlignmentsFile);
    
    foreach my $chain (@evidence_chains) {
        my $accession = $chain->{accession};
        my $ev_type = $chain->{ev_type};
        my $ev_class = $chain->{ev_class};
        my $gaps_aref = $chain->{gaps};
        foreach my $gap (@$gaps_aref) {
            my ($end5, $end3) = @$gap;
            my $gap_length = abs ($end3 - $end5) + 1;
            if ($gap_length <= $max_gap_length) {
                ## negative addition to coding coverage for '(imperfect) introns' inferred from protein alignments
                &add_match_coverage($end5, $end3, -1 * $EVIDENCE_WEIGHTING{$ev_type}, $accession, $ev_type);
            }
        }
    }
    
    return;
}


####
sub max {
    my @eles = @_;

    @eles = sort {$a<=>$b} @eles;
    
    my $max_val = pop @eles;
    return ($max_val);
}

####
sub min {
    my @eles = @_;
    @eles = sort {$a<=>$b} @eles;
    my $min_val = shift @eles;

    return ($min_val);

}

####
sub augment_intergenic_from_start_stop_peaks {
    
    ## make sure exons are in order.
    @EXONS = sort {$a->{end5}<=>$b->{end5}} @EXONS;
    
    &augment_intergenic_from_start_peaks();
    &augment_intergenic_from_stop_peaks();
    
    return;
}

####
sub augment_intergenic_from_start_peaks {
    
    print "\nAugmenting intergenic region scores from START peaks:\n" if $DEBUG;

    my @plus_strand_starts;
    my @minus_strand_starts;
    foreach my $start_peak (@START_PEAKS) {
        ## find the nearest upstream terminal exon on same strand, or initial exon from opposite strand
        my ($position, $strand) = ($start_peak->{position}, $start_peak->{strand});
        
        if ($strand eq '+') {
            push (@plus_strand_starts, $start_peak);
        }
        else {
            push (@minus_strand_starts, $start_peak);
        }
    }


    
    my $logfh;
    
    if ($DEBUG) {
        open ($logfh, ">augment_intergenic_from_start_peaks.dat") or die $!;    
    }

    
    {
        ## process the plus strand first:
        foreach my $plus_strand_start (@plus_strand_starts) {
            my $nearest_relevant_exon = undef;
            my $position = $plus_strand_start->{position};
            print "Got plus strand start at $position.\n" if $DEBUG;
            my $closest_exon = &find_closest_exon_within_range("initial|single", "+", $position, $START_STOP_RANGE);
            
            unless ($closest_exon) {
                print "-no initial+ exon in range.\n" if $DEBUG;
                next;
            }
            # reset position to the end of this terminal exon.
            $position = &min($closest_exon->get_coords());
            
            foreach my $exon (reverse @EXONS) {  # walking from last to first
                my $exon_type = $exon->getExonType();
                my $strand = $exon->get_orientation();
                my $exon_rend = &max($exon->get_coords());

                ## first, must be left of the start
                unless ($exon_rend < $position) { next; }
                
                if ($strand eq '+' && $exon_type =~ /terminal|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
                if ($strand eq '-' && $exon_type =~ /initial|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
            }
            
            my $exon_rend;
            
            if ($nearest_relevant_exon) {
                ## max out the intergenic!
                $exon_rend = &max($nearest_relevant_exon->get_coords());
             
                print "-found relevant exon at position $exon_rend\n" if $DEBUG;
            }
            else {
                print "-no relevant exon found. Using beginning of sequence.\n" if $DEBUG;
                $exon_rend = 1;
            }

            
            print $logfh "+strand\t$exon_rend - $position\n" if $DEBUG;
            print "-augmenting intergenic from $exon_rend to $position (start peaks, plus strand starts)\n" if $DEBUG;
            
            for (my $i = $exon_rend; $i <= $position; $i++) {
                
                unless ($MASK[$i]) {
                    $INTERGENIC_SCORES[$i] = $SUM_GENEPRED_WEIGHTS;
                }
            }
        }
    }

    {

        ## Now, do the same for the minus strand.
        
        ## process the plus strand first:
        foreach my $minus_strand_start (@minus_strand_starts) {
            my $nearest_relevant_exon = undef;
            my $position = $minus_strand_start->{position};
            
            print "Got minus strand start at $position.\n" if $DEBUG;

            my $closest_exon = &find_closest_exon_within_range("initial|single", "-", $position, $START_STOP_RANGE);
            
            unless ($closest_exon) {
                print "-no initial- exon in range.\n" if $DEBUG;
                next;
            }
            # reset position to the end of this terminal exon.
            $position = &max($closest_exon->get_coords());
            
            print "-found closest initial minus strand exon at position $position\n" if $DEBUG;
            
            foreach my $exon (@EXONS) {  # walking from first to last
                my $exon_type = $exon->getExonType();
                my $strand = $exon->get_orientation();
                my $exon_lend = &min($exon->get_coords());

                ## first, must be right of the start
                unless ($exon_lend > $position) { next; }
                
                if ($strand eq '+' && $exon_type =~ /initial|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
                if ($strand eq '-' && $exon_type =~ /terminal|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
            }
         

            my $exon_lend;
            if ($nearest_relevant_exon) {
                ## max out the intergenic!
                $exon_lend = &min($nearest_relevant_exon->get_coords());
                
                print "-found relevant exon at position $exon_lend\n" if $DEBUG;
                
            }
            else {
                print "-no relevant exon found. Using end of sequence.\n" if $DEBUG;
                $exon_lend = $genomic_seq_length;
            }
            
            print $logfh "-strand\t$position - $exon_lend\n" if $DEBUG;
            print "-augmenting intergenic from $position to $exon_lend (minus strand starts)\n" if $DEBUG;
            
            
            for (my $i = $position; $i <= $exon_lend; $i++) {
                
                unless ($MASK[$i]) {
                    $INTERGENIC_SCORES[$i] = $SUM_GENEPRED_WEIGHTS;
                }
            }
        }
    }

    close $logfh if $DEBUG;
    
    return;
}


####
sub find_closest_exon_within_range {
    my ($exon_type, $strand, $position, $range) = @_;
    
    ## search all exons, find one of the same type, strand, and closest to position within the range.

    my $closest_exon = undef;
    my $closest_distance = undef;

    foreach my $exon (@EXONS) {
        
        my $orient = $exon->get_orientation();
        my $type = $exon->getExonType();
        
        unless ($type =~ /$exon_type/ && $strand eq $orient) { next; }

        my ($exon_end5, $exon_end3) = $exon->get_coords();
        
        my $delta = &min (abs ($position - $exon_end5), abs ($position - $exon_end3));
        
        unless ($delta <= $range) { next; }

        if (defined($closest_distance)) {
            if ($closest_distance > $delta) {
                $closest_distance = $delta;
                $closest_exon = $exon;
            }
        }
        else {
            # init to first found in range
            $closest_distance = $delta;
            $closest_exon = $exon;
        }
    }
    
    return ($closest_exon);
}


####
sub augment_intergenic_from_stop_peaks {
    
    print "\nAugmenting intergenic region scores from STOP peaks:\n" if $DEBUG;

    my @plus_strand_ends;
    my @minus_strand_ends;
    foreach my $end_peak (@END_PEAKS) {
        ## find the nearest upstream initial exon on same strand, or terminal exon from opposite strand
        my ($position, $strand) = ($end_peak->{position}, $end_peak->{strand});
        
        if ($strand eq '+') {
            push (@plus_strand_ends, $end_peak);
        }
        else {
            push (@minus_strand_ends, $end_peak);
        }
    }

    my $logfh;
    
    if ($DEBUG) {
        open ($logfh, ">augment_intergenic_from_stop_peaks.dat") or die $!;    
    }
    
    

    {
        ## process the plus strand first:
        foreach my $plus_strand_end (@plus_strand_ends) {
            
            my $nearest_relevant_exon = undef;
            my $position = $plus_strand_end->{position};
            print "Got plus strand end at $position.\n" if $DEBUG;
            my $closest_exon = &find_closest_exon_within_range("terminal|single", "+", $position, $START_STOP_RANGE);
            
            unless ($closest_exon) {
                print "-no terminal+ exon in range.\n" if $DEBUG;
                next;
            }
            # reset position to the end of this terminal exon.
            $position = &max($closest_exon->get_coords());
            
            foreach my $exon (@EXONS) {  # walking from first to last
                my $exon_type = $exon->getExonType();
                my $strand = $exon->get_orientation();
                my $exon_lend = &min($exon->get_coords());

                ## first, must be left of the start
                unless ($exon_lend > $position) { next; }
                
                if ($strand eq '+' && $exon_type =~ /initial|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
                if ($strand eq '-' && $exon_type =~ /terminal|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
            }
            

            my $exon_lend;
            if ($nearest_relevant_exon) {
                ## max out the intergenic!
                $exon_lend = &min($nearest_relevant_exon->get_coords());
               
                print "-found relevant exon at position $exon_lend\n" if $DEBUG;
            }
            else {
                print "-no relevant exon found, using end of sequence.\n" if $DEBUG;
                $exon_lend = $genomic_seq_length;
            }


            print $logfh "+strand\t$position - $exon_lend\n" if $DEBUG;
            print "-augmenting intergenic from $position to $exon_lend (stop peaks, plus strand)\n" if $DEBUG;
            
            for (my $i = $position; $i <= $exon_lend; $i++) {
                
                unless ($MASK[$i]) {
                    $INTERGENIC_SCORES[$i] = $SUM_GENEPRED_WEIGHTS;
                }
            }
        }
    }

    {

        ## Now, do the same for the minus strand.
        
        ## process the plus strand first:
        foreach my $minus_strand_end (@minus_strand_ends) {
            my $nearest_relevant_exon = undef;
            my $position = $minus_strand_end->{position};
            
            print "Got minus strand end at position $position\n" if $DEBUG;

            my $closest_exon = &find_closest_exon_within_range("terminal|single", "-", $position, $START_STOP_RANGE);
            
            unless ($closest_exon) {
                print "-no terminal- exon in range.\n" if $DEBUG;
                next;
            }
            # reset position to the end of this terminal exon.
            $position = &min($closest_exon->get_coords());
            
            foreach my $exon (reverse @EXONS) {  # walking from first to last
                my $exon_type = $exon->getExonType();
                my $strand = $exon->get_orientation();
                my $exon_rend = &max($exon->get_coords());

                ## first, must be right of the start
                unless ($exon_rend < $position) { next; }
                
                if ($strand eq '+' && $exon_type =~ /terminal|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
                if ($strand eq '-' && $exon_type =~ /initial|single/) {
                    # got it!
                    $nearest_relevant_exon = $exon;
                    last;
                }
            }
    
            
            my $exon_rend;
            if ($nearest_relevant_exon) {
                ## max out the intergenic!
                $exon_rend = &max($nearest_relevant_exon->get_coords());
                
                print "-found relevant exon at position $exon_rend\n" if $DEBUG;
            }
            else {
                print "-no relevant exon found, using beginning of sequence.\n" if $DEBUG;
                $exon_rend = 1;
            }

            print $logfh "-strand\t$exon_rend - $position\n" if $DEBUG;
            print "-augmenting intergenic from $exon_rend to $position (stop peaks, minus strand)\n" if $DEBUG;
            
            for (my $i = $exon_rend; $i <= $position; $i++) {
                                        
                unless ($MASK[$i]) {
                    $INTERGENIC_SCORES[$i] = $SUM_GENEPRED_WEIGHTS;
                }
            }
        }
    }
    

    close $logfh if $DEBUG;

    
    return;
}


####
sub in_range_of_genomic_sequence {
    my ($lend, $rend) = @_;

    if ($lend < 1 || $rend < 1
        || $lend > $genomic_seq_length
        || $rend > $genomic_seq_length ) {

        return(0);
    }
    else {
        return(1);
    }
}

#############################################################
## Class exon ###############################################
#############################################################

package Exon;
use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my ($end5, $end3) = @_;
    
    ## check valid constructor call
    unless ($end5 >= 0 && $end3 >= 0) {
        confess "Exon::new()  ERROR, need valid end5 and end3 values in constructor.  Currently, end5: $end5, end3: $end3\n\n";
    }
    
    ## frame info is in frames 1, 2, or 3 / index directly via the frame attributes.
    
    my $self = {
        end5 => $end5,
        end3 => $end3,
        type => undef,             #  initial, internal, terminal, single, bound
        startFrame => undef,       #  allowable frame of exon
        endFrame => undef,         #  based on startFrame and exon length.
        
        orientation => '+',        # set to +|-   (+ default)
        
        evidence_tokens => {},             # accessions for evidence supporting entry 
        evidence_list => [],  ## contains [ [acc,ev_type], ...]
        
        leftSeqBoundary => undef,  ## need to check for stop codon creation across splice boundaries.
        rightSeqBoundary=> undef,
        
        
        ## scores for DP functions
        base_score => 0,    # base score = (num evidence) + splice_site_scores + start_score + stop_score
        sum_score => 0,     # score of the highest scoring path ending at this node
        link => 0,          # reference to highest scoring path to the left of this exon
        
        
        ## base score composition:
        coding_score => 0,

        ## masking info:
        masked_percent => 0, #init
        
        
    };
    
    bless ($self, $packagename);
    return ($self);
}


####
sub setLeftSeqBoundary {
    my $self = shift;
    my $seqBound = shift;
    
    $self->{leftSeqBoundary} = $seqBound;
}


####
sub setRightSeqBoundary {
    my $self = shift;
    my $seqBound = shift;
    
    $self->{rightSeqBoundary} = $seqBound;
}



#### public method, package Exon
sub setStartFrame {
    my $self = shift;
    my $start_frame = shift;
    unless ($start_frame > 0 && $start_frame <= 3) { 
        confess "Exon::setStartFrame ($start_frame) not allowed.\n"; 
    }
    
    my $length = $self->length();
    $self->{startFrame} = $start_frame;
    
    ## end frame calculation:
    my $endframe = ($length + $start_frame - 1) % 3;
    if ($endframe == 0) { $endframe = 3;}
    $self->{endFrame} = $endframe;
    
}

#### 
sub getStartFrame {
    my $self = shift;
    return ($self->{startFrame});
}

####
sub getEndFrame {
    my $self = shift;
    return ($self->{endFrame});
}


#### private method, package Exon
sub setEndFrame {
    my $self = shift;
    my $end_frame = shift;
    $self->{endFrame}->[$end_frame] = 1;
}

#### public method, package Exon
sub setExonType {
    my $self = shift;
    my $type = shift;
    
    unless ($type eq "initial" || $type eq "internal" || $type eq "terminal" || $type eq "single" || $type eq "bound") {
        confess "Exon::setExonType()  Error, type($type) isn't recognized.\n";
    }
    $self->{type} = $type;
}

sub getExonType {
    my $self = shift;
    return ($self->{type});
}

sub set_masked_percent {
    my $self = shift;
    my $percent = shift;

    $self->{masked_percent} = $percent;
}

sub get_masked_percent {
    my $self = shift;
    return ($self->{masked_percent});
}


#### public method, package Exon
sub length {
    my $self = shift;
    my $length = abs ($self->{end3} - $self->{end5}) + 1;
    return ($length);
}


sub appendEvidence {
    my $self = shift;
    my $evidence_acc = shift;
    my $ev_type = shift;
    
    my $token = "{$evidence_acc;$ev_type}";
    
    unless ($self->{evidence_tokens}->{$token}) {
        $self->{evidence_tokens}->{$token} = 1;
        # add it
        push (@{$self->{evidence_list}}, [ $evidence_acc, $ev_type ]);
    }
    return;
}


sub get_coords {
    my $self = shift;
    return ($self->{end5}, $self->{end3});
}


sub get_evidence {
    my $self = shift;
    return (@{$self->{evidence_list}});
}


sub get_orientation {
    my $self = shift;
    return ($self->{orientation});
}


sub get_exon_score {
    my $self = shift;
    
    my $score = $self->{coding_score};
    
    return ($score);
}
        


sub toString () {
    my $self = shift;
    my $text = $self->{end5} 
    . "\t" . $self->{end3} 
    . "\t" . $self->{type} . $self->{orientation} 
    . "\t" . $self->{startFrame} 
    . "\t" . $self->{endFrame} 
    . "\t";

    my @evidence = $self->get_evidence();
    foreach my $ev (@evidence) {
        my ($acc, $class) = @$ev;
        $text .= "{$acc;$class},";
    }
    chop $text; #remove trailing comma
    
    return ($text);
}



##########################################################################################
##########################################################################################


package EVM_prediction;
use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    
    my @exon_list = @_;
    
    
    my $self = { prediction_score => undef,
                 prediction_lend => undef,
                 prediction_rend => undef,
                 prediction_orient => undef,
                 prediction_exons => [], 
                 intron_coords => [],
                 eliminate => 0, # flag to indicate that this prediction should be eliminated.
                 noncoding_equivalent => -1,
                 
                 };

    bless ($self, $packagename);
    
    $self->_init(@exon_list);

    return ($self);
}


####
sub _init {
    my $self = shift;
    my @exon_list = @_;
    
    #print "Creating prediction object:\n";

    @exon_list = sort {$a->{end5}<=>$b->{end5}} @exon_list;

    # copy sorted list of exons to the gene prediction
    @{$self->{prediction_exons}} = @exon_list;

    my $orient = $exon_list[0]->get_orientation();
    $self->{prediction_orient} = $orient;

    ## get the intron scores:
    my $intron_score = 0;
    
    for (my $i = 1; $i <= $#exon_list; $i++) {
        my $indiv_intron_score = $self->_get_exon_pair_intron_score($exon_list[$i-1], $exon_list[$i]);
        $intron_score += $indiv_intron_score;
        #print "Adding intron score: $indiv_intron_score\n";
    }

    ## add up the rest of the prediction score by summing exon scores
    my $prediction_total_score = $intron_score;
    my @coords; # set up gene bounds too!

    foreach my $exon (@exon_list) {
        my $exon_score = $exon->get_exon_score();
        $prediction_total_score += $exon_score;
        #print "Adding exon score: $exon_score\n";
        push (@coords, $exon->get_coords());
    }

    #print "Total score: $prediction_total_score\n";
    
    $self->{prediction_score} = $prediction_total_score;
    
    @coords = sort {$a<=>$b} @coords;
    my $lend = shift @coords;
    my $rend = pop @coords;
    
    $self->{prediction_lend} = $lend;
    $self->{prediction_rend} = $rend;
}

sub replace_exons {
    my $self = shift;
    my @replacement_exons = @_;
    
    $self->_init(@replacement_exons);

    return;
}


sub eliminate {
    my $self = shift;
    $self->{eliminate} = 1;
}

sub is_eliminated {
    my $self = shift;
    return ($self->{eliminate});
}



sub get_score {
    my $self = shift;
    return ($self->{prediction_score});
}


sub get_orientation {
    my $self = shift;
    return ($self->{prediction_orient});
}

sub get_span {
    my $self = shift;
    return ($self->{prediction_lend}, $self->{prediction_rend});
}


sub set_noncoding_equivalent {
    my $self = shift;
    my $noncoding_value = shift;
    $self->{noncoding_equivalent} = $noncoding_value;
}


####
sub _get_exon_pair_intron_score {
    my $self = shift;
    my $exonA = shift;
    my $exonB = shift;

    my $orient = $self->get_orientation();
    
    my ($exonA_5prime, $exonA_3prime) = $exonA->get_coords();
    my ($exonB_5prime, $exonB_3prime) = $exonB->get_coords();

    my ($intron_5prime, $intron_3prime);
    if ($orient eq '+') {
        ($intron_5prime, $intron_3prime) = ($exonA_3prime+1, $exonB_5prime-2);
    } else {
        ($intron_5prime, $intron_3prime) = ($exonB_3prime-1, $exonA_5prime+2);
    }

    my $intron_key = "$intron_5prime" . "_" . "$intron_3prime";
    my $intron_score = $INTRONS_TO_SCORE{$intron_key};
    unless (defined $intron_score) {
        confess "Error, no intron available based on $intron_5prime-$intron_3prime\n" . $self->toString();
    }

    
    ## get simple intron coords:
    my ($exonA_lend, $exonA_rend) = sort {$a<=>$b} ($exonA_5prime, $exonA_3prime);
    my ($exonB_lend, $exonB_rend) = sort {$a<=>$b} ($exonB_5prime, $exonB_3prime);
    
    my ($intron_lend, $intron_rend) = ($exonA_rend + 1, $exonB_lend - 1);
    
    ## add intron if it's not already there:
    ## TODO: separate this into a separate function that's called only once...

    my $found_intron = 0;
    foreach my $intron_set (@{$self->{intron_coords}}) {
        my ($lend, $rend) = @$intron_set;
        if ($lend == $intron_lend && $rend == $intron_rend) {
            $found_intron=1;
            last;
        }
    }
    unless ($found_intron) {
        ## store it.
        push (@{$self->{intron_coords}}, [$intron_lend, $intron_rend]);
    }
    
    return ($intron_score);
}

sub get_exons {
    my $self = shift;
    return (@{$self->{prediction_exons}});
}

sub get_intron_coords {
    my $self = shift;
    return (@{$self->{intron_coords}});
}


sub is_complete_gene {
    my $self = shift;
    if ($self->is_5prime_partial() || $self->is_3prime_partial()) {
        return(0);
    }
    else {
        return (1);
    }
}


####
sub is_5prime_partial {
    my $self = shift;
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        if ($exon->getExonType() =~ /initial|single/) {
            return (0);
        }
    }

    return (1); # lacks initial exon
}

####
sub is_3prime_partial {
    my $self = shift;
    
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        if ($exon->getExonType() =~ /terminal|single/) {
            return(0);
        }
    }

    return (1); # lacks terminal exon
}
    

####
sub get_coding_length {
    my $self = shift;
    
    my $sum_exon_length = 0;
    
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        my $length = $exon->length();
        $sum_exon_length += $length;
    }

    return ($sum_exon_length);
}


####
sub toString {
    my $self = shift;
    
    my $text = "# EVM prediction: " . "Mode:" . $self->{mode} 
    . " S-ratio: " . $self->{score_ratio} . " "
    . $self->{prediction_lend} . "-" . $self->{prediction_rend} 
    . " orient(" . $self->{prediction_orient} . ") "
        . "score(" . sprintf ("%.2f", $self->{prediction_score}) . ") " 
        . "noncoding_equivalent(" . sprintf ("%.2f", $self->{noncoding_equivalent}) . ")" . 
        " raw_noncoding(" . sprintf ("%.2f", $self->{raw_noncoding}) . ") "
        . "offset(" . sprintf ("%.2f", $self->{offset_noncoding}) . ") " ;
    
    if ($self->is_eliminated()) {
        $text .= " *** ELIMINATED *** ";
    }
    
    $text .= "\n";
    
    my $orient = $self->{prediction_orient};
    
    my @components;
    ## order introns and exons.
    my @introns = $self->get_intron_coords();
    my @exons = $self->get_exons();

    foreach my $intron (@introns) {
        my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$intron;
        my ($intron_end5, $intron_end3) = ($intron_lend, $intron_rend);

        if ($orient eq '-') {
            ($intron_end5, $intron_end3) = ($intron_end3, $intron_end5);
        }
        
        # create intron key:
        my ($key_end5, $key_end3) = ($intron_end5, $intron_end3);
        if ($orient eq '+') {
            $key_end3--;
        }
        else {
            $key_end3++;
        }
        
        my $intron_key = "${key_end5}_${key_end3}";
        unless (exists $INTRONS_TO_SCORE{$intron_key}) {
            confess "Error, intron key $intron_key doesn't exist";
        }
                                                       
        push (@components, { type => "intron",
                             coords => [$intron_end5, $intron_end3],
                             key => $intron_key,
                         } );
    }
    
    foreach my $exon (@exons) {
        my ($end5, $end3) = $exon->get_coords();
        push (@components, { type => "exon",
                             coords => [$end5, $end3],
                             exon => $exon,
                         } );
    }
    
    foreach my $component (sort {$a->{coords}->[0]<=>$b->{coords}->[0]} @components) {
        my $type = $component->{type};
        my $coords = $component->{coords};
        if ($type eq 'exon') {
            my $exon = $component->{exon};
            $text .= $exon->toString() . "\n";
        }
        else {
            # intron
            my $key = $component->{key};
            my ($intron_end5, $intron_end3) = @$coords;
            $text .= "$intron_end5\t$intron_end3\tINTRON\t\t\t";
            
            my $evidence_aref = $INTRONS_TO_EVIDENCE{$key} or confess "Error, no evidence associated with intron $key";
            foreach my $ev_pair (@$evidence_aref) {
                my ($acc, $type) = @$ev_pair;
                $text .= "{$acc;$type},";
            }
            chop $text;
            $text .= "\n";
        }
    }
    
    return ($text);
}





