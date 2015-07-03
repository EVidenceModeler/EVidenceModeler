#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use Carp;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Data::Dumper;
use File::Basename;

umask(0000);

my $param_string = <<__PARAMS__;

###############################################################################################

The evidence for each contig is partitioned into a separate contig directory.

If the contig is shorter than the segmentSize value, the data will be further partitioned 
into segmentSize data sets.

The output file (specified by --partition_listing) will contain a list of all directories containing
the partitioned data and the partitioned data ranges.



###############################################################################################
#
# Options:
#  --genome                * :fasta file containing all genome sequences
#  --gene_predictions      * :file containing gene predictions
#  --protein_alignments      :file containing protein alignments
#  --transcript_alignments   :file containing transcript alignments
#  --pasaTerminalExons       :file containing terminal exons based on PASA long-orf data.
#  --repeats                 :file containing repeats to be masked
#
#
#  --segmentSize           * :length of a single sequence for running EVM
#  --overlapSize           * :length of sequence overlap between segmented sequences
#
#  --partition_listing     * :name of output file to be created that contains the list of partitions
#  -h                        : help message with more detailed information.
#
#   (*required options)
#
################################################################################################

__PARAMS__

    ;

my ($genome, $gene_predictions, $protein_alignments, $transcript_alignments,
    $pasaTerminalExons, $repeats, $segmentSize, $overlapSize,
    $partition_listing, $help
    );


&GetOptions ('genome=s' => \$genome,
             'gene_predictions=s' => \$gene_predictions,
             'protein_alignments=s' => \$protein_alignments,
             'transcript_alignments=s' => \$transcript_alignments,
             'pasaTerminalExons=s' => \$pasaTerminalExons,
             'repeats=s' => \$repeats,
             'segmentSize=i' => \$segmentSize,
             'overlapSize=i' => \$overlapSize,
             'partition_listing=s' => \$partition_listing,
             'h' => \$help,
             );

# process required options
if ($help) { die $param_string; }

unless ($genome && $gene_predictions && $segmentSize && defined($overlapSize) && $partition_listing) {
    die $param_string;
}


my $genome_basename = basename($genome);
my $gene_predictions_basename = basename($gene_predictions);
my $protein_alignments_basename = basename($protein_alignments) if $protein_alignments;
my $transcript_alignments_basename = basename($transcript_alignments) if $transcript_alignments;
my $pasaTerminalExons_basename = basename($pasaTerminalExons) if $pasaTerminalExons;
my $repeats_basename = basename($repeats) if $repeats;

my @files_to_process = ( { type => "gene_predictions",
                           basename => $gene_predictions_basename,
                           file => $gene_predictions,
                       } );
if ($protein_alignments) {
    push (@files_to_process, { type => "protein_alignments",
                               basename => $protein_alignments_basename,
                               file => $protein_alignments,
                           });
}
if ($transcript_alignments) {
    push (@files_to_process, { type => "transcript_alignments",
                               basename => $transcript_alignments_basename,
                               file => $transcript_alignments,
                           });
}
if ($pasaTerminalExons) {
    push (@files_to_process, { type => "pasaTerminalExons",
                               basename => $pasaTerminalExons_basename,
                               file => $pasaTerminalExons,
                           });
}

if ($repeats) {
    push (@files_to_process, { type => "repeats",
                               basename => $repeats_basename,
                               file => $repeats,
                           });
}


my $curr_dir = cwd();
my $util_dir = $FindBin::Bin;
my $gff_partitioner = "$util_dir/gff_range_retriever.pl";


&partition_files_based_on_contig(@files_to_process);


open (my $ofh, ">$partition_listing") or die "Error, cannot write $partition_listing file";

my $fasta_reader = new Fasta_reader($genome);
while (my $seq_obj = $fasta_reader->next()) {
    
    my $accession = $seq_obj->get_accession();
    my $sequence = $seq_obj->get_sequence();
    my $sequence_length = length($sequence);
    
    ## create accession partition area
    my $acc_dir = "$curr_dir/$accession";
    unless (-d $acc_dir) {
        mkdir ($acc_dir) or die "Error, cannot mkdir $accession";
    }
    
    ## write the genome sequence
    open (my $acc_ofh, ">$acc_dir/$genome_basename") or die $!;
    print $acc_ofh $seq_obj->get_FASTA_format();
    close $acc_ofh;
    
    ## pull out the relevant data:
    foreach my $file_struct (@files_to_process) {
        my $filename = $file_struct->{file};
        my $basename = $file_struct->{basename};
        my $acc_filename = "$acc_dir/$basename";
        $file_struct->{acc_file} = $acc_filename;

        unless (-e $acc_filename) {
            system "touch $acc_filename";
            if ($?) { die "Error, died running cmd: touch $acc_filename"; }
        }
        
        #my $cmd = "$gff_partitioner $accession 1 $sequence_length < $filename > $acc_filename";
        #&process_cmd($cmd);
    }
    
    my @ranges = &get_range_list($sequence_length);
    my $num_ranges = scalar (@ranges);
    print "$accession has $num_ranges partitions.\n";
    
    if ($num_ranges == 1) {
        # do not partition further. (as indicated below with 'N')
        print $ofh "$accession\t$acc_dir\tN\n";
    }
    else {
        ## partition the data into segments:
        foreach my $range (@ranges) {
            my ($range_lend, $range_rend) = @$range;
            my $range_length = $range_rend - $range_lend + 1;
            
            my $partition_dir = "$acc_dir/$accession" . "_$range_lend" . "-$range_rend";
            mkdir $partition_dir or die "Error, cannot mkdir $partition_dir";
            my $genome_subseq = substr($sequence, $range_lend - 1, $range_length);
            $genome_subseq =~ s/(\S{60})/$1\n/g; #make fasta format
            
            open (my $part_ofh, ">$partition_dir/$genome_basename") or die "Error, cannot write to $partition_dir/$genome_basename";
            print $part_ofh ">$accession\n";
            print $part_ofh $genome_subseq;
            close $part_ofh;
            
            foreach my $file_struct (@files_to_process) {
                my $filename = $file_struct->{acc_file};
                my $basename = $file_struct->{basename};
                my $cmd = "$gff_partitioner $accession $range_lend $range_rend ADJUST_TO_ONE < $filename > $partition_dir/$basename";
                &process_cmd($cmd);
            }
            print $ofh "$accession\t$acc_dir\tY\t$partition_dir\n";
        }
    }
}

close $ofh;


exit(0);



####
sub get_range_list {
    my ($sequence_length) = @_;

    my $range_lend = 1;
    my $range_rend = $segmentSize;
    my @ranges;
    
    while ($range_lend < $sequence_length - $overlapSize) {
        $range_rend = $range_lend + $segmentSize - 1;
        if ($range_rend > $sequence_length) {
            $range_rend = $sequence_length;
        }
        push (@ranges, [$range_lend, $range_rend]);

        $range_lend += $segmentSize - $overlapSize;
    }

    unless (@ranges) {
        ## must be a shorty such that the overlap length is greater than the sequence length.
        if ($sequence_length > $segmentSize) {
            die "Error, no ranges and sequence length ($sequence_length) does not exceed segmentSize ($segmentSize)";
        }
        push (@ranges, [1, $sequence_length]);
    }
    
    return (@ranges);
}



 

####
sub process_cmd {
    my $cmd = shift;
    print "CMD: $cmd\n";
    die "error, $cmd, $? " if system $cmd;
}

####
sub partition_files_based_on_contig {
    my @file_structs = @_;

    foreach my $file_struct (@file_structs) {

        my $file = $file_struct->{file};
        my $basename = $file_struct->{basename};
        
        my $ofh;
        my $curr_contig = "";
        open (my $fh, $file) or die "Error, cannot open $file";
        while (<$fh>) {
            my $line = $_;
            unless (/\w/) { next; }
            if (/^\#/) { next; } # comment lines
            my @x = split (/\t/);
            my $contig_id = $x[0];
            if ($contig_id ne $curr_contig) {
                $curr_contig = $contig_id;
                unless (-d "$curr_contig") {
                    print STDERR "writing $file for $curr_contig\n";
                    mkdir $curr_contig or die "Error, cannot mkdir $curr_contig";
                }
                close $ofh if $ofh;
                open ($ofh, ">>$curr_contig/$basename") or die "Error, cannot write to $curr_contig/$basename";
            }
            print $ofh $line;
        }
    }
    return;
}

