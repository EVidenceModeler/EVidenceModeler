#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
BEGIN {
    use lib "$FindBin::Bin/lib";
}
use TkAnnotationViewer;
use Getopt::Long qw(:config no_ignore_case bundling);


our $SEE = 0;

my $params = <<_EOPARAMS_;

################# Evidence Modeler ##############################
#
#
# Opts for evidence modeler:
#
#  --genome                      :genome sequence in fasta format
#  --gene_predictions_listing    :list of gene predictions files 
#  --alignments_listing          :list of alignment evidence files (proteins, transcripts, repeats, etc).
#
#################################################################

_EOPARAMS_

    ;


my ($genomicSeqFile, $gene_predictions_files, $alignments_files);

&GetOptions ( 
              ## evm opts
              "genome=s"=>\$genomicSeqFile,
              "gene_predictions_listing=s"=>\$gene_predictions_files,
              "alignments_listing=s"=>\$alignments_files,
              );

main: {
    
    unless ($genomicSeqFile && ($gene_predictions_files || $alignments_files) ) {
        die $params;
    }
    
    
    my %file_listings = (  genome_seq_file => $genomicSeqFile,
                           gene_predictions_files => [],
                           search_evidence_files => [],
                           );
    
    if ($gene_predictions_files) {
        $gene_predictions_files =~ s/\s//g;
        my @files = split (/,/, $gene_predictions_files);
        my $file_list_aref = $file_listings{gene_predictions_files};
        push (@$file_list_aref, @files);
    }

    if ($alignments_files) {
        $alignments_files =~ s/\s//g;
        my @files = split (/,/, $alignments_files);
        my $file_list_aref = $file_listings{search_evidence_files};
        push (@$file_list_aref, @files);
    }
    
    
    my $Annotation_viewer = new TkAnnotationViewer();
    $Annotation_viewer->set_up_gui();
    
    $Annotation_viewer->{clone_viewer}->display_gff3_annots( \%file_listings);
    
    Tk::MainLoop;
}

exit(0);

