#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib $ENV{EUK_MODULES};
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;

my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;


my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
  
  my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
  
  foreach my $gene_id (@gene_ids) {
    my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
    
    foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
      
      my $isoform_id = $isoform->{Model_feat_name};
      
      my $orient = $isoform->get_orientation();
      my $source = $isoform->{source} or die "Error, no source populated" . $isoform->toString();
      
      my ($model_lend, $model_rend) = sort {$a<=>$b} $isoform->get_model_span();
      
      print join ("\t", $asmbl_id, $source, "mRNA", $model_lend, $model_rend, ".", $orient, ".", "GenePrediction $isoform_id") . "\n";
      
      foreach my $exon ($isoform->get_exons()) {
	if (my $cds = $exon->get_CDS_exon_obj()) {
	  my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds->get_coords();
	  
	  print join ("\t", $asmbl_id, $source, "CDS", $cds_lend, $cds_rend, ".", $orient, ".", "GenePrediction $isoform_id") . "\n";
	}
      }
    }
  }
}


exit(0);

