#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../../PerlLib");
use Gene_obj;
use CdbTools;

my $usage = "usage: $0 gff2_file genome_fasta\n\n";

my $gff2_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;


my %data;
my %gene_id_to_source;

open (my $fh, $gff2_file) or die "Error, cannot open file $gff2_file";
while (<$fh>) {
  chomp;
  my @x = split (/\t/);
  my ($asmbl_id, $source, $feat_type, $lend, $rend, $score, $orient, $phase, $gene_info) = @x;
  
  unless ($feat_type eq 'CDS') {
    next;
  }
  
  my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
  
  $gene_info =~ /(Transcript|GenePrediction) (\S+)/ or die "Error, cannot parse the Transcript identifier";
  my $gene_id = $2;
  $gene_id =~ s/\"//g;

  $gene_id_to_source{$gene_id} = $source;

  $data{$asmbl_id}->{$gene_id}->{$end5} = $end3;
  
}
close $fh;


foreach my $asmbl_id (keys %data) {
  
  my $sequence = &cdbyank_linear($asmbl_id, $genome_fasta);
  
  foreach my $gene_id (keys %{$data{$asmbl_id}}) {
    
    my $coords_href = $data{$asmbl_id}->{$gene_id};
    
    my $gene_obj = new Gene_obj();
    $gene_obj->populate_gene_object($coords_href, $coords_href, \$sequence);
    
    $gene_obj->{asmbl_id} = $asmbl_id;
    $gene_obj->{TU_feat_name} = $gene_id;
    $gene_obj->{Model_feat_name} = "$gene_id.model";
    $gene_obj->{com_name} = "predicted gene";
    $gene_obj->{source} = $gene_id_to_source{$gene_id};

    
    print $gene_obj->to_GFF3_format();
    my $protein = $gene_obj->get_protein_sequence();
    print "#PROT $gene_id $protein\n\n";
    
  }
}

exit(0);
    
    
