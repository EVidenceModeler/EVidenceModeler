#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Data::Dumper;
use Gene_validator;

our $SEE = 0;

my $usage = "usage: $0 evm_output_file contig_accession \n\n";

my $evm_file = $ARGV[0] or die $usage;
my $contig_id = $ARGV[1] or die $usage;

# print STDERR "processing $evm_file, $contig_id\n";

my %model_num_to_coords;

my %model_id_to_ev_type;

my $model_id = 1;

my %end5_to_phase;

open (my $fh, $evm_file) or die "Error, cannot open $evm_file\n";
while (<$fh>) {
    if (/^\!/) { next; } ## comment line

    if (/^\#/) {
        my $ev_type;
        if (/ELIMINATED/) {
            $ev_type = "EVM_elm";
        } else {
            $ev_type = "EVM";
        }
        $model_id_to_ev_type{$model_id} = $ev_type;
        next;
    }
    
    chomp;
    
    unless (/\w/) {
        $model_id++;
        next;
    }
    
    
    my @x = split (/\t/);
    if (scalar (@x) == 6  && $x[0] =~ /^\d+$/ && $x[1] =~ /^\d+$/) {
        
        if ($x[2] eq 'INTRON') { next; }
        
        my $coords_ref = $model_num_to_coords{$model_id};
        unless (ref $coords_ref) {
            $coords_ref = $model_num_to_coords{$model_id} = {};
        }
        
        $coords_ref->{$x[0]} = $x[1];

        $end5_to_phase{$x[0]} = $x[3];
        
    }
}
close $fh;


my %phase_conversion = ( 1 => 0,
                         2 => 1,
                         3 => 2,
                         4 => 0,
                         5 => 1,
                         6 => 2
                         );



my @gene_objs;

foreach my $model_id (sort {$a<=>$b} keys %model_num_to_coords) {

    my $ev_type = $model_id_to_ev_type{$model_id};
    
    # if ($ev_type ne "EVM") { next; } #ignoring the EVM_elm for now.

    my $coords_ref = $model_num_to_coords{$model_id};
        
    my $gene_obj = new Gene_obj();
    $gene_obj->populate_gene_obj($coords_ref, $coords_ref);
    
    ## set phase:
    foreach my $exon ($gene_obj->get_exons()) {
        my $cds_obj = $exon->get_CDS_obj();
        my ($end5, $end3) = $cds_obj->get_coords();
        $cds_obj->set_phase( $phase_conversion{ $end5_to_phase{$end5} } );
    }
    
    $gene_obj->{Model_feat_name} = "evm.model.$contig_id.$model_id";
    $gene_obj->{TU_feat_name} = "evm.TU.$contig_id.$model_id";
    $gene_obj->{com_name} = "EVM prediction $contig_id.$model_id";
    
    $gene_obj->{asmbl_id} = "$contig_id";

    print $gene_obj->to_GFF3_format( source => $ev_type);
    print "\n"; # add spacer between genes
        
}

exit(0);

