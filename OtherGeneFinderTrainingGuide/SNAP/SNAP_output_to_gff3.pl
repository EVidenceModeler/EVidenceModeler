#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use CdbTools;

my $usage = "usage: $0 snap_output_file genome_fasta\n\n";

my $snap_out = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;

my %data;
my $asmbl_id = "";

open (my $fh, $snap_out) or die $!;
while (<$fh>) {
	chomp;
    if (/>(\S+)/) {
        $asmbl_id = $1;
    }
    else {
        my @x = split (/\t/);
        if (scalar (@x) != 9) {
            die "Error, line $_ cannot be parsed.\n";
        }
        my ($type, $end5, $end3, $orient, $score, $a, $b, $c, $id) = @x;
        if ($orient eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }

        $data{$asmbl_id}->{$id}->{$end5} = $end3;
    }
}
close $fh;

foreach my $asmbl_id (sort keys %data) {
    
    my $genome_seq = &cdbyank_linear ($asmbl_id, $genome_fasta);
    
    my @gene_objs;
    foreach my $gene_id (keys %{$data{$asmbl_id}}) {
        my $coord_href = $data{$asmbl_id}->{$gene_id};
		
        my $gene_obj = new Gene_obj ();
        $gene_obj->populate_gene_obj($coord_href, $coord_href, \$genome_seq);
		$gene_obj->{TU_feat_name} = $gene_id;
		$gene_obj->{Model_feat_name} = "mRNA.$gene_id";
		$gene_obj->{com_name} = "SNAP pred $gene_id";
		$gene_obj->{asmbl_id} = $asmbl_id;
		
        print $gene_obj->to_GFF3_format() . "\n";
		
        my $protein = $gene_obj->get_protein_sequence();
        print "## $gene_id: $protein\n";
        
	}
	
}


exit(0);

