#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use CdbTools;

my $usage = "\nusage: $0 jigsaw_gff genome.fasta\n\n";

my $jigsaw_gff = $ARGV[0] or die $usage;
my $genome = $ARGV[1] or die $usage;


main: {
	my %data;
	open (my $fh, $jigsaw_gff) or die "Error, cannot open fiel $jigsaw_gff";
	while (<$fh>) {
		chomp;
		unless (/\w/) { next;}
		if (/^\#/) { next; }
		
		chomp;
		my ($contig, $jigsaw, $exon_type, $lend, $rend, $score, $orient, $phase, $gene_info) = split (/\t/);
		
		$gene_info =~ /^([^;]+)/ or die "Error, cannot parse gene identifier from gene info: $gene_info";
		
		my $gene_id = $1;
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		$data{$contig}->{$gene_id}->{$end5} = $end3;
		
	}
	close $fh;

	foreach my $asmbl_id (keys %data) {
		my $genome_seq = &cdbyank_linear($asmbl_id, $genome);
		
		foreach my $gene_id (keys %{$data{$asmbl_id}}) {

			my $coords_href = $data{$asmbl_id}->{$gene_id};
			
			my $gene_obj = new Gene_obj();
			$gene_obj->populate_gene_object($coords_href, $coords_href, \$genome_seq);
			$gene_obj->{asmbl_id} = $asmbl_id;
			$gene_obj->{TU_feat_name} = "gene.$gene_id";
			$gene_obj->{Model_feat_name} = "model.$gene_id";
			$gene_obj->{source} = "Jigsaw";
			$gene_obj->{com_name} = "jigsaw prediction";

			my $protein = $gene_obj->get_protein_sequence();

			print $gene_obj->to_GFF3_format();
			print "#PROT $gene_id $protein\n\n";
			
		}
	}

	exit(0);
}


