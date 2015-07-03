#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use CdbTools;

my $usage = "usage: $0 glimmerhmm_output genome_fasta\n\n";

my $glimmerhmm_output = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;

my %data;


my $asmbl_id;
open (my $fh, $glimmerhmm_output) or die "Error, cannot open $glimmerhmm_output";
while (<$fh>) {
    chomp;
	s/^\s+//;
    
	if (/^Sequence name: (\S+)/) {
		$asmbl_id = $1;
		next;
	}
	
	
	my @x = split (/\s+/);
    unless (scalar(@x) == 7) { next; }

	
	my ($model_no, $exon_no, $orient, $type, $lend, $rend, @rest) = @x;
    
    unless ($orient =~ /^[\+\-]$/) {
        next; 
    }
    
    if ($orient eq '-') {
        ($lend, $rend) = ($rend, $lend);
    }
	
    $data{$asmbl_id}->{$model_no}->{$lend} = $rend;
}
close $fh;

print STDERR "-done parsing inputs.\n";

foreach my $asmbl_id (keys %data) {
	
	print STDERR "-processing results for $asmbl_id\n";
	my $genome_seq = &cdbyank_linear($asmbl_id, $genome_fasta);
	
	foreach my $model (keys %{$data{$asmbl_id}}) {

		print STDERR "processing $model\n";

		
		
		my $coords_href = $data{$asmbl_id}->{$model};
		
		my $gene_obj = new Gene_obj();
		$gene_obj->populate_gene_object($coords_href, $coords_href, \$genome_seq);
		
		$model = "$asmbl_id.$model";
		$gene_obj->{TU_feat_name} = "glimmerHMM.$model";
		$gene_obj->{Model_feat_name} = "mRNA.glimmerHMM.$model";
		$gene_obj->{asmbl_id} = $asmbl_id;
		$gene_obj->{com_name} = "GlimmerHMM model $model";
		
		print $gene_obj->to_GFF3_format()  . "\n\n";
		
		print "#PROT $model " . $gene_obj->get_protein_sequence() . "\n\n";
		
	}
	


}


exit(0);
