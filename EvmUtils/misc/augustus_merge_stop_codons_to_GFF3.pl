#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $model_type = "Augustus";

my $usage = "usage: $0 augustus.gff.output\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
	my %data;

	my $scaffold = "";

	## parse input file
	open (my $fh, $input_file) or die "Error, cannot open file $input_file";
	while (<$fh>) {
		if (/^\#/) { next; }
		chomp;
		unless (/\w/) { next; }
		
		my @x = split(/\t/);
		if ($x[2] eq 'CDS' || $x[2] eq "stop_codon") { 
			my $scaffold = $x[0];
			my $orient = $x[6];
			my $lend = $x[3];
			my $rend = $x[4];
			
			my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

			my $info = $x[8];
			$info =~ /Parent=([^;\s]+)/;
			my $model = $1 or die "Error, cannot parse model from $info";
			
            #print "MODEL: [$model]\n";
            
			$data{$scaffold}->{$model}->{$end5} = $end3;
		}
	}

	close $fh;


	## Generate gff3 output
	foreach my $scaffold (keys %data) {
		
		my $models_href = $data{$scaffold};

		foreach my $model (keys %$models_href) {

			my $coords_href = $models_href->{$model};

			my $gene_obj = new Gene_obj();
			$gene_obj->populate_gene_object($coords_href, $coords_href);
			$gene_obj->{asmbl_id} = $scaffold;
			$gene_obj->{TU_feat_name} = "gene.$model";
			$gene_obj->{Model_feat_name} = "model.$model";
			$gene_obj->{com_name} = "$model_type prediction";
            
            $gene_obj->join_adjacent_exons();
            
			print $gene_obj->to_GFF3_format(source => $model_type) . "\n";

            # print $gene_obj->toString() . "\n\n";
            
		}

	}


	exit(0);

}
		
		
