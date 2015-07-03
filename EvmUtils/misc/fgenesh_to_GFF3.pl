#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $model_type = "FGENESH";

my $usage = "usage: $0 fgenesh.output\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
	my %data;

	my $scaffold;

	## parse input file
	open (my $fh, $input_file) or die "Error, cannot open file $input_file";
	while (<$fh>) {
		chomp;
		if (/Seq name: (\S+)/) {
			$scaffold = $1;
			next;
		}
		
		s/^\s+//g;
		my @x = split(/\s+/);
		unless (scalar (@x) > 7) { next; } # not a line of interest
		if ($x[3] =~ /^CDS\w/ && $x[1] =~ /[\+\-]/) {
			my $model_num = $x[0];
			my $orient = $x[1];
			my $lend = $x[4];
			my $rend = $x[6];
			
			my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

			$data{$scaffold}->{$model_num}->{$end5} = $end3;
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
			$gene_obj->{TU_feat_name} = "gene.$scaffold.$model";
			$gene_obj->{Model_feat_name} = "model.$scaffold.$model";
			$gene_obj->{com_name} = "$model_type prediction";
			
			print $gene_obj->to_GFF3_format(source => $model_type) . "\n";
		}

	}


	exit(0);

}
		
		
