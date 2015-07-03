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
		if ($x[2] eq 'CDS') { 
			my $scaffold = $x[0];
			my $orient = $x[6];
			my $lend = $x[3];
			my $rend = $x[4];
			
			my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

			my $info = $x[8];

            my $gene_id;
            my $trans_id;
            if ($info =~ /transcript_id \"([^\"]+)"; gene_id \"([^\"]+)\"/) {
                $trans_id = $1;
                $gene_id = $2;

                $trans_id = "$scaffold-$trans_id";
                $gene_id = "$scaffold-$gene_id";

            }
            else {
                die "Error, cannot parse gene_id and transcript_id from $_";
            }
            
            my $model = join("$;", $gene_id, $trans_id);
			
			$data{$scaffold}->{$model}->{$end5} = $end3;
		}
	}

	close $fh;


	## Generate gff3 output
	foreach my $scaffold (keys %data) {
		
		my $models_href = $data{$scaffold};

		foreach my $model (keys %$models_href) {

            my ($gene_id, $trans_id) = split(/$;/, $model);

			my $coords_href = $models_href->{$model};

			my $gene_obj = new Gene_obj();
			$gene_obj->populate_gene_object($coords_href, $coords_href);
			$gene_obj->{asmbl_id} = $scaffold;
			$gene_obj->{TU_feat_name} = $gene_id;
			$gene_obj->{Model_feat_name} = $trans_id;
			$gene_obj->{com_name} = "$model_type prediction";
            
			print $gene_obj->to_GFF3_format(source => $model_type) . "\n";

        }

	}


	exit(0);

}
		
		
