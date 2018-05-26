#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $model_type = "GeMoMa";

my $usage = "usage: $0 GeMoMa.output.gff\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
	my %data;

	my $scaffold;

	## parse input file
	open (my $fh, $input_file) or die "Error, cannot open file $input_file";
	while (<$fh>) {
		chomp;
        if (/^\#/) { next; }

        my $line = $_;
                
		my @x = split(/\s+/);
        if ($x[2] eq "CDS") {
            my $scaffold = $x[0];
            my $info = $x[8];

            my $gene_id;
            if ($info =~ /Parent=([^;]+)/) {
                $gene_id = $1;
            }
            else {
                die "Error, cannot extract Parent info from $line";
            }
                        
            
			my $lend = $x[3];
			my $rend = $x[4];
            
            my $orient = $x[6];
            
			my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

			$data{$scaffold}->{$gene_id}->{$end5} = $end3;
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
		
		
