#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $model_type = "SNAP";

my $usage = "usage: $0 SNAP.output\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
	my %data;

    print STDERR "-parsing $input_file\n";
	## parse input file
	open (my $fh, $input_file) or die "Error, cannot open file $input_file";
	while (<$fh>) {
        
        chomp;
		if (/^\#/) { 
			next; 
		}
		
		my @x = split(/\s+/);
        my $scaffold = $x[0];
        my $orient = $x[6];
        my $lend = $x[3];
        my $rend = $x[4];
        
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        my $model_id = $x[8];
        
        $data{$scaffold}->{$model_id}->{$end5} = $end3;
		
	}

	close $fh;


    print STDERR "-done parsing $input_file\n";
    
	## Generate gff3 output
	foreach my $scaffold (keys %data) {
		
		my $models_href = $data{$scaffold};

		foreach my $model (keys %$models_href) {

            print STDERR "-writing gff3 for $model\n";
            
			my $coords_href = $models_href->{$model};

			my $gene_obj = new Gene_obj();
			$gene_obj->populate_gene_object($coords_href, $coords_href);
			$gene_obj->{asmbl_id} = $scaffold;
			$gene_obj->{TU_feat_name} = "gene.$model";
			$gene_obj->{Model_feat_name} = "$model";
			$gene_obj->{com_name} = "$model_type prediction";
			
			print $gene_obj->to_GFF3_format(source => $model_type) . "\n";
		}

	}


	exit(0);

}
		
		
