#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $model_type = "GeneID";

my $usage = "usage: $0 geneID.gff\n\n";
$usage = "Note, GeneID must be run like so:\n"
    . "geneid -GP param/paramfile.param genomic.fasta > output.gff\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
	my %data;
		
	## parse input file
	open (my $fh, $input_file) or die "Error, cannot open file $input_file";
	while (<$fh>) {
		if (/^\#/) { next; }
        chomp;
        
        my @x = split(/\t/);
        my $scaff = $x[0];
        my $feat_type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $model_id = $x[8];

        if ($feat_type eq "exon") {
            
            ## add feature
            
            my $model_token = join("-", $scaff, $model_id);
            
            my $exon_feature = { scaff => $scaff,
                                 lend => $lend,
                                 rend => $rend,
                                 orient => $orient,
                                 model_id => $model_id,
            };

            push (@{$data{$model_token}}, $exon_feature);
        }
    }
    close $fh;


    ## write gff3 output
    foreach my $model_id (sort keys %data) {

        my @exons = @{$data{$model_id}};

        my $scaff = $exons[0]->{scaff};
        my $orient = $exons[0]->{orient};
        
        my %coords;

        foreach my $exon (@exons) {
            my $lend = $exon->{lend};
            my $rend = $exon->{rend};
            
            my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
            
            $coords{$end5} = $end3;
        }

        my $gene_obj = new Gene_obj();
        $gene_obj->populate_gene_object(\%coords, \%coords);
        
        $gene_obj->{asmbl_id} = $scaff;
        $gene_obj->{TU_feat_name} = "t.$model_id";
        $gene_obj->{Model_feat_name} = $model_id;
        $gene_obj->{com_name} = "GeneID prediction $model_id";
        
        print $gene_obj->to_GFF3_format(source => "GeneID") . "\n";

        
    }
    

    exit(0);
}
		
