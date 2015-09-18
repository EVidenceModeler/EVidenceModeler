#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $usage = "\n\n\tusage: $0 GeneMark-ES_contig1.gff3\n\n\n";

my $genemark_file = $ARGV[0] or die $usage;


main: {

    my %gene_id_to_info;

    open (my $fh, $genemark_file) or die "Error, cannot open file $genemark_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        
        my @x = split(/\t/);
        my $scaffold = $x[0];
        
        my $feature_type = $x[2];
        if ($feature_type ne 'exon') {
            die "Error, the format of this file is expected to simply contain exon records only. You likely need a different converter";
        }
        my $lend = $x[3];
        my $rend = $x[4];

        my $orient = $x[6];
        
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

        my $info = $x[8];
        $info =~ /Parent=(\S+);/ or die "Error, cannot parse Parent tag from $info";

        my $gene_id = $1 or die "Error, no gene_id parsed";

        $gene_id = join("^", $scaffold, $gene_id); # just in case gene_ids aren't unique across scaffolds.

        $gene_id_to_info{$gene_id}->{chr} = $scaffold;
        push (@{$gene_id_to_info{$gene_id}->{CDS}}, [$end5, $end3]);

    }
    close $fh;

    foreach my $gene_id (keys %gene_id_to_info) {
        
        my $chr = $gene_id_to_info{$gene_id}->{chr};
        my @exons = @{$gene_id_to_info{$gene_id}->{CDS}};

        my $gene_obj = new Gene_obj();

        $gene_obj->populate_gene_object_via_CDS_coords(@exons);
        
        $gene_obj->{com_name} = "GeneMark-ES prediction";
        $gene_obj->{TU_feat_name} = "gene.$gene_id";
        $gene_obj->{Model_feat_name} = "model.$gene_id";
        $gene_obj->{asmbl_id} = $chr;

        print $gene_obj->to_GFF3_format();

    }

    exit(0);

}

        
        
        
