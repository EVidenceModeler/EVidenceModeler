#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $usage = "usage: $0 match.gff\n\n";

my $match_gff = $ARGV[0] or die $usage;

main: {


    my %genes;

    open (my $fh, $match_gff) or die "Error, cannot open file $match_gff";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        chomp;
       
        my @x = split(/\t/);
        my $scaff = $x[0];
        my $source = $x[1];
        my $feat_type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $info = $x[8];

        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);


        my %tags = &get_tags($info);
        
        if (my $parent = $tags{Parent}) {
            
            $genes{$parent}->{coords}->{$end5} = $end3;
        }
        elsif (my $id = $tags{ID}) {
            my $name = $tags{Name} || "No name";
            $genes{$id}->{name} = $name;
            $genes{$id}->{scaffold} = $scaff;
            $genes{$id}->{source} = $source;
        }
    }
    close $fh;

    foreach my $gene_id (keys %genes) {
        my $data_href = $genes{$gene_id};
        
        my $coords_href = $data_href->{coords};
      
        my $gene_obj = new Gene_obj();
        $gene_obj->populate_gene_object($coords_href, $coords_href);
        $gene_obj->{asmbl_id} = $data_href->{scaffold};
        $gene_obj->{com_name} = $data_href->{name};
        $gene_obj->{TU_feat_name} = $gene_id;
        $gene_obj->{Model_feat_name} = "m.$gene_id";
        
        print $gene_obj->to_GFF3_format(source => $data_href->{source});
        print "\n";
    }
        
    exit(0);

}

####
sub get_tags {
    my ($info) = @_;
    
    my %tags;
    
    my @fields = split(/;/, $info);
    foreach my $field (@fields) {
        my ($key, $val) = split(/=/, $field);
        
        $tags{$key} = $val;
    }

    return(%tags);
}

        
