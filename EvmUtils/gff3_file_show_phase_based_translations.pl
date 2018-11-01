#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Data::Dumper;

my $usage = "\n\nusage: $0 gff3_file genome_db\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;


my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = $genome{$asmbl_id} or die "Error, cannot find sequence for $asmbl_id"; #cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

    my @gene_objs;
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		$gene_obj_ref ->set_CDS_phases(\$genome_seq);
        push (@gene_objs, $gene_obj_ref);
    }
    @gene_objs = sort {$a->{gene_span}->[0]<=>$b->{gene_span}->[0]} @gene_objs;
    
    foreach my $gene_obj_ref (@gene_objs) {
        print $gene_obj_ref->to_GFF3_format() . "\n";

        print "# protein:\n" . $gene_obj_ref->get_protein_sequence() . "\n";
        print "# CDS\n" . $gene_obj_ref->get_CDS_sequence() . "\n";

        my $first_cds_exon = &get_first_cds_exon($gene_obj_ref);

        my $phase = $first_cds_exon->{phase};

        my $cds_seq = $gene_obj_ref->get_CDS_sequence();
        print "Start phase: $phase\n";
        my @chars = split(//, $cds_seq);

        my $prefix = "";
        if ($phase == 1) { 
            $prefix .= shift @chars;
            $prefix .= shift @chars;
        }
        elsif ($phase == 2) {
            $prefix .= shift @chars;
        }
        
        my @codon_strings;
        my @transl_strings;

        my $codon_string = "";
        my $transl_string = "";
        if ($prefix) {
            print "partial codon: $prefix\n";
        }
        while (@chars) {
            my $codon;
            for (1..3) {
                my $char = shift @chars;
                $codon .= $char if $char;
            }
            $codon_string .= "$codon ";
            my $transl = translate_sequence($codon, 1);
            $transl_string .= " $transl  ";

            if (length($codon_string) > 40) {
                push (@codon_strings, $codon_string);
                push (@transl_strings, $transl_string);
                $codon_string = "";
                $transl_string = "";
            }
        }

        if ($codon_string) {
            push (@codon_strings, $codon_string);
            push (@transl_strings, $transl_string);
        }

        while(@codon_strings) {
            my $codon_string = shift @codon_strings;
            my $transl_string = shift @transl_strings;

            print "$codon_string\n";
            print "$transl_string\n\n";
        }
    }
}


exit(0);



####
sub get_first_cds_exon {
    my ($gene_obj) = @_;

    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
        if (my $cds_exon = $exon->get_CDS_exon_obj()) {
            return($cds_exon);
        }
    }

    die "Error, no cds exon found for gene: " . Dumper($gene_obj);
}
