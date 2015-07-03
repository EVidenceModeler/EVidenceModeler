#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../..//PerlLib");
use File::Basename;

use Fasta_reader;
use GTF_to_geneobjs;
use GFF3_to_geneobjs;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;
use BHStats;

umask(0000);

my $usage = <<_EOUSAGE_;

####################################################################################################
#  
#  --annot_file              genes file in corresponding format 
#  --format                  GTF or GFF3
#  -v                        verbose
#  -X                        eXport features
#####################################################################################################

_EOUSAGE_

    ;

my ($annot_file, $format, $help, $EXPORT_FLAG);

our $SEE;

&GetOptions ( 
			  "annot_file=s" => \$annot_file,
			  "format=s" => \$format,
			  "help|h" => \$help,
			  "v" => \$SEE,
			  "X" => \$EXPORT_FLAG,
              );

if ($help) { die $usage; }
unless ($annot_file && $format =~ /^(GTF|GFF3)$/) { 
    die $usage;
}

my $core_annot_filename = basename($annot_file);
my ($exonfh, $intronfh, $genicfh, $intergenicfh);
if ($EXPORT_FLAG) {
	open ($exonfh, ">$core_annot_filename.exons") or die $!;
	open ($intronfh, ">$core_annot_filename.introns") or die $!;
	open ($genicfh, ">$core_annot_filename.genes") or die $!;
	open ($intergenicfh, ">$core_annot_filename.intergenic") or die $!;
}

# stats interested in:
my $gene_count = 0;
my $mRNA_count = 0;
my $alt_spliced_gene_count = 0;
my $intron_containing_gene_count = 0;
my $unique_exon_count = 0;
my $unique_cds_count = 0;
my $unique_intron_count = 0;

my $alt_splice_diff_CDSs_count = 0;
my $diff_splice_CDS_count = 0;

my $num_intergenic_regions = 0;
my $sum_intergenic_lengths = 0;
my $sum_gene_lengths = 0;
my $sum_intron_lengths = 0;
my $sum_exon_lengths = 0;

main: {
        
    ## get the genes:
    my %gene_id_to_gene;
    my $seqname_map_href;
	
    if ($format eq 'GTF') {
        print STDERR "-processing GTF files\n" if $SEE;
        $seqname_map_href = &GTF_to_geneobjs::parse_file($annot_file, \%gene_id_to_gene);
	}
    else {
        print STDERR "-processing GFF3 files\n" if $SEE;
        $seqname_map_href = &GFF3_to_geneobjs::parse_file($annot_file, \%gene_id_to_gene);
	}

	
	foreach my $contig (keys %$seqname_map_href) {
        
        print STDERR "// processing $contig\n" if $SEE;
		my $gene_ids_aref = $seqname_map_href->{$contig};
		
		my %exons;
		my %cdss;
		my %introns;

		my @gene_spans;

		foreach my $gene_id (@$gene_ids_aref) {
                
			my $gene_obj = $gene_id_to_gene{$gene_id} or die "Error, no gene retrieved from reference via $gene_id";
                
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj->get_coords();
			push (@gene_spans, [$gene_lend, $gene_rend]);
			
			$gene_count++;
		
			
			my %complete_CDS_tokens; # tracking the full CDS structure; ignores alt isoforms w/ splice vars in UTRs only
	
			if ($gene_obj->get_additional_isoforms()) {
				$alt_spliced_gene_count++;
			}


			my $got_intron_containing_gene_flag = 0;
			
			foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {
				$mRNA_count++;

				my @all_cds_coords;
				my @exons = $isoform->get_exons();
				foreach my $exon (@exons) {
					my $exon_token = join ("_", $exon->get_coords());
					$exons{$exon_token} = 1;

					if (my $cds = $exon->get_CDS_exon_obj()) {
						my $cds_token = join ("_", $cds->get_coords());
						$cdss{$cds_token} = 1;
					
						push (@all_cds_coords, $cds->get_coords());
					}
				}

				if (my @introns = $isoform->get_intron_coordinates()) {
					$got_intron_containing_gene_flag = 1;
					foreach my $intron_coordset (@introns) {
						my $intron_token = join ("_", @$intron_coordset);
						$introns{$intron_token}++;
					}
				}

				my $complete_cds_token = join ("_", sort @all_cds_coords);
				$complete_CDS_tokens{$complete_cds_token} = 1;
				
			}

			if ($got_intron_containing_gene_flag) {
				$intron_containing_gene_count++;
			}


			my $num_cds_tokens = scalar (keys %complete_CDS_tokens);
			
			$diff_splice_CDS_count += $num_cds_tokens;
			if ($num_cds_tokens > 1) {
				$alt_splice_diff_CDSs_count++; # count alt-spliced genes ignoring UTR diffs only.
			}
			
		} # end of foreach gene


		&analyze_intergenics($contig, \@gene_spans);
		
		# get exon and CDS stats:
		my @unique_exons = keys %exons;
		my @unique_cdss = keys %cdss;
		my @unique_introns = keys %introns;

		$unique_exon_count += scalar(@unique_exons);
		$unique_cds_count += scalar(@unique_cdss);
		$unique_intron_count += scalar(@unique_introns);

		## get sum lengths:
		foreach my $exon (@unique_exons) {
			my ($exon_lend, $exon_rend) = sort {$a<=>$b} split (/_/, $exon);
			my $exon_len =  ($exon_rend - $exon_lend) + 1;
			$sum_exon_lengths += $exon_len;
			print $exonfh "exon\t$contig\t$exon_lend\t$exon_rend\t$exon_len\n" if $EXPORT_FLAG;
		}

		foreach my $intron (@unique_introns) {
			my ($intron_lend, $intron_rend) = sort {$a<=>$b} split (/_/, $intron);
			my $intron_len = $intron_rend - $intron_lend + 1;
			$sum_intron_lengths += $intron_len;
			print $intronfh "intron\t$contig\t$intron_lend\t$intron_rend\t$intron_len\n" if $EXPORT_FLAG;
		}
		
		
	} # end of foreach contig

	
	## summarize statistics:
	
	print "\n\n";
	printf ("%d genes\n", $gene_count);
	printf ("%d mRNAs\n", $mRNA_count);
	printf ("%d exons\n", $unique_exon_count);
	printf ("%d cdss\n", $unique_cds_count);
	printf ("%d introns\n", $unique_intron_count);
	
	printf ("\n%.1f exons per gene\n", $unique_exon_count / $gene_count);
	printf ("%.1f cdss per gene\n", $unique_cds_count / $gene_count);
	printf ("%.1f introns per gene\n", $unique_intron_count / $gene_count);
	
	print "\n";
	printf ("%d genes alternatively spliced\n", $alt_spliced_gene_count);
	printf ("%.1f%% genes alternatively spliced\n", $alt_spliced_gene_count / $gene_count * 100);
	print "\n";
	printf ("%.2f transcripts per gene\n", $mRNA_count / $gene_count);
	
	if ($alt_spliced_gene_count) {
		my $genes_not_alt_spliced = $gene_count - $alt_spliced_gene_count;
		printf ("%.2f transcripts per alt-spliced gene\n", ( $mRNA_count - $genes_not_alt_spliced ) / $alt_spliced_gene_count );
		print "\n\n";
	}

	if ($alt_splice_diff_CDSs_count) {
		
		printf ("%d genes alt spliced w/ altsplicing in coding regions.\n", $alt_splice_diff_CDSs_count);
		my $genes_not_alt_spliced = $gene_count - $alt_splice_diff_CDSs_count;
		
		printf("%.2f CDS-structures per alt-spliced gene\n\n\n", ($diff_splice_CDS_count - $genes_not_alt_spliced) / $alt_splice_diff_CDSs_count);
		
	}



	## Lengths:
	print "\n";
	printf("%.2f avg gene length\n", ($sum_gene_lengths / $gene_count));
	printf("%.2f avg exon length\n", ($sum_exon_lengths / $unique_exon_count));
	printf("%.2f avg intron length\n", ($sum_intron_lengths / $unique_intron_count));
	
	printf("%.2f avg bp between genes\n", ($sum_intergenic_lengths / $num_intergenic_regions));
	
	
}


if ($EXPORT_FLAG) {
	close $exonfh;
	close $intronfh;
	close $genicfh;
	close $intergenicfh;
}



exit(0);


			
####
sub analyze_intergenics {
	my ($contig, $gene_spans_aref) = @_;
	
	my @gene_spans = sort {$a->[0]<=>$b->[0]} @$gene_spans_aref;
	
	my @intergenics;

	my $gene_span = shift @gene_spans;
	my $prev_rend = $gene_span->[1];
	
	while (@gene_spans) {
		$gene_span = shift @gene_spans;
		my ($gene_lend, $gene_rend) = sort {$a<=>$b} @$gene_span;

		my $gene_length = $gene_rend - $gene_lend + 1;
		$sum_gene_lengths += $gene_length;

		print $genicfh "GENE\t$contig\t$gene_lend\t$gene_rend\t$gene_length\n";
		
		push (@intergenics, [$prev_rend + 1, $gene_lend - 1] );
		
		$prev_rend = $gene_rend;
	}
	
	foreach my $intergenic (@intergenics) {
		my ($lend, $rend) = @$intergenic;
		my $len = $rend - $lend + 1;
		
		if ($len > 0) {
			$sum_intergenic_lengths += $len;
			$num_intergenic_regions++;
			print $intergenicfh "INTERGENIC\t$contig\t$lend\t$rend\t$len\n" if $EXPORT_FLAG;
		}
	}

	return;
}



		
		
