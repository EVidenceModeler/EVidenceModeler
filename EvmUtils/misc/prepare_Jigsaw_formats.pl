#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

use CdbTools;
use Gene_obj_indexer;
use GFF3_utils;
use Getopt::Long qw(:config no_ignore_case bundling);


umask(0000);

my $usage = <<_EOUSAGE_;

####################################################################################################
#  
#  --template_gff3_file      gold standard gene structures in gff3 format
#  --genome                  genome sequence in multi-fasta format
#  --alignment_gff3_files     comma-separated list of gff3 files containing alignments
#  --gene_gff3_files         comma-separated list of gff3 files containing gene predictions.
#
# -v   verbose setting.
#
#####################################################################################################

_EOUSAGE_

    ;


my  $gene_counter = 0;

my ($gold_standard_gff3_superset_file, $genome_fasta_file, $evidence_gff3_files, $gene_gff3_files, $help, $VERBOSE);

&GetOptions ( 
			  "template_gff3_file=s" => \$gold_standard_gff3_superset_file,
			  "genome=s" => \$genome_fasta_file,
			  "alignment_gff3_files=s" => \$evidence_gff3_files,
			  "gene_gff3_files=s" => \$gene_gff3_files,
			  "help|h" => \$help,
			  "v" => \$VERBOSE,
			  
			  );

if ($help) { die $usage; }

unless ($gold_standard_gff3_superset_file && $genome_fasta_file && $gene_gff3_files) {
	die $usage;
}

my @gene_files = split (/,/, $gene_gff3_files);
foreach my $gene_file (@gene_files) {
	$gene_file =~ s/\s//g;
}

my @evidence_files = split (/,/, $evidence_gff3_files) if $evidence_gff3_files;
foreach my $ev_file (@evidence_files) { 
	$ev_file =~ s/\s//g;
}

 main: {
	 
	 my %curation_asmbl_ids;
	 
	 {
		 ## write the curation data set:
		 my $gene_obj_indexer_href = {};
		 
		 ## associate gene identifiers with contig id's.
		 my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gold_standard_gff3_superset_file, $gene_obj_indexer_href);
		 
		 foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
			 $curation_asmbl_ids{$asmbl_id} = 1;
			 &ensure_dir($asmbl_id);
			 open (my $curation_file, ">$asmbl_id/$asmbl_id.curation.cd") or die $!;
			 print STDERR "-writing curation file for $asmbl_id\n";
			 my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
			 foreach my $gene_id (@gene_ids) {
				 my $gene_obj = $gene_obj_indexer_href->{$gene_id};
				 $gene_counter++;
				 foreach my $exon ($gene_obj->get_exons()) {
					 if (my $cds = $exon->get_CDS_exon_obj()) {
						 my ($cds_end5, $cds_end3) = $cds->get_coords();
						 print $curation_file "$gene_counter $cds_end5 $cds_end3\n";
					 }
				 }
			 }
			 close $curation_file;
		 }
		 
	 }
	 
	 my %all_asmbl_ids;
	 my %gene_ev_files;
	 
	 {
		 ## write other genefinding files
		 use File::Basename;
		 foreach my $gene_file (@gene_files) {
			 my  $basename = basename($gene_file);
			 $gene_ev_files{$basename } = 1;
			 
			 ## get gene objects
			 ## write the curation data set:                                                                                                                                                     
			 my $gene_obj_indexer_href = {};			
			 ## associate gene identifiers with contig id's.                                                                                                                                     
			 my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gene_file, $gene_obj_indexer_href);            
			 foreach my $asmbl_id (keys %$contig_to_gene_list_href) {
				 print STDERR "-processing $asmbl_id gene prediction $gene_file evidence\n";
				 $all_asmbl_ids{$asmbl_id} = 1;
				 &ensure_dir($asmbl_id);
				 
				 my $sequence = &cdbyank_linear($asmbl_id, $genome_fasta_file);
				 open (my $seq_fh, ">$asmbl_id/$asmbl_id.seq") or die $!;
				 print $seq_fh ">$asmbl_id\n$sequence\n";
				 close $seq_fh;
				 
				 open (my $gene_fh, ">$asmbl_id/$asmbl_id.$basename") or die $!;
				 
				 my @genes = @{$contig_to_gene_list_href->{$asmbl_id}};
				 foreach my $gene (@genes) {
					 my  $gene_obj = $gene_obj_indexer_href->{$gene};
					 $gene_counter++;
					 my @exons = $gene_obj->get_exons();	 
					 $gene_obj->create_all_sequence_types(\$sequence);
					 my $protein = $gene_obj->get_protein_sequence();
					 for (my $i = 0; $i <= $#exons; $i++) {
						 # determine exon type:
						 my ($exon_end5, $exon_end3) = $exons[$i]->get_coords();
						 
						 my $exon_type = "internal";
						 if ($i == 0 && scalar(@exons) != 1 && $protein =~ /^M/) {
							 $exon_type = "initial";
						 }
						 elsif ($i == $#exons && scalar(@exons) != 1 && $protein =~ /\*$/) {
							 $exon_type = "terminal";
						 }
						 elsif (scalar @exons == 1 && $protein =~ /^M\w*\*$/) {
							 $exon_type = "single";
						 }
						 print $gene_fh "$gene_counter $exon_type $exon_end5 $exon_end3 1\n";
					 }
				 }
				 close $gene_fh;
			 }
		 }
		 
	 }
	 
	 
	 ## write the alignment file data:
	 my %alignment_types;
	 {
		 foreach my $alignment_gff3 (@evidence_files) {
			 my %data = &parse_gff3_alignment_file($alignment_gff3);
			 my $basename = basename($alignment_gff3);
			 $alignment_types{$basename} = 1;
			 
			 foreach my $asmbl_id (keys %data) {
				 print STDERR "-process alignments $basename for $asmbl_id\n";
				 open (my $fh, ">$asmbl_id/$asmbl_id.$basename") or die $!;
				 
				 my $asmbl_href = $data{$asmbl_id};
				 foreach my $target_acc (keys %$asmbl_href) {
					 my $target_aref = $asmbl_href->{$target_acc};
					 my @segments = @{$target_aref};
					 foreach my $segment (@segments) {
						 my ($end5, $end3, $per_id, $mlend, $mrend) = ($segment->{end5},
																	   $segment->{end3},
																	   $segment->{per_id},
																	   $segment->{mlend},
																	   $segment->{mrend});
						 print $fh "$end5 $end3 $mlend $mrend $target_acc $per_id $per_id\n";
					 }
				 }
				 
				 close $fh;
				 
			 }
		 }
	 }
	 
	 ## write the training directory files:
	 {
		 mkdir("my_train_dir") or die $!;
		 {
			 open (my $fh, ">my_train_dir/my_dir_list.txt") or die $!;
			 foreach my $asmbl_id (keys %curation_asmbl_ids) {
				 print $fh "../$asmbl_id\n";
			 }
			 close $fh;
		 }
		 
		 {
			 
			 open (my $fh, ">my_train_dir/my_evidence_list.txt") or die $!;
			 open (my $weightings_fh, ">my_train_dir/my_evidence_list_wghts.txt") or die $!;
			 
			 print $fh "curation.cd default curation\n";
			 foreach my $gene_ev_file (keys %gene_ev_files) {
				 print $fh "$gene_ev_file default geneprediction acc don coding start stop intron\n";
				 print $weightings_fh "$gene_ev_file default geneprediction acc 1.0 don 1.0 coding 1.0 start 1.0 stop 1.0 intron 1.0\n";
			 }
			 
			 foreach my $alignment_type (keys %alignment_types) {
				 print $fh "$alignment_type default homology acc don coding start stop intron\n";
				 print $weightings_fh "$alignment_type  default homology acc 1.0 don 1.0 coding 1.0 intron 1.0\n";
			 }
			 
			 
			 close $fh;
			 close $weightings_fh;
			 
		 }
		 
		 {
			 # ensure that every asmbl_id has an ev file:
			 foreach my $asmbl_id (keys %all_asmbl_ids) {
				 foreach my $gene_ev_file (keys %gene_ev_files) {
					 my $filename = "$asmbl_id/$asmbl_id.$gene_ev_file";
					 unless (-e $filename) {
						 system "touch $filename";
					 }
				 }
			 }
		 }
	 }
	 
	 exit(0);
 }

####
sub ensure_dir {
	my ($dirname) = @_;
	
	if (! -d $dirname) {
		mkdir $dirname or die "Error, cannot mkdir $dirname";
	}
	
	return;
}

####
sub parse_gff3_alignment_file {
	my ($filename) = @_;
	
	my %data;
	my %acc_counter;
	
	open (my $fh, $filename) or die "Error, cannot open file $filename";
	while (<$fh>) {
		chomp;
		my @x = split (/\t/);
		my ($asmbl_id, $source, $type, $lend, $rend, $per_id, $orient, $dot, $match_info) = @x;
		
		$match_info =~ /ID=([^;]+)/ or die "Error, cannot parse match ID from $match_info";
		my $match_ID = $1;
		
		$match_info =~ /Target=(\S+)\s(\d+)\s+(\d+)/ or die "Error, can't parse the target match data from $match_info";
		my $target_acc = $1;
		my $match_lend = $2;
		my $match_rend = $3;
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
		
		my $match_chain_ID = "$match_ID$;$target_acc";
		
		## make new match ID
		if (! $acc_counter{$match_chain_ID}) {
			$acc_counter{$match_chain_ID} = ++$acc_counter{$target_acc};
			
		}
		
		$target_acc .= "." . $acc_counter{$target_acc};
		
		push (@{$data{$asmbl_id}->{$target_acc}}, { end5 => $end5,
													end3 => $end3,
													per_id => $per_id,
													mlend => $match_lend,
													mrend => $match_rend, } );
	}
	
	return (%data);
}
