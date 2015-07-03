#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

#####################################
#
#  -A    ab initio predictions
#
#  -P    protein alignments
#
#  -T    transcript alignments
#
#####################################

_EOUSAGE_

	;



my $help_flag;
my $ab_initios_file;
my $proteins_file;
my $transcripts_file;


&GetOptions ( 'h' => \$help_flag,
			  'A=s' => \$ab_initios_file,
			  'P=s' => \$proteins_file,
			  'T=s' => \$transcripts_file,
			  
	);

if ($help_flag) {
	die $usage;
}

unless ($ab_initios_file || $proteins_file || $transcripts_file) {
	die $usage;
}


main: {
	my %class_to_sources;
	if ($ab_initios_file) {
		my @ab_initio_sources = &get_source_types($ab_initios_file);
		push (@{$class_to_sources{"ABINITIO_PREDICTION"}}, @ab_initio_sources);
	}
	if ($proteins_file) {
		my @protein_sources = &get_source_types($proteins_file);
		push (@{$class_to_sources{"PROTEIN"}}, @protein_sources);
	}
	if ($transcripts_file) {
		my @transcript_sources = &get_source_types($transcripts_file);
		push (@{$class_to_sources{"TRANSCRIPT"}}, @transcript_sources);
	}

	foreach my $source_type (sort keys %class_to_sources) {
		my @sources = @{$class_to_sources{$source_type}};
		foreach my $source (@sources) {
			print "$source_type\t$source\t1\n";
		}
	}

	exit(0);
}



####
sub get_source_types {
	my ($file) = @_;

	my %source_types;

	open (my $fh, $file) or die "Error, cannot open file $file";
	while (<$fh>) {
		chomp;
		unless (/\w/) { next; }
		my @x = split (/\t/);
		my $source_type = $x[1];
		$source_types{$source_type} =1 ;
	}
	close $fh;

	return(keys %source_types);
}


