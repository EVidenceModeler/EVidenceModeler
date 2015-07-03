#!/usr/bin/env perl

use strict;
use warnings;
use URI::Escape;

my $usage = "usage: $0 btab_file ev_type=[N|P]\n\n";

my $btab_file = $ARGV[0] or die $usage;
my $ev_type = $ARGV[1] or die $usage;

my $match_counter = 0;

if ($ev_type eq "N") {
	$ev_type = "cDNA_match";
}
elsif ($ev_type eq "P") {
	$ev_type = "nucleotide_to_protein_match";
}
else {
	die $usage;
}

main: {
	my @matches;
	my $prev_match_ID = undef;
	open (my $fh, $btab_file) or die "Error, cannot open file $btab_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		if (/^\#/) { next; }
		my $line = $_;
		my @x = split (/\t/);
		my $match_id = $x[13];
		
		if ( (! defined($prev_match_ID)) || $prev_match_ID ne $match_id) {
			&process_matches(@matches) if @matches;
			@matches = ($line);
			$prev_match_ID = $match_id;
		}
		else {
			push (@matches, $line);
		}
	}
	
	&process_matches(@matches); # get last ones.
	
	
	exit(0);
}

####
sub process_matches {
	my @matches = @_;
	

	my $orient;
	
	# find orientation
	foreach my $match (@matches) {
		my @x = split (/\t/, $match);
		my ($genome_end5, $genome_end3) = ($x[6], $x[7]);
		unless ($orient) {
			if ($genome_end5 < $genome_end3) {
				$orient = '+';
			}
			elsif ($genome_end5 > $genome_end3) {
				$orient = '-';
			}
		}

		if ($orient) { last; }
		
	}
	
	$match_counter++;	
	foreach my $match (@matches) {
		my @x = split (/\t/, $match);
		my $genome_contig = $x[0];
		my $prot_name = $x[5] . " " . $x[15];
		$prot_name = uri_escape($prot_name);
		my ($genome_lend, $genome_rend) = sort {$a<=>$b} ($x[6], $x[7]);
		my ($hit_lend, $hit_rend) = ($x[8], $x[9]);
		my $per_id = $x[10];
		
		
		print join ("\t", $genome_contig, "AAT", $ev_type, $genome_lend, $genome_rend, $per_id, $orient, ".", "ID=match_$match_counter;Target=$prot_name $hit_lend $hit_rend") . "\n";

	}
	print "\n"; # spacer
	
	return;
}
