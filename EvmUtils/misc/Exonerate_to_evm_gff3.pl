#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 exonerate_output.txt\n\n";

my $exonerate_file = $ARGV[0] or die $usage;


my $ALIGNMENT_COUNTER = 0;

main: {
	
	my @alignments;
		
	open (my $fh, $exonerate_file) or die "Error, cannot open file $exonerate_file";
	while (<$fh>) {
		chomp;
		if (/END OF GFF DUMP/) { 
			if (@alignments) {
				&process_alignments(@alignments);
				@alignments = ();
			}
		}
				
		if (/^\#/) { next; }
		
		my @x = split(/\t/);
		unless (scalar (@x) >= 9) { next; }
		
		my ($contig, $prog_name, $feat_type, $lend, $rend, $score, $orient, $phase, $annots, @rest) = @x;
		
		if ($feat_type eq "gene") {
			
			my $query_name = "";
			if ($annots =~ /sequence (\S+)/) {
				$query_name = $1;
			}
			else {
				die "Error, cannot parse sequence identity from [$annots] in line $_";
			}

			my $alignment_struct = { query => $query_name,
									 contig => $contig,
									 orient => $orient,
									 segments => [],
			};
			
			push (@alignments, $alignment_struct);


		}
		elsif ($feat_type eq "exon") {
			
			my $curr_alignment_struct = $alignments[$#alignments];
			
			my $exon_seg = { lend => $lend,
							 rend => $rend,
			};

			my $seg_list = $curr_alignment_struct->{segments};
			push (@$seg_list, $exon_seg);
		}
	}


	exit(0);
}


####
sub process_alignments {
	my (@alignments) = @_;

	$ALIGNMENT_COUNTER++;

	foreach my $alignment (@alignments) {
		my $query = $alignment->{query};
		my $contig = $alignment->{contig};
		my $orient = $alignment->{orient};

		my $segs_aref = $alignment->{segments};

		foreach my $seg (@$segs_aref) {
			my $lend = $seg->{lend};
			my $rend = $seg->{rend};

			print join("\t", $contig, "exonerate", "match", $lend, $rend, ".", $orient, ".",
					   "ID=exonerate.$ALIGNMENT_COUNTER;Target=$query") . "\n";
		}
		print "\n";
	}
	
	return;
}
