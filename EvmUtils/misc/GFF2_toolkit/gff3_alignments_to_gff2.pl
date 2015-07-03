#!/usr/bin/env perl 

use strict;
use warnings;


my $usage = "\nusage: $0 gff3_alignments_file TYPE\n\nTYPe=EST|PROTEIN\n\n";
my $gff3_file = $ARGV[0] or die $usage;
my $type = $ARGV[1] or die $usage; # EST or Protein

main: {
	my $prev_ID = "";
	my $txt = "";
	open (my $fh, $gff3_file) or die "Error, cannot open file $gff3_file";
	while (<$fh>) {
		/ID=([^;]+)/ or die "Error, cannot parse ID info from $_";
		my $curr_ID = $1;
		if ($curr_ID ne $prev_ID) {
			&process_txt($txt) if $txt;
			$txt = $_;
		}
		else {
			$txt .= $_;
		}
		
		$prev_ID = $curr_ID;
	}
	
	&process_txt($txt); # get last one.
	
	exit(0);
}


####
sub process_txt {
	my ($text) = @_;

	my @coords;
	my $Orient;
	my @structs;
	my $Target;
	my $Asmbl_id;
	my $Source;

	my @lines = split (/\n/, $text);
	foreach my $line (@lines) {
		my @x = split (/\t/, $line);
		my ($asmbl_id, $source, $trash, $lend, $rend, $per_id, $orient, $trash2, $match_info) = @x;
		
		my (@y) = split (/\s+/, $match_info);
		my $match_rend = pop @y;
		my $match_lend = pop @y;
		
		push (@coords, $lend, $rend);
		
		$match_info =~ /ID=([^;]+);Target=(\S+)/ or die "Error, cannot parse match and target info from $match_info";
		
		my $match_ID = $1;
		my $target = $2;
		$target =~ s/:/_/g;

		unless ($Target) {
			$Target = $target;
			$Orient = $orient;
			$Asmbl_id = $asmbl_id;
			$Source = $source;
		}
				
		push (@structs, { 
			lend => $lend,
			rend => $rend,
			per_id => $per_id,
			match_lend => $match_lend,
			match_rend => $match_rend,
		} );
		
			 
	}
	
	@coords = sort {$a<=>$b} @coords;
	my $lend = shift @coords;
	my $rend = pop @coords;

	print join ("\t", $Asmbl_id, $Source, "match", $lend, $rend, ".", $Orient, ".", "Target $type:$Target") . "\n";
	foreach my $struct (@structs) {
		print join ("\t", $Asmbl_id, $Source, "HSP", $struct->{lend}, $struct->{rend}, $struct->{per_id}, $Orient, ".", "Target $type:$Target " 
					. $struct->{match_lend} . " " . $struct->{match_rend} . " +") . "\n";
	}
	

	return;
}

