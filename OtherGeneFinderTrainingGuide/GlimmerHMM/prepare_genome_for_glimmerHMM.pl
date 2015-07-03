#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Cwd;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 genomeMultiFasta GlimmerTrainingDirectory\n\n";

my $genome = $ARGV[0] or die $usage;
my $trainDir = $ARGV[1] or die $usage;

my $outdir = cwd() . "/glimmerHMM.$$.run";

mkdir ($outdir) or die "Error, cannot mkdir $outdir";

my $fasta_reader = new Fasta_reader($genome);
while (my $seq_obj = $fasta_reader->next()) {
	
	my $acc = $seq_obj->get_accession();
	my $seq = $seq_obj->get_FASTA_format();

	my $outfile = "$outDir/$acc";

	open (my $ofh, ">$outfile") or die "Error, cannot write to file [$outfile] ";
	print $ofh $seq;
	close $ofh;

	my $cmd = "/home/radon01/bhaas/bio_tools/GlimmerHMM/bin/glimmerhmm $outfile $trainDir > $outfile.glimmerHMM";

	print "$cmd\n";
}

exit(0);


