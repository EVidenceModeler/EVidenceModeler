#!/usr/bin/env perl

use strict;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;

$|++;

my $usage = "usage: $0 genome_multiFastaDB\n\n";
my $fastaFile = $ARGV[0] or die $usage;

my $baseDir = "genome_distributed";
mkdir ($baseDir);
chmod (0777, $baseDir);


my $fasta_reader = new Fasta_reader($fastaFile);
while (my $record = $fasta_reader->next()) {
    my $accession = $record->get_accession();
    my $fastaSeq = $record->get_FASTA_format();
    my $dir = "$baseDir/$accession";
    
    print "Processing $dir\n";

    if (! -d $dir) {
	mkdir ($dir) or die "Error, couldn't mkdir $dir\n";
    }
    
    open (FILE, ">$dir/$accession.seq") or die "Cannot write file $dir/$accession.seq";
    print FILE $fastaSeq;
    close FILE;
}

exit;
