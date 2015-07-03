#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 pasa_coordsfile genes_coordsfile genome_fastafile\n\n";

my $pasa = $ARGV[0] or die $usage;
my $genes = $ARGV[1] or die $usage;
my $genome = $ARGV[2] or die $usage; 

my %pasa_text;
my %gene_text;

&parse_file($pasa, \%pasa_text);
&parse_file($genes, \%gene_text);

open (my $fh, $genome) or die $!;
while (<$fh>) {
    if (/>(\S+)/) {
        my $asmbl = $1;
        unless (exists $pasa_text{$asmbl}) {
            $pasa_text{$asmbl} = "";
        }
        unless (exists $gene_text{$asmbl}) {
            $gene_text{$asmbl} = "";
        }
        
        print "$pasa_text{$asmbl}\n$gene_text{$asmbl}\n";
    }
}

exit(0);


####
sub parse_file {
    my ($filename, $data_href) = @_;

    my $current_asmbl = "";
    
    open (my $fh, $filename) or die $!;
    while (my $line = <$fh>) {
        if ($line =~ /\w/) {
            my @x = split (/\s+/, $line);
            my $asmbl = $x[0];
            $current_asmbl = $asmbl;
            
            $data_href->{$current_asmbl} .= $line;
        }
        else {
            # add blank line to current entry
            $data_href->{$current_asmbl} .= $line;
        }
    }
    close $fh;
}


