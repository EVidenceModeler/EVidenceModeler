#!/usr/bin/env perl

## Given the output from running pasa_asmbls_to_training_set.dbi
##   this pulls out the terminal CDS exons for use with
##   evidence_modeler

use strict;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use GFF_maker;

my $usage = "usage: $0 fastaFile gffFile\n\n";
my $fastaFile = $ARGV[0] or die $usage;
my $gffFile = $ARGV[1] or die $usage;


my %asmbls_with_termini; #every asmbl that has a terminal exon of interest
#  structure:   $asmbls{$pasa_acc}->{prime5} = 1;  or ->{prime3} = 1

## determine which asmbls are of interest:
open (my $FILE, $fastaFile) or die $!;
while (<$FILE>) {
    if (/^>/) {
        #print;
        my ($acc, $contig, $coords, $segs, $type, $len) = split (/\s+/);
        $acc =~ s/>//;
        $type =~ s/type://;
        my $have_5prime = 0;
        my $have_3prime = 0;
        
        if ($type eq "complete") {
            $have_5prime = 1;
            $have_3prime = 1;
        } 
        
        elsif ($type eq "3prime_partial") {
            $have_5prime = 1;
        }
        
        elsif ($type eq "5prime_partial") {
            $have_3prime = 1;
        }
        
        $asmbls_with_termini{$acc} = { prime5 => $have_5prime,
                                       prime3 => $have_3prime };
    }
}
close $FILE;

## Extract terminal exons from GFF file:

my %pasa_acc_to_asmbl_id;

my @pasa_coordsets;
my $EV_TYPE = undef;

open ($FILE, $gffFile) or die $!;
my $prev_pasa_acc = "";
while (<$FILE>) {
    #print;
    chomp;
    my @x = split (/\t/);
    if ($x[2] eq "CDS") {
        my ($asmbl_id, $ev_type, $end5, $end3, $strand, $phase, $acc_info) = ($x[0], $x[1], $x[3], $x[4], $x[6], $x[7], $x[8]);

        if (! defined($EV_TYPE)) {
            $EV_TYPE = $ev_type;
        }
        
        if ($strand eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }
        
        $acc_info =~ /LongOrf\.(asmbl_\d+)/;
        my $pasa_acc = $1 or die "Error, cannot get pasa_acc from $acc_info ";
        
        $pasa_acc_to_asmbl_id{$pasa_acc} = $asmbl_id;
        
        if ($pasa_acc ne $prev_pasa_acc) {
            if (@pasa_coordsets) {
                &process_pasa_coordsets($prev_pasa_acc, \@pasa_coordsets);
                @pasa_coordsets = ();
            }
        }
        if ($asmbls_with_termini{$pasa_acc}) {
            push (@pasa_coordsets, [$end5, $end3, $phase, $strand]);
        }
        $prev_pasa_acc = $pasa_acc;
    }
}

## process the last guy
if (@pasa_coordsets) {
    &process_pasa_coordsets($prev_pasa_acc, \@pasa_coordsets);
    @pasa_coordsets = ();
}

exit(0);




####
sub process_pasa_coordsets {
    my ($pasa_acc, $coordsets_aref) = @_;
    
    my $asmbl_id = $pasa_acc_to_asmbl_id{$pasa_acc};
    
    my @coordsets = sort {$a->[0]<=>$b->[0]} @$coordsets_aref;

    my $orient = $coordsets[0]->[3];
    
    ## check 5' exon:
    if ($asmbls_with_termini{$pasa_acc}->{prime5}) {
        
        my $exon;
        if ($orient eq '+') {
            $exon = shift @coordsets;
        } else {
            $exon = pop @coordsets;
        }
        
        my ($end5, $end3, $phase, $orient) = @$exon;
        my $type = "initial";
        if ( (! @coordsets) && $asmbls_with_termini{$pasa_acc}->{prime3}) {
            $type = "single";
        }

        my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
        
        print GFF_maker::get_GFF_line( { seq_id => $asmbl_id,
                                         source => $EV_TYPE,
                                         attributes => "ID=$pasa_acc-$end5-$end3",
                                         lend => $lend,
                                         rend => $rend,
                                         strand => $orient,
                                         type => $type,
                                         phase => $phase,
                                     }
                                       );
        
    }
    
    if ($asmbls_with_termini{$pasa_acc}->{prime3}) {
        my $exon;
        if ($orient eq '+') {
            $exon = pop @coordsets;
        } else {
            $exon = shift @coordsets;
        }
        
        if ($exon) {
            # if only a single exon and has both termini, then would only have one exon used earlier.
            my ($end5, $end3, $phase, $orient) = @$exon;
            
            my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
            
            print GFF_maker::get_GFF_line( { seq_id => $asmbl_id,
                                             source => $EV_TYPE,
                                             attributes => "ID=$pasa_acc-$end5-$end3",
                                             lend => $lend,
                                             rend => $rend,
                                             strand => $orient,
                                             type => "terminal",
                                             phase => $phase,
                                         }
                                           );
        }
    }
}




