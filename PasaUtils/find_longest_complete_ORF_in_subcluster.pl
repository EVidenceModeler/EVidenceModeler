#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 subcluster_listings trainingSetCandidates.fasta [min_ORF_length=300] \n\n";

my $subcluster_file = $ARGV[0] or die $usage;
my $training_candidates = $ARGV[1] or die $usage;
my $cutoff_size = $ARGV[2] || 300;

my %cdna_to_subcluster;
my %subcluster_to_cdna_list;

## populate the subcluster info:
{
    open (my $fh, $subcluster_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($subcluster_id, $cdna_acc) = split (/\s+/);
        $cdna_to_subcluster{$cdna_acc} = $subcluster_id;
        
        my $cdna_list_aref = $subcluster_to_cdna_list{$subcluster_id};
        unless (ref $cdna_list_aref) {
            $cdna_list_aref = $subcluster_to_cdna_list{$subcluster_id} = [];
        }

        push (@$cdna_list_aref, $cdna_acc);
    }
    close $fh;
}


my %complete_asmbls_meeting_criteria_to_length;
my %subcluster_containing_complete_entry;
## parse the training file, retrieve only those complete entries above cutoff
{
    open (my $fh, $training_candidates) or die $!;
    while (<$fh>) {
        if (/>/) {
            chomp;
            # ie. >asmbl_1 contig:1 coords:3631-5832 num_segs:6 type:5prime_partial len:472

            my ($acc, $contig, $coords, $num_segs, $type, $len) = split (/\s+/);
            $acc =~ s/>//;
            $len =~ s/len://;
            $type =~ s/type://;
            if ($type eq 'complete' && $len >= $cutoff_size) {
                $complete_asmbls_meeting_criteria_to_length{$acc} = $len;
                my $subcluster = $cdna_to_subcluster{$acc};
                $subcluster_containing_complete_entry{$subcluster} = 1;
            }
        }
    }
    close $fh;
}


## pull out the longest complete entry per subcluster:
foreach my $subcluster (keys %subcluster_containing_complete_entry) {
    
    my @accs = @{$subcluster_to_cdna_list{$subcluster}};

    my @entries;
    foreach my $acc (@accs) {
        
        if (my $length = $complete_asmbls_meeting_criteria_to_length{$acc}) {
            push (@entries, { acc => $acc,
                              length => $length,
                          } );
        }
    }
    
    @entries = reverse sort {$a->{length}<=>$b->{length}} @entries;
    
    # print Dumper (\@entries);

    my $longest_entry = shift @entries;
    
    print $longest_entry->{acc} . "\t" . $longest_entry->{length} . "\n";


}


exit(0);



    



                
