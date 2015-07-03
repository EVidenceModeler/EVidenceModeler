package EvidenceTierParser;

use strict;
use warnings;
use Carp;
use EvidenceTier;
use EvidenceFeature;


sub new {
    my $packagename = shift;
    my $params_href = shift;

    my $self = { max_rows_per_tier => undef,
                 range_start => undef,
                 range_end => undef,
    };
    
    if (ref $params_href eq "HASH") {
        foreach my $param (keys %$params_href) {
            if (exists $self->{$param}) {
                $self->{$param} = $params_href->{$param};
            }
            else {
                confess "Error, attribute $param isn't recognized";
            }
        }
    }

    bless ($self, $packagename);

    return ($self);
}





sub parse_evidence_tiers {
    my $self = shift;
    my ($filename) = @_;
    ## requires GFF-format (GFF3)

    srand();
    my %rand_IDs; # in case an ID isn't there, we'll create one.
    
    my ($range_start, $range_end) = ($self->{range_start},
                                     $self->{range_end} );
    
    
    my %ID_to_feature;  
    my %child_to_parent;
    
    open (my $fh, $filename) or confess "cannot open file: $filename";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        my @x = split (/\t/);
        
        if (scalar @x != 9) {
            print STDERR "error, something wrong with GFF entry: $_\nskipping it.\n";
            next;
        }

        my ($genome_contig, $source, $type, $lend, $rend, $score, $orient, $phase, $feature_info) = @x;
        
        if (defined ($range_start) && defined($range_end)) {
            unless ($lend <= $range_end && $rend >= $range_start) { next; } # no overlap
            if ($lend < $range_start) {
                $lend = $range_start; #truncate it.
            }
            if ($rend > $range_end) {
                $rend = $range_end;
            }
        }
        
        my $id = "";
        my $parent = "";
        if ($feature_info =~ /;?ID=([^\;]+)/) {
            $id = $1;
        }
        else {
            do {
                $id = int(rand(10000000));
            } while ($rand_IDs{$id});
            
        }
        if ($feature_info =~ /;?Parent=([^\;]+)/i) {
            $parent = $1;
        }

        if (my $feature = $ID_to_feature{$id}) {
            # got a feature already, add the segment:
            $feature->add_segment($lend, $rend, $orient);
        }
        else {
            ## new feature:
            my $feature = new EvidenceFeature($id, $type, $source);
            $feature->add_segment($lend, $rend, $orient);
            $ID_to_feature{$id} = $feature;
        }
        
        if ($parent) {
            $child_to_parent{$id} = $parent;
        }
        
    }

    

    ## process all the parent/child relationships:
    foreach my $child_id (keys %child_to_parent) {
        my $parent_id = $child_to_parent{$child_id};
        
        my $parent_obj = $ID_to_feature{$parent_id} or confess "error, no parent retrieved based on ID: $parent_id";
        my $child_obj = $ID_to_feature{$child_id} or confess "error, no child retrieved based on ID: $child_id" ;
        
        $parent_obj->add_child($child_obj);
    }

    ## Convert feature list into an EvidenceTier object
    
    # only those features that lack parents are top level features that should be tiered:
    my %feature_source_to_ID_list;
    foreach my $feature_obj (values %ID_to_feature) {
        unless ($feature_obj->{parent_ID}) {
            ## get the feature source
            my $feature_source = $feature_obj->{feature_source};
            my $id = $feature_obj->{ID};
            my $list_ref = $feature_source_to_ID_list{$feature_source};
            unless (ref $list_ref) {
                $list_ref = $feature_source_to_ID_list{$feature_source} = [];
            }
            push (@$list_ref, $id);
        }
    }


    my @tiers;
    
    foreach my $feature_source (keys %feature_source_to_ID_list) {
        my $id_list_aref = $feature_source_to_ID_list{$feature_source};
        my @feature_objs;
        foreach my $id (@$id_list_aref) {
            my $feature = $ID_to_feature{$id};
            push (@feature_objs, $feature);
        }

        my $evidence_tier_obj = new EvidenceTier( { max_rows_per_tier => $self->{max_rows_per_tier} });
        $evidence_tier_obj->tier_features(@feature_objs);
        push (@tiers, $evidence_tier_obj);
    }

    return (@tiers);
    
}
    
1; #EOM
