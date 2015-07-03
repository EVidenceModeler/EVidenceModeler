package EvidenceTier;

## EvidenceTier provides a collection of features, each feature tiered into a given row.

use strict;
use warnings;
use Carp;
use EvidenceFeature;

sub new {
    my $packagename = shift;
    my $params_href = shift;
    
    my $self = {
        top_rows => [], # plus strand features
        bottom_rows => [], # minus strand features
        max_rows_per_tier => undef,
    };
    
    if (ref $params_href eq 'HASH') {
        foreach my $param (keys %$params_href) {
            if (exists $self->{$param}) {
                $self->{$param} = $params_href->{$param};
            }
            else {
                confess "Error, don't recognize param $param";
            }
        }
    }
    
    bless ($self, $packagename);
    
    return ($self);
}

####
sub tier_features {
    my $self = shift;

    my @feature_objs = @_;

    ## sort by feature length:
    #@feature_objs = reverse sort { $a->{span_length} <=> $b->{span_length} } @feature_objs;

    my $top_rows_aref = $self->{top_rows};
    my $bottom_rows_aref = $self->{bottom_rows};
    
    foreach my $feature_obj (@feature_objs) {

        my $orient = $feature_obj->{orient};
        
        if ($orient eq '+') {
            $self->_place_feature($top_rows_aref, $feature_obj);
        }
        else {
            $self->_place_feature($bottom_rows_aref, $feature_obj);
        }
    }

    return;
}

####
sub get_top_strand_rows {
    my $self = shift;
    return (@{$self->{top_rows}});
}

####
sub get_bottom_strand_rows {
    my $self = shift;
    return (@{$self->{bottom_rows}});
}

####
sub toString {
    my $self = shift;

    my @top_strand_rows = $self->get_top_strand_rows();
    my @bottom_strand_rows = $self->get_bottom_strand_rows();

    my $text = "Features:\n";
    my $row_count = 0;
    foreach my $row (@top_strand_rows, @bottom_strand_rows) {
        $row_count++;
        foreach my $feature (@$row) {
            $text .= "row($row_count)" .  $feature->toString();
        }
    }
    return ($text);
}

######### Private methods


####
sub _place_feature {
    my $self = shift;
    my ($rows_aref, $feature) = @_;

    my $feature_placed_flag = 0;
    
    foreach my $row (@$rows_aref) {
        unless (&_has_overlapping_feature($row, $feature)) {
            push (@$row, $feature);
            $feature_placed_flag = 1;
            last;
        }
    }

    unless ($feature_placed_flag) {
        # add this as a new row:
        if (my $max_row_count = $self->{max_rows_per_tier}) {
            my $num_current_rows = scalar (@$rows_aref);
            if ($num_current_rows < $max_row_count) {
                push (@$rows_aref, [ $feature ] );
            }
        } else {
            ## row numbers unlimited
            push (@$rows_aref, [ $feature ] );
        } 
    }
    
    return;
}


####
sub _has_overlapping_feature {
    my ($row, $feature) = @_;
    foreach my $other_feature (@$row) {
        my ($lend, $rend) = $other_feature->get_coords();
        if ($feature->overlaps_coords($lend, $rend)) {
            return (1);
        }
    }
    return (0); #no overlapping feature found
}



1; #EOM



    
